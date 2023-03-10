import io
import os
import sys
import logging
from contextlib import ExitStack
from typing import (
    List,
    Optional,
    BinaryIO,
    TextIO,
    Any,
    Tuple,
    Dict,
    Sequence,
    Union,
    TYPE_CHECKING,
    Iterator,
)
from abc import ABC, abstractmethod
import multiprocessing
from multiprocessing.connection import Connection
from pathlib import Path
import traceback

import dnaio

from .utils import (
    Progress,
    DummyProgress,
    open_raise_limit,
    xopen_rb_raise_limit,
)
from .modifiers import (
    SingleEndModifier,
    PairedEndModifier,
    PairedEndModifierWrapper,
    ModificationInfo,
)
from .report import Statistics
from .predicates import (
    TooShort,
    TooLong,
    TooManyN,
    TooManyExpectedErrors,
    CasavaFiltered,
    DiscardTrimmed,
    DiscardUntrimmed,
    Predicate,
)
from .steps import (
    SingleEndSink,
    PairedEndSink,
    SingleEndFilter,
    PairedEndFilter,
    Demultiplexer,
    PairedDemultiplexer,
    CombinatorialDemultiplexer,
    RestFileWriter,
    WildcardFileWriter,
    InfoFileWriter,
    SingleEndStep,
    PairedSingleEndStep,
)

logger = logging.getLogger()

mpctx = multiprocessing.get_context()

# See https://github.com/python/typeshed/issues/9860
if TYPE_CHECKING:
    mpctx_Process = multiprocessing.Process
else:
    mpctx_Process = mpctx.Process


class InputFiles:
    def __init__(
        self,
        *files: BinaryIO,
        interleaved: bool = False,
    ):
        self._files = files
        self.interleaved = interleaved
        for f in self._files:
            assert f is not None

    def open(self):
        return dnaio.open(*self._files, interleaved=self.interleaved, mode="r")

    def close(self) -> None:
        for file in self._files:
            file.close()


class InputPaths:
    def __init__(self, *paths: str, interleaved: bool = False):
        self.paths = paths
        self.interleaved = interleaved

    def open(self) -> InputFiles:
        files = [xopen_rb_raise_limit(path) for path in self.paths]
        return InputFiles(*files, interleaved=self.interleaved)


class OutputFiles:
    """
    The attributes are either None or open file-like objects except for demultiplex_out
    and demultiplex_out2, which are dictionaries that map an adapter name
    to a file-like object.
    """

    def __init__(
        self,
        out: Optional[BinaryIO] = None,
        out2: Optional[BinaryIO] = None,
        untrimmed: Optional[BinaryIO] = None,
        untrimmed2: Optional[BinaryIO] = None,
        too_short: Optional[BinaryIO] = None,
        too_short2: Optional[BinaryIO] = None,
        too_long: Optional[BinaryIO] = None,
        too_long2: Optional[BinaryIO] = None,
        info: Optional[BinaryIO] = None,
        rest: Optional[BinaryIO] = None,
        wildcard: Optional[BinaryIO] = None,
        demultiplex_out: Optional[Dict[str, BinaryIO]] = None,
        demultiplex_out2: Optional[Dict[str, BinaryIO]] = None,
        combinatorial_out: Optional[Dict[Tuple[str, str], BinaryIO]] = None,
        combinatorial_out2: Optional[Dict[Tuple[str, str], BinaryIO]] = None,
        force_fasta: Optional[bool] = None,
    ):
        self.out = out
        self.out2 = out2
        self.untrimmed = untrimmed
        self.untrimmed2 = untrimmed2
        self.too_short = too_short
        self.too_short2 = too_short2
        self.too_long = too_long
        self.too_long2 = too_long2
        self.info = info
        self.rest = rest
        self.wildcard = wildcard
        self.demultiplex_out = demultiplex_out
        self.demultiplex_out2 = demultiplex_out2
        self.combinatorial_out = combinatorial_out
        self.combinatorial_out2 = combinatorial_out2
        self.force_fasta = force_fasta

    def __iter__(self):
        for f in [
            self.out,
            self.out2,
            self.untrimmed,
            self.untrimmed2,
            self.too_short,
            self.too_short2,
            self.too_long,
            self.too_long2,
            self.info,
            self.rest,
            self.wildcard,
        ]:
            if f is not None:
                yield f
        for outs in (
            self.demultiplex_out,
            self.demultiplex_out2,
            self.combinatorial_out,
            self.combinatorial_out2,
        ):
            if outs is not None:
                for f in outs.values():
                    assert f is not None
                    yield f

    def as_bytesio(self) -> "OutputFiles":
        """
        Create a new OutputFiles instance that has BytesIO instances for each non-None output file
        """
        result = OutputFiles(force_fasta=self.force_fasta)
        for attr in (
            "out",
            "out2",
            "untrimmed",
            "untrimmed2",
            "too_short",
            "too_short2",
            "too_long",
            "too_long2",
            "info",
            "rest",
            "wildcard",
        ):
            if getattr(self, attr) is not None:
                setattr(result, attr, io.BytesIO())
        for attr in (
            "demultiplex_out",
            "demultiplex_out2",
            "combinatorial_out",
            "combinatorial_out2",
        ):
            if getattr(self, attr) is not None:
                setattr(result, attr, dict())
                for k, v in getattr(self, attr).items():
                    getattr(result, attr)[k] = io.BytesIO()
        return result

    def close(self) -> None:
        """Close all output files that are not stdout"""
        for f in self:
            if f is sys.stdout or f is sys.stdout.buffer:
                continue
            f.close()


class Pipeline(ABC):
    """
    Processing pipeline that loops over reads and applies modifiers and filters
    """

    n_adapters = 0
    paired = False

    def __init__(self) -> None:
        self._reader: Any = None
        self._steps: List[Any] = []
        self._infiles: Optional[InputFiles] = None
        self._outfiles: Optional[OutputFiles] = None
        self._demultiplexer = None
        self._textiowrappers: List[TextIO] = []

        # Filter settings
        self._minimum_length = None
        self._maximum_length = None
        self.max_n = None
        self.max_expected_errors = None
        self.discard_casava = False
        self.discard_trimmed = False
        self.discard_untrimmed = False

    def connect_io(self, infiles: InputFiles, outfiles: OutputFiles) -> None:
        self._infiles = infiles
        self._reader = infiles.open()
        self._set_output(outfiles)

    @abstractmethod
    def _open_writer(
        self,
        file: BinaryIO,
        file2: Optional[BinaryIO] = None,
        force_fasta: Optional[bool] = None,
    ):
        pass

    def _set_output(self, outfiles: OutputFiles) -> None:
        self._steps = []
        self._textiowrappers = []
        self._outfiles = outfiles

        for step_class, outfile in (
            (RestFileWriter, outfiles.rest),
            (InfoFileWriter, outfiles.info),
            (WildcardFileWriter, outfiles.wildcard),
        ):
            if outfile:
                textiowrapper = io.TextIOWrapper(outfile)
                self._textiowrappers.append(textiowrapper)
                self._steps.append(
                    self._wrap_single_end_step(step_class(textiowrapper))
                )

        # minimum length and maximum length
        for lengths, file1, file2, predicate_class in (
            (self._minimum_length, outfiles.too_short, outfiles.too_short2, TooShort),
            (self._maximum_length, outfiles.too_long, outfiles.too_long2, TooLong),
        ):
            if lengths is None:
                continue
            writer = self._open_writer(file1, file2) if file1 else None
            f1 = predicate_class(lengths[0]) if lengths[0] is not None else None
            if len(lengths) == 2 and lengths[1] is not None:
                f2 = predicate_class(lengths[1])
            else:
                f2 = None
            self._steps.append(
                self._make_filter(predicate1=f1, predicate2=f2, writer=writer)
            )

        if self.max_n is not None:
            f1 = f2 = TooManyN(self.max_n)
            self._steps.append(self._make_filter(f1, f2, None))

        if self.max_expected_errors is not None:
            if not self._reader.delivers_qualities:
                logger.warning(
                    "Ignoring option --max-ee because input does not contain quality values"
                )
            else:
                f1 = f2 = TooManyExpectedErrors(self.max_expected_errors)
                self._steps.append(self._make_filter(f1, f2, None))

        if self.discard_casava:
            f1 = f2 = CasavaFiltered()
            self._steps.append(self._make_filter(f1, f2, None))

        if (
            int(self.discard_trimmed)
            + int(self.discard_untrimmed)
            + int(outfiles.untrimmed is not None)
            > 1
        ):
            raise ValueError(
                "discard_trimmed, discard_untrimmed and outfiles.untrimmed must not "
                "be set simultaneously"
            )

        if (
            outfiles.demultiplex_out is not None
            or outfiles.combinatorial_out is not None
        ):
            self._demultiplexer = self._create_demultiplexer(outfiles)
            self._steps.append(self._demultiplexer)
        else:
            # Some special handling to allow overriding the wrapper for
            # --discard-untrimmed/--untrimmed-(paired-)output

            # Set up the remaining filters to deal with --discard-trimmed,
            # --discard-untrimmed and --untrimmed-output. These options
            # are mutually exclusive in order to avoid brain damage.
            if self.discard_trimmed:
                self._steps.append(
                    self._make_filter(DiscardTrimmed(), DiscardTrimmed(), None)
                )
            elif self.discard_untrimmed:
                self._steps.append(self._make_untrimmed_filter(None))
            elif outfiles.untrimmed:
                untrimmed_writer = self._open_writer(
                    outfiles.untrimmed, outfiles.untrimmed2
                )
                self._steps.append(self._make_untrimmed_filter(untrimmed_writer))

            self._steps.append(self._final_filter(outfiles))
        for i, step in enumerate(self._steps, 1):
            logger.debug("Pipeline step %d: %s", i, step)

    def flush(self) -> None:
        for f in self._textiowrappers:
            f.flush()
        assert self._outfiles is not None
        for f in self._outfiles:
            f.flush()

    def close(self) -> None:
        self._close_input()
        self._close_output()

    def _close_input(self) -> None:
        self._reader.close()
        if self._infiles is not None:
            self._infiles.close()

    def _close_output(self) -> None:
        for f in self._textiowrappers:
            f.close()
        # Closing a TextIOWrapper also closes the underlying file, so
        # this closes some files a second time.
        assert self._outfiles is not None
        self._outfiles.close()

    @property
    def uses_qualities(self) -> bool:
        assert self._reader is not None
        return self._reader.delivers_qualities

    @abstractmethod
    def process_reads(
        self, progress: Optional[Progress] = None
    ) -> Tuple[int, int, Optional[int]]:
        pass

    @abstractmethod
    def _make_filter(
        self, predicate1: Optional[Predicate], predicate2: Optional[Predicate], writer
    ):
        pass

    @abstractmethod
    def _make_untrimmed_filter(self, writer):
        pass

    @abstractmethod
    def _final_filter(self, outfiles: OutputFiles):
        pass

    @abstractmethod
    def _create_demultiplexer(self, outfiles: OutputFiles):
        pass

    @abstractmethod
    def _wrap_single_end_step(self, step: SingleEndStep):
        pass

    def run(
        self,
        inpaths: InputPaths,
        outfiles: OutputFiles,
        cores: int,
        progress: Union[bool, Progress, None] = None,
        buffer_size: Optional[int] = None,
    ) -> Statistics:
        """
        Run this pipeline.

        This uses a SerialPipelineRunner if cores is 1 and a ParallelPipelineRunner otherwise.

        Args:
            inpaths:
            outfiles:
            cores: number of cores to run the pipeline on (this is actually the number of worker
                processes, there will be one extra process for reading the input file(s))
            progress: Set to False for no progress bar, True for Cutadapt’s default progress bar,
                or use an object that supports .update() and .close() (e.g. a tqdm instance)
            buffer_size: Forwarded to `ParallelPipelineRunner()`. Ignored if cores is 1.

        Returns:
            A Statistics object
        """
        with self.make_runner(
            inpaths, outfiles, cores, progress, buffer_size
        ) as runner:
            statistics = runner.run()
        return statistics

    def make_runner(
        self,
        inpaths: InputPaths,
        outfiles: OutputFiles,
        cores: int,
        progress: Union[bool, Progress, None] = None,
        buffer_size: Optional[int] = None,
    ):
        if progress is None or progress is False:
            progress = DummyProgress()
        elif progress is True:
            progress = Progress()
        if cores > 1:
            return ParallelPipelineRunner(
                self,
                inpaths,
                outfiles,
                progress,
                n_workers=cores,
                buffer_size=buffer_size,
            )
        else:
            return SerialPipelineRunner(self, inpaths.open(), outfiles, progress)


class SingleEndPipeline(Pipeline):
    """
    Processing pipeline for single-end reads
    """

    def __init__(self, modifiers: List[SingleEndModifier]):
        super().__init__()
        self._modifiers: List[SingleEndModifier] = modifiers

    def process_reads(
        self, progress: Optional[Progress] = None
    ) -> Tuple[int, int, Optional[int]]:
        """Run the pipeline. Return statistics"""
        n = 0  # no. of processed reads
        total_bp = 0
        for read in self._reader:
            n += 1
            if n % 10000 == 0 and progress is not None:
                progress.update(10000)
            total_bp += len(read)
            info = ModificationInfo(read)
            for modifier in self._modifiers:
                read = modifier(read, info)
            for filter_ in self._steps:
                if filter_(read, info):
                    break
        if progress is not None:
            progress.update(n % 10000)
        return (n, total_bp, None)

    def _open_writer(
        self,
        file: BinaryIO,
        file2: Optional[BinaryIO] = None,
        force_fasta: Optional[bool] = None,
    ):
        assert file2 is None
        assert not isinstance(file, (str, bytes, Path))
        return open_raise_limit(
            dnaio.open,
            file,
            mode="w",
            qualities=self.uses_qualities,
            fileformat="fasta" if force_fasta else None,
        )

    def _make_filter(
        self, predicate1: Optional[Predicate], predicate2: Optional[Predicate], writer
    ):
        _ = predicate2
        assert predicate1 is not None
        return SingleEndFilter(predicate1, writer)

    def _make_untrimmed_filter(self, writer):
        return SingleEndFilter(DiscardUntrimmed(), writer)

    def _final_filter(self, outfiles: OutputFiles):
        assert outfiles.out2 is None and outfiles.out is not None
        writer = self._open_writer(outfiles.out, force_fasta=outfiles.force_fasta)
        return SingleEndSink(writer)

    def _create_demultiplexer(self, outfiles: OutputFiles) -> Demultiplexer:
        writers: Dict[Optional[str], Any] = dict()
        if outfiles.untrimmed is not None:
            writers[None] = self._open_writer(
                outfiles.untrimmed, force_fasta=outfiles.force_fasta
            )
        assert outfiles.demultiplex_out is not None
        for name, file in outfiles.demultiplex_out.items():
            writers[name] = self._open_writer(file, force_fasta=outfiles.force_fasta)
        return Demultiplexer(writers)

    def _wrap_single_end_step(self, step: SingleEndStep):
        return step

    @property
    def minimum_length(self):
        return self._minimum_length

    @minimum_length.setter
    def minimum_length(self, value):
        assert value is None or len(value) == 1
        self._minimum_length = value

    @property
    def maximum_length(self):
        return self._maximum_length

    @maximum_length.setter
    def maximum_length(self, value):
        assert value is None or len(value) == 1
        self._maximum_length = value


class PairedEndPipeline(Pipeline):
    """
    Processing pipeline for paired-end reads.
    """

    paired = True

    def __init__(
        self,
        modifiers: List[
            Union[
                PairedEndModifier,
                Tuple[Optional[SingleEndModifier], Optional[SingleEndModifier]],
            ]
        ],
        pair_filter_mode: str,
    ):
        super().__init__()
        self._modifiers: List[PairedEndModifier] = []
        self._pair_filter_mode = pair_filter_mode
        self._reader = None
        # Whether to ignore pair_filter mode for discard-untrimmed filter
        self.override_untrimmed_pair_filter = False
        self._add_modifiers(modifiers)

    def _add_modifiers(self, modifiers):
        for modifier in modifiers:
            if isinstance(modifier, tuple):
                self._add_two_single_modifiers(*modifier)
            else:
                self._add_modifier(modifier)

    def _add_two_single_modifiers(
        self,
        modifier1: Optional[SingleEndModifier],
        modifier2: Optional[SingleEndModifier],
    ) -> None:
        """
        Add two single-end modifiers that modify R1 and R2, respectively.
        One of them can be None, in which case the modifier
        is only applied to the respective other read.
        """
        if modifier1 is None and modifier2 is None:
            raise ValueError("Not both modifiers can be None")
        self._modifiers.append(PairedEndModifierWrapper(modifier1, modifier2))

    def _add_modifier(self, modifier: PairedEndModifier) -> None:
        """Add a Modifier (without wrapping it in a PairedEndModifierWrapper)"""
        self._modifiers.append(modifier)

    def process_reads(
        self, progress: Optional[Progress] = None
    ) -> Tuple[int, int, Optional[int]]:
        n = 0  # no. of processed reads
        total1_bp = 0
        total2_bp = 0
        assert self._reader is not None
        for read1, read2 in self._reader:
            n += 1
            if n % 10000 == 0 and progress is not None:
                progress.update(10000)
            total1_bp += len(read1)
            total2_bp += len(read2)
            info1 = ModificationInfo(read1)
            info2 = ModificationInfo(read2)
            for modifier in self._modifiers:
                read1, read2 = modifier(read1, read2, info1, info2)
            for filter_ in self._steps:
                # Stop writing as soon as one of the filters was successful.
                if filter_(read1, read2, info1, info2):
                    break
        if progress is not None:
            progress.update(n % 10000)
        return (n, total1_bp, total2_bp)

    def _open_writer(
        self,
        file: BinaryIO,
        file2: Optional[BinaryIO] = None,
        force_fasta: Optional[bool] = None,
    ):
        # file and file2 must already be file-like objects because we don’t want to
        # take care of threads and compression levels here.
        for f in (file, file2):
            assert not isinstance(f, (str, bytes, Path))
        return open_raise_limit(
            dnaio.open,
            file,
            file2=file2,
            mode="w",
            qualities=self.uses_qualities,
            fileformat="fasta" if force_fasta else None,
            interleaved=file2 is None,
        )

    def _make_filter(
        self,
        predicate1: Optional[Predicate],
        predicate2: Optional[Predicate],
        writer,
        pair_filter_mode=None,
    ):
        if pair_filter_mode is None:
            pair_filter_mode = self._pair_filter_mode
        return PairedEndFilter(
            predicate1, predicate2, writer, pair_filter_mode=pair_filter_mode
        )

    def _make_untrimmed_filter(self, writer):
        """
        Return a different filter wrapper when adapters were given only for R1
        or only for R2 (then override_untrimmed_pair_filter will be set)
        """
        return self._make_filter(
            DiscardUntrimmed(),
            DiscardUntrimmed(),
            writer,
            pair_filter_mode="both" if self.override_untrimmed_pair_filter else None,
        )

    def _final_filter(self, outfiles):
        writer = self._open_writer(
            outfiles.out, outfiles.out2, force_fasta=outfiles.force_fasta
        )
        return PairedEndSink(writer)

    def _create_demultiplexer(self, outfiles):
        def open_writer(file, file2):
            return self._open_writer(file, file2, force_fasta=outfiles.force_fasta)

        if outfiles.combinatorial_out is not None:
            assert outfiles.untrimmed is None and outfiles.untrimmed2 is None
            writers = dict()
            for key, out in outfiles.combinatorial_out.items():
                writers[key] = open_writer(out, outfiles.combinatorial_out2[key])
            return CombinatorialDemultiplexer(writers)
        else:
            writers = dict()
            if outfiles.untrimmed is not None:
                writers[None] = open_writer(outfiles.untrimmed, outfiles.untrimmed2)
            for name, file in outfiles.demultiplex_out.items():
                writers[name] = open_writer(file, outfiles.demultiplex_out2[name])
            return PairedDemultiplexer(writers)

    def _wrap_single_end_step(self, step: SingleEndStep):
        return PairedSingleEndStep(step)

    @property
    def minimum_length(self):
        return self._minimum_length

    @minimum_length.setter
    def minimum_length(self, value):
        assert value is None or len(value) == 2
        self._minimum_length = value

    @property
    def maximum_length(self):
        return self._maximum_length

    @maximum_length.setter
    def maximum_length(self, value):
        assert value is None or len(value) == 2
        self._maximum_length = value


class ReaderProcess(mpctx_Process):
    """
    Read chunks of FASTA or FASTQ data (single-end or paired) and send them to a worker.

    The reader repeatedly

    - reads a chunk from the file(s)
    - reads a worker index from the Queue
    - sends the chunk to connections[index]

    and finally sends the stop token -1 ("poison pills") to all connections.
    """

    def __init__(
        self,
        *paths: str,
        connections: Sequence[Connection],
        queue: multiprocessing.Queue,
        buffer_size: int,
        stdin_fd,
    ):
        """
        Args:
            paths: path to input files
            connections: a list of Connection objects, one for each worker.
            queue: a Queue of worker indices. A worker writes its own index into this
                queue to notify the reader that it is ready to receive more data.
            buffer_size:
            stdin_fd:

        Note:
            This expects the paths to the input files as strings because these can be pickled
            while file-like objects such as BufferedReader cannot. When using multiprocessing with
            the "spawn" method, which is the default method on macOS, function arguments must be
            picklable.
        """
        super().__init__()
        if len(paths) > 2:
            raise ValueError("Reading from more than two files currently not supported")
        if not paths:
            raise ValueError("Must provide at least one file")
        self._paths = paths
        self.connections = connections
        self.queue = queue
        self.buffer_size = buffer_size
        self.stdin_fd = stdin_fd

    def run(self):
        if self.stdin_fd != -1:
            sys.stdin.close()
            sys.stdin = os.fdopen(self.stdin_fd)
        try:
            with ExitStack() as stack:
                files = [
                    stack.enter_context(xopen_rb_raise_limit(path))
                    for path in self._paths
                ]
                for index, chunks in enumerate(self._read_chunks(*files)):
                    self.send_to_worker(index, *chunks)
            self.shutdown()
        except Exception as e:
            # TODO better send this to a common "something went wrong" Queue
            for connection in self.connections:
                connection.send(-2)
                connection.send((e, traceback.format_exc()))

    def _read_chunks(self, *files) -> Iterator[Tuple[memoryview, ...]]:
        if len(files) == 1:
            for chunk in dnaio.read_chunks(files[0], self.buffer_size):
                yield (chunk,)
        elif len(files) == 2:
            for chunks in dnaio.read_paired_chunks(
                files[0], files[1], self.buffer_size
            ):
                yield chunks
        else:
            raise NotImplementedError

    def send_to_worker(self, chunk_index, chunk1, chunk2=None):
        worker_index = self.queue.get()
        connection = self.connections[worker_index]
        connection.send(chunk_index)
        connection.send_bytes(chunk1)
        if chunk2 is not None:
            connection.send_bytes(chunk2)

    def shutdown(self):
        # Send poison pills to all workers
        for _ in range(len(self.connections)):
            worker_index = self.queue.get()
            self.connections[worker_index].send(-1)


class WorkerProcess(mpctx_Process):
    """
    The worker repeatedly reads chunks of data from the read_pipe, runs the pipeline on it
    and sends the processed chunks to the write_pipe.

    To notify the reader process that it wants data, it puts its own identifier into the
    need_work_queue before attempting to read data from the read_pipe.
    """

    def __init__(
        self,
        id_: int,
        pipeline: Pipeline,
        n_input_files: int,
        interleaved_input: bool,
        orig_outfiles: OutputFiles,
        read_pipe: Connection,
        write_pipe: Connection,
        need_work_queue: multiprocessing.Queue,
    ):
        super().__init__()
        self._id = id_
        self._pipeline = pipeline
        self._n_input_files = n_input_files
        self._interleaved_input = interleaved_input
        self._read_pipe = read_pipe
        self._write_pipe = write_pipe
        self._need_work_queue = need_work_queue
        # Do not store orig_outfiles directly because it contains
        # _io.BufferedWriter attributes, which cannot be pickled.
        self._original_outfiles = orig_outfiles.as_bytesio()

    def run(self):
        try:
            stats = Statistics()
            while True:
                # Notify reader that we need data
                self._need_work_queue.put(self._id)
                chunk_index = self._read_pipe.recv()
                if chunk_index == -1:
                    # reader is done
                    break
                elif chunk_index == -2:
                    # An exception has occurred in the reader
                    e, tb_str = self._read_pipe.recv()
                    logger.error("%s", tb_str)
                    raise e

                infiles = self._make_input_files()
                outfiles = self._original_outfiles.as_bytesio()
                self._pipeline.connect_io(infiles, outfiles)
                (n, bp1, bp2) = self._pipeline.process_reads()
                self._pipeline.flush()
                cur_stats = Statistics().collect(n, bp1, bp2, [], self._pipeline._steps)
                stats += cur_stats
                self._send_outfiles(outfiles, chunk_index, n)
                self._pipeline.close()

            m = self._pipeline._modifiers
            modifier_stats = Statistics().collect(
                0, 0, 0 if self._pipeline.paired else None, m, []
            )
            stats += modifier_stats
            self._write_pipe.send(-1)
            self._write_pipe.send(stats)
        except Exception as e:
            self._write_pipe.send(-2)
            self._write_pipe.send((e, traceback.format_exc()))

    def _make_input_files(self) -> InputFiles:
        files = [
            io.BytesIO(self._read_pipe.recv_bytes()) for _ in range(self._n_input_files)
        ]
        return InputFiles(*files, interleaved=self._interleaved_input)

    def _send_outfiles(self, outfiles: OutputFiles, chunk_index: int, n_reads: int):
        self._write_pipe.send(chunk_index)
        self._write_pipe.send(n_reads)

        for f in outfiles:
            f.flush()
            assert isinstance(f, io.BytesIO)
            processed_chunk = f.getvalue()
            self._write_pipe.send_bytes(processed_chunk)


class OrderedChunkWriter:
    """
    We may receive chunks of processed data from worker processes
    in any order. This class writes them to an output file in
    the correct order.
    """

    def __init__(self, outfile):
        self._chunks = dict()
        self._current_index = 0
        self._outfile = outfile

    def write(self, data, index):
        """ """
        self._chunks[index] = data
        while self._current_index in self._chunks:
            self._outfile.write(self._chunks[self._current_index])
            del self._chunks[self._current_index]
            self._current_index += 1

    def wrote_everything(self):
        return not self._chunks


class PipelineRunner(ABC):
    """
    A read processing pipeline
    """

    def __init__(self, pipeline: Pipeline, progress: Progress):
        self._pipeline = pipeline
        self._progress = progress

    @abstractmethod
    def run(self) -> Statistics:
        pass

    @abstractmethod
    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class ParallelPipelineRunner(PipelineRunner):
    """
    Run a Pipeline in parallel

    - When connect_io() is called, a reader process is spawned.
    - When run() is called, as many worker processes as requested are spawned.
    - In the main process, results are written to the output files in the correct
      order, and statistics are aggregated.

    If a worker needs work, it puts its own index into a Queue() (_need_work_queue).
    The reader process listens on this queue and sends the raw data to the
    worker that has requested work. For sending the data from reader to worker,
    a Connection() is used. There is one such connection for each worker (self._pipes).

    For sending the processed data from the worker to the main process, there
    is a second set of connections, again one for each worker.

    When the reader is finished, it sends 'poison pills' to all workers.
    When a worker receives this, it sends a poison pill to the main process,
    followed by a Statistics object that contains statistics about all the reads
    processed by that worker.
    """

    def __init__(
        self,
        pipeline: Pipeline,
        inpaths: InputPaths,
        outfiles: OutputFiles,
        progress: Progress,
        n_workers: int,
        buffer_size: Optional[int] = None,
    ):
        super().__init__(pipeline, progress)
        self._n_workers = n_workers
        self._need_work_queue: multiprocessing.Queue = mpctx.Queue()
        self._buffer_size = 4 * 1024**2 if buffer_size is None else buffer_size
        self._outfiles = outfiles

        self._assign_input(*inpaths.paths, interleaved=inpaths.interleaved)

    def _assign_input(
        self,
        *paths: str,
        interleaved: bool = False,
    ) -> None:
        self._n_input_files = len(paths)
        self._interleaved_input = interleaved
        # the workers read from these connections
        connections = [mpctx.Pipe(duplex=False) for _ in range(self._n_workers)]
        self._connections, connw = zip(*connections)
        try:
            fileno = sys.stdin.fileno()
        except io.UnsupportedOperation:
            # This happens during tests: pytest sets sys.stdin to an object
            # that does not have a file descriptor.
            fileno = -1
        self._reader_process = ReaderProcess(
            *paths,
            connections=connw,
            queue=self._need_work_queue,
            buffer_size=self._buffer_size,
            stdin_fd=fileno,
        )
        self._reader_process.daemon = True
        self._reader_process.start()

    def _start_workers(self) -> Tuple[List[WorkerProcess], List[Connection]]:
        workers = []
        connections = []
        for index in range(self._n_workers):
            conn_r, conn_w = mpctx.Pipe(duplex=False)
            connections.append(conn_r)
            worker = WorkerProcess(
                index,
                self._pipeline,
                self._n_input_files,
                self._interleaved_input,
                self._outfiles,
                self._connections[index],
                conn_w,
                self._need_work_queue,
            )
            worker.daemon = True
            worker.start()
            workers.append(worker)
        return workers, connections

    def run(self) -> Statistics:
        workers, connections = self._start_workers()
        writers = []
        for f in self._outfiles:
            writers.append(OrderedChunkWriter(f))
        stats = Statistics()
        while connections:
            ready_connections: List[Any] = multiprocessing.connection.wait(connections)
            for connection in ready_connections:
                chunk_index = self._try_receive(connection)
                if chunk_index == -1:
                    # the worker is done
                    cur_stats = self._try_receive(connection)
                    stats += cur_stats
                    connections.remove(connection)
                    continue

                number_of_reads = self._try_receive(connection)
                self._progress.update(number_of_reads)
                for writer in writers:
                    data = connection.recv_bytes()
                    writer.write(data, chunk_index)
        for writer in writers:
            assert writer.wrote_everything()
        for w in workers:
            w.join()
        self._reader_process.join()
        self._progress.close()
        return stats

    @staticmethod
    def _try_receive(connection):
        """
        Try to receive data over `self.connection` and return it.
        If an exception was received, raise it.
        """
        result = connection.recv()
        if result == -2:
            # An exception has occurred on the other end
            e, tb_str = connection.recv()
            # The other end does not send an actual traceback object because these are
            # not picklable, but a string representation.
            logger.debug("%s", tb_str)
            for child in multiprocessing.active_children():
                child.terminate()
            raise e
        return result

    def close(self) -> None:
        self._outfiles.close()


class SerialPipelineRunner(PipelineRunner):
    """
    Run a Pipeline on a single core
    """

    def __init__(
        self,
        pipeline: Pipeline,
        infiles: InputFiles,
        outfiles: OutputFiles,
        progress: Progress,
    ):
        super().__init__(pipeline, progress)
        self._pipeline.connect_io(infiles, outfiles)

    def run(self) -> Statistics:
        (n, total1_bp, total2_bp) = self._pipeline.process_reads(
            progress=self._progress
        )
        if self._progress is not None:
            self._progress.close()
        # TODO
        modifiers = getattr(self._pipeline, "_modifiers", None)
        assert modifiers is not None
        return Statistics().collect(
            n, total1_bp, total2_bp, modifiers, self._pipeline._steps
        )

    def close(self) -> None:
        self._pipeline.close()
