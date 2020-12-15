import io
import os
import sys
import copy
import logging
import functools
from typing import List, Optional, BinaryIO, TextIO, Any, Tuple, Dict
from abc import ABC, abstractmethod
from multiprocessing import Process, Pipe, Queue
from pathlib import Path
import multiprocessing.connection
from multiprocessing.connection import Connection
import traceback

from xopen import xopen
import dnaio

from .utils import Progress, FileOpener
from .modifiers import SingleEndModifier, PairedEndModifier, PairedEndModifierWrapper, ModificationInfo
from .report import Statistics
from .filters import (Redirector, PairedRedirector, NoFilter, PairedNoFilter, InfoFileWriter,
    RestFileWriter, WildcardFileWriter, TooShortReadFilter, TooLongReadFilter, NContentFilter,
    MaximumExpectedErrorsFilter,
    CasavaFilter, DiscardTrimmedFilter, DiscardUntrimmedFilter, Demultiplexer,
    PairedDemultiplexer, CombinatorialDemultiplexer)

logger = logging.getLogger()


class InputFiles:
    def __init__(self, file1: BinaryIO, file2: Optional[BinaryIO] = None, interleaved: bool = False):
        self.file1 = file1
        self.file2 = file2
        self.interleaved = interleaved

    def open(self):
        return dnaio.open(self.file1, file2=self.file2, interleaved=self.interleaved, mode="r")

    def close(self) -> None:
        self.file1.close()
        if self.file2 is not None:
            self.file2.close()


class InputPaths:
    def __init__(self, path1: str, path2: Optional[str] = None, interleaved: bool = False):
        self.path1 = path1
        self.path2 = path2
        self.interleaved = interleaved

    def open(self, file_opener: FileOpener) -> InputFiles:
        file1, file2 = file_opener.xopen_pair(self.path1, self.path2, "rb")
        return InputFiles(file1, file2, self.interleaved)


class OutputFiles:
    """
    The attributes are either None or open file-like objects except for demultiplex_out
    and demultiplex_out2, which are dictionaries that map an adapter name
    to file-like objects.
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
            self.demultiplex_out, self.demultiplex_out2,
            self.combinatorial_out, self.combinatorial_out2,
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
            "out", "out2", "untrimmed", "untrimmed2", "too_short", "too_short2", "too_long",
            "too_long2", "info", "rest", "wildcard"
        ):
            if getattr(self, attr) is not None:
                setattr(result, attr, io.BytesIO())
        for attr in "demultiplex_out", "demultiplex_out2", "combinatorial_out", "combinatorial_out2":
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

    def __init__(self, file_opener: FileOpener):
        self._reader = None  # type: Any
        self._filters = []  # type: List[Any]
        self._infiles = None  # type: Optional[InputFiles]
        self._outfiles = None  # type: Optional[OutputFiles]
        self._demultiplexer = None
        self._textiowrappers = []  # type: List[TextIO]

        # Filter settings
        self._minimum_length = None
        self._maximum_length = None
        self.max_n = None
        self.max_expected_errors = None
        self.discard_casava = False
        self.discard_trimmed = False
        self.discard_untrimmed = False
        self.file_opener = file_opener

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
        self._filters = []
        self._outfiles = outfiles
        filter_wrapper = self._filter_wrapper()

        for filter_class, outfile in (
            (RestFileWriter, outfiles.rest),
            (InfoFileWriter, outfiles.info),
            (WildcardFileWriter, outfiles.wildcard),
        ):
            if outfile:
                textiowrapper = io.TextIOWrapper(outfile)
                self._textiowrappers.append(textiowrapper)
                self._filters.append(filter_wrapper(None, filter_class(textiowrapper), None))

        # minimum length and maximum length
        for lengths, file1, file2, filter_class in (
                (self._minimum_length, outfiles.too_short, outfiles.too_short2, TooShortReadFilter),
                (self._maximum_length, outfiles.too_long, outfiles.too_long2, TooLongReadFilter)
        ):
            if lengths is None:
                continue
            writer = self._open_writer(file1, file2) if file1 else None
            f1 = filter_class(lengths[0]) if lengths[0] is not None else None
            if len(lengths) == 2 and lengths[1] is not None:
                f2 = filter_class(lengths[1])
            else:
                f2 = None
            self._filters.append(filter_wrapper(writer, filter=f1, filter2=f2))

        if self.max_n is not None:
            f1 = f2 = NContentFilter(self.max_n)
            self._filters.append(filter_wrapper(None, f1, f2))

        if self.max_expected_errors is not None:
            if not self._reader.delivers_qualities:
                logger.warning("Ignoring option --max-ee as input does not contain quality values")
            else:
                f1 = f2 = MaximumExpectedErrorsFilter(self.max_expected_errors)
                self._filters.append(filter_wrapper(None, f1, f2))

        if self.discard_casava:
            f1 = f2 = CasavaFilter()
            self._filters.append(filter_wrapper(None, f1, f2))

        if int(self.discard_trimmed) + int(self.discard_untrimmed) + int(outfiles.untrimmed is not None) > 1:
            raise ValueError('discard_trimmed, discard_untrimmed and outfiles.untrimmed must not '
                'be set simultaneously')

        if outfiles.demultiplex_out is not None or outfiles.combinatorial_out is not None:
            self._demultiplexer = self._create_demultiplexer(outfiles)
            self._filters.append(self._demultiplexer)
        else:
            # Allow overriding the wrapper for --discard-untrimmed/--untrimmed-(paired-)output
            untrimmed_filter_wrapper = self._untrimmed_filter_wrapper()

            # Set up the remaining filters to deal with --discard-trimmed,
            # --discard-untrimmed and --untrimmed-output. These options
            # are mutually exclusive in order to avoid brain damage.
            if self.discard_trimmed:
                self._filters.append(
                    filter_wrapper(None, DiscardTrimmedFilter(), DiscardTrimmedFilter()))
            elif self.discard_untrimmed:
                self._filters.append(
                    untrimmed_filter_wrapper(None, DiscardUntrimmedFilter(), DiscardUntrimmedFilter()))
            elif outfiles.untrimmed:
                untrimmed_writer = self._open_writer(outfiles.untrimmed, outfiles.untrimmed2)
                self._filters.append(
                    untrimmed_filter_wrapper(untrimmed_writer, DiscardUntrimmedFilter(), DiscardUntrimmedFilter()))
            self._filters.append(self._final_filter(outfiles))
        logger.debug("Filters: %s", self._filters)

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
    def process_reads(self, progress: Progress = None) -> Tuple[int, int, Optional[int]]:
        pass

    @abstractmethod
    def _filter_wrapper(self):
        pass

    @abstractmethod
    def _untrimmed_filter_wrapper(self):
        pass

    @abstractmethod
    def _final_filter(self, outfiles: OutputFiles):
        pass

    @abstractmethod
    def _create_demultiplexer(self, outfiles: OutputFiles):
        pass


class SingleEndPipeline(Pipeline):
    """
    Processing pipeline for single-end reads
    """
    def __init__(self, file_opener: FileOpener):
        super().__init__(file_opener)
        self._modifiers = []  # type: List[SingleEndModifier]

    def add(self, modifier: SingleEndModifier):
        if modifier is None:
            raise ValueError("Modifier must not be None")
        self._modifiers.append(modifier)

    def process_reads(self, progress: Progress = None) -> Tuple[int, int, Optional[int]]:
        """Run the pipeline. Return statistics"""
        n = 0  # no. of processed reads
        total_bp = 0
        for read in self._reader:
            n += 1
            if n % 10000 == 0 and progress:
                progress.update(n)
            total_bp += len(read)
            info = ModificationInfo(read)
            for modifier in self._modifiers:
                read = modifier(read, info)
            for filter_ in self._filters:
                if filter_(read, info):
                    break
        return (n, total_bp, None)

    def _open_writer(
        self,
        file: BinaryIO,
        file2: Optional[BinaryIO] = None,
        force_fasta: Optional[bool] = None,
    ):
        assert file2 is None
        assert not isinstance(file, (str, bytes, Path))
        return self.file_opener.dnaio_open_raise_limit(
            file, mode="w", qualities=self.uses_qualities, fileformat="fasta" if force_fasta else None)

    def _filter_wrapper(self):
        return Redirector

    def _untrimmed_filter_wrapper(self):
        return Redirector

    def _final_filter(self, outfiles: OutputFiles):
        assert outfiles.out2 is None and outfiles.out is not None
        writer = self._open_writer(outfiles.out, force_fasta=outfiles.force_fasta)
        return NoFilter(writer)

    def _create_demultiplexer(self, outfiles: OutputFiles):
        writers = dict()  # type: Dict[Optional[str], Any]
        if outfiles.untrimmed is not None:
            writers[None] = self._open_writer(outfiles.untrimmed, force_fasta=outfiles.force_fasta)
        assert outfiles.demultiplex_out is not None
        for name, file in outfiles.demultiplex_out.items():
            writers[name] = self._open_writer(file, force_fasta=outfiles.force_fasta)
        return Demultiplexer(writers)

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

    def __init__(self, pair_filter_mode, file_opener: FileOpener):
        super().__init__(file_opener)
        self._modifiers = []  # type: List[PairedEndModifier]
        self._pair_filter_mode = pair_filter_mode
        self._reader = None
        # Whether to ignore pair_filter mode for discard-untrimmed filter
        self.override_untrimmed_pair_filter = False

    def add(self, modifier1: Optional[SingleEndModifier], modifier2: Optional[SingleEndModifier]) -> None:
        """
        Add a modifier for R1 and R2. One of them can be None, in which case the modifier
        will only be added for the respective read.
        """
        if modifier1 is None and modifier2 is None:
            raise ValueError("Not both modifiers can be None")
        self._modifiers.append(PairedEndModifierWrapper(modifier1, modifier2))

    def add_both(self, modifier: SingleEndModifier) -> None:
        """
        Add one modifier for both R1 and R2
        """
        assert modifier is not None
        self._modifiers.append(PairedEndModifierWrapper(modifier, copy.copy(modifier)))

    def add_paired_modifier(self, modifier: PairedEndModifier) -> None:
        """Add a Modifier (without wrapping it in a PairedEndModifierWrapper)"""
        self._modifiers.append(modifier)

    def process_reads(self, progress: Progress = None) -> Tuple[int, int, Optional[int]]:
        n = 0  # no. of processed reads
        total1_bp = 0
        total2_bp = 0
        assert self._reader is not None
        for read1, read2 in self._reader:
            n += 1
            if n % 10000 == 0 and progress:
                progress.update(n)
            total1_bp += len(read1)
            total2_bp += len(read2)
            info1 = ModificationInfo(read1)
            info2 = ModificationInfo(read2)
            for modifier in self._modifiers:
                read1, read2 = modifier(read1, read2, info1, info2)
            for filter_ in self._filters:
                # Stop writing as soon as one of the filters was successful.
                if filter_(read1, read2, info1, info2):
                    break
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
        return self.file_opener.dnaio_open_raise_limit(
            file,
            file2=file2,
            mode="w",
            qualities=self.uses_qualities,
            fileformat="fasta" if force_fasta else None,
            interleaved=file2 is None,
        )

    def _filter_wrapper(self, pair_filter_mode=None):
        if pair_filter_mode is None:
            pair_filter_mode = self._pair_filter_mode
        return functools.partial(PairedRedirector, pair_filter_mode=pair_filter_mode)

    def _untrimmed_filter_wrapper(self):
        """
        Return a different filter wrapper when adapters were given only for R1
        or only for R2 (then override_untrimmed_pair_filter will be set)
        """
        if self.override_untrimmed_pair_filter:
            return self._filter_wrapper(pair_filter_mode='both')
        else:
            return self._filter_wrapper()

    def _final_filter(self, outfiles):
        writer = self._open_writer(outfiles.out, outfiles.out2, force_fasta=outfiles.force_fasta)
        return PairedNoFilter(writer)

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


class ReaderProcess(Process):
    """
    Read chunks of FASTA or FASTQ data (single-end or paired) and send to a worker.

    The reader repeatedly

    - reads a chunk from the file(s)
    - reads a worker index from the Queue
    - sends the chunk to connections[index]

    and finally sends the stop token -1 ("poison pills") to all connections.
    """

    def __init__(self, path: str, path2: Optional[str], connections, queue, buffer_size, stdin_fd):
        """
        queue -- a Queue of worker indices. A worker writes its own index into this
            queue to notify the reader that it is ready to receive more data.
        connections -- a list of Connection objects, one for each worker.
        """
        super().__init__()
        self.path = path
        self.path2 = path2
        self.connections = connections
        self.queue = queue
        self.buffer_size = buffer_size
        self.stdin_fd = stdin_fd

    def run(self):
        if self.stdin_fd != -1:
            sys.stdin.close()
            sys.stdin = os.fdopen(self.stdin_fd)
        try:
            with xopen(self.path, 'rb') as f:
                if self.path2:
                    with xopen(self.path2, 'rb') as f2:
                        for chunk_index, (chunk1, chunk2) in enumerate(
                                dnaio.read_paired_chunks(f, f2, self.buffer_size)):
                            self.send_to_worker(chunk_index, chunk1, chunk2)
                else:
                    for chunk_index, chunk in enumerate(dnaio.read_chunks(f, self.buffer_size)):
                        self.send_to_worker(chunk_index, chunk)

            # Send poison pills to all workers
            for _ in range(len(self.connections)):
                worker_index = self.queue.get()
                self.connections[worker_index].send(-1)
        except Exception as e:
            # TODO better send this to a common "something went wrong" Queue
            for connection in self.connections:
                connection.send(-2)
                connection.send((e, traceback.format_exc()))

    def send_to_worker(self, chunk_index, chunk1, chunk2=None):
        worker_index = self.queue.get()
        connection = self.connections[worker_index]
        connection.send(chunk_index)
        connection.send_bytes(chunk1)
        if chunk2 is not None:
            connection.send_bytes(chunk2)


class WorkerProcess(Process):
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
        two_input_files: bool,
        interleaved_input: bool,
        orig_outfiles: OutputFiles,
        read_pipe: Connection,
        write_pipe: Connection,
        need_work_queue: Queue,
    ):
        super().__init__()
        self._id = id_
        self._pipeline = pipeline
        self._two_input_files = two_input_files
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
                    logger.error('%s', tb_str)
                    raise e

                infiles = self._make_input_files()
                outfiles = self._original_outfiles.as_bytesio()
                self._pipeline.connect_io(infiles, outfiles)
                (n, bp1, bp2) = self._pipeline.process_reads()
                self._pipeline.flush()
                cur_stats = Statistics().collect(n, bp1, bp2, [], self._pipeline._filters)
                stats += cur_stats
                self._send_outfiles(outfiles, chunk_index, n)

            m = self._pipeline._modifiers
            modifier_stats = Statistics().collect(0, 0, 0 if self._pipeline.paired else None, m, [])
            stats += modifier_stats
            self._write_pipe.send(-1)
            self._write_pipe.send(stats)
        except Exception as e:
            self._write_pipe.send(-2)
            self._write_pipe.send((e, traceback.format_exc()))

    def _make_input_files(self) -> InputFiles:
        data = self._read_pipe.recv_bytes()
        input = io.BytesIO(data)

        if self._two_input_files:
            data = self._read_pipe.recv_bytes()
            input2 = io.BytesIO(data)  # type: Optional[BinaryIO]
        else:
            input2 = None
        return InputFiles(input, input2, interleaved=self._interleaved_input)

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
        """
        """
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
    def run(self):
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
        infiles: InputPaths,
        outfiles: OutputFiles,
        progress: Progress,
        n_workers: int,
        buffer_size: int = 4 * 1024**2,
    ):
        super().__init__(pipeline, progress)
        self._n_workers = n_workers
        self._need_work_queue = Queue()  # type: Queue
        self._buffer_size = buffer_size
        self._assign_input(infiles.path1, infiles.path2, infiles.interleaved)
        self._outfiles = outfiles

    def _assign_input(
        self,
        path1: str,
        path2: Optional[str] = None,
        interleaved: bool = False,
    ) -> None:
        self._two_input_files = path2 is not None
        self._interleaved_input = interleaved
        # the workers read from these connections
        connections = [Pipe(duplex=False) for _ in range(self._n_workers)]
        self._connections, connw = zip(*connections)
        try:
            fileno = sys.stdin.fileno()
        except io.UnsupportedOperation:
            # This happens during tests: pytest sets sys.stdin to an object
            # that does not have a file descriptor.
            fileno = -1
        self._reader_process = ReaderProcess(path1, path2, connw,
            self._need_work_queue, self._buffer_size, fileno)
        self._reader_process.daemon = True
        self._reader_process.start()

    def _start_workers(self) -> Tuple[List[WorkerProcess], List[Connection]]:
        workers = []
        connections = []
        for index in range(self._n_workers):
            conn_r, conn_w = Pipe(duplex=False)
            connections.append(conn_r)
            worker = WorkerProcess(
                index, self._pipeline,
                self._two_input_files,
                self._interleaved_input, self._outfiles,
                self._connections[index], conn_w, self._need_work_queue)
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
        n = 0  # A running total of the number of processed reads (for progress indicator)
        while connections:
            ready_connections = multiprocessing.connection.wait(connections)
            for connection in ready_connections:
                assert isinstance(connection, Connection)
                chunk_index = connection.recv()
                if chunk_index == -1:
                    # the worker is done
                    cur_stats = connection.recv()
                    if stats == -2:
                        # An exception has occurred in the worker (see below,
                        # this happens only when there is an exception sending
                        # the statistics)
                        e, tb_str = connection.recv()
                        logger.error('%s', tb_str)
                        raise e
                    stats += cur_stats
                    connections.remove(connection)
                    continue
                elif chunk_index == -2:
                    # An exception has occurred in the worker
                    e, tb_str = connection.recv()

                    # We should use the worker's actual traceback object
                    # here, but traceback objects are not picklable.
                    logger.error('%s', tb_str)
                    raise e

                # No. of reads processed in this chunk
                chunk_n = connection.recv()
                if chunk_n == -2:
                    e, tb_str = connection.recv()
                    logger.error('%s', tb_str)
                    raise e
                n += chunk_n
                self._progress.update(n)
                for writer in writers:
                    data = connection.recv_bytes()
                    writer.write(data, chunk_index)
        for writer in writers:
            assert writer.wrote_everything()
        for w in workers:
            w.join()
        self._reader_process.join()
        self._progress.stop(n)
        return stats

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
        (n, total1_bp, total2_bp) = self._pipeline.process_reads(progress=self._progress)
        if self._progress:
            self._progress.stop(n)
        # TODO
        modifiers = getattr(self._pipeline, "_modifiers", None)
        assert modifiers is not None
        return Statistics().collect(n, total1_bp, total2_bp, modifiers, self._pipeline._filters)

    def close(self) -> None:
        self._pipeline.close()
