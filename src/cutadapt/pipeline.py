import io
import os
import sys
import copy
import logging
import functools
from typing import List, IO, Optional, BinaryIO, TextIO, Any, Tuple
from abc import ABC, abstractmethod
from multiprocessing import Process, Pipe, Queue
from pathlib import Path
import multiprocessing.connection
from multiprocessing.connection import Connection
import traceback

from xopen import xopen
import dnaio

from .utils import Progress, FileOpener
from .modifiers import SingleEndModifier, PairedModifier, PairedModifierWrapper, ModificationInfo
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
        return dnaio.open(self.file1, file2=self.file2,
            interleaved=self.interleaved, mode="r")


class OutputFiles:
    """
    The attributes are open file-like objects except when demultiplex is True. In that case,
    untrimmed, untrimmed2, out and out2 are file names or templates
    as required by the used demultiplexer ('{name}' etc.).

    Files may also be None.
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
        demultiplex: bool = False,
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
        self.demultiplex = demultiplex
        self.force_fasta = force_fasta

    def __iter__(self):
        yield self.out
        yield self.out2
        yield self.untrimmed
        yield self.untrimmed2
        yield self.too_short
        yield self.too_short2
        yield self.too_long
        yield self.too_long2
        yield self.info
        yield self.rest
        yield self.wildcard


class Pipeline(ABC):
    """
    Processing pipeline that loops over reads and applies modifiers and filters
    """
    n_adapters = 0

    def __init__(self, file_opener: FileOpener):
        self._close_files = []  # type: List[IO]
        self._reader = None  # type: Any
        self._filters = []  # type: List[Any]
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

    def connect_io(self, infiles: InputFiles, outfiles: OutputFiles):
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

    def _set_output(self, outfiles: OutputFiles):
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

        if outfiles.demultiplex:
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

    def flush(self) -> None:
        for f in self._textiowrappers:
            f.flush()
        assert self._outfiles is not None
        for f in self._outfiles:
            if f is not None:
                f.flush()

    def close(self) -> None:
        self._reader.close()
        for f in self._textiowrappers:
            f.close()  # This also closes the underlying files; a second close occurs below
        assert self._outfiles is not None
        for f in self._outfiles:
            # TODO do not use hasattr
            if f is not None and f is not sys.stdin and f is not sys.stdout and hasattr(f, 'close'):
                f.close()
        if self._demultiplexer is not None:
            self._demultiplexer.close()

    @property
    def uses_qualities(self) -> bool:
        assert self._reader is not None
        return self._reader.delivers_qualities

    @abstractmethod
    def process_reads(self, progress: Progress = None):
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
    paired = False

    def __init__(self, file_opener: FileOpener):
        super().__init__(file_opener)
        self._modifiers = []  # type: List[SingleEndModifier]

    def add(self, modifier: SingleEndModifier):
        if modifier is None:
            raise ValueError("Modifier must not be None")
        self._modifiers.append(modifier)

    def process_reads(self, progress: Progress = None):
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
        return dnaio.open(
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
        return Demultiplexer(outfiles.out, outfiles.untrimmed, qualities=self.uses_qualities,
            file_opener=self.file_opener)

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
        self._modifiers = []  # type: List[PairedModifier]
        self._pair_filter_mode = pair_filter_mode
        self._reader = None
        # Whether to ignore pair_filter mode for discard-untrimmed filter
        self.override_untrimmed_pair_filter = False

    def add(self, modifier1: Optional[SingleEndModifier], modifier2: Optional[SingleEndModifier]):
        """
        Add a modifier for R1 and R2. One of them can be None, in which case the modifier
        will only be added for the respective read.
        """
        if modifier1 is None and modifier2 is None:
            raise ValueError("Not both modifiers can be None")
        self._modifiers.append(PairedModifierWrapper(modifier1, modifier2))

    def add_both(self, modifier: SingleEndModifier):
        """
        Add one modifier for both R1 and R2
        """
        assert modifier is not None
        self._modifiers.append(PairedModifierWrapper(modifier, copy.copy(modifier)))

    def add_paired_modifier(self, paired_modifier: PairedModifier):
        """Add a Modifier (without wrapping it in a PairedModifierWrapper)"""
        self._modifiers.append(paired_modifier)

    def process_reads(self, progress: Progress = None):
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
        return dnaio.open(
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
        if '{name1}' in outfiles.out and '{name2}' in outfiles.out:
            return CombinatorialDemultiplexer(outfiles.out, outfiles.out2,
                outfiles.untrimmed, qualities=self.uses_qualities, file_opener=self.file_opener)
        else:
            return PairedDemultiplexer(outfiles.out, outfiles.out2,
                outfiles.untrimmed, outfiles.untrimmed2, qualities=self.uses_qualities,
                file_opener=self.file_opener)

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

    def __init__(self, file, file2, connections, queue, buffer_size, stdin_fd):
        """
        queue -- a Queue of worker indices. A worker writes its own index into this
            queue to notify the reader that it is ready to receive more data.
        connections -- a list of Connection objects, one for each worker.
        """
        super().__init__()
        self.file = file
        self.file2 = file2
        self.connections = connections
        self.queue = queue
        self.buffer_size = buffer_size
        self.stdin_fd = stdin_fd

    def run(self):
        if self.stdin_fd != -1:
            sys.stdin.close()
            sys.stdin = os.fdopen(self.stdin_fd)
        try:
            with xopen(self.file, 'rb') as f:
                if self.file2:
                    with xopen(self.file2, 'rb') as f2:
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
    def __init__(self, id_, pipeline, two_input_files,
            interleaved_input, orig_outfiles, read_pipe, write_pipe, need_work_queue):
        super().__init__()
        self._id = id_
        self._pipeline = pipeline
        self._two_input_files = two_input_files
        self._interleaved_input = interleaved_input
        self._orig_outfiles = orig_outfiles
        self._read_pipe = read_pipe
        self._write_pipe = write_pipe
        self._need_work_queue = need_work_queue

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
                outfiles = self._make_output_files()
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

    def _make_input_files(self):
        data = self._read_pipe.recv_bytes()
        input = io.BytesIO(data)

        if self._two_input_files:
            data = self._read_pipe.recv_bytes()
            input2 = io.BytesIO(data)
        else:
            input2 = None
        return InputFiles(input, input2, interleaved=self._interleaved_input)

    def _make_output_files(self):
        """
        Using self._orig_outfiles as a template, make a new OutputFiles instance
        that has BytesIO instances for each non-None output file
        """
        output_files = copy.copy(self._orig_outfiles)
        for attr in (
            "out", "out2", "untrimmed", "untrimmed2", "too_short", "too_short2", "too_long",
            "too_long2", "info", "rest", "wildcard"
        ):
            orig_outfile = getattr(self._orig_outfiles, attr)
            if orig_outfile is not None:
                output = io.BytesIO()
                setattr(output_files, attr, output)

        return output_files

    def _send_outfiles(self, outfiles: OutputFiles, chunk_index: int, n_reads: int):
        self._write_pipe.send(chunk_index)
        self._write_pipe.send(n_reads)

        for f in outfiles:
            if f is None:
                continue
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
    def __init__(self, pipeline: Pipeline, progress: Progress, *args, **kwargs):
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
        infiles: InputFiles,
        outfiles: OutputFiles,
        progress: Progress,
        n_workers: int,
        buffer_size: int = 4 * 1024**2,
    ):
        super().__init__(pipeline, progress)
        self._n_workers = n_workers
        self._need_work_queue = Queue()  # type: Queue
        self._buffer_size = buffer_size
        self._assign_input(infiles.file1, infiles.file2, infiles.interleaved)
        self._assign_output(outfiles)

    def _assign_input(
        self,
        file1: BinaryIO,
        file2: Optional[BinaryIO] = None,
        interleaved: bool = False,
    ):
        self._two_input_files = file2 is not None
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
        self._reader_process = ReaderProcess(file1, file2, connw,
            self._need_work_queue, self._buffer_size, fileno)
        self._reader_process.daemon = True
        self._reader_process.start()

    @staticmethod
    def can_output_to(outfiles: OutputFiles) -> bool:
        return outfiles.out is not None and not outfiles.demultiplex

    def _assign_output(self, outfiles: OutputFiles):
        if not self.can_output_to(outfiles):
            raise ValueError()
        self._outfiles = outfiles

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

    def run(self):
        workers, connections = self._start_workers()
        writers = []
        for outfile in self._outfiles:
            if outfile is None:
                continue
            writers.append(OrderedChunkWriter(outfile))
        stats = None
        n = 0  # A running total of the number of processed reads (for progress indicator)
        while connections:
            ready_connections = multiprocessing.connection.wait(connections)
            for connection in ready_connections:
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
                    if stats is None:
                        stats = cur_stats
                    else:
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

    def close(self):
        for f in self._outfiles:
            # TODO do not use hasattr
            if f is not None and f is not sys.stdin and f is not sys.stdout and hasattr(f, 'close'):
                f.close()


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

    def run(self):
        (n, total1_bp, total2_bp) = self._pipeline.process_reads(progress=self._progress)
        if self._progress:
            self._progress.stop(n)
        # TODO
        return Statistics().collect(n, total1_bp, total2_bp, self._pipeline._modifiers, self._pipeline._filters)

    def close(self):
        self._pipeline.close()
