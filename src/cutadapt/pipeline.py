import io
import os
import sys
import copy
import logging
import functools
from abc import ABC, abstractmethod
from multiprocessing import Process, Pipe, Queue
import multiprocessing.connection
import traceback

from xopen import xopen
import dnaio

from .utils import Progress
from .modifiers import PairedModifier, Modifier
from .report import Statistics
from .filters import (Redirector, PairedRedirector, NoFilter, PairedNoFilter, InfoFileWriter,
    RestFileWriter, WildcardFileWriter, TooShortReadFilter, TooLongReadFilter, NContentFilter,
    CasavaFilter, DiscardTrimmedFilter, DiscardUntrimmedFilter, Demultiplexer,
    PairedDemultiplexer)

logger = logging.getLogger()


class InputFiles:
    def __init__(self, file1, file2=None, interleaved=False):
        self.file1 = file1
        self.file2 = file2
        self.interleaved = interleaved


class OutputFiles:
    """
    The attributes are open file-like objects except when demultiplex is True. In that case,
    untrimmed, untrimmed2 are file names, and out and out2 are file name templates
    containing '{name}'.
    If interleaved is True, then out is written interleaved.
    Files may also be None.
    """
    # TODO interleaving for the other file pairs (too_short, too_long, untrimmed)?
    def __init__(
            self,
            out=None,
            out2=None,
            untrimmed=None,
            untrimmed2=None,
            too_short=None,
            too_short2=None,
            too_long=None,
            too_long2=None,
            info=None,
            rest=None,
            wildcard=None,
            demultiplex=False,
            interleaved=False,
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
        self.interleaved = interleaved

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

    def __init__(self):
        self._close_files = []
        self._reader = None
        self._filters = []
        self._modifiers = []
        self._outfiles = None
        self._demultiplexer = None

        # Filter settings
        self._minimum_length = None
        self._maximum_length = None
        self.max_n = None
        self.discard_casava = False
        self.discard_trimmed = False
        self.discard_untrimmed = False

    def set_input(self, infiles: InputFiles):
        self._reader = dnaio.open(infiles.file1, file2=infiles.file2,
            interleaved=infiles.interleaved, mode='r')

    def _open_writer(self, file, file2, **kwargs):
        # TODO backwards-incompatible change (?) would be to use outfiles.interleaved
        # for all outputs
        return dnaio.open(file, file2=file2, mode='w', qualities=self.uses_qualities,
            **kwargs)

    def set_output(self, outfiles):
        self._filters = []
        self._outfiles = outfiles
        filter_wrapper = self._filter_wrapper()

        if outfiles.rest:
            self._filters.append(RestFileWriter(outfiles.rest))
        if outfiles.info:
            self._filters.append(InfoFileWriter(outfiles.info))
        if outfiles.wildcard:
            self._filters.append(WildcardFileWriter(outfiles.wildcard))

        # minimum length and maximum length
        for lengths, file1, file2, filter_class in (
                (self._minimum_length, outfiles.too_short, outfiles.too_short2, TooShortReadFilter),
                (self._maximum_length, outfiles.too_long, outfiles.too_long2, TooLongReadFilter)
        ):
            writer = None
            if lengths is not None:
                if file1:
                    writer = self._open_writer(file1, file2)
                f1 = filter_class(lengths[0]) if lengths[0] is not None else None
                if len(lengths) == 2 and lengths[1] is not None:
                    f2 = filter_class(lengths[1])
                else:
                    f2 = None
                self._filters.append(filter_wrapper(writer, filter=f1, filter2=f2))

        if self.max_n is not None:
            f1 = f2 = NContentFilter(self.max_n)
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
                self._filters.append(filter_wrapper(None, DiscardTrimmedFilter(), DiscardTrimmedFilter()))
            elif self.discard_untrimmed:
                self._filters.append(untrimmed_filter_wrapper(None, DiscardUntrimmedFilter(), DiscardUntrimmedFilter()))
            elif outfiles.untrimmed:
                untrimmed_writer = self._open_writer(outfiles.untrimmed, outfiles.untrimmed2)
                self._filters.append(untrimmed_filter_wrapper(untrimmed_writer, DiscardUntrimmedFilter(), DiscardUntrimmedFilter()))
            self._filters.append(self._final_filter(outfiles))

    def close(self):
        for f in self._outfiles:
            # TODO do not use hasattr
            if f is not None and f is not sys.stdin and f is not sys.stdout and hasattr(f, 'close'):
                f.close()
        if self._demultiplexer is not None:
            self._demultiplexer.close()

    @property
    def uses_qualities(self):
        return self._reader.delivers_qualities

    @abstractmethod
    # TODO progress shouldn’t be a parameter
    def process_reads(self, progress: Progress=None):
        pass

    @abstractmethod
    def _filter_wrapper(self):
        pass

    @abstractmethod
    def _untrimmed_filter_wrapper(self):
        pass

    @abstractmethod
    def _final_filter(self, outfiles):
        pass

    @abstractmethod
    def _create_demultiplexer(self, outfiles):
        pass


class SingleEndPipeline(Pipeline):
    """
    Processing pipeline for single-end reads
    """
    paired = False

    def __init__(self):
        super().__init__()
        self._modifiers = []

    def add(self, modifier):
        if modifier is None:
            raise ValueError("Modifier must not be None")
        self._modifiers.append(modifier)

    def process_reads(self, progress: Progress = None):
        """Run the pipeline. Return statistics"""
        n = 0  # no. of processed reads  # TODO turn into attribute
        total_bp = 0
        for read in self._reader:
            n += 1
            if n % 10000 == 0 and progress:
                progress.update(n)
            total_bp += len(read.sequence)
            matches = []
            for modifier in self._modifiers:
                read = modifier(read, matches)
            for filter_ in self._filters:
                if filter_(read, matches):
                    break
        return (n, total_bp, None)

    def _filter_wrapper(self):
        return Redirector

    def _untrimmed_filter_wrapper(self):
        return Redirector

    def _final_filter(self, outfiles):
        writer = self._open_writer(outfiles.out, outfiles.out2)
        return NoFilter(writer)

    def _create_demultiplexer(self, outfiles):
        return Demultiplexer(outfiles.out, outfiles.untrimmed, qualities=self.uses_qualities)

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

    def __init__(self, pair_filter_mode):
        """
        Setting modify_first_read_only to True enables "legacy mode"
        """
        super().__init__()
        self._pair_filter_mode = pair_filter_mode
        self._reader = None
        # Whether to gnore pair_filter mode for discard-untrimmed filter
        self.override_untrimmed_pair_filter = False

    def add(self, modifier1, modifier2):
        """
        Add a modifier for R1 and R2. One of them can be None, in which case the modifier
        will only be added for the respective read.
        """
        if modifier1 is None and modifier2 is None:
            raise ValueError("Not both modifiers can be None")
        self._modifiers.append(PairedModifier(modifier1, modifier2))

    def add_both(self, modifier):
        """
        Add one modifier for both R1 and R2
        """
        assert modifier is not None
        self._modifiers.append(PairedModifier(modifier, copy.copy(modifier)))

    def add_paired_modifier(self, paired_modifier):
        """Add a Modifier without wrapping it in a PairedModifier"""
        self._modifiers.append(paired_modifier)

    def process_reads(self, progress: Progress = None):
        n = 0  # no. of processed reads
        total1_bp = 0
        total2_bp = 0
        for read1, read2 in self._reader:
            n += 1
            if n % 10000 == 0 and progress:
                progress.update(n)
            total1_bp += len(read1.sequence)
            total2_bp += len(read2.sequence)
            matches1 = []
            matches2 = []
            for modifier in self._modifiers:
                read1, read2 = modifier(read1, read2, matches1, matches2)
            for filter_ in self._filters:
                # Stop writing as soon as one of the filters was successful.
                if filter_(read1, read2, matches1, matches2):
                    break
        return (n, total1_bp, total2_bp)

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
        writer = self._open_writer(outfiles.out, outfiles.out2, interleaved=outfiles.interleaved)
        return PairedNoFilter(writer)

    def _create_demultiplexer(self, outfiles):
        return PairedDemultiplexer(outfiles.out, outfiles.out2,
            outfiles.untrimmed, outfiles.untrimmed2, qualities=self.uses_qualities)

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


def reader_process(file, file2, connections, queue, buffer_size, stdin_fd):
    """
    Read chunks of FASTA or FASTQ data from *file* and send to a worker.

    queue -- a Queue of worker indices. A worker writes its own index into this
        queue to notify the reader that it is ready to receive more data.
    connections -- a list of Connection objects, one for each worker.

    The function repeatedly

    - reads a chunk from the file
    - reads a worker index from the Queue
    - sends the chunk to connections[index]

    and finally sends "poison pills" (the value -1) to all connections.
    """
    if stdin_fd != -1:
        sys.stdin.close()
        sys.stdin = os.fdopen(stdin_fd)
    try:
        with xopen(file, 'rb') as f:
            if file2:
                with xopen(file2, 'rb') as f2:
                    for chunk_index, (chunk1, chunk2) in enumerate(dnaio.read_paired_chunks(f, f2, buffer_size)):
                        # Determine the worker that should get this chunk
                        worker_index = queue.get()
                        pipe = connections[worker_index]
                        pipe.send(chunk_index)
                        pipe.send_bytes(chunk1)
                        pipe.send_bytes(chunk2)
            else:
                for chunk_index, chunk in enumerate(dnaio.read_chunks(f, buffer_size)):
                    # Determine the worker that should get this chunk
                    worker_index = queue.get()
                    pipe = connections[worker_index]
                    pipe.send(chunk_index)
                    pipe.send_bytes(chunk)

        # Send poison pills to all workers
        for _ in range(len(connections)):
            worker_index = queue.get()
            connections[worker_index].send(-1)
    except Exception as e:
        # TODO better send this to a common "something went wrong" Queue
        for worker_index in range(len(connections)):
            connections[worker_index].send(-2)
            connections[worker_index].send((e, traceback.format_exc()))


class WorkerProcess(Process):
    """
    The worker repeatedly reads chunks of data from the read_pipe, runs the pipeline on it
    and sends the processed chunks to the write_pipe.

    To notify the reader process that it wants data, it puts its own identifier into the
    need_work_queue before attempting to read data from the read_pipe.
    """
    def __init__(self, id_, pipeline, input_path1, input_path2,
            interleaved_input, orig_outfiles, read_pipe, write_pipe, need_work_queue):
        super().__init__()
        self._id = id_
        self._pipeline = pipeline
        self._input_path1 = input_path1
        self._input_path2 = input_path2
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

                # Setting the .buffer.name attributess below is necessary because
                # file format detection uses the file name
                data = self._read_pipe.recv_bytes()
                input = io.BytesIO(data)
                input.name = self._input_path1

                if self._input_path2:
                    data = self._read_pipe.recv_bytes()
                    input2 = io.BytesIO(data)
                    input2.name = self._input_path2
                else:
                    input2 = None
                output = io.BytesIO()
                output.name = self._orig_outfiles.out.name

                if self._orig_outfiles.out2 is not None:
                    output2 = io.BytesIO()
                    output2.name = self._orig_outfiles.out2.name
                else:
                    output2 = None

                infiles = InputFiles(input, input2, interleaved=self._interleaved_input)
                outfiles = OutputFiles(out=output, out2=output2, interleaved=self._orig_outfiles.interleaved)
                self._pipeline.set_input(infiles)
                self._pipeline.set_output(outfiles)
                (n, bp1, bp2) = self._pipeline.process_reads()
                cur_stats = Statistics().collect(n, bp1, bp2, [], self._pipeline._filters)
                stats += cur_stats

                output.flush()
                processed_chunk = output.getvalue()

                self._write_pipe.send(chunk_index)
                self._write_pipe.send(n)  # no. of reads processed in this chunk
                self._write_pipe.send_bytes(processed_chunk)
                if self._orig_outfiles.out2 is not None:
                    output2.flush()
                    processed_chunk2 = output2.getvalue()
                    self._write_pipe.send_bytes(processed_chunk2)

            m = self._pipeline._modifiers
            modifier_stats = Statistics().collect(0, 0, 0 if self._pipeline.paired else None, m, [])
            stats += modifier_stats
            self._write_pipe.send(-1)
            self._write_pipe.send(stats)
        except Exception as e:
            self._write_pipe.send(-2)
            self._write_pipe.send((e, traceback.format_exc()))


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

    def write(self, data, chunk_index):
        """
        """
        self._chunks[chunk_index] = data
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
    def __init__(self, pipeline):
        self._pipeline = pipeline
        self._progress = Progress()

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def close(self):
        pass


class ParallelPipelineRunner(PipelineRunner):
    """
    Run a Pipeline in parallel

    - When set_input() is called, a reader process is spawned.
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
        n_workers: int,
        buffer_size=4*1024**2
    ):
        super().__init__(pipeline)
        self._pipes = []  # the workers read from these
        self._reader_process = None
        self._outfiles = None
        self._input_path1 = None
        self._input_path2 = None
        self._interleaved_input = None
        self._n_workers = n_workers
        self._need_work_queue = Queue()
        self._buffer_size = buffer_size
        self._assign_input(infiles.file1, infiles.file2, infiles.interleaved)
        self._assign_output(outfiles)

    def _assign_input(self, file1, file2=None, interleaved=False):
        if self._reader_process is not None:
            raise RuntimeError('Do not call set_input more than once')
        self._input_path1 = file1 if type(file1) is str else file1.name
        self._input_path2 = file2 if type(file2) is str or file2 is None else file2.name
        self._interleaved_input = interleaved
        connections = [Pipe(duplex=False) for _ in range(self._n_workers)]
        self._pipes, connw = zip(*connections)
        try:
            fileno = sys.stdin.fileno()
        except io.UnsupportedOperation:
            # This happens during tests: pytest sets sys.stdin to an object
            # that does not have a file descriptor.
            fileno = -1
        self._reader_process = Process(target=reader_process, args=(file1, file2, connw,
            self._need_work_queue, self._buffer_size, fileno))
        self._reader_process.daemon = True
        self._reader_process.start()

    @staticmethod
    def can_output_to(outfiles):
        return (
            outfiles.out is not None
            and outfiles.rest is None
            and outfiles.info is None
            and outfiles.wildcard is None
            and outfiles.too_short is None
            and outfiles.too_short2 is None
            and outfiles.too_long is None
            and outfiles.too_long2 is None
            and outfiles.untrimmed is None
            and outfiles.untrimmed2 is None
            and not outfiles.demultiplex
        )

    def _assign_output(self, outfiles):
        if not self.can_output_to(outfiles):
            raise ValueError()
        self._outfiles = outfiles

    def _start_workers(self):
        workers = []
        connections = []
        for index in range(self._n_workers):
            conn_r, conn_w = Pipe(duplex=False)
            connections.append(conn_r)
            worker = WorkerProcess(
                index, self._pipeline,
                self._input_path1, self._input_path2,
                self._interleaved_input, self._outfiles,
                self._pipes[index], conn_w, self._need_work_queue)
            worker.daemon = True
            worker.start()
            workers.append(worker)
        return workers, connections

    def run(self):
        workers, connections = self._start_workers()
        writers = []
        for outfile in [self._outfiles.out, self._outfiles.out2]:
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
                        # TODO traceback should only be printed in development
                        logger.debug('%s', tb_str)
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

                    # TODO traceback should only be printed in development
                    # We should use the worker's actual traceback object
                    # here, but traceback objects are not picklable.
                    logger.debug('%s', tb_str)
                    raise e

                # No. of reads processed in this chunk
                chunk_n = connection.recv()
                if chunk_n == -2:
                    e, tb_str = connection.recv()
                    logger.debug('%s', tb_str)
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

    def __init__(self, pipeline: Pipeline, infiles: InputFiles, outfiles: OutputFiles):
        super().__init__(pipeline)
        self._pipeline.set_input(infiles)
        self._pipeline.set_output(outfiles)

    def run(self):
        (n, total1_bp, total2_bp) = self._pipeline.process_reads(progress=self._progress)
        if self._progress:
            self._progress.stop(n)
        # TODO
        return Statistics().collect(n, total1_bp, total2_bp, self._pipeline._modifiers, self._pipeline._filters)

    def close(self):
        self._pipeline.close()
