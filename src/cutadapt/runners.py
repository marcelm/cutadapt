import io
import logging
import multiprocessing
import os
import sys
import traceback
from abc import ABC, abstractmethod
from contextlib import ExitStack
from multiprocessing.connection import Connection
from typing import Any, List, Optional, Tuple, Sequence, Iterator, TYPE_CHECKING

import dnaio

from cutadapt.files import (
    InputFiles,
    OutputFiles,
    InputPaths,
    xopen_rb_raise_limit,
    detect_file_format,
    FileFormat,
    ProxyWriter,
)
from cutadapt.pipeline import Pipeline
from cutadapt.report import Statistics
from cutadapt.utils import Progress

logger = logging.getLogger()

mpctx = multiprocessing.get_context("spawn")

# See https://github.com/python/typeshed/issues/9860
if TYPE_CHECKING:
    mpctx_Process = multiprocessing.Process
else:
    mpctx_Process = mpctx.Process


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
        file_format_connection: Connection,
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
        self._file_format_connection = file_format_connection
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
                file_format = detect_file_format(files[0])
                self._file_format_connection.send(file_format)
                if file_format is not None:
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
        inpaths: InputPaths,
        proxy_files: List[ProxyWriter],
        read_pipe: Connection,
        write_pipe: Connection,
        need_work_queue: multiprocessing.Queue,
    ):
        super().__init__()
        self._id = id_
        self._pipeline = pipeline
        self._n_input_files = len(inpaths.paths)
        self._interleaved_input = inpaths.interleaved
        self._read_pipe = read_pipe
        self._write_pipe = write_pipe
        self._need_work_queue = need_work_queue
        self._proxy_files = proxy_files

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

                files = [
                    io.BytesIO(self._read_pipe.recv_bytes())
                    for _ in range(self._n_input_files)
                ]
                infiles = InputFiles(*files, interleaved=self._interleaved_input)
                (n, bp1, bp2) = self._pipeline.process_reads(infiles)
                stats += Statistics().collect(n, bp1, bp2, [], [])
                self._send_outfiles(chunk_index, n)

            stats += Statistics().collect(
                0,
                0,
                0 if self._pipeline.paired else None,
                self._pipeline._modifiers,
                self._pipeline._steps,
            )
            self._write_pipe.send(-1)
            self._write_pipe.send(stats)
        except Exception as e:
            self._write_pipe.send(-2)
            self._write_pipe.send((e, traceback.format_exc()))

    def _send_outfiles(self, chunk_index: int, n_reads: int):
        self._write_pipe.send(chunk_index)
        self._write_pipe.send(n_reads)
        for pf in self._proxy_files:
            for chunk in pf.drain():
                self._write_pipe.send_bytes(chunk)


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

    def write(self, data: bytes, index: int):
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

    @abstractmethod
    def run(self, pipeline, progress: Progress, outfiles: OutputFiles) -> Statistics:
        """
        progress: Use an object that supports .update() and .close() such
        as DummyProgress, cutadapt.utils.Progress or a tqdm instance
        """

    @abstractmethod
    def close(self):
        pass

    @abstractmethod
    def input_file_format(self) -> FileFormat:
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class ParallelPipelineRunner(PipelineRunner):
    """
    Run a Pipeline in parallel

    - At construction, a reader process is spawned.
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
        inpaths: InputPaths,
        n_workers: int,
        buffer_size: Optional[int] = None,
    ):
        self._n_workers = n_workers
        self._need_work_queue: multiprocessing.Queue = mpctx.Queue()
        self._buffer_size = 4 * 1024**2 if buffer_size is None else buffer_size
        self._inpaths = inpaths
        # the workers read from these connections
        connections = [mpctx.Pipe(duplex=False) for _ in range(self._n_workers)]
        self._connections, connw = zip(*connections)
        try:
            fileno = sys.stdin.fileno()
        except io.UnsupportedOperation:
            # This happens during tests: pytest sets sys.stdin to an object
            # that does not have a file descriptor.
            fileno = -1

        file_format_connection_r, file_format_connection_w = mpctx.Pipe(duplex=False)
        self._reader_process = ReaderProcess(
            *inpaths.paths,
            file_format_connection=file_format_connection_w,
            connections=connw,
            queue=self._need_work_queue,
            buffer_size=self._buffer_size,
            stdin_fd=fileno,
        )
        self._reader_process.daemon = True
        self._reader_process.start()
        file_format: Optional[FileFormat] = file_format_connection_r.recv()
        if file_format is None:
            raise dnaio.exceptions.UnknownFileFormat(
                f"Format of input file '{self._inpaths.paths[0]}' not recognized."
            )
        self._input_file_format = file_format

    def _start_workers(
        self, pipeline, proxy_files
    ) -> Tuple[List[WorkerProcess], List[Connection]]:
        workers = []
        connections = []
        for index in range(self._n_workers):
            conn_r, conn_w = mpctx.Pipe(duplex=False)
            connections.append(conn_r)
            worker = WorkerProcess(
                index,
                pipeline,
                self._inpaths,
                proxy_files,
                self._connections[index],
                conn_w,
                self._need_work_queue,
            )
            worker.daemon = True
            worker.start()
            workers.append(worker)
        return workers, connections

    def run(self, pipeline, progress, outfiles: OutputFiles) -> Statistics:
        workers, connections = self._start_workers(pipeline, outfiles.proxy_files())
        chunk_writers = []
        for f in outfiles.binary_files():
            chunk_writers.append(OrderedChunkWriter(f))
        stats = Statistics()
        while connections:
            ready_connections: List[Any] = multiprocessing.connection.wait(connections)
            for connection in ready_connections:
                chunk_index: int = self._try_receive(connection)
                if chunk_index == -1:
                    # the worker is done
                    cur_stats = self._try_receive(connection)
                    stats += cur_stats
                    connections.remove(connection)
                    continue

                number_of_reads: int = self._try_receive(connection)
                progress.update(number_of_reads)
                for writer in chunk_writers:
                    data = connection.recv_bytes()
                    writer.write(data, chunk_index)
        for writer in chunk_writers:
            assert writer.wrote_everything()
        for w in workers:
            w.join()
        self._reader_process.join()
        progress.close()
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
        pass

    def input_file_format(self) -> FileFormat:
        return self._input_file_format


class SerialPipelineRunner(PipelineRunner):
    """
    Run a Pipeline on a single core
    """

    def __init__(
        self,
        infiles: InputFiles,
    ):
        self._infiles = infiles
        self._input_file_format = infiles

    def run(
        self, pipeline: Pipeline, progress: Progress, outfiles: OutputFiles
    ) -> Statistics:
        (n, total1_bp, total2_bp) = pipeline.process_reads(
            self._infiles, progress=progress
        )
        if progress is not None:
            progress.close()
        # TODO
        modifiers = getattr(pipeline, "_modifiers", None)
        assert modifiers is not None
        return Statistics().collect(n, total1_bp, total2_bp, modifiers, pipeline._steps)  # type: ignore[attr-defined]

    def close(self):
        self._infiles.close()

    def input_file_format(self) -> FileFormat:
        detected = detect_file_format(self._infiles._files[0])
        if detected is None:
            raise dnaio.exceptions.UnknownFileFormat(
                f"Format of input file '{self._infiles._files[0].name}' not recognized."
            )
        return detected


def make_runner(
    inpaths: InputPaths,
    cores: int,
    buffer_size: Optional[int] = None,
) -> PipelineRunner:
    """
    Run a pipeline.

    This uses a SerialPipelineRunner if cores is 1 and a ParallelPipelineRunner otherwise.

    Args:
        inpaths:
        cores: number of cores to run the pipeline on (this is actually the number of worker
            processes, there will be one extra process for reading the input file(s))
        buffer_size: Forwarded to `ParallelPipelineRunner()`. Ignored if cores is 1.
    """
    runner: PipelineRunner
    if cores > 1:
        runner = ParallelPipelineRunner(
            inpaths,
            n_workers=cores,
            buffer_size=buffer_size,
        )
    else:
        runner = SerialPipelineRunner(inpaths.open())

    return runner
