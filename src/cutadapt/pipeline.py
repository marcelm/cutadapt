from __future__ import print_function, division, absolute_import

import io
import os
import re
import sys
import copy
import logging
import functools
from multiprocessing import Process, Pipe, Queue
import multiprocessing.connection
import traceback

from xopen import xopen

from . import seqio
from .modifiers import ZeroCapper
from .report import Statistics
from .filters import (Redirector, PairedRedirector, NoFilter, PairedNoFilter, InfoFileWriter,
	RestFileWriter, WildcardFileWriter, TooShortReadFilter, TooLongReadFilter, NContentFilter,
	DiscardTrimmedFilter, DiscardUntrimmedFilter, Demultiplexer, PairedEndDemultiplexer)
from .seqio import read_chunks_from_file, read_paired_chunks

logger = logging.getLogger()


class OutputFiles(object):
	"""
	The attributes are open file-like objects except when demultiplex is True. In that case,
	untrimmed, untrimmed2 are file names, and out and out2 are file name templates
	containing '{name}'.
	If interleaved is True, then out is written interleaved.
	Files may also be None.
	"""
	# TODO interleaving for the other file pairs (too_short, too_long, untrimmed)?
	def __init__(self,
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


class Pipeline(object):
	"""
	Processing pipeline that loops over reads and applies modifiers and filters
	"""
	should_warn_legacy = False
	n_adapters = 0

	def __init__(self, ):
		self._close_files = []
		self._reader = None
		self._filters = []
		self._modifiers = []
		self._colorspace = None
		self._outfiles = None
		self._demultiplexer = None

	def set_input(self, file1, file2=None, qualfile=None, colorspace=False, fileformat=None,
			interleaved=False):
		self._reader = seqio.open(file1, file2, qualfile, colorspace, fileformat,
			interleaved, mode='r')
		self._colorspace = colorspace
		# Special treatment: Disable zero-capping if no qualities are available
		if not self._reader.delivers_qualities:
			self._modifiers = [m for m in self._modifiers if not isinstance(m, ZeroCapper)]

	def _open_writer(self, file, file2, **kwargs):
		# TODO backwards-incompatible change (?) would be to use outfiles.interleaved
		# for all outputs
		return seqio.open(file, file2, mode='w', qualities=self.uses_qualities,
			colorspace=self._colorspace, **kwargs)

	# TODO set max_n default to None
	def set_output(self, outfiles, minimum_length=0, maximum_length=sys.maxsize,
			max_n=-1, discard_trimmed=False, discard_untrimmed=False):
		self._filters = []
		self._outfiles = outfiles
		filter_wrapper = self._filter_wrapper()

		if outfiles.rest:
			self._filters.append(RestFileWriter(outfiles.rest))
		if outfiles.info:
			self._filters.append(InfoFileWriter(outfiles.info))
		if outfiles.wildcard:
			self._filters.append(WildcardFileWriter(outfiles.wildcard))

		too_short_writer = None
		if minimum_length > 0:
			if outfiles.too_short:
				too_short_writer = self._open_writer(outfiles.too_short, outfiles.too_short2)
			self._filters.append(
				filter_wrapper(too_short_writer, TooShortReadFilter(minimum_length)))

		too_long_writer = None
		if maximum_length < sys.maxsize:
			if outfiles.too_long:
				too_long_writer = self._open_writer(outfiles.too_long, outfiles.too_long2)
			self._filters.append(
				filter_wrapper(too_long_writer, TooLongReadFilter(maximum_length)))

		if max_n != -1:
			self._filters.append(filter_wrapper(None, NContentFilter(max_n)))

		if int(discard_trimmed) + int(discard_untrimmed) + int(outfiles.untrimmed is not None) > 1:
			raise ValueError('discard_trimmed, discard_untrimmed and outfiles.untrimmed must not '
				'be set simultaneously')

		if outfiles.demultiplex:
			self._demultiplexer = self._create_demultiplexer(outfiles)
			self._filters.append(self._demultiplexer)
		else:
			# Set up the remaining filters to deal with --discard-trimmed,
			# --discard-untrimmed and --untrimmed-output. These options
			# are mutually exclusive in order to avoid brain damage.
			if discard_trimmed:
				self._filters.append(filter_wrapper(None, DiscardTrimmedFilter()))
			elif discard_untrimmed:
				self._filters.append(filter_wrapper(None, DiscardUntrimmedFilter()))
			elif outfiles.untrimmed:
				untrimmed_writer = self._open_writer(outfiles.untrimmed, outfiles.untrimmed2)
				self._filters.append(filter_wrapper(untrimmed_writer, DiscardUntrimmedFilter()))
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

	def run(self):
		(n, total1_bp, total2_bp) = self.process_reads()
		# TODO
		m = self._modifiers
		m2 = getattr(self, '_modifiers2', [])
		stats = Statistics()
		stats.collect(n, total1_bp, total2_bp, m, m2, self._filters)
		return stats

	def process_reads(self):
		raise NotImplementedError()

	def _filter_wrapper(self):
		raise NotImplementedError()

	def _final_filter(self, outfiles):
		raise NotImplementedError()

	def _create_demultiplexer(self, outfiles):
		raise NotImplementedError()


class SingleEndPipeline(Pipeline):
	"""
	Processing pipeline for single-end reads
	"""
	paired = False

	def __init__(self):
		super(SingleEndPipeline, self).__init__()
		self._modifiers = []

	def add(self, modifier):
		self._modifiers.append(modifier)

	def add1(self, modifier):
		"""An alias for the add() function. Makes the interface similar to PairedEndPipeline"""
		self.add(modifier)

	def process_reads(self):
		"""Run the pipeline. Return statistics"""
		n = 0  # no. of processed reads  # TODO turn into attribute
		total_bp = 0
		for read in self._reader:
			n += 1
			total_bp += len(read.sequence)
			for modifier in self._modifiers:
				read = modifier(read)
			for filter in self._filters:
				if filter(read):
					break
		return (n, total_bp, None)

	def _filter_wrapper(self):
		return Redirector

	def _final_filter(self, outfiles):
		writer = self._open_writer(outfiles.out, outfiles.out2)
		return NoFilter(writer)

	def _create_demultiplexer(self, outfiles):
		return Demultiplexer(outfiles.out, outfiles.untrimmed, qualities=self.uses_qualities,
			colorspace=self._colorspace)


class PairedEndPipeline(Pipeline):
	"""
	Processing pipeline for paired-end reads.
	"""

	def __init__(self, pair_filter_mode, modify_first_read_only=False):
		"""Setting modify_first_read_only to True enables "legacy mode"
		"""
		super(PairedEndPipeline, self).__init__()
		self._modifiers2 = []
		self._pair_filter_mode = pair_filter_mode
		self._modify_first_read_only = modify_first_read_only
		self._add_both_called = False
		self._should_warn_legacy = False
		self._reader = None

	def set_input(self, *args, **kwargs):
		super(PairedEndPipeline, self).set_input(*args, **kwargs)
		if not self._reader.delivers_qualities:
			self._modifiers2 = [m for m in self._modifiers2 if not isinstance(m, ZeroCapper)]

	def add(self, modifier):
		"""
		Add a modifier for R1 and R2. If modify_first_read_only is True,
		the modifier is *not* added for R2.
		"""
		self._modifiers.append(modifier)
		if not self._modify_first_read_only:
			modifier2 = copy.copy(modifier)
			self._modifiers2.append(modifier2)
		else:
			self._should_warn_legacy = True

	def add1(self, modifier):
		"""Add a modifier for R1 only"""
		self._modifiers.append(modifier)

	def add2(self, modifier):
		"""Add a modifier for R2 only"""
		assert not self._modify_first_read_only
		self._modifiers2.append(modifier)

	def process_reads(self):
		n = 0  # no. of processed reads
		total1_bp = 0
		total2_bp = 0
		for read1, read2 in self._reader:
			n += 1
			total1_bp += len(read1.sequence)
			total2_bp += len(read2.sequence)
			for modifier in self._modifiers:
				read1 = modifier(read1)
			for modifier in self._modifiers2:
				read2 = modifier(read2)
			for filter in self._filters:
				# Stop writing as soon as one of the filters was successful.
				if filter(read1, read2):
					break
		return (n, total1_bp, total2_bp)

	@property
	def should_warn_legacy(self):
		return self._should_warn_legacy

	@should_warn_legacy.setter
	def should_warn_legacy(self, value):
		self._should_warn_legacy = bool(value)

	@property
	def paired(self):
		return 'first' if self._modify_first_read_only else 'both'

	def _filter_wrapper(self):
		return functools.partial(PairedRedirector, pair_filter_mode=self._pair_filter_mode)

	def _final_filter(self, outfiles):
		writer = self._open_writer(outfiles.out, outfiles.out2, interleaved=outfiles.interleaved)
		return PairedNoFilter(writer)

	def _create_demultiplexer(self, outfiles):
		return PairedEndDemultiplexer(outfiles.out, outfiles.out2,
			outfiles.untrimmed, outfiles.untrimmed2, qualities=self.uses_qualities,
			colorspace=self._colorspace)


def available_cpu_count():
	"""
	Return the number of available virtual or physical CPUs on this system.
	The number of available CPUs can be smaller than the total number of CPUs
	when the cpuset(7) mechanism is in use, as is the case on some cluster
	systems.

	Adapted from http://stackoverflow.com/a/1006301/715090
	"""
	try:
		with open('/proc/self/status') as f:
			status = f.read()
		m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', status)
		if m:
			res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
			if res > 0:
				return min(res, multiprocessing.cpu_count())
	except IOError:
		pass

	return multiprocessing.cpu_count()


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
		sys.stdin = os.fdopen(stdin_fd)
	try:
		with xopen(file, 'rb') as f:
			if file2:
				with xopen(file2, 'rb') as f2:
					for chunk_index, (chunk1, chunk2) in enumerate(read_paired_chunks(f, f2, buffer_size)):
						# Determine the worker that should get this chunk
						worker_index = queue.get()
						pipe = connections[worker_index]
						pipe.send(chunk_index)
						pipe.send_bytes(chunk1)
						pipe.send_bytes(chunk2)
			else:
				for chunk_index, chunk in enumerate(read_chunks_from_file(f, buffer_size)):
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
	def __init__(self, id_, pipeline, filtering_options, input_path1, input_path2,
			interleaved_input, orig_outfiles, read_pipe, write_pipe, need_work_queue):
		super(WorkerProcess, self).__init__()
		self._id = id_
		self._pipeline = pipeline
		self._filtering_options = filtering_options
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
				input = io.TextIOWrapper(io.BytesIO(data), encoding='ascii')
				input.buffer.name = self._input_path1

				if self._input_path2:
					data = self._read_pipe.recv_bytes()
					input2 = io.TextIOWrapper(io.BytesIO(data), encoding='ascii')
					input2.buffer.name = self._input_path2
				else:
					input2 = None
				output = io.TextIOWrapper(io.BytesIO(), encoding='ascii')
				output.buffer.name = self._orig_outfiles.out.name

				if self._orig_outfiles.out2 is not None:
					output2 = io.TextIOWrapper(io.BytesIO(), encoding='ascii')
					output2.buffer.name = self._orig_outfiles.out2.name
				else:
					output2 = None

				outfiles = OutputFiles(out=output, out2=output2, interleaved=self._orig_outfiles.interleaved)
				self._pipeline.set_input(input, input2, interleaved=self._interleaved_input)
				self._pipeline.set_output(outfiles, *self._filtering_options[0], **self._filtering_options[1])
				(n, bp1, bp2) = self._pipeline.process_reads()
				cur_stats = Statistics()
				cur_stats.collect(n, bp1, bp2, [], [], self._pipeline._filters)
				stats += cur_stats

				output.flush()
				processed_chunk = output.buffer.getvalue()

				self._write_pipe.send(chunk_index)
				self._write_pipe.send_bytes(processed_chunk)
				if self._orig_outfiles.out2 is not None:
					output2.flush()
					processed_chunk2 = output2.buffer.getvalue()
					self._write_pipe.send_bytes(processed_chunk2)

			m = self._pipeline._modifiers
			m2 = getattr(self._pipeline, '_modifiers2', [])
			modifier_stats = Statistics()
			modifier_stats.collect(0, 0, 0 if self._pipeline.paired else None, m, m2, [])
			stats += modifier_stats
			self._write_pipe.send(-1)
			self._write_pipe.send(stats)
		except Exception as e:
			self._write_pipe.send(-2)
			self._write_pipe.send((e, traceback.format_exc()))


class OrderedChunkWriter(object):
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
			# TODO 1) do not decode 2) use .buffer.write
			self._outfile.write(self._chunks[self._current_index].decode('utf-8'))
			del self._chunks[self._current_index]
			self._current_index += 1

	def wrote_everything(self):
		return not self._chunks


class ParallelPipelineRunner(object):
	"""
	Wrap a SingleEndPipeline, running it in parallel

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

	def __init__(self, pipeline, n_workers, buffer_size=4*1024**2):
		self._pipeline = pipeline
		self._pipes = []  # the workers read from these
		self._reader_process = None
		self._filtering_options = None
		self._outfiles = None
		self._input_path1 = None
		self._input_path2 = None
		self._interleaved_input = None
		self._n_workers = n_workers
		self._need_work_queue = Queue()
		self._buffer_size = buffer_size

	def set_input(self, file1, file2=None, qualfile=None, colorspace=False, fileformat=None,
			interleaved=False):
		if self._reader_process is not None:
			raise RuntimeError('Do not call set_input more than once')
		assert qualfile is None and colorspace is False and fileformat is None
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

	def set_output(self, outfiles, *args, **kwargs):
		if not self.can_output_to(outfiles):
			raise ValueError()
		self._filtering_options = args, kwargs
		self._outfiles = outfiles

	def _start_workers(self):
		workers = []
		connections = []
		for index in range(self._n_workers):
			conn_r, conn_w = Pipe(duplex=False)
			connections.append(conn_r)
			worker = WorkerProcess(index, self._pipeline,
				self._filtering_options, self._input_path1, self._input_path2,
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
					logger.error('%s', tb_str)
					# We should use the worker's actual traceback object
					# here, but traceback objects are not picklable.
					raise e

				for writer in writers:
					data = connection.recv_bytes()
					writer.write(data, chunk_index)
		for writer in writers:
			assert writer.wrote_everything()
		for w in workers:
			w.join()
		self._reader_process.join()
		return stats

	def close(self):
		for f in self._outfiles:
			# TODO do not use hasattr
			if f is not None and f is not sys.stdin and f is not sys.stdout and hasattr(f, 'close'):
				f.close()
