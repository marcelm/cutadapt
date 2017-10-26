from __future__ import print_function, division, absolute_import

import io
import sys
import time
import logging
import functools
from multiprocessing import Process, Pipe, Queue
import multiprocessing.connection

from xopen import xopen

from . import seqio
from .modifiers import ZeroCapper
from .report import Statistics
from .filters import (Redirector, PairedRedirector, NoFilter, PairedNoFilter, InfoFileWriter,
	RestFileWriter, WildcardFileWriter, TooShortReadFilter, TooLongReadFilter, NContentFilter,
	DiscardTrimmedFilter, DiscardUntrimmedFilter, Demultiplexer, PairedEndDemultiplexer)

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
		start_time = time.clock()
		(n, total1_bp, total2_bp) = self.process_reads()
		#self.close_files()
		elapsed_time = time.clock() - start_time
		# TODO
		m = self._modifiers
		m2 = getattr(self, '_modifiers2', [])
		stats = Statistics()
		stats.collect(n, total1_bp, total2_bp, elapsed_time, m, m2, self._filters)
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
			self._modifiers2.append(modifier)
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
		return self._should_warn_legacy or self._modify_first_read_only and len(self._filters) > 1

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


BUFFER_SIZE = 4000000


def find_fasta_record_end(buf, end):
	"""
	Search for the end of the last complete FASTA record within buf[start:end]
	"""
	pos = buf.rfind(b'\n>', 0, end)
	if pos != -1:
		return pos + 1
	if buf[0:1] == b'>':
		return 0
	# TODO
	assert False


def find_fastq_record_end(buf, end, _newline=ord('\n')):
	"""
	Search for the end of the last complete FASTQ record in buf[:end]
	"""
	linebreaks = buf.count(_newline, 0, end)
	right = end
	for _ in range(linebreaks % 4 + 1):
		right = buf.rfind(_newline, 0, right)
		assert right != -1  # TODO
	return right + 1


def read_chunks_from_file(f, buffer_size=4000000):
	"""
	f needs to be a file opened in binary mode
	"""
	# This buffer is re-used in each iteration.
	buf = bytearray(buffer_size)

	# Read one byte to determine file format
	# TODO if there is a comment char, we assume FASTA
	start = f.readinto(memoryview(buf)[0:1])
	if start == 1 and buf[0:1] == b'@':
		find_record_end = find_fastq_record_end
	elif start == 1 and buf[0:1] == b'#' or buf[0:1] == b'>':
		find_record_end = find_fasta_record_end
	elif start > 0:
		raise ValueError('input file format unknown')

	while True:
		bufend = f.readinto(memoryview(buf)[start:]) + start
		if start == bufend:
			# End of file
			break
		end = find_record_end(buf, bufend)
		assert end <= bufend
		if end > 0:
			yield memoryview(buf)[0:end]
		start = bufend - end
		assert start >= 0
		buf[0:start] = buf[end:bufend]

	if start > 0:
		yield memoryview(buf)[0:start]


def reader_process(file, connections, queue):
	"""
	Read chunks of FASTA or FASTQ data from *file* and send to a worker.

	queue -- a Queue of worker indices. A worker writes its own index into this
		queue to notify us that it is ready to receive more data.
	connections -- a list of Connection objects, one for each worker.

	The function repeatedly

	- reads a chunk from the file
	- reads a worker index from the Queue
	- sends the chunk to connections[index]

	and finally sends "poison pills" (the value -1) to all connections.
	"""
	with xopen(file, 'rb') as f:
		for chunk_index, chunk in enumerate(read_chunks_from_file(f)):
			# Determine the worker that should get this chunk
			worker_index = queue.get()
			pipe = connections[worker_index]
			pipe.send(chunk_index)
			pipe.send_bytes(chunk)

	# Send poison pills to all workers
	for _ in range(len(connections)):
		worker_index = queue.get()
		connections[worker_index].send(-1)

	# try:
	#   ...
	# except Exception as e:
	# 	traceb = traceback.format_exc()
	# 	reads_queue.put((e, traceb))


class WorkerProcess(Process):
	"""
	The worker repeatedly reads chunks of data from the read_pipe, runs the pipeline on it
	and sends the processed chunks to the write_pipe.

	To notify the reader process that it wants data, it puts its own identifier into the queue
	before attempting to read data from the read_pipe.
	"""
	def __init__(self, id_, pipeline, filtering_options, orig_outfiles, read_pipe, write_pipe,
			need_work_queue):
		super(WorkerProcess, self).__init__()
		self._id = id_
		self._pipeline = pipeline
		self._filtering_options = filtering_options
		self._orig_outfiles = orig_outfiles
		self._read_pipe = read_pipe
		self._write_pipe = write_pipe
		self._need_work_queue = need_work_queue

	def run(self):
		n = 0  # no. of processed reads  # TODO turn into attribute
		total_bp = 0

		while True:
			# Notify reader that we need data
			self._need_work_queue.put(self._id)
			chunk_index = self._read_pipe.recv()
			if chunk_index == -1:
				# reader is done
				break
			data = self._read_pipe.recv_bytes()

			# logger.info('WORKER: Read %d bytes', len(data))
			input = io.TextIOWrapper(io.BytesIO(data), encoding='ascii')
			output = io.TextIOWrapper(io.BytesIO(), encoding='ascii')
			# Output format depends on file name, so make output file name available
			output.buffer.name = self._orig_outfiles.out.name

			outfiles = OutputFiles(out=output)
			self._pipeline.set_input(input)
			self._pipeline.set_output(outfiles, *self._filtering_options[0], **self._filtering_options[1])
			stats = self._pipeline.run()  # TODO accumulate stats
			output.flush()
			processed_chunk = output.buffer.getvalue()

			self._write_pipe.send((chunk_index, stats))  # TODO does not work
			self._write_pipe.send_bytes(processed_chunk)


		self._write_pipe.send(-1)


class ParallelPipelineRunner(object):
	"""
	Wrap a SingleEndPipeline, running it in parallel


	"""

	def __init__(self, pipeline, n_workers):
		self._pipeline = pipeline
		self._pipes = []  # the workers read from these
		self._reader_process = None
		self._filtering_options = None
		self._outfiles = None
		self._n_workers = n_workers
		self._need_work_queue = Queue()

	def set_input(self, file1, file2=None, qualfile=None, colorspace=False, fileformat=None,
			interleaved=False):
		if self._reader_process is not None:
			raise RuntimeError('Do not call set_input more than once')
		assert file2 is None and qualfile is None and colorspace is False and fileformat is None and interleaved is False
		connections = [Pipe(duplex=False) for _ in range(self._n_workers)]
		self._pipes, connw = zip(*connections)
		p = Process(target=reader_process, args=(file1, connw, self._need_work_queue))
		p.daemon = True
		p.start()
		self._reader_process = p

	@staticmethod
	def can_output_to(outfiles):
		return (
			outfiles.out is not None
			and outfiles.out2 is None
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
			and not outfiles.interleaved
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
				self._filtering_options, self._outfiles, self._pipes[index], conn_w,
				self._need_work_queue)
			worker.daemon = True
			worker.start()
			workers.append(worker)
		return workers, connections

	def run(self):
		start_time = time.clock()
		workers, connections = self._start_workers()
		chunks = dict()
		current_chunk = 0
		while connections:
			ready_connections = multiprocessing.connection.wait(connections)
			for connection in ready_connections:
				chunk_index = connection.recv()
				if chunk_index == -1:
					# the worker is done
					connections.remove(connection)
					continue
				data = connection.recv_bytes()
				chunks[chunk_index] = data
				while current_chunk in chunks:
					self._outfiles.out.write(chunks[current_chunk].decode('utf-8'))  # TODO could be written as binary
					del chunks[current_chunk]
					current_chunk += 1

		for w in workers:
			w.join()
		self._reader_process.join()

		(n, total1_bp, total2_bp) = (0, 0, 0)  # TODO
		#self._pipeline.close()
		elapsed_time = time.clock() - start_time
		# TODO
		m = self._pipeline._modifiers
		m2 = getattr(self, '_modifiers2', [])
		stats = Statistics()
		stats.collect(n, total1_bp, total2_bp, elapsed_time, m, m2, self._pipeline._filters)
		return stats

	def close(self):
		for f in self._outfiles:
			# TODO do not use hasattr
			if f is not None and f is not sys.stdin and f is not sys.stdout and hasattr(f, 'close'):
				f.close()
