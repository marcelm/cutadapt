from __future__ import print_function, division, absolute_import

import sys
import time
import logging
import functools
from collections import namedtuple

from . import seqio
from .modifiers import ZeroCapper
from .report import Statistics
from .filters import (Redirector, PairedRedirector, NoFilter, PairedNoFilter, InfoFileWriter,
	RestFileWriter, WildcardFileWriter, TooShortReadFilter, TooLongReadFilter, NContentFilter,
	DiscardTrimmedFilter, DiscardUntrimmedFilter, Demultiplexer, PairedEndDemultiplexer)


logger = logging.getLogger()


# The attributes are open file-like objects except when demultiplex is True. In that case,
# untrimmed, untrimmed2 are file names, and out and out2 are file name templates
# containing '{name}'.
# If interleaved is True, then out is written interleaved
# TODO interleaving for the other file pairs (too_short, too_long, untrimmed)?
# Files may also be set to None if not specified on the command-line.
# TODO move to pipeline
OutputFiles = namedtuple('OutputFiles',
	'rest info wildcard too_short too_short2 too_long too_long2 untrimmed untrimmed2 out out2 '
	'demultiplex interleaved')


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

	def register_file_to_close(self, file):
		if file is not None and file is not sys.stdin and file is not sys.stdout:
			self._close_files.append(file)

	def close_files(self):
		for f in self._close_files:
			f.close()

	def run(self):
		start_time = time.clock()
		(n, total1_bp, total2_bp) = self.process_reads()
		self.close_files()
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
