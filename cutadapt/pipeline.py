from __future__ import print_function, division, absolute_import

import sys
import time
import logging

from . import seqio
from .modifiers import ZeroCapper
from .report import Statistics

logger = logging.getLogger()


class Pipeline(object):
	"""
	Processing pipeline that loops over reads and applies modifiers and filters
	"""
	should_warn_legacy = False
	n_adapters = 0

	def __init__(self):
		self._close_files = []
		self._reader = None
		self._filters = []
		self._modifiers = []

	def set_filters(self, filters):
		self._filters = filters

	def open_input(self, file1, file2=None, qualfile=None, colorspace=False, fileformat=None,
			interleaved=False):
		self._reader = seqio.open(file1, file2, qualfile, colorspace, fileformat,
			interleaved, mode='r')
		# Special treatment: Disable zero-capping if no qualities are available
		if not self._reader.delivers_qualities:
			self._modifiers = [m for m in self._modifiers if not isinstance(m, ZeroCapper)]

	@property
	def uses_qualities(self):
		return self._reader.delivers_qualities

	def register_file_to_close(self, file):
		if file is not None and file is not sys.stdin and file is not sys.stdout:
			self._close_files.append(file)

	def close_files(self):
		for f in self._close_files:
			f.close()

	def process_reads(self):
		raise NotImplementedError()

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


class PairedEndPipeline(Pipeline):
	"""
	Processing pipeline for paired-end reads.
	"""
	def __init__(self, modify_first_read_only):
		"""Setting modify_first_read_only to True enables "legacy mode"
		"""
		super(PairedEndPipeline, self).__init__()
		self._modifiers2 = []
		self._modify_first_read_only = modify_first_read_only
		self._add_both_called = False
		self._should_warn_legacy = False
		self._reader = None

	def open_input(self, *args, **kwargs):
		super(PairedEndPipeline, self).open_input(*args, **kwargs)
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
