# coding: utf-8
"""
Classes for writing and filtering of processed reads.

A Filter is a callable that has the read as its only argument. If it is called,
it returns True if the read should be filtered (discarded), and False if not. A
filter also has a code used to identify the filter that rejected a read. To be 
used, a filter needs to be wrapped in one of the Wrapper classes. To determine 
what happens to a read, a list of Wrappers with different filters is created and 
each Wrapper is called in turn until one returns True. The main program will
determine whether and where to write the read(s) based on whether it was rejected
by a filter, and which one.
"""
from __future__ import print_function, division, absolute_import
from .xopen import xopen
from . import seqio
from enum import Enum

DISCARD = True
KEEP = False

class SingleWrapper(object):
	"""
	This is for single-end reads and for paired-end reads, using the 'legacy' filtering mode 
	(backwards compatibility). That is, if the first read matches the filtering criteria, the 
	pair is discarded. The second read is not inspected.
	"""
	def __init__(self, f):
		self.filtered = 0
		self.filter = f
		
	def __call__(self, read1, read2=None):
		if self.filter(read1):
			self.filtered += 1
			return DISCARD
		return KEEP

class PairedWrapper(object):
	"""
	This is for paired-end reads, using the 'new-style' filtering where both reads are inspected. 
	That is, the entire pair is discarded if at least 1 or 2 of the reads match the
	filtering criteria.
	"""
	def __init__(self, f, min_affected=1):
		"""
		min_affected -- values 1 and 2 are allowed.
			1 means: the pair is discarded if any read matches
			2 means: the pair is discarded if both reads match
		"""
		if not min_affected in (1, 2):
			raise ValueError("min_affected must be 1 or 2")
		self.filtered = 0
		self.filter = f
		self.min_affected = min_affected

	def __call__(self, read1, read2):
		failures = 0
		if self.filter(read1):
			failures += 1
		if (self.min_affected - failures == 1) and self.filter(read2):
			failures += 1
		if failures >= self.min_affected:
			self.filtered += 1
			return DISCARD
		return KEEP

class FilterFactory(object):
	def __init__(self, paired, min_affected):
		self.paired = paired
		self.min_affected = min_affected
	
	def __call__(self, filter_type, *args, **kwargs):
		f = filter_type.value(*args, **kwargs)
		if self.paired == "both":
			return PairedWrapper(f, self.min_affected)
		else:
			return SingleWrapper(f)

class TooShortReadFilter(object):
	def __init__(self, minimum_length):
		self.minimum_length = minimum_length
	
	def __call__(self, read):
		return len(read) < self.minimum_length

class TooLongReadFilter(object):
	def __init__(self, maximum_length):
		self.maximum_length = maximum_length
	
	def __call__(self, read):
		return len(read) > self.maximum_length

class NContentFilter(object):
	"""
	Discards a reads that has a number of 'N's over a given threshold. It handles both raw 
	counts of Ns as well as proportions. Note, for raw counts, it is a greater than comparison, 
	so a cutoff of '1' will keep reads with a single N in it.
	"""
	def __init__(self, count):
		"""
		Count -- if it is below 1.0, it will be considered a proportion, and above and equal to
		1 will be considered as discarding reads with a number of N's greater than this cutoff.
		"""
		assert count >= 0
		self.is_proportion = count < 1.0
		self.cutoff = count

	def __call__(self, read):
		"""Return True when the read should be discarded"""
		n_count = read.sequence.lower().count('n')
		if self.is_proportion:
			if len(read) == 0:
				return False
			return n_count / len(read) > self.cutoff
		else:
			return n_count > self.cutoff

class UntrimmedFilter(object):
	"""
	Return True if read is untrimmed.
	"""
	def __call__(self, read):
		return read.match is None

class TrimmedFilter(object):
	"""
	Return True if read is trimmed.
	"""
	def __call__(self, read):
		return read.match is not None

class NoFilter(object):
	def __call__(self, read):
		return False

# Enum of filter types
class FilterType(Enum):
	TOO_SHORT	= TooShortReadFilter
	TOO_LONG	= TooLongReadFilter
	N_CONTENT	= NContentFilter
	TRIMMED		= TrimmedFilter
	UNTRIMMED	= UntrimmedFilter
	NONE		= NoFilter
