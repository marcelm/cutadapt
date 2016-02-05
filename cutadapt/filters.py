# coding: utf-8
"""
Classes for writing and filtering of processed reads.

A Filter is a callable that has the read as its only argument. If it is called,
it returns True if the read should be filtered (discarded), and False if not.

To be used, a filter needs to be wrapped in one of the redirector classes.
They are called so because they can redirect filtered reads to a file if so
desired. They also keep statistics.

To determine what happens to a read, a list of redirectors with different
filters is created and each redirector is called in turn until one returns True.
The read is then assumed to have been "consumed", that is, either written
somewhere or filtered (should be discarded).
"""
from __future__ import print_function, division, absolute_import
from .xopen import xopen
from . import seqio

# Constants used when returning from a Filter’s __call__ method to improve
# readability (it is unintuitive that "return True" means "discard the read").
DISCARD = True
KEEP = False


class NoFilter(object):
	"""
	No filtering, just send each read to the given writer.
	"""
	def __init__(self, writer):
		self.filtered = 0
		self.writer = writer
		self.filter = filter
		self.written = 0  # no of written reads  TODO move to writer
		self.written_bp = [0, 0]

	def __call__(self, read):
		self.writer.write(read)
		self.written += 1
		self.written_bp[0] += len(read)
		return DISCARD


class PairedNoFilter(object):
	"""
	No filtering, just send each paired-end read to the given writer.
	"""
	def __init__(self, writer):
		self.filtered = 0
		self.writer = writer
		self.written = 0  # no of written reads or read pairs  TODO move to writer
		self.written_bp = [0, 0]

	def __call__(self, read1, read2):
		self.writer.write(read1, read2)
		self.written += 1
		self.written_bp[0] += len(read1)
		self.written_bp[1] += len(read2)
		return DISCARD


class Redirector(object):
	"""
	Redirect discarded reads to the given writer. This is for single-end reads.
	"""
	def __init__(self, writer, filter):
		self.filtered = 0
		self.writer = writer
		self.filter = filter
		self.written = 0  # no of written reads  TODO move to writer
		self.written_bp = [0, 0]

	def __call__(self, read):
		if self.filter(read):
			self.filtered += 1
			if self.writer is not None:
				self.writer.write(read)
				self.written += 1
				self.written_bp[0] += len(read)
			return DISCARD
		return KEEP


class PairedRedirector(object):
	"""
	Redirect discarded reads to the given writer. This is for paired-end reads,
	using the 'new-style' filtering where both reads are inspected. That is,
	the entire pair is discarded if at least 1 or 2 of the reads match the
	filtering criteria.
	"""
	def __init__(self, writer, filter, min_affected=1):
		"""
		min_affected -- values 1 and 2 are allowed.
			1 means: the pair is discarded if any read matches
			2 means: the pair is discarded if both reads match
		"""
		if not min_affected in (1, 2):
			raise ValueError("min_affected must be 1 or 2")
		self.filtered = 0
		self.writer = writer
		self.filter = filter
		self._min_affected = min_affected
		self.written = 0  # no of written reads or read pairs  TODO move to writer
		self.written_bp = [0, 0]

	def __call__(self, read1, read2):
		if self.filter(read1) + self.filter(read2) >= self._min_affected:
			self.filtered += 1
			# discard read
			if self.writer is not None:
				self.writer.write(read1, read2)
				self.written += 1
				self.written_bp[0] += len(read1)
				self.written_bp[1] += len(read2)
			return DISCARD
		return KEEP


class LegacyPairedRedirector(object):
	"""
	Redirect discarded reads to the given writer. This is for paired-end reads,
	using the 'legacy' filtering mode (backwards compatibility). That is, if
	the first read matches the filtering criteria, the pair is discarded. The
	second read is not inspected.
	"""
	def __init__(self, writer, filter):
		self.filtered = 0
		self.writer = writer
		self.filter = filter
		self.written = 0  # no of written reads or read pairs  TODO move to writer
		self.written_bp = [0, 0]

	def __call__(self, read1, read2):
		if self.filter(read1):
			self.filtered += 1
			# discard read
			if self.writer is not None:
				self.writer.write(read1, read2)
				self.written += 1
				self.written_bp[0] += len(read1)
				self.written_bp[1] += len(read2)
			return DISCARD
		return KEEP


class TooShortReadFilter(object):
	# TODO paired_outfile is left at its default value None (read2 is silently discarded)
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
	Discards a reads that has a number of 'N's over a given threshold. It handles both raw counts of Ns as well
	as proportions. Note, for raw counts, it is a greater than comparison, so a cutoff
	of '1' will keep reads with a single N in it.
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


class DiscardUntrimmedFilter(object):
	"""
	Return True if read is untrimmed.
	"""
	def __call__(self, read):
		return read.match is None


class DiscardTrimmedFilter(object):
	"""
	Return True if read is trimmed.
	"""
	def __call__(self, read):
		return read.match is not None


class Demultiplexer(object):
	"""
	Demultiplex trimmed reads. Reads are written to different output files
	depending on which adapter matches. Files are created when the first read
	is written to them.
	"""
	def __init__(self, path_template, untrimmed_path, colorspace, qualities):
		"""
		path_template must contain the string '{name}', which will be replaced
		with the name of the adapter to form the final output path.
		Reads without an adapter match are written to the file named by
		untrimmed_path.
		"""
		assert '{name}' in path_template
		self.template = path_template
		self.untrimmed_path = untrimmed_path
		self.untrimmed_writer = None
		self.writers = dict()
		self.written = 0
		self.written_bp = [0, 0]
		self.colorspace = colorspace
		self.qualities = qualities

	def __call__(self, read1, read2=None):
		if read2 is None:
			# single-end read
			if read1.match is None:
				if self.untrimmed_writer is None and self.untrimmed_path is not None:
					self.untrimmed_writer = seqio.open(self.untrimmed_path,
						mode='w', colorspace=self.colorspace, qualities=self.qualities)
				if self.untrimmed_writer is not None:
					self.written += 1
					self.written_bp[0] += len(read1)
					self.untrimmed_writer.write(read1)
			else:
				name = read1.match.adapter.name
				if name not in self.writers:
					self.writers[name] = seqio.open(self.template.replace('{name}', name),
						mode='w', colorspace=self.colorspace, qualities=self.qualities)
				self.written += 1
				self.written_bp[0] += len(read1)
				self.writers[name].write(read1)
			return DISCARD
		else:
			assert False, "Not supported"  # pragma: no cover

	def close(self):
		for w in self.writers.values():
			w.close()
		if self.untrimmed_writer is not None:
			self.untrimmed_writer.close()
