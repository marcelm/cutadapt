# coding: utf-8
"""
Classes for writing and filtering of processed reads.

To determine what happens to a read, a list of filters is created and each
one is called in turn (via its __call__ method) until one returns True.
The read is then assumed to have been "consumed", that is, either written
somewhere or filtered (should be discarded). Filters and writers are currently
not distinguished: The idea is that at least one of the filters will apply.
"""
from __future__ import print_function, division, absolute_import
from .xopen import xopen
from . import seqio

# Constants used when returning from a Filter’s __call__ method to improve
# readability (it is unintuitive that "return True" means "discard the read").
DISCARD = True
KEEP = False


class Filter(object):
	"""Abstract base class for filters"""
	def __init__(self, check_second=True):
		"""
		check_second -- whether the second read in a pair is also checked for
		its length. If True, the read is discarded if *any* of the two reads
		fulfills the criteria for discarding a read.
		"""
		self.check_second = check_second
		self.filtered = 0  # statistics

	def discard(self, read):
		"""
		Return True if read should be discarded (implement this in a derived class).
		"""
		raise NotImplementedError()

	def __call__(self, read1, read2=None):
		if self.discard(read1) or (
				self.check_second and read2 is not None and self.discard(read2)):
			self.filtered += 1
			return DISCARD
		return KEEP


class RedirectingFilter(Filter):
	"""
	Abstract base class for a filter that writes the reads it discards to a
	separate output file.
	"""
	def __init__(self, writer=None, check_second=True):
		super(RedirectingFilter, self).__init__(check_second)
		self.writer = writer
		self.written = 0  # no of written reads or read pairs  TODO move to writer
		self.written_bp = [0, 0]

	def __call__(self, read1, read2=None):
		if super(RedirectingFilter, self).__call__(read1, read2) == DISCARD:
			if self.writer is not None:
				self.writer.write(read1, read2)
				self.written += 1
				self.written_bp[0] += len(read1)
				if read2 is not None:
					self.written_bp[1] += len(read2)
			return DISCARD
		return KEEP


class TooShortReadFilter(RedirectingFilter):
	def __init__(self, minimum_length, too_short_writer, check_second=True):
		# TODO paired_outfile is left at its default value None (read2 is silently discarded)
		super(TooShortReadFilter, self).__init__(writer=too_short_writer, check_second=check_second)
		self.minimum_length = minimum_length

	def discard(self, read):
		return len(read) < self.minimum_length


class TooLongReadFilter(RedirectingFilter):
	def __init__(self, maximum_length, too_long_writer, check_second=True):
		super(TooLongReadFilter, self).__init__(writer=too_long_writer, check_second=check_second)
		self.maximum_length = maximum_length

	def discard(self, read):
		return len(read) > self.maximum_length


class NContentFilter(Filter):
	"""
	Discards reads over a given threshold of N's. It handles both raw counts of Ns as well
	as proportions. Note, for raw counts, it is a greater than comparison, so a cutoff
	of '1' will keep reads with a single N in it.
	"""
	def __init__(self, count, check_second=True):
		"""
		Count -- if it is below 1.0, it will be considered a proportion, and above and equal to
		1 will be considered as discarding reads with a number of N's greater than this cutoff.
		"""
		super(NContentFilter, self).__init__(check_second)
		assert count >= 0
		self.is_proportion = count < 1.0
		self.cutoff = count

	def discard(self, read):
		n_count = read.sequence.lower().count('n')
		if self.is_proportion:
			if len(read) == 0:
				return False
			return n_count / len(read) > self.cutoff
		else:
			return n_count > self.cutoff
		return False


class DiscardUntrimmedFilter(RedirectingFilter):
	"""
	A Filter that discards untrimmed reads.
	"""
	def __init__(self, untrimmed_writer, check_second=True):
		super(DiscardUntrimmedFilter, self).__init__(
			writer=untrimmed_writer,
			check_second=check_second)

	def discard(self, read):
		return read.match is None


class DiscardTrimmedFilter(RedirectingFilter):
	"""
	A filter that discards trimmed reads.
	"""
	def __init__(self, trimmed_writer, check_second=True):
		super(DiscardTrimmedFilter, self).__init__(
			writer=trimmed_writer,
			check_second=check_second)

	def discard(self, read):
		return read.match is not None


class Demultiplexer(object):
	"""
	Demultiplex trimmed reads. Reads are written to different output files
	depending on which adapter matches. Files are created when the first read
	is written to them.
	"""
	def __init__(self, path_template, untrimmed_path, fileformat, colorspace):
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
		self.fileformat = fileformat
		self.colorspace = colorspace

	def __call__(self, read1, read2=None):
		if read2 is None:
			# single-end read
			if read1.match is None:
				if self.untrimmed_writer is None and self.untrimmed_path is not None:
					self.untrimmed_writer = seqio.open(self.untrimmed_path,
						mode='w', fileformat=self.fileformat, colorspace=self.colorspace)
				if self.untrimmed_writer is not None:
					self.written += 1
					self.written_bp[0] += len(read1)
					self.untrimmed_writer.write(read1)
			else:
				name = read1.match.adapter.name
				if name not in self.writers:
					self.writers[name] = seqio.open(self.template.replace('{name}', name),
						mode='w', fileformat=self.fileformat, colorspace=self.colorspace)
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
