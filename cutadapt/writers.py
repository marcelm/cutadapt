# coding: utf-8
"""
Classes for writing of processed reads
"""
from __future__ import print_function, division, absolute_import
from cutadapt.xopen import xopen


class TooShortReadFilter(object):
	def __init__(self, minimum_length, too_short_outfile, check_second=True):
		"""
		check_second -- whether the second read in a pair is also checked for
		its length. If True, the read is discarded if *any* of the two reads are
		too short.
		"""
		self.too_short_outfile = too_short_outfile
		self.minimum_length = minimum_length
		self.too_short = 0
		self.check_second = check_second

	def __call__(self, read1, read2=None):
		"""
		Return whether the read was written somewhere.
		"""
		if len(read1.sequence) < self.minimum_length or (read2 is not None and
				self.check_second and len(read2.sequence) < self.minimum_length):
			self.too_short += 1
			if self.too_short_outfile is not None:
				read1.write(self.too_short_outfile)
			# TODO read2 is silently discarded
			return True
		return False


class TooLongReadFilter(object):
	def __init__(self, maximum_length, too_long_outfile, check_second=True):
		"""
		check_second -- whether the second read in a pair is also checked for
		its length. If True, the read is discarded if *any* of the two reads are
		too long.
		"""
		self.too_long_outfile = too_long_outfile
		self.maximum_length = maximum_length
		self.too_long = 0
		self.check_second = check_second

	def __call__(self, read1, read2=None):
		if len(read1.sequence) > self.maximum_length or (read2 is not None and
				self.check_second and len(read2.sequence) > self.minimum_length):
			self.too_long += 1
			if self.too_long_outfile is not None:
				read1.write(self.too_long_outfile)
			# TODO read2 is silently discarded
			return True
		return False


class ProcessedReadWriter(object):
	"""
	Write trimmed and untrimmed reads to the proper output file(s).
	"""
	def __init__(self,
			trimmed_outfile,
			trimmed_paired_outfile,
			untrimmed_outfile,
			untrimmed_paired_outfile):
		self.trimmed_outfile = trimmed_outfile
		self.untrimmed_outfile = untrimmed_outfile
		self.trimmed_paired_outfile = trimmed_paired_outfile
		self.untrimmed_paired_outfile = untrimmed_paired_outfile

	def __call__(self, read1, read2=None):
		"""
		Write this read to the proper file.

		If read2 is not None, this is a paired-end read.
		"""
		if read2 is None:
			# single end
			if read1.match is not None and self.trimmed_outfile:
				read1.write(self.trimmed_outfile)
			if read1.match is None and self.untrimmed_outfile:
				read1.write(self.untrimmed_outfile)
		else:
			# paired end
			if read1.match is not None:  # or (not self.legacy and read2.match is not None):
				if self.trimmed_outfile:
					read1.write(self.trimmed_outfile)
				if self.trimmed_paired_outfile:
					read2.write(self.trimmed_paired_outfile)
			else:
				if self.untrimmed_outfile:
					read1.write(self.untrimmed_outfile)
				if self.untrimmed_paired_outfile:
					read2.write(self.untrimmed_paired_outfile)
		return True


class Demultiplexer(object):
	"""
	Demultiplexed trimmed reads. Reads are written to different output files
	depending on which adapter matches. Files are created when the first read
	is written to them.
	"""
	def __init__(self, path_template, untrimmed_path):
		"""
		path_template must contain the string '{name}', which will be replaced
		with the name of the adapter to form the final output path.
		Reads without an adapter match are written to the file named by
		untrimmed_path.
		"""
		assert '{name}' in path_template
		self.template = path_template
		self.untrimmed_path = untrimmed_path
		self.untrimmed_outfile = None
		self.files = dict()

	def __call__(self, read1, read2=None):
		if read2 is None:
			if read1.match is None:
				if self.untrimmed_outfile is None and self.untrimmed_path is not None:
					self.untrimmed_outfile = xopen(self.untrimmed_path, 'w')
				if self.untrimmed_outfile is not None:
					read1.write(self.untrimmed_outfile)
			else:
				name = read1.match.adapter.name
				if name not in self.files:
					self.files[name] = xopen(self.template.format(name=name), 'w')
				read1.write(self.files[name])
		else:
			assert False, "Not supported"  # pragma: no cover

	def close(self):
		for f in self.files.values():
			f.close()
		if self.untrimmed_outfile is not None:
			self.untrimmed_outfile.close()


class NContentTrimmer(object):
	"""
	Discards reads over a given threshold of N's. It handles both raw counts of Ns as well
	as proportions. Note, for raw counts, it is a greater than comparison, so a cutoff
	of '1' will keep reads with a single N in it.
	"""
	def __init__(self, count):
		"""
		Count -- if it is below 1.0, it will be considered a proportion, and above and equal to
		1 will be considered as discarding reads with a number of N's greater than this cutoff.
		"""
		self.proportion = count < 1.0
		self.cutoff = count
		self.too_many_n = 0

	def __call__(self, read1, read2=None):
		if self.proportion:
			if len(read1.sequence) and read1.sequence.lower().count('n') / len(read1.sequence) > self.cutoff or \
					(read2 is not None and len(read2.sequence) and read2.sequence.lower().count('n')/len(read2.sequence) > self.cutoff):
				self.too_many_n += 1
				return True
		else:
			if read1.sequence.lower().count('n') > self.cutoff or \
					(read2 is not None and read2.sequence.lower().count('n') > self.cutoff):
				self.too_many_n += 1
				return True
		return False
