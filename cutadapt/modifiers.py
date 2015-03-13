# coding: utf-8
from __future__ import print_function, division, absolute_import
import re
from cutadapt.qualtrim import quality_trim_index
from cutadapt.compat import maketrans


class UnconditionalCutter(object):
	"""
	A modifier that unconditionally removes the first n or the last n bases from a read.

	If the length is positive, the bases are removed from the beginning of the read.
	If the length is negative, the bases are removed from the end of the read.
	"""
	def __init__(self, length):
		self.length = length

	def __call__(self, read):
		if self.length > 0:
			return read[self.length:]
		elif self.length < 0:
			return read[:self.length]


class LengthTagModifier(object):
	"""
	Replace "length=..." strings in read names.
	"""
	def __init__(self, length_tag):
		self.regex = re.compile(r"\b" + length_tag + r"[0-9]*\b")
		self.length_tag = length_tag

	def __call__(self, read):
		read = read[:]
		if read.name.find(self.length_tag) >= 0:
			read.name = self.regex.sub(self.length_tag + str(len(read.sequence)), read.name)
		return read


class SuffixRemover(object):
	"""
	Remove a given suffix from read names.
	"""
	def __init__(self, suffix):
		self.suffix = suffix

	def __call__(self, read):
		read = read[:]
		if read.name.endswith(self.suffix):
			read.name = read.name[:-len(self.suffix)]
		return read


class PrefixSuffixAdder(object):
	"""
	Add a suffix and a prefix to read names
	"""
	def __init__(self, prefix, suffix):
		self.prefix = prefix
		self.suffix = suffix

	def __call__(self, read):
		read = read[:]
		read.name = self.prefix + read.name + self.suffix
		return read


class DoubleEncoder(object):
	"""
	Double-encode colorspace reads, using characters ACGTN to represent colors.
	"""
	def __init__(self):
		self.double_encode_trans = maketrans('0123.', 'ACGTN')

	def __call__(self, read):
		read = read[:]
		read.sequence = read.sequence.translate(self.double_encode_trans)
		return read


class ZeroCapper(object):
	"""
	Change negative quality values of a read to zero
	"""
	def __init__(self, quality_base=33):
		qb = quality_base
		self.zero_cap_trans = maketrans(''.join(map(chr, range(qb))), chr(qb) * qb)

	def __call__(self, read):
		read = read[:]
		read.qualities = read.qualities.translate(self.zero_cap_trans)
		return read


def PrimerTrimmer(read):
	"""Trim primer base from colorspace reads"""
	read = read[1:]
	read.primer = ''
	return read


class QualityTrimmer(object):
	def __init__(self, cutoff_front, cutoff_back, base):
		self.cutoff_front = cutoff_front
		self.cutoff_back = cutoff_back
		self.base = base
		self.trimmed_bases = 0

	def __call__(self, read):
		start, stop = quality_trim_index(read.qualities, self.cutoff_front, self.cutoff_back, self.base)
		self.trimmed_bases += len(read) - (stop - start)
		return read[start:stop]


class NEndTrimmer(object):
	"""Trims Ns from the 3' and 5' end of reads"""
	def __init__(self):
		self.start_trim = re.compile(r'^N+')
		self.end_trim = re.compile(r'N+$')

	def __call__(self, read):
		sequence = read.sequence
		start_cut = self.start_trim.match(sequence)
		end_cut = self.end_trim.search(sequence)
		start_cut = start_cut.end() if start_cut else 0
		end_cut = end_cut.start() if end_cut else len(read)
		return read[start_cut:end_cut]
