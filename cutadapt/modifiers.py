# coding: utf-8
"""
This module implements all the read modifications that cutadapt supports.
A modifier must be callable. It is implemented as a function if no parameters
need to be stored, and as a class with a __call__ method if there are parameters
(or statistics).
"""
from __future__ import print_function, division, absolute_import
import re
from cutadapt.qualtrim import quality_trim_index, nextseq_trim_index
from cutadapt.compat import maketrans
from collections import OrderedDict

class AdapterCutter(object):
	"""
	Repeatedly find one of multiple adapters in reads.
	The number of times the search is repeated is specified by the
	times parameter.
	"""

	def __init__(self, adapters, times=1, wildcard_file=None, info_file=None,
			rest_writer=None, action='trim'):
		"""
		adapters -- list of Adapter objects

		action -- What to do with a found adapter: None, 'trim', or 'mask'
		"""
		self.adapters = adapters
		self.times = times
		self.wildcard_file = wildcard_file
		self.info_file = info_file
		self.rest_writer = rest_writer
		self.action = action
		self.with_adapters = 0
		self.keep_match_info = self.info_file is not None

	def _best_match(self, read):
		"""
		Find the best matching adapter in the given read.

		Return either a Match instance or None if there are no matches.
		"""
		best = None
		for adapter in self.adapters:
			match = adapter.match_to(read)
			if match is None:
				continue

			# the no. of matches determines which adapter fits best
			if best is None or match.matches > best.matches:
				best = match
		return best

	def _write_info(self, read):
		"""
		Write to the info, wildcard and rest files.
		# TODO
		# This design with a read having a .match attribute and
		# a match having a .read attribute is really confusing.
		"""
		match = read.match
		if self.rest_writer and match:
			self.rest_writer.write(match)

		if self.wildcard_file and match:
			print(match.wildcards(), read.name, file=self.wildcard_file)

		if self.info_file:
			if read.match_info:
				for m in read.match_info:
					print(*m, sep='\t', file=self.info_file)
			else:
				seq = read.sequence
				qualities = read.qualities if read.qualities is not None else ''
				print(read.name, -1, seq, qualities, sep='\t', file=self.info_file)

	def __call__(self, read):
		"""
		Determine the adapter that best matches the given read.
		Since the best adapter is searched repeatedly, a list
		of Match instances is returned, which
		need to be applied consecutively to the read.
		The list is empty if there are no adapter matches.

		The read is converted to uppercase before it is compared to the adapter
		sequences.

		Cut found adapters from a single read. Return modified read.
		"""
		matches = []

		# try at most self.times times to remove an adapter
		trimmed_read = read
		for t in range(self.times):
			match = self._best_match(trimmed_read)
			if match is None:
				# nothing found
				break
			matches.append(match)
			trimmed_read = match.adapter.trimmed(match)
		
		if not matches:
			trimmed_read.match = None
			trimmed_read.match_info = None
			self._write_info(trimmed_read)
			return trimmed_read
		
		if __debug__:
			assert len(trimmed_read) < len(read), "Trimmed read isn't shorter than original"

		if self.action == 'trim':
			# read is already trimmed, nothing to do
			pass
		elif self.action == 'mask':
			# add N from last modification
			masked_sequence = trimmed_read.sequence
			for match in sorted(matches, reverse=True, key=lambda m: m.astart):
				ns = 'N' * (len(match.read.sequence) -
							len(match.adapter.trimmed(match).sequence))
				# add N depending on match position
				if match.front:
					masked_sequence = ns + masked_sequence
				else:
					masked_sequence += ns
			# set masked sequence as sequence with original quality
			trimmed_read.sequence = masked_sequence
			trimmed_read.qualities = matches[0].read.qualities

			assert len(trimmed_read.sequence) == len(read)
		elif self.action is None:
			trimmed_read = read
		
		trimmed_read.match = matches[-1]
		if self.keep_match_info:
			trimmed_read.match_info = [match.get_info_record() for match in matches]
		self._write_info(trimmed_read)
		
		self.with_adapters += 1
		return trimmed_read


class UnconditionalCutter(object):
	"""
	A modifier that unconditionally removes the first n or the last n bases from a read.

	If the length is positive, the bases are removed from the beginning of the read.
	If the length is negative, the bases are removed from the end of the read.
	"""
	def __init__(self, lengths):
		self.beg_length = sum(l for l in lengths if l > 0)
		self.end_length = sum(l for l in lengths if l < 0)

	def __call__(self, read):
		if self.end_length < 0:
			return read[self.beg_length:self.end_length]
		else:
			return read[self.beg_length:]


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
	def __init__(self, suffixes):
		self.suffixes = suffixes

	def __call__(self, read):
		name = read.name
		for s in self.suffixes:
			if name.endswith(s):
				name = name[:-len(s)]
		read = read[:]
		read.name = name
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
		adapter_name = 'no_adapter' if read.match is None else read.match.adapter.name
		read.name = self.prefix.replace('{name}', adapter_name) + read.name + \
			self.suffix.replace('{name}', adapter_name)
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

class PrimerTrimmer(object):
	def __call__(self, read):
		"""Trim primer base from colorspace reads"""
		read = read[1:]
		read.primer = ''
		return read


class NextseqQualityTrimmer(object):
	def __init__(self, cutoff, base):
		self.cutoff = cutoff
		self.base = base
		self.trimmed_bases = 0

	def __call__(self, read):
		stop = nextseq_trim_index(read, self.cutoff, self.base)
		self.trimmed_bases += len(read) - stop
		return read[:stop]


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


class Modifiers(object):
	def __init__(self):
		self.mod1 = OrderedDict()
		self.mod2 = OrderedDict()
	
	def add_modifier(self, mod_type, read=1, **kwargs):
		mod = ModType.create_modifier(mod_type, **kwargs)
		if read & 1 > 0:
			self.mod1[mod_type] = mod
		if read & 2 > 0:
			self.mod2[mod_type] = mod
	
	def modify(self, read1, read2=None):
		for mod in self.mod1.values():
			read1 = mod(read1)
		if read2:
			for mod in self.mod2.values():
				read2 = mod(read2)
			return (read1, read2)
		else:
			return read1


class _ModType(dict):
	"""
	Enumeration of modifier types.
	
	This is an enum-like class that is compatible with python 2x. 
	This class is only meant to be instantiated once, below (this 
	would be an anonymous class, if python allowed such a thing).
	"""
	def __init__(self):
		super(_ModType, self).__init__(
			ADAPTER					= AdapterCutter,
			CUT						= UnconditionalCutter,
			LENGTH_TAG				= LengthTagModifier,
			REMOVE_SUFFIX			= SuffixRemover,
			ADD_PREFIX_SUFFIX		= PrefixSuffixAdder,
			ZERO_CAP				= ZeroCapper,
			TRIM_QUAL				= QualityTrimmer,
			TRIM_NEXTSEQ_QUAL		= NextseqQualityTrimmer,
			CS_DOUBLE_ENCODE		= DoubleEncoder,
			CS_TRIM_PRIMER			= PrimerTrimmer,
			TRIM_END_N				= NEndTrimmer
		)
	
	def create_modifier(self, mod_type, *args, **kwargs):
		return self[mod_type](*args, **kwargs)
	
	def __getattr__(self, name):
		assert name in self
		return name

ModType = _ModType()