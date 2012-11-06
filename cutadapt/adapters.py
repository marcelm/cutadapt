#!/usr/bin/env python
from __future__ import print_function, division
from collections import namedtuple, defaultdict
from . import align, colorspace
from .seqio import ColorspaceSequence

# Constants for the find_best_alignment function.
# The function is called with SEQ1 as the adapter, SEQ2 as the read.
BACK = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ1
FRONT = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.START_WITHIN_SEQ1
PREFIX = align.STOP_WITHIN_SEQ2
ANYWHERE = align.SEMIGLOBAL


class AdapterMatch:
	def __init__(self, astart, astop, rstart, rstop, matches, errors, front, adapter, read):
		self.astart, self.astop, self.rstart, self.rstop = astart, astop, rstart, rstop
		self.matches = matches
		self.errors = errors
		self.front = self._guess_is_front() if front is None else front
		self.adapter = adapter
		self.read = read

	def _guess_is_front(self):
		"""
		Return whether this is guessed to be a front adapter.

		The match is assumed to be a front adapter when the first base of
		the read is involved in the alignment to the adapter.
		"""
		# TODO remove
		# if match.rstart != 0  ==>  match.astart == 0
		assert self.rstart == 0 or self.astart == 0
		#return not (match.rstart > 0 and match.astart == 0)
		return self.rstart == 0

	def wildcards(self, wildcard_char='N'):
		"""
		Return a string that contains, for each wildcard character,
		the character that it matches. For example, if the adapter
		ATNGNA matches ATCGTA, then the string 'CT' is returned.

		If there are indels, this is not reliable as the full alignment
		is not available.
		"""
		wildcards = [ self.read.sequence[self.rstart + i] for i in range(self.length)
			if self.adapter.sequence[self.astart + i] == wildcard_char and self.rstart + i < len(self.read.sequence) ]
		return ''.join(wildcards)

	@property
	def length(self):
		"""
		Number of aligned characters in the adapter. If there are
		indels, this may be different from the number of characters
		in the read.
		"""
		return self.astop - self.astart

	def rest(self):
		"""
		Return the part of the read before this match if this is a
		'front' (5') adapter,
		return the part after the match if this is not a 'front' adapter (3').
		This can be an empty string.
		"""
		if self.front:
			return self.read.sequence[:self.rstart]
		else:
			return self.read.sequence[self.rstop:]


class Adapter(object):
	"""
	An adapter knows how to match itself to a read.
	In particular, it knows where it should be within the read and how to interpret
	wildcard characters.

	where --  One of the BACK, FRONT, PREFIX or ANYWHERE constants.
		If the adapter is located in the middle of the read,
		the constant influences which part of the read gets removed.

	sequence -- The adapter sequence as string. Will be converted to uppercase.

	max_error_rate -- Maximum allowed error rate. The error rate is
		the number of errors in the alignment divided by the length
		of the part of the alignment that matches the adapter.

	minimum_overlap -- Minimum length of the part of the alignment
		that matches the adapter.

	match_read_wildcards -- Whether wildcards ('N' characters) in the read
		are matches (at zero cost).

	match_adapter_wildcards -- Whether wildcards in the adapter are allowed
		to match any character in the read (at zero cost).

	name -- optional name of the adapter
	"""
	def __init__(self, sequence, where, max_error_rate, min_overlap=3,
			match_read_wildcards=False, match_adapter_wildcards=False,
			name=None):
		self.name = name
		self.sequence = sequence.upper()
		self.where = where
		self.max_error_rate = max_error_rate
		self.min_overlap = min_overlap
		self.wildcard_flags = 0
		self.match_adapter_wildcards = match_adapter_wildcards and 'N' in self.sequence
		if match_read_wildcards:
			self.wildcard_flags |= align.ALLOW_WILDCARD_SEQ2
		if self.match_adapter_wildcards:
			self.wildcard_flags |= align.ALLOW_WILDCARD_SEQ1
		# redirect to appropriate trimmed() function depending on
		# adapter type
		trimmers = {
			FRONT: self._trimmed_front,
			PREFIX: self._trimmed_front,
			BACK: self._trimmed_back,
			ANYWHERE: self._trimmed_anywhere
		}
		self.trimmed = trimmers[where]
		if where == ANYWHERE:
			self._front_flag = None # means: guess
		else:
			self._front_flag = where != BACK
		# statistics about length of removed sequences
		self.lengths_front = defaultdict(int)
		self.lengths_back = defaultdict(int)

	def __repr__(self):
		match_read_wildcards = bool(align.ALLOW_WILDCARD_SEQ2 & self.wildcard_flags)
		match_adapter_wildcards = bool(align.ALLOW_WILDCARD_SEQ1 & self.wildcard_flags)
		return '<Adapter(name="{name}", sequence="{sequence}", where={where}, '\
			'max_error_rate={max_error_rate}, min_overlap={min_overlap}, '\
			'match_read_wildcards={match_read_wildcards}, '\
			'match_adapter_wildcards={match_adapter_wildcards})>'.format(
				match_read_wildcards=match_read_wildcards,
				match_adapter_wildcards=match_adapter_wildcards,
				**vars(self))

	def match(self, read):
		"""
		Try to match this adapter to the given read and return an AdapterMatch instance.

		Return None if the minimum overlap length is not met or the error rate is too high.
		"""
		read_seq = read.sequence.upper()
		pos = -1
		# try to find an exact match first unless wildcards are allowed
		if not self.match_adapter_wildcards:
			if self.where == PREFIX:
				pos = 0 if read_seq.startswith(self.sequence) else -1
			else:
				pos = read_seq.find(self.sequence)
		if pos >= 0:
			match = AdapterMatch(
				0, len(self.sequence), pos, pos + len(self.sequence),
				len(self.sequence), 0, self._front_flag, self, read)
		else:
			# try approximate matching
			alignment = align.globalalign_locate(self.sequence, read_seq,
				self.max_error_rate, self.where, self.wildcard_flags)
			match = AdapterMatch(*(alignment + (self._front_flag, self, read)))

		# TODO globalalign_locate should be modified to allow the following
		# assertion.
		# assert length == 0 or match.errors / length <= self.max_error_rate
		if match.length < self.min_overlap or match.errors / match.length > self.max_error_rate:
			return None
		return match

	def _trimmed_anywhere(self, match):
		"""Return a trimmed read"""
		if match.front:
			return self._trimmed_front(match)
		else:
			return self._trimmed_back(match)

	def _trimmed_front(self, match):
		"""Return a trimmed read"""
		# TODO move away
		self.lengths_front[match.rstop] += 1
		return match.read[match.rstop:]

	def _trimmed_back(self, match):
		"""Return a trimmed read without the 3' (back) adapter"""
		# TODO move away
		self.lengths_back[len(match.read) - match.rstart] += 1
		return match.read[:match.rstart]

	def __len__(self):
		return len(self.sequence)


class ColorspaceAdapter(Adapter):
	def __init__(self, *args, **kwargs):
		super(ColorspaceAdapter, self).__init__(*args, **kwargs)
		has_nucleotide_seq = False
		if set(self.sequence) <= set('ACGT'):
			# adapter was given in basespace
			self.nucleotide_sequence = self.sequence
			has_nucleotide_seq = True
			self.sequence = colorspace.encode(self.sequence)[1:]
		if self.where in (PREFIX, FRONT) and not has_nucleotide_seq:
			raise ValueError("A 5' colorspace adapter needs to be given in nucleotide space")

	def match(self, read):
		"""Return AdapterMatch instance"""
		if self.where != PREFIX:
			return super(ColorspaceAdapter, self).match(read)
		# create artificial adapter that includes a first color that encodes the
		# transition from primer base into adapter
		asequence = colorspace.ENCODE[read.primer + self.nucleotide_sequence[0]] + self.sequence
		pos = 0 if read.sequence.startswith(asequence) else -1
		if pos >= 0:
			match = AdapterMatch(
				0, len(asequence), pos, pos + len(asequence),
				len(asequence), 0, self._front_flag, self, read)
		else:
			# try approximate matching
			alignment = align.globalalign_locate(asequence, read.sequence,
				self.max_error_rate, self.where, self.wildcard_flags)
			match = AdapterMatch(*(alignment + (self._front_flag, self, read)))

		# TODO globalalign_locate should be modified to allow the following
		# assertion.
		# assert length == 0 or match.errors / length <= self.max_error_rate
		if match.length < self.min_overlap or match.errors / match.length > self.max_error_rate:
			return None
		return match

	def _trimmed_front(self, match):
		"""Return a trimmed read"""
		read = match.read
		self.lengths_front[match.rstop] += 1
		# to remove a front adapter, we need to re-encode the first color following the adapter match
		color_after_adapter = read.sequence[match.rstop:match.rstop + 1]
		if not color_after_adapter:
			# the read is empty
			return read[match.rstop:]
		base_after_adapter = colorspace.DECODE[self.nucleotide_sequence[-1] + color_after_adapter]
		new_first_color = colorspace.ENCODE[read.primer + base_after_adapter]
		seq = new_first_color + read.sequence[(match.rstop + 1):]
		qual = read.qualities[match.rstop:] if read.qualities else None
		return ColorspaceSequence(read.name, seq, qual, read.primer)

	def _trimmed_back(self, match):
		"""Return a trimmed read"""
		read = Adapter._trimmed_back(self, match)
		# TODO avoid the copy (previously, this was just an index operation)
		# trim one more color if long enough
		#rstart = max(0, match.rstart - 1)
		read = read[:-1]
		return read

	def __repr__(self):
		return '<ColorspaceAdapter(sequence="{0}", where={1})>'.format(self.sequence, self.where)
