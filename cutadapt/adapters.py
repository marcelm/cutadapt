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


_AdapterMatchBase = namedtuple('AdapterMatch', ['astart', 'astop', 'rstart', 'rstop', 'matches', 'errors', 'adapter', 'read'])

class AdapterMatch(_AdapterMatchBase):
	def wildcards(self, wildcard_char='N'):
		"""
		TODO doc
		"""
		wildcards = [ self.read.sequence[self.rstart + i] for i in range(self.length)
			if self.adapter.sequence[self.astart + i] == wildcard_char ]
		return ''.join(wildcards)

	@property
	def length(self):
		return self.astop - self.astart


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
	"""
	def __init__(self, sequence, where, max_error_rate, min_overlap=3,
			match_read_wildcards=False, match_adapter_wildcards=False,
			rest_file=None):
		self.sequence = sequence.upper()
		self.where = where
		self.max_error_rate = max_error_rate
		self.min_overlap = min_overlap
		self.wildcard_flags = 0
		if match_read_wildcards:
			self.wildcard_flags |= align.ALLOW_WILDCARD_SEQ2
		if match_adapter_wildcards and 'N' in self.sequence:
			self.wildcard_flags |= align.ALLOW_WILDCARD_SEQ1
		removers = {
			FRONT: self.remove_front,
			PREFIX: self.remove_front,
			BACK: self.remove_back,
			ANYWHERE: self.remove_anywhere
		}
		self.remove = removers[where]
		self.rest_file = rest_file
		# statistics about length of removed sequences
		self.lengths_front = defaultdict(int)
		self.lengths_back = defaultdict(int)

	def __repr__(self):
		return '<Adapter(sequence="{0}", where={1})>'.format(self.sequence, self.where)

	def match(self, read):
		"""
		Try to match this adapter to the given read and return an AdapterMatch instance.

		Return None if the minimum overlap length is not met or the error rate is too high.
		"""
		# try to find an exact match first
		# TODO do not do this when wildcards are allowed!
		read_seq = read.sequence.upper()
		if self.where == PREFIX:
			pos = 0 if read_seq.startswith(self.sequence) else -1
		else:
			pos = read_seq.find(self.sequence)
		if pos >= 0:
			match = AdapterMatch(0, len(self.sequence), pos, pos + len(self.sequence), len(self.sequence), 0, self, read)
		else:
			# try approximate matching
			alignment = align.globalalign_locate(self.sequence, read_seq,
				self.max_error_rate, self.where, self.wildcard_flags)
			match = AdapterMatch(*(alignment + (self, read)))
		# TODO globalalign_locate should be modified to allow the following
		# assertion.
		# assert length == 0 or match.errors / length <= self.max_error_rate
		if match.length < self.min_overlap or match.errors / match.length > self.max_error_rate:
			return None
		return match

	def remove_anywhere(self, read, match):
		if match.astart == 0 and match.rstart > 0:
			return self.remove_back(read, match)
		else:
			return self.remove_front(read, match)

	def remove_front(self, read, match):
		self.lengths_front[match.rstop] += 1
		self._write_rest(read.sequence[:match.rstart], read)
		return read[match.rstop:]

	def remove_back(self, read, match):
		# The adapter is at the end of the read or within the read
		self._write_rest(read.sequence[match.rstop:], read)
		self.lengths_back[len(read) - match.rstart] += 1
		return read[:match.rstart]

	def _write_rest(self, rest, read):
		if len(rest) > 0 and self.rest_file:
			# The adapter is within the read
			print(rest, read.name, file=self.rest_file)

	def __len__(self):
		return len(self.sequence)


class ColorspaceAdapter(Adapter):
	def __init__(self, *args):
		super(ColorspaceAdapter, self).__init__(*args)
		has_nucleotide_seq = False
		if set(self.sequence) <= set('ACGT'):
			# adapter was given in basespace
			self.nucleotide_sequence = self.sequence
			has_nucleotide_seq = True
			self.sequence = colorspace.encode(self.sequence)[1:]
		if self.where in (PREFIX, FRONT) and not has_nucleotide_seq:
			raise ValueError("A 5' colorspace adapter needs to be given in nucleotide space")

	def match(self, read):
		if self.where != PREFIX:
			return super(ColorspaceAdapter, self).match(read)
		# create artificial adapter that includes a first color that encodes the
		# transition from primer base into adapter
		asequence = colorspace.ENCODE[read.primer + self.nucleotide_sequence[0]] + self.sequence
		pos = 0 if read.sequence.startswith(asequence) else -1
		if pos >= 0:
			match = AdapterMatch(0, len(asequence), pos, pos + len(asequence), len(asequence), 0, self, read)
		else:
			# try approximate matching
			alignment = align.globalalign_locate(asequence, read.sequence,
				self.max_error_rate, self.where, self.wildcard_flags)
			match = AdapterMatch(*(alignment) + (self, read))

		# TODO globalalign_locate should be modified to allow the following
		# assertion.
		# assert length == 0 or match.errors / length <= self.max_error_rate
		if match.length < self.min_overlap or match.errors / match.length > self.max_error_rate:
			return None
		return match

	def remove_front(self, read, match):
		self.lengths_front[match.rstop] += 1
		self._write_rest(read.sequence[:match.rstart], read)
		# to remove a front adapter, we need to re-encode the first color following the adapter match
		color_after_adapter = read.sequence[match.rstop:match.rstop + 1]
		if not color_after_adapter:
			# the read is empty
			return read[match.rstop:]
		base_after_adapter = colorspace.DECODE[self.nucleotide_sequence[-1] + color_after_adapter]
		new_first_color = colorspace.ENCODE[read.primer + base_after_adapter]
		seq = new_first_color + read.sequence[match.rstop + 1:]
		qual = read.qualities[match.rstop:] if read.qualities else None
		return ColorspaceSequence(read.name, seq, qual, read.primer)

	def remove_back(self, read, match):
		read = Adapter.remove_back(self, read, match)
		# TODO avoid the copy (previously, this was just an index operation)
		# trim one more color if long enough
		#rstart = max(0, match.rstart - 1)
		read = read[:-1]
		return read

	def __repr__(self):
		return '<ColorspaceAdapter(sequence="{0}", where={1})>'.format(self.sequence, self.where)
