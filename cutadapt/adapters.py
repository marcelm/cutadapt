from __future__ import print_function, division, absolute_import
import sys
from collections import defaultdict
from cutadapt import align, colorspace
from cutadapt.seqio import ColorspaceSequence, FastaReader

# Constants for the find_best_alignment function.
# The function is called with SEQ1 as the adapter, SEQ2 as the read.
BACK = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ1
FRONT = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.START_WITHIN_SEQ1
PREFIX = align.STOP_WITHIN_SEQ2
SUFFIX = align.START_WITHIN_SEQ2
ANYWHERE = align.SEMIGLOBAL


def parse_adapter_name(seq):
	"""
	Parse an adapter given as 'name=adapt' into 'name' and 'adapt'.
	"""
	fields = seq.split('=', 1)
	if len(fields) > 1:
		name, seq = fields
		name = name.strip()
	else:
		name = None
	seq = seq.strip()
	return name, seq


def parse_adapter(sequence, where):
	"""
	Recognize anchored adapter sequences and return a corrected tuple
	(sequence, where).
	"""
	if where == FRONT and sequence.startswith('^'):
		return (sequence[1:],  PREFIX)
	if where == BACK and sequence.endswith('$'):
		return (sequence[:-1], SUFFIX)
	return (sequence, where)


def gather_adapters(back, anywhere, front):
	"""
	Yield (name, seq, where) tuples from which Adapter instances can be built.
	This generator deals with the notation for anchored 5' adapters and also
	understands the ``file:`` syntax for reading adapters from an external FASTA
	file.
	"""
	for adapter_list, where in ((back, BACK), (anywhere, ANYWHERE), (front, FRONT)):
		for seq in adapter_list:
			if seq.startswith('file:'):
				# read adapter sequences from a file
				path = seq[5:]
				with FastaReader(path) as fasta:
					for record in fasta:
						name = record.name.split(None, 1)[0]
						seq, w = parse_adapter(record.sequence, where)
						yield (name, seq, w)
			else:
				name, seq = parse_adapter_name(seq)
				seq, w = parse_adapter(seq, where)
				yield (name, seq, w)


class AdapterMatch(object):
	"""
	TODO creating instances of this class is relatively slow and responsible for quite some runtime.
	"""
	__slots__ = ['astart', 'astop', 'rstart', 'rstop', 'matches', 'errors', 'front', 'adapter', 'read', 'length']
	def __init__(self, astart, astop, rstart, rstop, matches, errors, front, adapter, read):
		self.astart = astart
		self.astop = astop
		self.rstart = rstart
		self.rstop = rstop
		self.matches = matches
		self.errors = errors
		self.front = self._guess_is_front() if front is None else front
		self.adapter = adapter
		self.read = read
		# Number of aligned characters in the adapter. If there are
		# indels, this may be different from the number of characters
		# in the read.
		self.length = self.astop - self.astart

	def __str__(self):
		return 'AdapterMatch(astart={}, astop={}, rstart={}, rstop={}, matches={}, errors={})'.format(
			self.astart, self.astop, self.rstart, self.rstop, self.matches, self.errors)

	def _guess_is_front(self):
		"""
		Return whether this is guessed to be a front adapter.

		The match is assumed to be a front adapter when the first base of
		the read is involved in the alignment to the adapter.
		"""
		return self.rstart == 0

	def wildcards(self, wildcard_char='N'):
		"""
		Return a string that contains, for each wildcard character,
		the character that it matches. For example, if the adapter
		ATNGNA matches ATCGTA, then the string 'CT' is returned.

		If there are indels, this is not reliable as the full alignment
		is not available.
		"""
		wildcards = [ self.read.sequence[self.rstart + i:self.rstart + i + 1] for i in range(self.length)
			if self.adapter.sequence[self.astart + i] == wildcard_char and self.rstart + i < len(self.read.sequence) ]
		return ''.join(wildcards)

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

	where --  One of the BACK, FRONT, PREFIX, SUFFIX or ANYWHERE constants.
		This influences where the adapter is allowed to appear within in the
		read and also which part of the read is removed.

	sequence -- The adapter sequence as string. Will be converted to uppercase.
		Also, Us will be converted to Ts.

	max_error_rate -- Maximum allowed error rate. The error rate is
		the number of errors in the alignment divided by the length
		of the part of the alignment that matches the adapter.

	minimum_overlap -- Minimum length of the part of the alignment
		that matches the adapter.

	read_wildcards -- Whether IUPAC wildcards in the read are allowed.

	adapter_wildcards -- Whether IUPAC wildcards in the adapter are
		allowed.

	name -- optional name of the adapter. If not provided, the name is set to a
		unique number.
	"""
	automatic_name = 1

	def __init__(self, sequence, where, max_error_rate, min_overlap=3,
			read_wildcards=False, adapter_wildcards=True,
			name=None, indels=True):
		if name is None:
			self.name = str(self.__class__.automatic_name)
			self.__class__.automatic_name += 1
			self.name_is_generated = True
		else:
			self.name = name
			self.name_is_generated = False

		self.sequence = sequence.upper().replace('U', 'T')
		self.where = where
		self.max_error_rate = max_error_rate
		self.min_overlap = min_overlap
		self.indels = indels
		assert where in (PREFIX, SUFFIX) or self.indels
		self.wildcard_flags = 0
		self.adapter_wildcards = adapter_wildcards and not set(self.sequence) <= set('ACGT')
		if read_wildcards:
			self.wildcard_flags |= align.ALLOW_WILDCARD_SEQ2
		if self.adapter_wildcards:
			self.wildcard_flags |= align.ALLOW_WILDCARD_SEQ1
		# redirect to appropriate trimmed() function depending on
		# adapter type
		trimmers = {
			FRONT: self._trimmed_front,
			PREFIX: self._trimmed_front,
			BACK: self._trimmed_back,
			SUFFIX: self._trimmed_back,
			ANYWHERE: self._trimmed_anywhere
		}
		self.trimmed = trimmers[where]
		if where == ANYWHERE:
			self._front_flag = None  # means: guess
		else:
			self._front_flag = where not in (BACK, SUFFIX)
		# statistics about length of removed sequences
		self.lengths_front = defaultdict(int)
		self.lengths_back = defaultdict(int)
		self.errors_front = defaultdict(lambda: defaultdict(int))
		self.errors_back = defaultdict(lambda: defaultdict(int))
		self.adjacent_bases = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, '': 0 }

		self.aligner = align.Aligner(self.sequence, self.max_error_rate,
			flags=self.where, degenerate=self.wildcard_flags,
			min_overlap=self.min_overlap)

	def __repr__(self):
		read_wildcards = bool(align.ALLOW_WILDCARD_SEQ2 & self.wildcard_flags)
		return '<Adapter(name="{name}", sequence="{sequence}", where={where}, '\
			'max_error_rate={max_error_rate}, min_overlap={min_overlap}, '\
			'read_wildcards={read_wildcards}, '\
			'adapter_wildcards={adapter_wildcards}, '\
			'indels={indels})>'.format(
				read_wildcards=read_wildcards,
				**vars(self))

	def match_to(self, read):
		"""
		Try to match this adapter to the given read and return an AdapterMatch instance.

		Return None if the minimum overlap length is not met or the error rate is too high.
		"""
		read_seq = read.sequence.upper()
		pos = -1
		# try to find an exact match first unless wildcards are allowed
		if not self.adapter_wildcards:
			if self.where == PREFIX:
				pos = 0 if read_seq.startswith(self.sequence) else -1
			elif self.where == SUFFIX:
				pos = (len(read_seq) - len(self.sequence)) if read_seq.endswith(self.sequence) else -1
			else:
				pos = read_seq.find(self.sequence)
		if pos >= 0:
			match = AdapterMatch(
				0, len(self.sequence), pos, pos + len(self.sequence),
				len(self.sequence), 0, self._front_flag, self, read)
		else:
			# try approximate matching
			if not self.indels:
				assert self.where in (PREFIX, SUFFIX)
				if self.where == PREFIX:
					alignment = align.compare_prefixes(self.sequence, read_seq, self.wildcard_flags)
				else:
					alignment = align.compare_suffixes(self.sequence, read_seq, self.wildcard_flags)
				astart, astop, rstart, rstop, matches, errors = alignment
				if astop - astart >= self.min_overlap and errors / (astop - astart) <= self.max_error_rate:
					match = AdapterMatch(*(alignment + (self._front_flag, self, read)))
				else:
					match = None
			else:
				alignment = self.aligner.locate(read_seq)
				if alignment is None:
					match = None
				else:
					astart, astop, rstart, rstop, matches, errors = alignment
					match = AdapterMatch(astart, astop, rstart, rstop, matches, errors, self._front_flag, self, read)

		if match is None:
			return None
		assert match.length > 0 and match.errors / match.length <= self.max_error_rate, match
		assert match.length >= self.min_overlap
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
		self.errors_front[match.rstop][match.errors] += 1
		return match.read[match.rstop:]

	def _trimmed_back(self, match):
		"""Return a trimmed read without the 3' (back) adapter"""
		# TODO move away
		self.lengths_back[len(match.read) - match.rstart] += 1
		self.errors_back[len(match.read) - match.rstart][match.errors] += 1
		adjacent_base = match.read.sequence[match.rstart-1:match.rstart]
		if adjacent_base not in 'ACGT':
			adjacent_base = ''
		self.adjacent_bases[adjacent_base] += 1
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
		self.aligner.reference = self.sequence

	def match_to(self, read):
		"""Return AdapterMatch instance"""
		if self.where != PREFIX:
			return super(ColorspaceAdapter, self).match_to(read)
		# create artificial adapter that includes a first color that encodes the
		# transition from primer base into adapter
		asequence = colorspace.ENCODE[read.primer + self.nucleotide_sequence[0:1]] + self.sequence

		pos = 0 if read.sequence.startswith(asequence) else -1
		if pos >= 0:
			match = AdapterMatch(
				0, len(asequence), pos, pos + len(asequence),
				len(asequence), 0, self._front_flag, self, read)
		else:
			# try approximate matching
			self.aligner.reference = asequence
			alignment = self.aligner.locate(read.sequence)
			if alignment is not None:
				match = AdapterMatch(*(alignment + (self._front_flag, self, read)))
			else:
				match = None

		if match is None:
			return None
		assert match.length > 0 and match.errors / match.length <= self.max_error_rate
		assert match.length >= self.min_overlap
		return match

	def _trimmed_front(self, match):
		"""Return a trimmed read"""
		read = match.read
		self.lengths_front[match.rstop] += 1
		self.errors_front[match.rstop][match.errors] += 1
		# to remove a front adapter, we need to re-encode the first color following the adapter match
		color_after_adapter = read.sequence[match.rstop:match.rstop + 1]
		if not color_after_adapter:
			# the read is empty
			return read[match.rstop:]
		base_after_adapter = colorspace.DECODE[self.nucleotide_sequence[-1:] + color_after_adapter]
		new_first_color = colorspace.ENCODE[read.primer + base_after_adapter]
		read.sequence = new_first_color + read.sequence[(match.rstop + 1):]
		read.qualities = read.qualities[match.rstop:] if read.qualities else None
		return read

	def _trimmed_back(self, match):
		"""Return a trimmed read"""
		# trim one more color if long enough
		adjusted_rstart = max(match.rstart - 1, 0)
		self.lengths_back[len(match.read) - adjusted_rstart] += 1
		self.errors_back[len(match.read) - adjusted_rstart][match.errors] += 1
		return match.read[:adjusted_rstart]

	def __repr__(self):
		return '<ColorspaceAdapter(sequence={0!r}, where={1})>'.format(self.sequence, self.where)
