# coding: utf-8
"""
Adapter finding and trimming classes

The ...Adapter classes are responsible for finding adapters.
The ...Match classes trim the reads.
"""
from __future__ import print_function, division, absolute_import
import sys
import re
from collections import defaultdict
from cutadapt import align, colorspace
from cutadapt.seqio import ColorspaceSequence, FastaReader

# Constants for the find_best_alignment function.
# The function is called with SEQ1 as the adapter, SEQ2 as the read.
# TODO get rid of those constants, use strings instead
BACK = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ1
FRONT = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.START_WITHIN_SEQ1
PREFIX = align.STOP_WITHIN_SEQ2
SUFFIX = align.START_WITHIN_SEQ2
ANYWHERE = align.SEMIGLOBAL
LINKED = 'linked'


def parse_braces(sequence):
	"""
	Replace all occurrences of ``x{n}`` (where x is any character) with n
	occurrences of x. Raise ValueError if the expression cannot be parsed.

	>>> parse_braces('TGA{5}CT')
	TGAAAAACT
	"""
	# Simple DFA with four states, encoded in prev
	result = ''
	prev = None
	for s in re.split('(\{|\})', sequence):
		if s == '':
			continue
		if prev is None:
			if s == '{':
				raise ValueError('"{" must be used after a character')
			if s == '}':
				raise ValueError('"}" cannot be used here')
			prev = s
			result += s
		elif prev == '{':
			prev = int(s)
			if not 0 <= prev <= 10000:
				raise ValueError('Value {} invalid'.format(prev))
		elif isinstance(prev, int):
			if s != '}':
				raise ValueError('"}" expected')
			result = result[:-1] + result[-1] * prev
			prev = None
		else:
			if s != '{':
				raise ValueError('Expected "{"')
			prev = '{'
	# Check if we are in a non-terminating state
	if isinstance(prev, int) or prev == '{':
		raise ValueError("Unterminated expression")
	return result


class AdapterParser(object):
	"""
	Factory for Adapter classes that all use the same parameters (error rate,
	indels etc.). The given **kwargs will be passed to the Adapter constructors.
	"""
	def __init__(self, colorspace=False, **kwargs):
		self.colorspace = colorspace
		self.constructor_args = kwargs
		self.adapter_class = ColorspaceAdapter if colorspace else Adapter

	def _parse_no_file(self, spec, name=None, cmdline_type='back'):
		"""
		Parse an adapter specification not using ``file:`` notation and return
		an object of an appropriate Adapter class. The notation for anchored
		5' and 3' adapters is supported. If the name parameter is None, then
		an attempt is made to extract the name from the specification
		(If spec is 'name=ADAPTER', name will be 'name'.)

		cmdline_type -- describes which commandline parameter was used (``-a``
		is 'back', ``-b`` is 'anywhere', and ``-g`` is 'front').
		"""
		if name is None:
			name, spec = self._extract_name(spec)
		orig_spec = spec
		types = dict(back=BACK, front=FRONT, anywhere=ANYWHERE)
		if cmdline_type not in types:
			raise ValueError('cmdline_type cannot be {0!r}'.format(cmdline_type))
		where = types[cmdline_type]

		front_anchored, back_anchored = False, False
		if spec.startswith('^'):
			spec = spec[1:]
			front_anchored = True
		if spec.endswith('$'):
			spec = spec[:-1]
			back_anchored = True

		sequence1, middle, sequence2 = spec.partition('...')
		if where == ANYWHERE:
			if front_anchored or back_anchored:
				raise ValueError("'anywhere' (-b) adapters may not be anchored")
			if middle == '...':
				raise ValueError("'anywhere' (-b) adapters may not be linked")
			return self.adapter_class(sequence=spec, where=where, name=name, **self.constructor_args)

		assert where == FRONT or where == BACK
		if middle == '...':
			if not sequence1:
				if where == BACK:  # -a ...ADAPTER
					spec = sequence2
				else:  # -g ...ADAPTER
					raise ValueError('Invalid adapter specification')
			elif not sequence2:
				if where == BACK:  # -a ADAPTER...
					spec = sequence1
					where = FRONT
					front_anchored = True
				else:  # -g ADAPTER...
					spec = sequence1
			else:
				# linked adapter
				if self.colorspace:
					raise NotImplementedError(
						'Using linked adapters in colorspace is not supported')
				# automatically anchor 5' adapter if -a is used
				if where == BACK:
					front_anchored = True

				return LinkedAdapter(sequence1, sequence2, name=name,
					front_anchored=front_anchored, back_anchored=back_anchored,
					**self.constructor_args)
		if front_anchored and back_anchored:
			raise ValueError('Trying to use both "^" and "$" in adapter specification {!r}'.format(orig_spec))
		if front_anchored:
			if where == BACK:
				raise ValueError("Cannot anchor the 3' adapter at its 5' end")
			where = PREFIX
		elif back_anchored:
			if where == FRONT:
				raise ValueError("Cannot anchor 5' adapter at 3' end")
			where = SUFFIX

		return self.adapter_class(sequence=spec, where=where, name=name, **self.constructor_args)

	def parse(self, spec, cmdline_type='back'):
		"""
		Parse an adapter specification and yield appropriate Adapter classes.
		This works like the _parse_no_file() function above, but also supports the
		``file:`` notation for reading adapters from an external FASTA
		file. Since a file can contain multiple adapters, this
		function is a generator.
		"""
		if spec.startswith('file:'):
			# read adapter sequences from a file
			with FastaReader(spec[5:]) as fasta:
				for record in fasta:
					name = record.name.split(None, 1)[0]
					yield self._parse_no_file(record.sequence, name, cmdline_type)
		else:
			name, spec = self._extract_name(spec)
			yield self._parse_no_file(spec, name, cmdline_type)

	def _extract_name(self, spec):
		"""
		Parse an adapter specification given as 'name=adapt' into 'name' and 'adapt'.
		"""
		fields = spec.split('=', 1)
		if len(fields) > 1:
			name, spec = fields
			name = name.strip()
		else:
			name = None
		spec = spec.strip()
		return name, spec

	def parse_multi(self, back, anywhere, front):
		"""
		Parse all three types of commandline options that can be used to
		specify adapters. back, anywhere and front are lists of strings,
		corresponding to the respective commandline types (-a, -b, -g).

		Return a list of appropriate Adapter classes.
		"""
		adapters = []
		for specs, cmdline_type in (back, 'back'), (anywhere, 'anywhere'), (front, 'front'):
			for spec in specs:
				adapters.extend(self.parse(spec, cmdline_type))
		return adapters


class Match(object):
	"""
	Representation of a single adapter matched to a single read.

	TODO creating instances of this class is relatively slow and responsible for quite some runtime.
	"""
	__slots__ = ['astart', 'astop', 'rstart', 'rstop', 'matches', 'errors', 'remove_before',
		'adapter', 'read', 'length', '_trimmed_read', 'adjacent_base']

	def __init__(self, astart, astop, rstart, rstop, matches, errors, remove_before, adapter, read):
		"""
		remove_before -- True: remove bases before adapter. False: remove after
		"""
		self.astart = astart
		self.astop = astop
		self.rstart = rstart
		self.rstop = rstop
		self.matches = matches
		self.errors = errors
		self.adapter = adapter
		self.read = read
		if remove_before:
			self._trim_front()
		else:
			self._trim_back()
		self.remove_before = remove_before
		# Number of aligned characters in the adapter. If there are
		# indels, this may be different from the number of characters
		# in the read.
		self.length = self.astop - self.astart
		assert self.length > 0
		assert self.errors / self.length <= self.adapter.max_error_rate
		assert self.length - self.errors > 0

	def __repr__(self):
		return 'Match(astart={0}, astop={1}, rstart={2}, rstop={3}, matches={4}, errors={5})'.format(
			self.astart, self.astop, self.rstart, self.rstop, self.matches, self.errors)

	def wildcards(self, wildcard_char='N'):
		"""
		Return a string that contains, for each wildcard character,
		the character that it matches. For example, if the adapter
		ATNGNA matches ATCGTA, then the string 'CT' is returned.

		If there are indels, this is not reliable as the full alignment
		is not available.
		"""
		wildcards = [ self.read.sequence[self.rstart + i] for i in range(self.length)
			if self.adapter.sequence[self.astart + i] == wildcard_char and
				self.rstart + i < len(self.read.sequence) ]
		return ''.join(wildcards)

	def rest(self):
		"""
		Return the part of the read before this match if this is a
		'front' (5') adapter,
		return the part after the match if this is not a 'front' adapter (3').
		This can be an empty string.
		"""
		if self.remove_before:
			return self.read.sequence[:self.rstart]
		else:
			return self.read.sequence[self.rstop:]

	def get_info_record(self):
		seq = self.read.sequence
		qualities = self.read.qualities
		info = (
			self.read.name,
			self.errors,
			self.rstart,
			self.rstop,
			seq[0:self.rstart],
			seq[self.rstart:self.rstop],
			seq[self.rstop:],
			self.adapter.name
		)
		if qualities:
			info += (
				qualities[0:self.rstart],
				qualities[self.rstart:self.rstop],
				qualities[self.rstop:]
			)
		else:
			info += ('', '', '')

		return info

	def trimmed(self):
		return self._trimmed_read

	def _trim_front(self):
		"""Compute the trimmed read, assuming it’s a 'front' adapter"""
		self._trimmed_read = self.read[self.rstop:]
		self.adjacent_base = ''

	def _trim_back(self):
		"""Compute the trimmed read, assuming it’s a 'back' adapter"""
		adjacent_base = self.read.sequence[self.rstart-1:self.rstart]
		if adjacent_base not in 'ACGT':
			adjacent_base = ''
		self.adjacent_base = adjacent_base
		self._trimmed_read = self.read[:self.rstart]

	def update_statistics(self, statistics):
		"""Update AdapterStatistics in place"""
		if self.remove_before:
			statistics.errors_front[self.rstop][self.errors] += 1
		else:
			statistics.errors_back[len(self.read) - len(self._trimmed_read)][self.errors] += 1
			statistics.adjacent_bases[self.adjacent_base] += 1


class ColorspaceMatch(Match):
	adjacent_base = ''

	def _trim_front(self):
		"""Return a trimmed read"""
		read = self.read
		# to remove a front adapter, we need to re-encode the first color following the adapter match
		color_after_adapter = read.sequence[self.rstop:self.rstop + 1]
		if not color_after_adapter:
			# the read is empty
			new_read = read[self.rstop:]
		else:
			base_after_adapter = colorspace.DECODE[self.adapter.nucleotide_sequence[-1:] + color_after_adapter]
			new_first_color = colorspace.ENCODE[read.primer + base_after_adapter]
			new_read = read[:]
			new_read.sequence = new_first_color + read.sequence[(self.rstop + 1):]
			new_read.qualities = read.qualities[self.rstop:] if read.qualities else None
		self._trimmed_read = new_read

	def _trim_back(self):
		"""Return a trimmed read"""
		# trim one more color if long enough
		adjusted_rstart = max(self.rstart - 1, 0)
		self._trimmed_read = self.read[:adjusted_rstart]

	def update_statistics(self, statistics):
		"""Update AdapterStatistics in place"""
		if self.remove_before:
			statistics.errors_front[self.rstop][self.errors] += 1
		else:
			statistics.errors_back[len(self.read) - len(self._trimmed_read)][self.errors] += 1


def _generate_adapter_name(_start=[1]):
	name = str(_start[0])
	_start[0] += 1
	return name


class Adapter(object):
	"""
	This class can find a single adapter characterized by sequence, error rate,
	type etc. within reads.

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

	def __init__(self, sequence, where, max_error_rate=0.1, min_overlap=3,
			read_wildcards=False, adapter_wildcards=True, name=None, indels=True):
		self.debug = False
		self.name = _generate_adapter_name() if name is None else name
		self.sequence = parse_braces(sequence.upper().replace('U', 'T'))  # TODO move away
		if not self.sequence:
			raise ValueError('Sequence is empty')
		self.where = where
		self.max_error_rate = max_error_rate
		self.min_overlap = min(min_overlap, len(self.sequence))
		self.indels = indels
		iupac = frozenset('XACGTURYSWKMBDHVN')
		if adapter_wildcards and not set(self.sequence) <= iupac:
			for c in self.sequence:
				if c not in iupac:
					raise ValueError('Character {!r} in adapter sequence {!r} is '
						'not a valid IUPAC code. Use only characters '
						'XACGTURYSWKMBDHVN.'.format(c, self.sequence))
		# Optimization: Use non-wildcard matching if only ACGT is used
		self.adapter_wildcards = adapter_wildcards and not set(self.sequence) <= set('ACGT')
		self.read_wildcards = read_wildcards
		self.remove_before = where not in (BACK, SUFFIX)

		self.aligner = align.Aligner(self.sequence, self.max_error_rate,
			flags=self.where, wildcard_ref=self.adapter_wildcards, wildcard_query=self.read_wildcards)
		self.aligner.min_overlap = self.min_overlap
		if not self.indels:
			# TODO
			# When indels are disallowed, an entirely different algorithm
			# should be used.
			self.aligner.indel_cost = 100000

	def __repr__(self):
		return '<Adapter(name={name!r}, sequence={sequence!r}, where={where}, '\
			'max_error_rate={max_error_rate}, min_overlap={min_overlap}, '\
			'read_wildcards={read_wildcards}, '\
			'adapter_wildcards={adapter_wildcards}, '\
			'indels={indels})>'.format(**vars(self))

	def enable_debug(self):
		"""
		Print out the dynamic programming matrix after matching a read to an
		adapter.
		"""
		self.debug = True
		self.aligner.enable_debug()

	def match_to(self, read, match_class=Match):
		"""
		Attempt to match this adapter to the given read.

		Return a Match instance if a match was found;
		return None if no match was found given the matching criteria (minimum
		overlap length, maximum error rate).
		"""
		read_seq = read.sequence.upper()  # temporary copy
		remove_before = self.remove_before
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
			if self.where == ANYWHERE:
				# guess: if alignment starts at pos 0, it’s a 5' adapter
				remove_before = pos == 0
			match = match_class(
				0, len(self.sequence), pos, pos + len(self.sequence),
				len(self.sequence), 0, remove_before, self, read)
		else:
			# try approximate matching
			if not self.indels and self.where in (PREFIX, SUFFIX):
				if self.where == PREFIX:
					alignment = align.compare_prefixes(self.sequence, read_seq,
						wildcard_ref=self.adapter_wildcards, wildcard_query=self.read_wildcards)
				else:
					alignment = align.compare_suffixes(self.sequence, read_seq,
						wildcard_ref=self.adapter_wildcards, wildcard_query=self.read_wildcards)
				astart, astop, rstart, rstop, matches, errors = alignment
				if astop - astart >= self.min_overlap and errors / (astop - astart) <= self.max_error_rate:
					match = match_class(*(alignment + (remove_before, self, read)))
				else:
					match = None
			else:
				alignment = self.aligner.locate(read_seq)
				if self.debug:
					print(self.aligner.dpmatrix)  # pragma: no cover
				if alignment is None:
					match = None
				else:
					astart, astop, rstart, rstop, matches, errors = alignment
					if self.where == ANYWHERE:
						# guess: if alignment starts at pos 0, it’s a 5' adapter
						remove_before = rstart == 0
					match = match_class(astart, astop, rstart, rstop, matches, errors, remove_before, self, read)

		if match is None:
			return None
		assert match.length > 0 and match.errors / match.length <= self.max_error_rate, match
		assert match.length >= self.min_overlap
		return match

	def __len__(self):
		return len(self.sequence)

	def random_match_probabilities(self, gc_content):
		"""
		Estimate probabilities that this adapter matches a
		random sequence. Indels are not taken into account.

		Returns a list p, where p[i] is the probability that
		i bases of this adapter match a random sequence with
		GC content gc_content.
		"""
		if self.remove_before:
			seq = self.sequence[::-1]
		else:
			seq = self.sequence
		allowed_bases = 'CGRYSKMBDHVN' if self.adapter_wildcards else 'GC'
		p = 1
		probabilities = [p]
		for i, c in enumerate(seq):
			if c in allowed_bases:
				p *= gc_content / 2.
			else:
				p *= (1 - gc_content) / 2
			probabilities.append(p)
		return probabilities


class ColorspaceAdapter(Adapter):
	"""
	An Adapter, but in color space. It does not support all adapter types
	(see the 'where' parameter).
	"""

	def __init__(self, *args, **kwargs):
		"""
		sequence -- the adapter sequence as a str, can be given in nucleotide space or in color space
		where -- PREFIX, FRONT, BACK
		"""
		if kwargs.get('adapter_wildcards', False):
			raise ValueError('Wildcards not supported for colorspace adapters')
		kwargs['adapter_wildcards'] = False
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

	def __repr__(self):
		return '<ColorspaceAdapter(sequence={0!r}, where={1})>'.format(self.sequence, self.where)

	def match_to(self, read, match_class=ColorspaceMatch):
		"""
		Match the adapter to the given read

		Return a ColorspaceMatch instance or None if the adapter was not found
		"""
		if self.where != PREFIX:
			return super(ColorspaceAdapter, self).match_to(read, match_class=match_class)
		# create artificial adapter that includes a first color that encodes the
		# transition from primer base into adapter
		asequence = colorspace.ENCODE[read.primer + self.nucleotide_sequence[0:1]] + self.sequence

		pos = 0 if read.sequence.startswith(asequence) else -1
		if pos >= 0:
			match = ColorspaceMatch(
				0, len(asequence), pos, pos + len(asequence),
				len(asequence), 0, self.remove_before, self, read)
		else:
			# try approximate matching
			self.aligner.reference = asequence
			alignment = self.aligner.locate(read.sequence)
			if self.debug:
				print(self.aligner.dpmatrix)  # pragma: no cover
			if alignment is not None:
				match = ColorspaceMatch(*(alignment + (self.remove_before, self, read)))
			else:
				match = None

		if match is None:
			return None
		assert match.length > 0 and match.errors / match.length <= self.max_error_rate
		assert match.length >= self.min_overlap
		return match


class LinkedMatch(object):
	"""
	Represent a match of a LinkedAdapter
	"""
	def __init__(self, front_match, back_match, adapter):
		"""
		One of front_match and back_match must be not None!
		"""
		self.front_match = front_match
		self.back_match = back_match
		self.adapter = adapter
		assert not adapter.front_anchored or front_match is not None
		assert not adapter.back_anchored or back_match is not None

	def __repr__(self):
		return '<LinkedMatch(front_match={0!r}, back_match={1}, adapter={2})>'.format(
			self.front_match, self.back_match, self.adapter)

	@property
	def matches(self):
		"""Number of matching bases"""
		m = getattr(self.front_match, 'matches', [])
		if self.back_match is not None:
			m += self.back_match.matches
		return m

	def trimmed(self):
		if self.back_match:
			# back match is relative to front match, so even if a front match exists,
			# this is correct
			return self.back_match.trimmed()
		else:
			assert self.front_match
			return self.front_match.trimmed()

	@property
	def adjacent_base(self):
		return self.back_match.adjacent_base

	def update_statistics(self, statistics):
		"""Update AdapterStatistics in place"""
		if self.front_match:
			statistics.errors_front[self.front_match.rstop][self.front_match.errors] += 1
		if self.back_match:
			statistics.errors_back[len(self.back_match.read) - self.back_match.rstart][self.back_match.errors] += 1


class LinkedAdapter(object):
	"""
	"""
	def __init__(self, front_sequence, back_sequence,
			front_anchored=True, back_anchored=False, name=None, **kwargs):
		"""
		kwargs are passed on to individual Adapter constructors
		"""
		where1 = PREFIX if front_anchored else FRONT
		where2 = SUFFIX if back_anchored else BACK
		self.front_anchored = front_anchored
		self.back_anchored = back_anchored

		# The following attributes are needed for the report
		self.where = LINKED
		self.name = _generate_adapter_name() if name is None else name
		self.front_adapter = Adapter(front_sequence, where=where1, name=None, **kwargs)
		self.back_adapter = Adapter(back_sequence, where=where2, name=None, **kwargs)

	def enable_debug(self):
		self.front_adapter.enable_debug()
		self.back_adapter.enable_debug()

	def match_to(self, read):
		"""
		Match the linked adapters against the given read. Any anchored adapters are
		required to exist for a successful match.
		"""
		front_match = self.front_adapter.match_to(read)
		if self.front_anchored and front_match is None:
			return None

		if front_match is not None:
			# TODO statistics
			read = front_match.trimmed()
		back_match = self.back_adapter.match_to(read)
		if back_match is None and (self.back_anchored or front_match is None):
			return None
		return LinkedMatch(front_match, back_match, self)
