"""
Adapter finding and trimming classes

The ...Adapter classes are responsible for finding adapters.
The ...Match classes trim the reads.
"""
import re
from collections import defaultdict
from cutadapt import align
from dnaio.readers import FastaReader


# Constants for the Aligner.locate() function.
# The function is called with SEQ1 as the adapter, SEQ2 as the read.
# TODO get rid of those constants, use strings instead
BACK = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ1
FRONT = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.START_WITHIN_SEQ1
PREFIX = align.STOP_WITHIN_SEQ2
SUFFIX = align.START_WITHIN_SEQ2
# Just like FRONT/BACK, but without internal matches
FRONT_NOT_INTERNAL = align.START_WITHIN_SEQ1 | align.STOP_WITHIN_SEQ2
BACK_NOT_INTERNAL = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ1
ANYWHERE = align.SEMIGLOBAL
LINKED = 'linked'

# TODO put this in some kind of "list of pre-defined adapter types" along with the info above
WHERE_TO_REMOVE_MAP = {
	PREFIX: 'prefix',
	FRONT_NOT_INTERNAL: 'prefix',
	FRONT: 'prefix',
	BACK: 'suffix',
	SUFFIX: 'suffix',
	BACK_NOT_INTERNAL: 'suffix',
	ANYWHERE: 'auto',
}


def parse_braces(sequence):
	"""
	Replace all occurrences of ``x{n}`` (where x is any character) with n
	occurrences of x. Raise ValueError if the expression cannot be parsed.

	>>> parse_braces('TGA{5}CT')
	'TGAAAAACT'
	"""
	# Simple DFA with four states, encoded in prev
	result = ''
	prev = None
	for s in re.split('([{}])', sequence):
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


class AdapterParser:
	"""
	Factory for Adapter classes that all use the same parameters (error rate,
	indels etc.). The given **kwargs will be passed to the Adapter constructors.
	"""
	def __init__(self, **kwargs):
		# kwargs: max_error_rate, min_overlap, read_wildcards, adapter_wildcards, indels
		self.default_parameters = kwargs

	@staticmethod
	def _extract_name(spec):
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

	parameters = {
		# abbreviations
		'e': 'max_error_rate',
		'error_rate': 'max_error_rate',
		'o': 'min_overlap',

		# allowed parameters
		'max_error_rate': None,
		'min_overlap': None,
		'anywhere': None,
	}

	@staticmethod
	def _parse_parameters(spec):
		"""Parse key=value;key=value;key=value into a dict"""

		fields = spec.split(';')
		result = dict()
		for field in fields:
			field = field.strip()
			if not field:
				continue
			key, equals, value = field.partition('=')
			if equals == '=' and value == '':
				raise ValueError('No value given')
			key = key.strip()
			if key not in AdapterParser.parameters:
				raise KeyError('Unknown parameter {}'.format(key))
			# unabbreviate
			while AdapterParser.parameters[key] is not None:
				key = AdapterParser.parameters[key]
			value = value.strip()
			if value == '':
				value = True
			else:
				try:
					value = int(value)
				except ValueError:
					value = float(value)
			if key in result:
				raise KeyError('Key {} specified twice'.format(key))
			result[key] = value
		return result

	@staticmethod
	def _parse_not_linked(spec, cmdline_type):
		"""
		Parse an adapter specification for a non-linked adapter (without '...')

		Allow:
		'back' and ADAPTER
		'back' and ADAPTERX
		'back' and ADAPTER$
		'front' and ADAPTER
		'front' and XADAPTER
		'front' and ^ADAPTER
		'anywhere' and ADAPTER
		"""
		error = ValueError(
				"You cannot use multiple placement restrictions for an adapter at the same time. "
				"Choose one of ^ADAPTER, ADAPTER$, XADAPTER or ADAPTERX")
		spec, middle, parameters_spec = spec.partition(';')
		name, spec = AdapterParser._extract_name(spec)
		spec = spec.strip()

		parameters = AdapterParser._parse_parameters(parameters_spec)

		# Special case for adapters consisting of only X characters:
		# This needs to be supported for backwards-compatibilitity
		if len(spec.strip('X')) == 0:
			return name, None, spec, None, {}

		front_restriction = None
		if spec.startswith('^'):
			front_restriction = 'anchored'
			spec = spec[1:]
		if spec.upper().startswith('X'):
			if front_restriction is not None:
				raise error
			front_restriction = 'noninternal'
			spec = spec.lstrip('xX')

		back_restriction = None
		if spec.endswith('$'):
			back_restriction = 'anchored'
			spec = spec[:-1]
		if spec.upper().endswith('X'):
			if back_restriction is not None:
				raise error
			back_restriction = 'noninternal'
			spec = spec.rstrip('xX')

		n_placement_restrictions = int(bool(front_restriction)) + int(bool(back_restriction))
		if n_placement_restrictions > 1:
			raise error

		if cmdline_type == 'front' and back_restriction:
			raise ValueError(
				"Allowed placement restrictions for a 5' adapter are XADAPTER and ^ADAPTER")
		if cmdline_type == 'back' and front_restriction:
			raise ValueError(
				"Allowed placement restrictions for a 3' adapter are ADAPTERX and ADAPTER$")

		if cmdline_type == 'anywhere' and n_placement_restrictions > 0:
			raise ValueError(
				"Placement restrictions (with X, ^, $) not supported for 'anywhere' (-b) adapters")
		assert front_restriction is None or back_restriction is None
		return name, front_restriction, spec, back_restriction, parameters

	def _parse(self, spec, cmdline_type='back', name=None):
		"""
		Parse an adapter specification not using ``file:`` notation and return
		an object of an appropriate Adapter class.

		name -- Adapter name if not included as part of the spec. (If spec is
		'name=ADAPTER', name will be 'name'.)

		cmdline_type -- describes which commandline parameter was used (``-a``
		is 'back', ``-b`` is 'anywhere', and ``-g`` is 'front').
		"""
		if cmdline_type not in ('front', 'back', 'anywhere'):
			raise ValueError('cmdline_type cannot be {!r}'.format(cmdline_type))
		spec1, middle, spec2 = spec.partition('...')
		del spec

		# Handle linked adapter
		if middle == '...' and spec1 and spec2:
			if cmdline_type == 'anywhere':
				raise ValueError("'anywhere' (-b) adapters may not be linked")
			name1, front1, sequence1, back1, parameters1 = self._parse_not_linked(spec1, 'front')
			assert back1 is None
			name2, front2, sequence2, back2, parameters2 = self._parse_not_linked(spec2, 'back')
			assert front2 is None
			if not name:
				name = name1

			# Automatically anchor the 5' adapter if -a is used
			if cmdline_type == 'back' and front1 is None:
				front1 = 'anchored'

			front_anchored = front1 == 'anchored'
			back_anchored = back2 == 'anchored'
			require_both = True if not front_anchored and not back_anchored else None
			front_parameters = self.default_parameters.copy()
			front_parameters.update(parameters1)
			back_parameters = self.default_parameters.copy()
			back_parameters.update(parameters2)
			return LinkedAdapter(
				sequence1, sequence2, name=name,
				front_restriction=front1,
				back_restriction=back2,
				require_both=require_both,
				front_parameters=front_parameters,
				back_parameters=back_parameters)

		if middle == '...':
			if not spec1:
				if cmdline_type == 'back':  # -a ...ADAPTER
					spec = spec2
				else:  # -g ...ADAPTER
					raise ValueError('Invalid adapter specification')
			elif not spec2:
				if cmdline_type == 'back':  # -a ADAPTER...
					cmdline_type = 'front'
					spec = '^' + spec1
				else:  # -g ADAPTER...
					spec = spec1
			else:
				assert False, 'This should not happen'
		else:
			spec = spec1

		specname, front_restriction, sequence, back_restriction, parameters = self._parse_not_linked(
			spec, cmdline_type)
		del spec
		if front_restriction == 'anchored':
			where = PREFIX
		elif front_restriction == 'noninternal':
			where = FRONT_NOT_INTERNAL
		elif back_restriction == 'anchored':
			where = SUFFIX
		elif back_restriction == 'noninternal':
			where = BACK_NOT_INTERNAL
		elif cmdline_type == 'front':
			where = FRONT
		elif cmdline_type == 'back':
			where = BACK
		else:
			assert cmdline_type == 'anywhere'
			where = ANYWHERE

		if not name:
			name = specname
		if parameters.get('anywhere', False):
			parameters['remove'] = WHERE_TO_REMOVE_MAP[where]
			where = ANYWHERE
			del parameters['anywhere']
		params = self.default_parameters.copy()
		params.update(parameters)
		return Adapter(sequence=sequence, where=where, name=name, **params)

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
					yield self._parse(record.sequence, cmdline_type, name=name)
		else:
			yield self._parse(spec, cmdline_type, name=None)

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


def returns_defaultdict_int():
	# We need this function to make EndStatistics picklable.
	# Even a @staticmethod of EndStatistics is not sufficient
	# as that is not picklable before Python 3.5.
	return defaultdict(int)


class EndStatistics:
	"""Statistics about the 5' or 3' end"""

	def __init__(self, adapter):
		self.where = adapter.where
		self.max_error_rate = adapter.max_error_rate
		self.sequence = adapter.sequence
		self.has_wildcards = adapter.adapter_wildcards
		# self.errors[l][e] == n iff n times a sequence of length l matching at e errors was removed
		self.errors = defaultdict(returns_defaultdict_int)
		self._remove = adapter.remove
		self.adjacent_bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '': 0}

	def __iadd__(self, other):
		if (self.where != other.where or self._remove != other._remove or
				self.max_error_rate != other.max_error_rate or self.sequence != other.sequence):
			raise RuntimeError('Incompatible EndStatistics, cannot be added')
		for base in ('A', 'C', 'G', 'T', ''):
			self.adjacent_bases[base] += other.adjacent_bases[base]
		for length, error_dict in other.errors.items():
			for errors in error_dict:
				self.errors[length][errors] += other.errors[length][errors]
		return self

	@property
	def lengths(self):
		# Python 2.6 has no dict comprehension
		d = dict((length, sum(errors.values())) for length, errors in self.errors.items())
		return d

	def random_match_probabilities(self, gc_content):
		"""
		Estimate probabilities that this adapter end matches a
		random sequence. Indels are not taken into account.

		Returns a list p, where p[i] is the probability that
		i bases of this adapter match a random sequence with
		GC content gc_content.

		The where parameter is necessary for linked adapters to
		specify which (front or back) of the two adapters is meant.
		"""
		seq = self.sequence
		# FIXME this is broken for self._remove == 'auto'
		if self._remove == 'prefix':
			seq = seq[::-1]
		allowed_bases = 'CGRYSKMBDHVN' if self.has_wildcards else 'GC'
		p = 1
		probabilities = [p]
		for i, c in enumerate(seq):
			if c in allowed_bases:
				p *= gc_content / 2.
			else:
				p *= (1 - gc_content) / 2
			probabilities.append(p)
		return probabilities


class AdapterStatistics:
	"""
	Statistics about an adapter. An adapter can work on the 5' end (front)
	or 3' end (back) of a read, and statistics for that are captured
	separately.
	"""

	def __init__(self, adapter, adapter2=None, where=None):
		self.name = adapter.name
		self.where = where if where is not None else adapter.where
		self.front = EndStatistics(adapter)
		if adapter2 is None:
			self.back = EndStatistics(adapter)
		else:
			self.back = EndStatistics(adapter2)

	def __iadd__(self, other):
		if self.where != other.where:  # TODO self.name != other.name or
			raise ValueError('incompatible objects')
		self.front += other.front
		self.back += other.back
		return self


class Match:
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
		return 'Match(astart={}, astop={}, rstart={}, rstop={}, matches={}, errors={})'.format(
			self.astart, self.astop, self.rstart, self.rstop, self.matches, self.errors)

	def wildcards(self, wildcard_char='N'):
		"""
		Return a string that contains, for each wildcard character,
		the character that it matches. For example, if the adapter
		ATNGNA matches ATCGTA, then the string 'CT' is returned.

		If there are indels, this is not reliable as the full alignment
		is not available.
		"""
		wildcards = [self.read.sequence[self.rstart + i] for i in range(self.length)
			if self.adapter.sequence[self.astart + i] == wildcard_char and
				self.rstart + i < len(self.read.sequence)]
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
			statistics.front.errors[self.rstop][self.errors] += 1
		else:
			statistics.back.errors[len(self.read) - len(self._trimmed_read)][self.errors] += 1
			statistics.back.adjacent_bases[self.adjacent_base] += 1


def _generate_adapter_name(_start=[1]):
	name = str(_start[0])
	_start[0] += 1
	return name


class Adapter:
	"""
	This class can find a single adapter characterized by sequence, error rate,
	type etc. within reads.

	where --  One of the BACK, FRONT, PREFIX, SUFFIX or ANYWHERE constants.
		This influences where the adapter is allowed to appear within in the
		read.

	remove -- describes which part of the read to remove if the adapter was found:
		  * "prefix" (for a 3' adapter)
		  * "suffix" (for a 5' adapter)
		  * "auto" for a 5'/3' mixed adapter (if the match involves the first base of the read, it
		    is assumed to be a 5' adapter and a 3' otherwise)
		  * None: One of the above is chosen depending on the 'where' parameter

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

	def __init__(self, sequence, where, remove=None, max_error_rate=0.1, min_overlap=3,
			read_wildcards=False, adapter_wildcards=True, name=None, indels=True):
		self._debug = False
		self.name = _generate_adapter_name() if name is None else name
		self.sequence = parse_braces(sequence.upper().replace('U', 'T'))  # TODO move away
		if not self.sequence:
			raise ValueError('Sequence is empty')
		self.where = where
		if remove not in (None, 'prefix', 'suffix'):
			raise ValueError('remove parameter must be "prefix", "suffix", "auto" or None')
		self.remove = WHERE_TO_REMOVE_MAP[where] if remove is None else remove
		self.max_error_rate = max_error_rate
		self.min_overlap = min(min_overlap, len(self.sequence))
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

		self.aligner = align.Aligner(self.sequence, self.max_error_rate,
			flags=self.where, wildcard_ref=self.adapter_wildcards, wildcard_query=self.read_wildcards)
		self.aligner.min_overlap = self.min_overlap
		self.indels = indels
		if not self.indels:
			# TODO
			# When indels are disallowed, an entirely different algorithm
			# should be used.
			self.aligner.indel_cost = 100000

	def __repr__(self):
		return '<Adapter(name={name!r}, sequence={sequence!r}, where={where}, '\
			'remove={remove}, max_error_rate={max_error_rate}, min_overlap={min_overlap}, '\
			'read_wildcards={read_wildcards}, '\
			'adapter_wildcards={adapter_wildcards}, '\
			'indels={indels})>'.format(**vars(self))

	def enable_debug(self):
		"""
		Print out the dynamic programming matrix after matching a read to an
		adapter.
		"""
		self._debug = True
		self.aligner.enable_debug()

	def match_to(self, read, match_class=Match):
		"""
		Attempt to match this adapter to the given read.

		Return a Match instance if a match was found;
		return None if no match was found given the matching criteria (minimum
		overlap length, maximum error rate).
		"""
		read_seq = read.sequence.upper()  # temporary copy
		pos = -1

		# try to find an exact match first unless wildcards are allowed
		if not self.adapter_wildcards:
			if self.where == PREFIX:
				pos = 0 if read_seq.startswith(self.sequence) else -1
			elif self.where == SUFFIX:
				pos = (len(read_seq) - len(self.sequence)) if read_seq.endswith(self.sequence) else -1
			elif self.where == BACK or self.where == FRONT:
				pos = read_seq.find(self.sequence)
			# TODO BACK_NOT_INTERNAL, FRONT_NOT_INTERNAL
		if pos >= 0:
			match_args = (
				0, len(self.sequence), pos, pos + len(self.sequence),
				len(self.sequence), 0)
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
					match_args = alignment
				else:
					match_args = None
			else:
				alignment = self.aligner.locate(read_seq)
				if self._debug:
					print(self.aligner.dpmatrix)  # pragma: no cover
				if alignment is None:
					match_args = None
				else:
					astart, astop, rstart, rstop, matches, errors = alignment
					match_args = (astart, astop, rstart, rstop, matches, errors)

		if match_args is None:
			return None
		if self.remove == 'auto':
			# guess: if alignment starts at pos 0, it’s a 5' adapter
			remove_before = match_args[2] == 0  # index 2 is rstart
		else:
			remove_before = self.remove == 'prefix'
		match = match_class(*match_args, remove_before=remove_before, adapter=self, read=read)

		assert match.length > 0 and match.errors / match.length <= self.max_error_rate, match
		assert match.length >= self.min_overlap
		return match

	def __len__(self):
		return len(self.sequence)

	def create_statistics(self):
		return AdapterStatistics(self)


class LinkedMatch:
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

	def __repr__(self):
		return '<LinkedMatch(front_match={!r}, back_match={}, adapter={})>'.format(
			self.front_match, self.back_match, self.adapter)

	@property
	def matches(self):
		"""Number of matching bases"""
		m = getattr(self.front_match, 'matches', 0)
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
			statistics.front.errors[self.front_match.rstop][self.front_match.errors] += 1
		if self.back_match:
			statistics.back.errors[len(self.back_match.read) - self.back_match.rstart][self.back_match.errors] += 1


class LinkedAdapter:
	"""
	"""
	def __init__(self, front_sequence, back_sequence, front_restriction='anchored',
			back_restriction=None, require_both=None, name=None,
			front_parameters=dict(), back_parameters=dict()):
		"""
		require_both -- require both adapters to match. If not specified, the default is to
			require only anchored adapters to match.
		kwargs are passed on to individual Adapter constructors
		"""
		if front_restriction is None:
			where1 = FRONT
		elif front_restriction == 'anchored':
			where1 = PREFIX
		elif front_restriction == 'noninternal':
			where1 = FRONT_NOT_INTERNAL
		else:
			raise ValueError('Value {} for front_restriction not allowed'.format(front_restriction))

		if back_restriction is None:
			where2 = BACK
		elif back_restriction == 'anchored':
			where2 = SUFFIX
		elif back_restriction == 'noninternal':
			where2 = BACK_NOT_INTERNAL
		else:
			raise ValueError(
				'Value {} for back_restriction not allowed'.format(back_restriction))
		if require_both:
			self._require_back_match = True
			self._require_front_match = True
		else:
			self._require_front_match = front_restriction == 'anchored'
			self._require_back_match = back_restriction == 'anchored'

		# The following attributes are needed for the report
		self.where = LINKED
		self.name = _generate_adapter_name() if name is None else name
		self.front_adapter = Adapter(front_sequence, where=where1, name=None, **front_parameters)
		self.back_adapter = Adapter(back_sequence, where=where2, name=None, **back_parameters)

	def enable_debug(self):
		self.front_adapter.enable_debug()
		self.back_adapter.enable_debug()

	def match_to(self, read):
		"""
		Match the linked adapters against the given read. Any anchored adapters are
		required to exist for a successful match. If both adapters are unanchored,
		both need to match.
		"""
		front_match = self.front_adapter.match_to(read)
		if self._require_front_match and front_match is None:
			return None

		if front_match is not None:
			# TODO statistics
			read = front_match.trimmed()
		back_match = self.back_adapter.match_to(read)
		if back_match is None and (self._require_back_match or front_match is None):
			return None
		return LinkedMatch(front_match, back_match, self)

	def create_statistics(self):
		return AdapterStatistics(self.front_adapter, self.back_adapter, where=LINKED)
