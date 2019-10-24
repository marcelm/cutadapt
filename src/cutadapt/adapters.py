"""
Adapter finding and trimming classes

The ...Adapter classes are responsible for finding adapters.
The ...Match classes trim the reads.
"""
import logging
from enum import Enum
from collections import defaultdict

from . import align

logger = logging.getLogger()


class Where(Enum):
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
    Where.PREFIX: 'prefix',
    Where.FRONT_NOT_INTERNAL: 'prefix',
    Where.FRONT: 'prefix',
    Where.BACK: 'suffix',
    Where.SUFFIX: 'suffix',
    Where.BACK_NOT_INTERNAL: 'suffix',
    Where.ANYWHERE: 'auto',
}


ADAPTER_TYPE_NAMES = {
    Where.BACK: "regular 3'",
    Where.BACK_NOT_INTERNAL: "non-internal 3'",
    Where.FRONT: "regular 5'",
    Where.FRONT_NOT_INTERNAL: "non-internal 5'",
    Where.PREFIX: "anchored 5'",
    Where.SUFFIX: "anchored 3'",
    Where.ANYWHERE: "variable 5'/3'",
    Where.LINKED: "linked",
}


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
        self.effective_length = adapter.effective_length
        self.has_wildcards = adapter.adapter_wildcards
        # self.errors[l][e] == n iff n times a sequence of length l matching at e errors was removed
        self.errors = defaultdict(returns_defaultdict_int)
        self._remove = adapter.remove
        self.adjacent_bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '': 0}

    def __iadd__(self, other):
        if not isinstance(other, self.__class__):
            raise ValueError("Cannot compare")
        if (
            self.where != other.where
            or self._remove != other._remove
            or self.max_error_rate != other.max_error_rate
            or self.sequence != other.sequence
            or self.effective_length != other.effective_length
        ):
            raise RuntimeError('Incompatible EndStatistics, cannot be added')
        for base in ('A', 'C', 'G', 'T', ''):
            self.adjacent_bases[base] += other.adjacent_bases[base]
        for length, error_dict in other.errors.items():
            for errors in error_dict:
                self.errors[length][errors] += other.errors[length][errors]
        return self

    @property
    def lengths(self):
        d = {length: sum(errors.values()) for length, errors in self.errors.items()}
        return d

    def random_match_probabilities(self, gc_content: float):
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
    """
    __slots__ = ['astart', 'astop', 'rstart', 'rstop', 'matches', 'errors', 'remove_before',
        'adapter', 'read', 'length', '_trimmed_read', 'adjacent_base']

    # TODO Can remove_before be removed from the constructor parameters?
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
            # Compute the trimmed read, assuming it’s a 'front' adapter
            self._trimmed_read = read[rstop:]
            self.adjacent_base = ''
        else:
            # Compute the trimmed read, assuming it’s a 'back' adapter
            self.adjacent_base = read.sequence[rstart - 1:rstart]
            self._trimmed_read = read[:rstart]
        self.remove_before = remove_before
        # Number of aligned characters in the adapter. If there are
        # indels, this may be different from the number of characters
        # in the read.
        self.length = astop - astart

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

    def update_statistics(self, statistics):
        """Update AdapterStatistics in place"""
        if self.remove_before:
            statistics.front.errors[self.rstop][self.errors] += 1
        else:
            statistics.back.errors[len(self.read) - len(self._trimmed_read)][self.errors] += 1
            try:
                statistics.back.adjacent_bases[self.adjacent_base] += 1
            except KeyError:
                statistics.back.adjacent_bases[''] = 1


def _generate_adapter_name(_start=[1]):
    name = str(_start[0])
    _start[0] += 1
    return name


class Adapter:
    """
    This class can find a single adapter characterized by sequence, error rate,
    type etc. within reads.

    where --  A Where enum value. This influences where the adapter is allowed to appear within the
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
        self.sequence = sequence.upper().replace('U', 'T')
        if not self.sequence:
            raise ValueError('Sequence is empty')
        self.where = where
        if remove not in (None, 'prefix', 'suffix', 'auto'):
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
        self.indels = indels
        if self.is_anchored and not self.indels:
            aligner_class = align.PrefixComparer if self.where is Where.PREFIX else align.SuffixComparer
            self.aligner = aligner_class(
                self.sequence,
                self.max_error_rate,
                wildcard_ref=self.adapter_wildcards,
                wildcard_query=self.read_wildcards,
                min_overlap=self.min_overlap
            )
        else:
            # TODO
            # Indels are suppressed by setting their cost very high, but a different algorithm
            # should be used instead.
            indel_cost = 1 if self.indels else 100000
            self.aligner = align.Aligner(
                self.sequence,
                self.max_error_rate,
                flags=self.where.value,
                wildcard_ref=self.adapter_wildcards,
                wildcard_query=self.read_wildcards,
                indel_cost=indel_cost,
                min_overlap=self.min_overlap,
            )

    def __repr__(self):
        return '<Adapter(name={name!r}, sequence={sequence!r}, where={where}, '\
            'remove={remove}, max_error_rate={max_error_rate}, min_overlap={min_overlap}, '\
            'read_wildcards={read_wildcards}, '\
            'adapter_wildcards={adapter_wildcards}, '\
            'indels={indels})>'.format(**vars(self))

    @property
    def is_anchored(self):
        """Return whether this adapter is anchored"""
        return self.where in {Where.PREFIX, Where.SUFFIX}

    @property
    def effective_length(self):
        return self.aligner.effective_length

    def enable_debug(self):
        """
        Print out the dynamic programming matrix after matching a read to an
        adapter.
        """
        self._debug = True
        self.aligner.enable_debug()

    def match_to(self, read):
        """
        Attempt to match this adapter to the given read.

        Return a Match instance if a match was found;
        return None if no match was found given the matching criteria (minimum
        overlap length, maximum error rate).
        """
        read_seq = read.sequence
        pos = -1

        # try to find an exact match first unless wildcards are allowed
        if not self.adapter_wildcards:
            if self.where is Where.PREFIX:
                pos = 0 if read_seq.startswith(self.sequence) else -1
            elif self.where is Where.SUFFIX:
                pos = (len(read_seq) - len(self.sequence)) if read_seq.endswith(self.sequence) else -1
            elif self.where is Where.BACK or self.where is Where.FRONT:
                pos = read_seq.find(self.sequence)
            # TODO BACK_NOT_INTERNAL, FRONT_NOT_INTERNAL
        if pos >= 0:
            match_args = (
                0, len(self.sequence), pos, pos + len(self.sequence),
                len(self.sequence), 0)
        else:
            # try approximate matching
            alignment = self.aligner.locate(read_seq)
            if self._debug:
                try:
                    print(self.aligner.dpmatrix)  # pragma: no cover
                except AttributeError:
                    pass
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
        match = Match(*match_args, remove_before=remove_before, adapter=self, read=read)

        assert match.length >= self.min_overlap
        return match

    def __len__(self):
        return len(self.sequence)

    def create_statistics(self):
        return AdapterStatistics(self)


class BackOrFrontAdapter(Adapter):
    """A 5' or 3' adapter.

    This is separate from the Adapter class so that a specialized match_to
    method can be implemented that reduces some of the runtime checks.

    TODO The generic Adapter class should become abstract, and the other
    adapter types should also get their own classes.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.where is Where.BACK or self.where is Where.FRONT
        self._remove_before = self.remove == 'prefix'

    def match_to(self, read):
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
            pos = read_seq.find(self.sequence)
        if pos >= 0:
            alignment = (
                0, len(self.sequence), pos, pos + len(self.sequence),
                len(self.sequence), 0)
        else:
            alignment = self.aligner.locate(read_seq)
        if self._debug:
            print(self.aligner.dpmatrix)  # pragma: no cover
        if alignment is None:
            return None

        match = Match(*alignment, remove_before=self._remove_before, adapter=self, read=read)
        return match


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
        m = 0
        if self.front_match is not None:
            m += self.front_match.matches
        if self.back_match is not None:
            m += self.back_match.matches
        return m

    @property
    def errors(self):
        e = 0
        if self.front_match is not None:
            e += self.front_match.errors
        if self.back_match is not None:
            e += self.back_match.errors
        return e

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
    def __init__(
        self,
        front_adapter,
        back_adapter,
        front_required,
        back_required,
        name,
    ):
        self.front_required = front_required
        self.back_required = back_required

        # The following attributes are needed for the report
        self.where = Where.LINKED
        self.name = _generate_adapter_name() if name is None else name
        self.front_adapter = front_adapter
        self.front_adapter.name = self.name
        self.back_adapter = back_adapter

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
        if self.front_required and front_match is None:
            return None

        if front_match is not None:
            # TODO statistics
            read = front_match.trimmed()
        back_match = self.back_adapter.match_to(read)
        if back_match is None and (self.back_required or front_match is None):
            return None
        return LinkedMatch(front_match, back_match, self)

    def create_statistics(self):
        return AdapterStatistics(self.front_adapter, self.back_adapter, where=Where.LINKED)

    @property
    def sequence(self):
        return self.front_adapter.sequence + "..." + self.back_adapter.sequence

    @property
    def remove(self):
        return None


class MultiAdapter:
    """
    Represent multiple adapters of the same type at once and use an index data structure
    to speed up matching. This acts like a "normal" Adapter as it provides a match_to
    method, but is faster with lots of adapters.

    There are quite a few restrictions:
    - the adapters need to be either all PREFIX or all SUFFIX adapters
    - no indels are allowed
    - the error rate allows at most 2 mismatches
    - wildcards in the adapter are not allowed
    - wildcards in the read are not allowed

    Use the is_acceptable() method to check individual adapters.
    """

    def __init__(self, adapters):
        """All given adapters must be of the same type, either Where.PREFIX or Where.SUFFIX"""
        if not adapters:
            raise ValueError("Adapter list is empty")
        self._where = adapters[0].where
        for adapter in adapters:
            self._accept(adapter)
            if adapter.where is not self._where:
                raise ValueError("All adapters must have identical 'where' attributes")
        self._adapters = adapters
        self._longest, self._index = self._make_index()

    def __repr__(self):
        return 'MultiAdapter(adapters={!r}, where={})'.format(self._adapters, self._where)

    @staticmethod
    def _accept(adapter):
        """Raise a ValueError if the adapter is not acceptable"""
        if adapter.where is not Where.PREFIX and adapter.where is not Where.SUFFIX:
            raise ValueError("Only anchored adapter types are allowed")
        if adapter.read_wildcards:
            raise ValueError("Wildcards in the read not supported")
        if adapter.adapter_wildcards:
            raise ValueError("Wildcards in the adapter not supported")
        if adapter.indels:
            raise ValueError("Indels not allowed")
        k = int(len(adapter) * adapter.max_error_rate)
        if k > 2:
            raise ValueError("Error rate too high")

    @staticmethod
    def is_acceptable(adapter):
        """
        Return whether this adapter is acceptable for being used by MultiAdapter

        Adapters are not acceptable if they allow wildcards, allow too many errors,
        or would lead to a very large index.
        """
        try:
            MultiAdapter._accept(adapter)
        except ValueError:
            return False
        return True

    def _make_index(self):
        logger.info('Building index of %s adapters ...', len(self._adapters))
        index = dict()
        longest = 0
        has_warned = False
        for adapter in self._adapters:
            sequence = adapter.sequence
            k = int(adapter.max_error_rate * len(sequence))
            for s, errors, matches in align.hamming_environment(sequence, k):
                if s in index:
                    other_adapter, other_errors, other_matches = index[s]
                    if matches < other_matches:
                        continue
                    if other_matches == matches and not has_warned:
                        logger.warning(
                            "Adapters %s %r and %s %r are very similar. At %s allowed errors, "
                            "the sequence %r cannot be assigned uniquely because the number of "
                            "matches is %s compared to both adapters.",
                            other_adapter.name, other_adapter.sequence, adapter.name,
                            adapter.sequence, k, s, matches
                        )
                        has_warned = True
                else:
                    index[s] = (adapter, errors, matches)
                longest = max(longest, len(s))
        logger.info('Built an index containing %s strings.', len(index))

        return longest, index

    def match_to(self, read):
        """
        Match the adapters against the read and return a Match that represents
        the best match or None if no match was found
        """
        if self._where is Where.PREFIX:
            def make_affix(n):
                return read.sequence[:n]
        else:
            def make_affix(n):
                return read.sequence[-n:]

        # Check all the prefixes of the read that could match
        best_adapter = None
        best_length = 0
        best_m = -1
        best_e = 1000
        # TODO do not go through all the lengths, only those that actually exist in the index
        for length in range(self._longest, -1, -1):
            if length < best_m:
                # No chance of getting the same or a higher number of matches, so we can stop early
                break

            affix = make_affix(length)
            try:
                adapter, e, m = self._index[affix]
            except KeyError:
                continue
            if m > best_m or (m == best_m and e < best_e):
                best_adapter = adapter
                best_e = e
                best_m = m
                best_length = length

        if best_m == -1:
            return None
        else:
            if self._where is Where.PREFIX:
                rstart, rstop = 0, best_length
            else:
                assert self._where is Where.SUFFIX
                rstart, rstop = len(read) - best_length, len(read)
            return Match(
                astart=0,
                astop=len(best_adapter.sequence),
                rstart=rstart,
                rstop=rstop,
                matches=best_m,
                errors=best_e,
                remove_before=best_adapter.remove == 'prefix',
                adapter=best_adapter,
                read=read
            )


def warn_duplicate_adapters(adapters):
    d = dict()
    for adapter in adapters:
        key = (adapter.sequence, adapter.where, adapter.remove)
        if key in d:
            logger.warning("Adapter %r (%s) was specified multiple times! "
                "Please make sure that this is what you want.",
                adapter.sequence, ADAPTER_TYPE_NAMES[adapter.where])
        d[key] = adapter.name
