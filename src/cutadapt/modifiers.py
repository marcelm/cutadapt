"""
This module implements all the read modifications that cutadapt supports.
A modifier must be callable and typically implemented as a class with a
__call__ method.
"""
import re
from abc import ABC, abstractmethod
from collections import OrderedDict
from .qualtrim import quality_trim_index, nextseq_trim_index
from .adapters import Where, MultiAdapter


class PairedModifier:
    """
    Wrapper for modifiers that work on both reads in a paired-end read
    """
    paired = True

    def __init__(self, modifier1, modifier2):
        """Set one of the modifiers to None to work on R1 or R2 only"""
        self._modifier1 = modifier1
        self._modifier2 = modifier2

    def __repr__(self):
        return 'PairedModifier(modifier1={!r}, modifier2={!r})'.format(
            self._modifier1, self._modifier2)

    def __call__(self, read1, read2, matches1, matches2):
        if self._modifier1 is None:
            return read1, self._modifier2(read2, matches2)
        if self._modifier2 is None:
            return self._modifier1(read1, matches1), read2
        return self._modifier1(read1, matches1), self._modifier2(read2, matches2)


class Modifier(ABC):
    @abstractmethod
    def __call__(self, read, matches):
        pass


class AdapterCutter(Modifier):
    """
    Repeatedly find one of multiple adapters in reads.
    The number of times the search is repeated is specified by the
    times parameter.
    """

    def __init__(self, adapters, times=1, action='trim'):
        """
        adapters -- list of Adapter objects

        action -- What to do with a found adapter: None, 'trim', or 'mask'
        """
        self.times = times
        assert action in ('trim', 'mask', None)
        self.action = action
        self.with_adapters = 0
        self.adapter_statistics = OrderedDict((a, a.create_statistics()) for a in adapters)

        prefix, suffix, other = self._split_adapters(adapters)
        # For somewhat better backwards compatibility, avoid re-ordering
        # the adapters when we don’t need to
        if len(prefix) > 1 or len(suffix) > 1:
            adapters = other
            for affix in (prefix, suffix):
                if len(affix) > 1:
                    adapters.append(MultiAdapter(affix))
                else:
                    adapters.extend(affix)
        self.adapters = adapters

    def __repr__(self):
        return 'AdapterCutter(adapters={!r}, times={}, action={!r})'.format(
            self.adapters, self.times, self.action)

    @staticmethod
    def _split_adapters(adapters):
        """
        Split adapters for MultiAdapter:
        - acceptable anchored 5'
        - acceptable anchored 3'
        - other
        """
        prefix, suffix, other = [], [], []
        for a in adapters:

            if MultiAdapter.is_acceptable(a):
                if a.where == Where.PREFIX:
                    lst = prefix
                else:
                    assert a.where == Where.SUFFIX
                    lst = suffix
            else:
                lst = other
            lst.append(a)
        return prefix, suffix, other

    def _best_match(self, read):
        """
        Find the best matching adapter in the given read.

        Return either a Match instance or None if there are no matches.
        """
        best_match = None
        for adapter in self.adapters:
            match = adapter.match_to(read)
            if match is None:
                continue

            # the no. of matches determines which adapter fits best
            if best_match is None or match.matches > best_match.matches:
                best_match = match
        return best_match

    @staticmethod
    def masked_read(trimmed_read, matches):
        # add N from last modification
        masked_sequence = trimmed_read.sequence
        for match in sorted(matches, reverse=True, key=lambda m: m.astart):
            ns = 'N' * (len(match.read.sequence) -
                        len(match.trimmed().sequence))
            # add N depending on match position
            if match.remove_before:
                masked_sequence = ns + masked_sequence
            else:
                masked_sequence += ns
        # set masked sequence as sequence with original quality
        # TODO modification in place
        trimmed_read.sequence = masked_sequence
        trimmed_read.qualities = matches[0].read.qualities
        return trimmed_read

    def __call__(self, read, matches):
        """
        Search for the best-matching adapter in a read, perform the requested action
        ('trim', 'mask', or None as determined by self.action) it and return the
        (possibly) modified read.

        *self.times* adapter removal rounds are done. During each round,
        only the best-matching adapter is trimmed. If no adapter was found in a round,
        no further rounds are attempted.

        The 'matches' parameter needs to be a list. Every time an adapter is found,
        a Match object describing the match will be appended to it.

        The read is converted to uppercase before it is compared to the adapter
        sequences.
        """
        trimmed_read = read
        for _ in range(self.times):
            match = self._best_match(trimmed_read)
            if match is None:
                # if nothing found, attempt no further rounds
                break
            matches.append(match)
            trimmed_read = match.trimmed()
            match.update_statistics(self.adapter_statistics[match.adapter])

        if not matches:
            return trimmed_read

        if self.action == 'trim':
            # read is already trimmed, nothing to do
            pass
        elif self.action == 'mask':
            trimmed_read = self.masked_read(trimmed_read, matches)
            assert len(trimmed_read.sequence) == len(read)
        elif self.action is None:  # --no-trim
            trimmed_read = read[:]

        self.with_adapters += 1
        return trimmed_read


class UnconditionalCutter(Modifier):
    """
    A modifier that unconditionally removes the first n or the last n bases from a read.

    If the length is positive, the bases are removed from the beginning of the read.
    If the length is negative, the bases are removed from the end of the read.
    """
    def __init__(self, length):
        self.length = length

    def __call__(self, read, matches):
        if self.length > 0:
            return read[self.length:]
        elif self.length < 0:
            return read[:self.length]


class LengthTagModifier(Modifier):
    """
    Replace "length=..." strings in read names.
    """
    def __init__(self, length_tag):
        self.regex = re.compile(r"\b" + length_tag + r"[0-9]*\b")
        self.length_tag = length_tag

    def __call__(self, read, matches):
        read = read[:]
        if read.name.find(self.length_tag) >= 0:
            read.name = self.regex.sub(self.length_tag + str(len(read.sequence)), read.name)
        return read


class SuffixRemover(Modifier):
    """
    Remove a given suffix from read names.
    """
    def __init__(self, suffix):
        self.suffix = suffix

    def __call__(self, read, matches):
        read = read[:]
        if read.name.endswith(self.suffix):
            read.name = read.name[:-len(self.suffix)]
        return read


class PrefixSuffixAdder(Modifier):
    """
    Add a suffix and a prefix to read names
    """
    def __init__(self, prefix, suffix):
        self.prefix = prefix
        self.suffix = suffix

    def __call__(self, read, matches):
        read = read[:]
        adapter_name = matches[-1].adapter.name if matches else 'no_adapter'
        read.name = self.prefix.replace('{name}', adapter_name) + read.name + \
            self.suffix.replace('{name}', adapter_name)
        return read


class ZeroCapper(Modifier):
    """
    Change negative quality values of a read to zero
    """
    def __init__(self, quality_base=33):
        qb = quality_base
        self.zero_cap_trans = str.maketrans(''.join(map(chr, range(qb))), chr(qb) * qb)

    def __call__(self, read, matches):
        read = read[:]
        read.qualities = read.qualities.translate(self.zero_cap_trans)
        return read


class NextseqQualityTrimmer(Modifier):
    def __init__(self, cutoff, base):
        self.cutoff = cutoff
        self.base = base
        self.trimmed_bases = 0

    def __call__(self, read, matches):
        stop = nextseq_trim_index(read, self.cutoff, self.base)
        self.trimmed_bases += len(read) - stop
        return read[:stop]


class QualityTrimmer(Modifier):
    def __init__(self, cutoff_front, cutoff_back, base):
        self.cutoff_front = cutoff_front
        self.cutoff_back = cutoff_back
        self.base = base
        self.trimmed_bases = 0

    def __call__(self, read, matches):
        start, stop = quality_trim_index(read.qualities, self.cutoff_front, self.cutoff_back, self.base)
        self.trimmed_bases += len(read) - (stop - start)
        return read[start:stop]


class Shortener(Modifier):
    """Unconditionally shorten a read to the given length

    If the length is positive, the bases are removed from the end of the read.
    If the length is negative, the bases are removed from the beginning of the read.
    """
    def __init__(self, length):
        self.length = length

    def __call__(self, read, matches):
        if self.length >= 0:
            return read[:self.length]
        else:
            return read[self.length:]


class NEndTrimmer(Modifier):
    """Trims Ns from the 3' and 5' end of reads"""
    def __init__(self):
        self.start_trim = re.compile(r'^N+')
        self.end_trim = re.compile(r'N+$')

    def __call__(self, read, matches):
        sequence = read.sequence
        start_cut = self.start_trim.match(sequence)
        end_cut = self.end_trim.search(sequence)
        start_cut = start_cut.end() if start_cut else 0
        end_cut = end_cut.start() if end_cut else len(read)
        return read[start_cut:end_cut]
