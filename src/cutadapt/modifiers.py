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
        assert action in ('trim', 'mask', 'lowercase', None)
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
        Split adapters into three different categories so that they can possibly be used
        with a MultiAdapter. Return a tuple (prefix, suffix, other), where
        - prefix is a list of all anchored 5' adapters that MultiAdapter would accept
        - suffix is a list of all anchored 3' adapters that MultiAdapter would accept
        - other is a list of all remaining adapters.
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

    @staticmethod
    def best_match(adapters, read):
        """
        Find the best matching adapter in the given read.

        Return either a Match instance or None if there are no matches.
        """
        best_match = None
        for adapter in adapters:
            match = adapter.match_to(read)
            if match is None:
                continue

            # the no. of matches determines which adapter fits best
            if best_match is None or match.matches > best_match.matches:
                best_match = match
        return best_match

    @staticmethod
    def remainder(matches):
        """
        Determine which part of the read was not trimmed. Return a tuple (start, stop)
        that gives the interval of the untrimmed part relative to the original read.

        matches is a list of Match objects. The original read is assumed to be
        matches[0].read
        """
        # Start with the full read
        read = matches[0].read
        start, stop = 0, len(read)
        for match in matches:
            if match.remove_before:
                # Length of the prefix that was removed
                start += match.rstop
            else:
                # Length of the suffix that was removed
                stop -= len(match.read) - match.rstart
        return (start, stop)

    @staticmethod
    def masked_read(trimmed_read, matches):
        start, stop = AdapterCutter.remainder(matches)
        read = matches[0].read
        # TODO modification in place
        trimmed_read.sequence = (
            'N' * start
            + read.sequence[start:stop]
            + 'N' * (len(read) - stop))
        trimmed_read.qualities = read.qualities
        return trimmed_read

    @staticmethod
    def lowercased_read(trimmed_read, matches):
        start, stop = AdapterCutter.remainder(matches)
        read_sequence = matches[0].read.sequence
        # TODO modification in place
        trimmed_read.sequence = (
            read_sequence[:start].lower()
            + read_sequence[start:stop].upper()
            + read_sequence[stop:].lower()
        )
        trimmed_read.qualities = matches[0].read.qualities
        return trimmed_read

    def __call__(self, read, matches):
        """
        Search for the best-matching adapter in a read, perform the requested action
        ('trim', 'mask', 'lowercase' or None as determined by self.action) and return the
        (possibly) modified read.

        *self.times* adapter removal rounds are done. During each round,
        only the best-matching adapter is trimmed. If no adapter was found in a round,
        no further rounds are attempted.

        The 'matches' parameter needs to be a list. Every time an adapter is found,
        a Match object describing the match will be appended to it.
        """
        trimmed_read = read
        if self.action == 'lowercase':
            trimmed_read.sequence = trimmed_read.sequence.upper()
        for _ in range(self.times):
            match = AdapterCutter.best_match(self.adapters, trimmed_read)
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
        elif self.action == 'lowercase':
            trimmed_read = self.lowercased_read(trimmed_read, matches)
            assert len(trimmed_read.sequence) == len(read)
        elif self.action is None:  # --no-trim
            trimmed_read = read[:]

        self.with_adapters += 1
        return trimmed_read


class PairedAdapterCutterError(Exception):
    pass


class PairedAdapterCutter:
    """
    A Modifier that trims adapter pairs from R1 and R2.
    """

    def __init__(self, adapters1, adapters2, action='trim'):
        """
        adapters1 -- list of Adapters to be removed from R1
        adapters2 -- list of Adapters to be removed from R1

        Both lists must have the same, non-zero length.
         read pair is trimmed if adapters1[i] is found in R1 and adapters2[i] in R2.

        action -- What to do with a found adapter: None, 'trim', or 'mask'
        """
        if len(adapters1) != len(adapters2):
            raise PairedAdapterCutterError(
                "The number of reads to trim from R1 and R2 must be the same. "
                "Given: {} for R1, {} for R2".format(len(adapters1), len(adapters2)))
        if not adapters1:
            raise PairedAdapterCutterError("No adapters given")
        self._adapters1 = adapters1
        self._adapter_indices = {a: i for i, a in enumerate(adapters1)}
        self._adapters2 = adapters2
        self.action = action
        self.with_adapters = 0
        self.adapter_statistics = [None, None]
        self.adapter_statistics[0] = OrderedDict((a, a.create_statistics()) for a in adapters1)
        self.adapter_statistics[1] = OrderedDict((a, a.create_statistics()) for a in adapters2)

    def __repr__(self):
        return 'PairedAdapterCutter(adapters1={!r}, adapters2={!r})'.format(
            self._adapters1, self._adapters2)

    def __call__(self, read1, read2, matches1, matches2):
        """
        """
        match1 = AdapterCutter.best_match(self._adapters1, read1)
        if match1 is None:
            return read1, read2
        adapter1 = match1.adapter
        adapter2 = self._adapters2[self._adapter_indices[adapter1]]
        match2 = adapter2.match_to(read2)
        if match2 is None:
            return read1, read2

        self.with_adapters += 1
        result = []
        for i, match, read in zip([0, 1], [match1, match2], [read1, read2]):
            trimmed_read = read
            if self.action == 'lowercase':
                trimmed_read.sequence = trimmed_read.sequence.upper()

            trimmed_read = match.trimmed()
            match.update_statistics(self.adapter_statistics[i][match.adapter])

            if self.action == 'trim':
                # read is already trimmed, nothing to do
                pass
            elif self.action == 'mask':
                trimmed_read = AdapterCutter.masked_read(trimmed_read, [match])
            elif self.action == 'lowercase':
                trimmed_read = AdapterCutter.lowercased_read(trimmed_read, [match])
                assert len(trimmed_read.sequence) == len(read)
            elif self.action is None:  # --no-trim
                trimmed_read = read[:]
            result.append(trimmed_read)
        matches1.append(match1)
        matches2.append(match2)
        return result


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
