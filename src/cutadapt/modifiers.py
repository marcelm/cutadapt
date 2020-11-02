"""
This module implements all the read modifications that cutadapt supports.
A modifier must be callable and typically implemented as a class with a
__call__ method.
"""
import re
from typing import Sequence, List, Tuple, Optional
from abc import ABC, abstractmethod
from collections import OrderedDict

from .qualtrim import quality_trim_index, nextseq_trim_index
from .adapters import MultipleAdapters, SingleAdapter, IndexedPrefixAdapters, IndexedSuffixAdapters, Match, remainder
from .utils import reverse_complemented_sequence


class ModificationInfo:
    """
    An object of this class is created for each read that passes through the pipeline.
    Any information (except the read itself) that needs to be passed from one modifier
    to one later in the pipeline or from one modifier to the filters is recorded here.
    """
    __slots__ = ["matches", "original_read"]

    def __init__(self, read):
        self.matches = []  # type: List[Match]
        self.original_read = read


class SingleEndModifier(ABC):
    @abstractmethod
    def __call__(self, read, info: ModificationInfo):
        pass


class PairedModifier(ABC):
    @abstractmethod
    def __call__(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        pass


class PairedModifierWrapper(PairedModifier):
    """
    Wrapper for modifiers that work on both reads in a paired-end read
    """
    paired = True

    def __init__(self, modifier1: Optional[SingleEndModifier], modifier2: Optional[SingleEndModifier]):
        """Set one of the modifiers to None to work on R1 or R2 only"""
        self._modifier1 = modifier1
        self._modifier2 = modifier2
        if self._modifier1 is None and self._modifier2 is None:
            raise ValueError("Not both modifiers may be None")

    def __repr__(self):
        return 'PairedModifier(modifier1={!r}, modifier2={!r})'.format(
            self._modifier1, self._modifier2)

    def __call__(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        if self._modifier1 is None:
            return read1, self._modifier2(read2, info2)  # type: ignore
        if self._modifier2 is None:
            return self._modifier1(read1, info1), read2
        return self._modifier1(read1, info1), self._modifier2(read2, info2)


class AdapterCutter(SingleEndModifier):
    """
    Repeatedly find one of multiple adapters in reads.
    The number of times the search is repeated is specified by the
    times parameter.
    """

    def __init__(
        self,
        adapters: List[SingleAdapter],
        times: int = 1,
        action: str = "trim",
        index: bool = True,
    ):
        """
        action -- What to do with a found adapter: None, 'trim', 'mask' or 'lowercase'

        index -- if True, an adapter index (for multiple adapters) is created if possible
        """
        self.times = times
        assert action in ('trim', 'mask', 'lowercase', None)
        self.action = action
        self.with_adapters = 0
        self.adapter_statistics = OrderedDict((a, a.create_statistics()) for a in adapters)
        if index:
            self.adapters = MultipleAdapters(self._regroup_into_indexed_adapters(adapters))
        else:
            self.adapters = MultipleAdapters(adapters)

    def __repr__(self):
        return 'AdapterCutter(adapters={!r}, times={}, action={!r})'.format(
            self.adapters, self.times, self.action)

    def _regroup_into_indexed_adapters(self, adapters):
        prefix, suffix, single = self._split_adapters(adapters)
        # For somewhat better backwards compatibility, avoid re-ordering
        # the adapters when we don’t need to
        if len(prefix) > 1 or len(suffix) > 1:
            result = single
            if len(prefix) > 1:
                result.append(IndexedPrefixAdapters(prefix))
            else:
                result.extend(prefix)
            if len(suffix) > 1:
                result.append(IndexedSuffixAdapters(suffix))
            else:
                result.extend(suffix)
            return result
        else:
            return adapters

    @staticmethod
    def _split_adapters(
        adapters: Sequence[SingleAdapter]
    ) -> Tuple[Sequence[SingleAdapter], Sequence[SingleAdapter], Sequence[SingleAdapter]]:
        """
        Split adapters into three different categories so that they can possibly be used
        with a MultiAdapter. Return a tuple (prefix, suffix, other), where
        - prefix is a list of all anchored 5' adapters that MultiAdapter would accept
        - suffix is a list of all anchored 3' adapters that MultiAdapter would accept
        - other is a list of all remaining adapters.
        """
        prefix = []  # type: List[SingleAdapter]
        suffix = []  # type: List[SingleAdapter]
        other = []  # type: List[SingleAdapter]
        for a in adapters:
            if IndexedPrefixAdapters.is_acceptable(a):
                prefix.append(a)
            elif IndexedSuffixAdapters.is_acceptable(a):
                suffix.append(a)
            else:
                other.append(a)
        return prefix, suffix, other

    @staticmethod
    def masked_read(read, matches: Sequence[Match]):
        start, stop = remainder(matches)
        result = read[:]
        result.sequence = (
            'N' * start
            + read.sequence[start:stop]
            + 'N' * (len(read) - stop))
        return result

    @staticmethod
    def lowercased_read(read, matches: Sequence[Match]):
        start, stop = remainder(matches)
        result = read[:]
        result.sequence = (
            read.sequence[:start].lower()
            + read.sequence[start:stop].upper()
            + read.sequence[stop:].lower()
        )
        return result

    def __call__(self, read, info: ModificationInfo):
        trimmed_read, matches = self.match_and_trim(read)
        if matches:
            self.with_adapters += 1
            for match in matches:
                match.update_statistics(self.adapter_statistics[match.adapter])
        info.matches.extend(matches)  # TODO extend or overwrite?
        return trimmed_read

    def match_and_trim(self, read):
        """
        Search for the best-matching adapter in a read, perform the requested action
        ('trim', 'mask', 'lowercase' or None as determined by self.action) and return the
        (possibly) modified read.

        *self.times* adapter removal rounds are done. During each round,
        only the best-matching adapter is trimmed. If no adapter was found in a round,
        no further rounds are attempted.

        Return a pair (trimmed_read, matches), where matches is a list of Match instances.
        """
        matches = []
        if self.action == 'lowercase':
            read.sequence = read.sequence.upper()
        trimmed_read = read
        for _ in range(self.times):
            match = self.adapters.match_to(trimmed_read.sequence)
            if match is None:
                # if nothing found, attempt no further rounds
                break
            matches.append(match)
            trimmed_read = match.trimmed(trimmed_read)

        if not matches:
            return trimmed_read, []

        if self.action == 'trim':
            # read is already trimmed, nothing to do
            pass
        elif self.action == 'mask':
            trimmed_read = self.masked_read(read, matches)
        elif self.action == 'lowercase':
            trimmed_read = self.lowercased_read(read, matches)
            assert len(trimmed_read.sequence) == len(read)
        elif self.action is None:  # --no-trim
            trimmed_read = read[:]

        return trimmed_read, matches


class ReverseComplementer(SingleEndModifier):
    """Trim adapters from a read and its reverse complement"""

    def __init__(self, adapter_cutter: AdapterCutter, rc_suffix: Optional[str] = " rc"):
        """
        rc_suffix -- suffix to add to the read name if sequence was reverse-complemented
        """
        self.adapter_cutter = adapter_cutter
        self.reverse_complemented = 0
        self._suffix = rc_suffix

    def __call__(self, read, info: ModificationInfo):
        reverse_read = reverse_complemented_sequence(read)

        forward_trimmed_read, forward_matches = self.adapter_cutter.match_and_trim(read)
        reverse_trimmed_read, reverse_matches = self.adapter_cutter.match_and_trim(reverse_read)

        forward_match_count = sum(m.matches for m in forward_matches)
        reverse_match_count = sum(m.matches for m in reverse_matches)
        use_reverse_complement = reverse_match_count > forward_match_count

        if use_reverse_complement:
            self.reverse_complemented += 1
            assert reverse_matches
            trimmed_read, matches = reverse_trimmed_read, reverse_matches
            if self._suffix:
                trimmed_read.name += self._suffix
        else:
            trimmed_read, matches = forward_trimmed_read, forward_matches

        if matches:
            self.adapter_cutter.with_adapters += 1
            for match in matches:
                stats = self.adapter_cutter.adapter_statistics[match.adapter]
                match.update_statistics(stats)
                stats.reverse_complemented += bool(use_reverse_complement)
            info.matches.extend(matches)  # TODO extend or overwrite?
        return trimmed_read


class PairedAdapterCutterError(Exception):
    pass


class PairedAdapterCutter(PairedModifier):
    """
    A Modifier that trims adapter pairs from R1 and R2.
    """

    def __init__(self, adapters1, adapters2, action='trim'):
        """
        adapters1 -- list of Adapters to be removed from R1
        adapters2 -- list of Adapters to be removed from R1

        Both lists must have the same, non-zero length.
         read pair is trimmed if adapters1[i] is found in R1 and adapters2[i] in R2.

        action -- What to do with a found adapter: None, 'trim', 'lowercase' or 'mask'
        """
        super().__init__()
        if len(adapters1) != len(adapters2):
            raise PairedAdapterCutterError(
                "The number of reads to trim from R1 and R2 must be the same. "
                "Given: {} for R1, {} for R2".format(len(adapters1), len(adapters2)))
        if not adapters1:
            raise PairedAdapterCutterError("No adapters given")
        self._adapters1 = MultipleAdapters(adapters1)
        self._adapter_indices = {a: i for i, a in enumerate(adapters1)}
        self._adapters2 = MultipleAdapters(adapters2)
        self.action = action
        self.with_adapters = 0
        self.adapter_statistics = [None, None]
        self.adapter_statistics[0] = OrderedDict((a, a.create_statistics()) for a in adapters1)
        self.adapter_statistics[1] = OrderedDict((a, a.create_statistics()) for a in adapters2)

    def __repr__(self):
        return 'PairedAdapterCutter(adapters1={!r}, adapters2={!r})'.format(
            self._adapters1, self._adapters2)

    def __call__(self, read1, read2, info1, info2):
        """
        """
        match1 = self._adapters1.match_to(read1.sequence)
        if match1 is None:
            return read1, read2
        adapter1 = match1.adapter
        adapter2 = self._adapters2[self._adapter_indices[adapter1]]
        match2 = adapter2.match_to(read2.sequence)
        if match2 is None:
            return read1, read2

        self.with_adapters += 1
        result = []
        for i, match, read in zip([0, 1], [match1, match2], [read1, read2]):
            trimmed_read = read
            if self.action == 'lowercase':
                trimmed_read.sequence = trimmed_read.sequence.upper()

            trimmed_read = match.trimmed(trimmed_read)
            match.update_statistics(self.adapter_statistics[i][match.adapter])

            if self.action == 'trim':
                # read is already trimmed, nothing to do
                pass
            elif self.action == 'mask':
                trimmed_read = AdapterCutter.masked_read(read, [match])
            elif self.action == 'lowercase':
                trimmed_read = AdapterCutter.lowercased_read(read, [match])
                assert len(trimmed_read.sequence) == len(read)
            elif self.action is None:  # --no-trim
                trimmed_read = read[:]
            result.append(trimmed_read)
        info1.matches.append(match1)
        info2.matches.append(match2)
        return result


class UnconditionalCutter(SingleEndModifier):
    """
    A modifier that unconditionally removes the first n or the last n bases from a read.

    If the length is positive, the bases are removed from the beginning of the read.
    If the length is negative, the bases are removed from the end of the read.
    """
    def __init__(self, length: int):
        self.length = length

    def __call__(self, read, info: ModificationInfo):
        if self.length > 0:
            return read[self.length:]
        elif self.length < 0:
            return read[:self.length]


class LengthTagModifier(SingleEndModifier):
    """
    Replace "length=..." strings in read names.
    """
    def __init__(self, length_tag):
        self.regex = re.compile(r"\b" + length_tag + r"[0-9]*\b")
        self.length_tag = length_tag

    def __call__(self, read, info: ModificationInfo):
        read = read[:]
        if read.name.find(self.length_tag) >= 0:
            read.name = self.regex.sub(self.length_tag + str(len(read.sequence)), read.name)
        return read


class SuffixRemover(SingleEndModifier):
    """
    Remove a given suffix from read names.
    """
    def __init__(self, suffix):
        self.suffix = suffix

    def __call__(self, read, info: ModificationInfo):
        read = read[:]
        if read.name.endswith(self.suffix):
            read.name = read.name[:-len(self.suffix)]
        return read


class PrefixSuffixAdder(SingleEndModifier):
    """
    Add a suffix and a prefix to read names
    """
    def __init__(self, prefix, suffix):
        self.prefix = prefix
        self.suffix = suffix

    def __call__(self, read, info):
        read = read[:]
        adapter_name = info.matches[-1].adapter.name if info.matches else 'no_adapter'
        read.name = self.prefix.replace('{name}', adapter_name) + read.name + \
            self.suffix.replace('{name}', adapter_name)
        return read


class ZeroCapper(SingleEndModifier):
    """
    Change negative quality values of a read to zero
    """
    def __init__(self, quality_base=33):
        qb = quality_base
        self.zero_cap_trans = str.maketrans(''.join(map(chr, range(qb))), chr(qb) * qb)

    def __call__(self, read, info: ModificationInfo):
        read = read[:]
        read.qualities = read.qualities.translate(self.zero_cap_trans)
        return read


class NextseqQualityTrimmer(SingleEndModifier):
    def __init__(self, cutoff, base):
        self.cutoff = cutoff
        self.base = base
        self.trimmed_bases = 0

    def __call__(self, read, info: ModificationInfo):
        stop = nextseq_trim_index(read, self.cutoff, self.base)
        self.trimmed_bases += len(read) - stop
        return read[:stop]


class QualityTrimmer(SingleEndModifier):
    def __init__(self, cutoff_front, cutoff_back, base):
        self.cutoff_front = cutoff_front
        self.cutoff_back = cutoff_back
        self.base = base
        self.trimmed_bases = 0

    def __call__(self, read, info: ModificationInfo):
        start, stop = quality_trim_index(read.qualities, self.cutoff_front, self.cutoff_back, self.base)
        self.trimmed_bases += len(read) - (stop - start)
        return read[start:stop]


class Shortener(SingleEndModifier):
    """Unconditionally shorten a read to the given length

    If the length is positive, the bases are removed from the end of the read.
    If the length is negative, the bases are removed from the beginning of the read.
    """
    def __init__(self, length):
        self.length = length

    def __call__(self, read, info: ModificationInfo):
        if self.length >= 0:
            return read[:self.length]
        else:
            return read[self.length:]


class NEndTrimmer(SingleEndModifier):
    """Trims Ns from the 3' and 5' end of reads"""
    def __init__(self):
        self.start_trim = re.compile(r'^N+')
        self.end_trim = re.compile(r'N+$')

    def __call__(self, read, info: ModificationInfo):
        sequence = read.sequence
        start_cut = self.start_trim.match(sequence)
        end_cut = self.end_trim.search(sequence)
        start_cut = start_cut.end() if start_cut else 0
        end_cut = end_cut.start() if end_cut else len(read)
        return read[start_cut:end_cut]
