"""
This module implements all the read modifications that cutadapt supports.
A modifier must be callable and typically implemented as a class with a
__call__ method.
"""
import re
import logging
from collections import defaultdict
from types import SimpleNamespace
from typing import Sequence, List, Tuple, Optional, Set
from abc import ABC, abstractmethod

from dnaio import record_names_match, SequenceRecord

from .qualtrim import quality_trim_index, nextseq_trim_index, poly_a_trim_index
from .adapters import (
    MultipleAdapters,
    SingleAdapter,
    IndexedPrefixAdapters,
    IndexedSuffixAdapters,
    Match,
    remainder,
    Adapter,
)
from .tokenizer import tokenize_braces, TokenizeError, Token, BraceToken
from .info import ModificationInfo

logger = logging.getLogger()


# If the number of prefix or suffix adapters is higher than this, switch to using an index
INDEXING_THRESHOLD = 5


class SingleEndModifier(ABC):
    @abstractmethod
    def __call__(self, read: SequenceRecord, info: ModificationInfo):
        pass


class PairedEndModifier(ABC):
    @abstractmethod
    def __call__(
        self,
        read1: SequenceRecord,
        read2: SequenceRecord,
        info1: ModificationInfo,
        info2: ModificationInfo,
    ) -> Tuple[SequenceRecord, SequenceRecord]:
        pass


class PairedEndModifierWrapper(PairedEndModifier):
    """
    Wrap two SingleEndModifiers that work on both reads in a paired-end read
    """

    paired = True

    def __init__(
        self,
        modifier1: Optional[SingleEndModifier],
        modifier2: Optional[SingleEndModifier],
    ):
        """Set one of the modifiers to None to work on R1 or R2 only"""
        self._modifier1 = modifier1
        self._modifier2 = modifier2
        if self._modifier1 is None and self._modifier2 is None:
            raise ValueError("Not both modifiers may be None")

    def __repr__(self):
        return (
            "PairedEndModifierWrapper("
            f"modifier1={self._modifier1!r}, modifier2={self._modifier2!r})"
        )

    def __call__(self, read1, read2, info1: ModificationInfo, info2: ModificationInfo):
        if self._modifier1 is None:
            return read1, self._modifier2(read2, info2)  # type: ignore
        if self._modifier2 is None:
            return self._modifier1(read1, info1), read2
        return self._modifier1(read1, info1), self._modifier2(read2, info2)


class AdapterCutter(SingleEndModifier):
    """
    Repeatedly find one of multiple adapters in reads.

    Arguments:
        adapters: Adapters to be searched
        times: Repeat the search this number of times.
        action: What to do with a found adapter.
            - *None*: Do nothing, only update the ModificationInfo appropriately
            - "trim": Remove the adapter and down- or upstream sequence depending on adapter type
            - "mask": Replace the part of the sequence that would have been removed with "N" bases
            - "lowercase": Convert the part of the sequence that would have been removed to lowercase
            - "retain": Like "trim", but leave the adapter sequence itself in the read
        index: If True, attempt to create an index to speed up the search (if possible)
    """

    def __init__(
        self,
        adapters: Sequence[Adapter],
        times: int = 1,
        action: Optional[str] = "trim",
        index: bool = True,
    ):
        self.times = times
        assert action in ("trim", "mask", "lowercase", "retain", None)
        self.action = action
        self.with_adapters = 0
        self.adapter_statistics = {a: a.create_statistics() for a in adapters}
        if index:
            self.adapters = MultipleAdapters(
                self._regroup_into_indexed_adapters(adapters)
            )
        else:
            self.adapters = MultipleAdapters(adapters)
        if action == "retain" and times > 1:
            raise ValueError("'retain' cannot be combined with times > 1")
        if self.times == 1 and self.action == "trim":
            self.match_and_trim = self._match_and_trim_once_action_trim  # type: ignore

    def __repr__(self):
        return (
            "AdapterCutter("
            f"adapters={self.adapters!r}, times={self.times}, action='{self.action}')"
        )

    def _regroup_into_indexed_adapters(self, adapters):
        prefix, suffix, single = self._split_adapters(adapters)
        # For somewhat better backwards compatibility, avoid re-ordering
        # the adapters when we don’t need to
        if len(prefix) > INDEXING_THRESHOLD or len(suffix) > INDEXING_THRESHOLD:
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
        adapters: Sequence[SingleAdapter],
    ) -> Tuple[
        Sequence[SingleAdapter], Sequence[SingleAdapter], Sequence[SingleAdapter]
    ]:
        """
        Split adapters into three different categories so that they can possibly be used
        with a MultiAdapter. Return a tuple (prefix, suffix, other), where
        - prefix is a list of all anchored 5' adapters that MultiAdapter would accept
        - suffix is a list of all anchored 3' adapters that MultiAdapter would accept
        - other is a list of all remaining adapters.
        """
        prefix: List[SingleAdapter] = []
        suffix: List[SingleAdapter] = []
        other: List[SingleAdapter] = []
        for a in adapters:
            if IndexedPrefixAdapters.is_acceptable(a):
                prefix.append(a)
            elif IndexedSuffixAdapters.is_acceptable(a):
                suffix.append(a)
            else:
                other.append(a)
        return prefix, suffix, other

    @staticmethod
    def trim_but_retain_adapter(read, matches: Sequence[Match]):
        start, stop = matches[-1].retained_adapter_interval()
        return read[start:stop]

    @staticmethod
    def masked_read(read, matches: Sequence[Match]):
        start, stop = remainder(matches)
        result = read[:]
        result.sequence = (
            "N" * start + read.sequence[start:stop] + "N" * (len(read) - stop)
        )
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
                self.adapter_statistics[match.adapter].add_match(match)
        info.matches.extend(matches)  # TODO extend or overwrite?
        return trimmed_read

    def match_and_trim(self, read):
        """
        Search for the best-matching adapter in a read, perform the requested action
        ('trim', 'mask' etc. as determined by self.action) and return the
        (possibly) modified read.

        *self.times* adapter removal rounds are done. During each round,
        only the best-matching adapter is trimmed. If no adapter was found in a round,
        no further rounds are attempted.

        Return a pair (trimmed_read, matches), where matches is a list of Match instances.
        """
        matches = []
        if self.action == "lowercase":  # TODO this should not be needed
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

        if self.action == "trim":
            # read is already trimmed, nothing to do
            pass
        elif self.action == "retain":
            trimmed_read = self.trim_but_retain_adapter(read, matches)
        elif self.action == "mask":
            trimmed_read = self.masked_read(read, matches)
        elif self.action == "lowercase":
            trimmed_read = self.lowercased_read(read, matches)
            assert len(trimmed_read.sequence) == len(read)
        elif self.action is None:
            trimmed_read = read[:]

        return trimmed_read, matches

    def _match_and_trim_once_action_trim(self, read):
        """
        Specalization of match_and_trim for the case that self.times == 1 and self.action == 'trim'
        """
        match = self.adapters.match_to(read.sequence)
        if match is not None:
            return match.trimmed(read), [match]
        else:
            return read, []


class ReverseComplementer(SingleEndModifier):
    """Trim adapters from a read and its reverse complement"""

    def __init__(self, adapter_cutter: AdapterCutter, rc_suffix: Optional[str] = " rc"):
        """
        rc_suffix -- suffix to add to the read name if sequence was reverse-complemented
        """
        self.adapter_cutter = adapter_cutter
        self.reverse_complemented = 0
        self._suffix = rc_suffix

    def __repr__(self):
        return f"ReverseComplementer(adapter_cutter={self.adapter_cutter})"

    def __call__(self, read: SequenceRecord, info: ModificationInfo):
        reverse_read = read.reverse_complement()

        forward_trimmed_read, forward_matches = self.adapter_cutter.match_and_trim(read)
        reverse_trimmed_read, reverse_matches = self.adapter_cutter.match_and_trim(
            reverse_read
        )

        forward_score = sum(m.score for m in forward_matches)
        reverse_score = sum(m.score for m in reverse_matches)
        use_reverse_complement = reverse_score > forward_score

        if use_reverse_complement:
            self.reverse_complemented += 1
            assert reverse_matches
            trimmed_read, matches = reverse_trimmed_read, reverse_matches
            info.is_rc = True
            if self._suffix:
                trimmed_read.name += self._suffix
        else:
            info.is_rc = False
            trimmed_read, matches = forward_trimmed_read, forward_matches

        if matches:
            self.adapter_cutter.with_adapters += 1
            for match in matches:
                stats = self.adapter_cutter.adapter_statistics[match.adapter]
                stats.add_match(match)
                stats.reverse_complemented += bool(use_reverse_complement)
            info.matches.extend(matches)  # TODO extend or overwrite?
        return trimmed_read


class PairedAdapterCutterError(Exception):
    pass


class PairedAdapterCutter(PairedEndModifier):
    """
    Trim adapters in pairs from R1 and R2.
    """

    def __init__(self, adapters1, adapters2, action="trim"):
        """
        adapters1 -- list of Adapters to be removed from R1
        adapters2 -- list of Adapters to be removed from R2

        Both lists must have the same, non-zero length.
         read pair is trimmed if adapters1[i] is found in R1 and adapters2[i] in R2.

        action -- What to do with a found adapter: None, 'trim', 'lowercase' or 'mask'
        """
        super().__init__()
        if len(adapters1) != len(adapters2):
            raise PairedAdapterCutterError(
                "The number of adapters to trim from R1 and R2 must be the same. "
                "Given: {} for R1, {} for R2".format(len(adapters1), len(adapters2))
            )
        if not adapters1:
            raise PairedAdapterCutterError("No adapters given")
        self._adapter_pairs = list(zip(adapters1, adapters2))
        logger.debug("Adapter pairs:")
        for a1, a2 in self._adapter_pairs:
            logger.debug(" • %s=%s -- %s=%s", a1.name, a1.spec(), a2.name, a2.spec())
        self.action = action
        self.with_adapters = 0
        self.adapter_statistics = [None, None]
        self.adapter_statistics[0] = {a: a.create_statistics() for a in adapters1}
        self.adapter_statistics[1] = {a: a.create_statistics() for a in adapters2}

    def __repr__(self):
        return f"PairedAdapterCutter(adapter_pairs={self._adapter_pairs!r})"

    def __call__(self, read1, read2, info1, info2):
        """ """
        best_matches = self._find_best_match_pair(read1.sequence, read2.sequence)
        if best_matches is None:
            return read1, read2
        match1, match2 = best_matches
        self.with_adapters += 1
        result = []
        for i, match, read in zip([0, 1], [match1, match2], [read1, read2]):
            trimmed_read = read
            if self.action == "lowercase":
                trimmed_read.sequence = trimmed_read.sequence.upper()

            trimmed_read = match.trimmed(trimmed_read)
            self.adapter_statistics[i][match.adapter].add_match(match)

            if self.action == "trim":
                # read is already trimmed, nothing to do
                pass
            elif self.action == "mask":
                trimmed_read = AdapterCutter.masked_read(read, [match])
            elif self.action == "lowercase":
                trimmed_read = AdapterCutter.lowercased_read(read, [match])
                assert len(trimmed_read.sequence) == len(read)
            elif self.action == "retain":
                trimmed_read = AdapterCutter.trim_but_retain_adapter(read, [match])
            elif self.action is None:  # --no-trim
                trimmed_read = read[:]
            result.append(trimmed_read)
        info1.matches.append(match1)
        info2.matches.append(match2)
        return result

    def _find_best_match_pair(
        self, sequence1: str, sequence2: str
    ) -> Optional[Tuple[Match, Match]]:
        best = None
        best_score = None
        best_errors = None
        for adapter1, adapter2 in self._adapter_pairs:
            match1 = adapter1.match_to(sequence1)
            if match1 is None:
                continue
            match2 = adapter2.match_to(sequence2)
            if match2 is None:
                continue
            total_score = match1.score + match2.score
            total_errors = match1.errors + match2.errors
            if (
                best is None
                or total_score > best_score
                or (total_score == best_score and total_errors < best_errors)
            ):
                best = match1, match2
                best_score = total_score
                best_errors = total_errors
        return best


class UnconditionalCutter(SingleEndModifier):
    """
    A modifier that unconditionally removes the first n or the last n bases from a read.

    If the length is positive, the bases are removed from the beginning of the read.
    If the length is negative, the bases are removed from the end of the read.
    """

    def __init__(self, length: int):
        self.length = length

    def __repr__(self):
        return f"UnconditionalCutter(length={self.length})"

    def __call__(self, read, info: ModificationInfo):
        if self.length > 0:
            info.cut_prefix = read.sequence[: self.length]
            return read[self.length :]
        elif self.length < 0:
            info.cut_suffix = read.sequence[self.length :]
            return read[: self.length]


class LengthTagModifier(SingleEndModifier):
    """
    Replace "length=..." strings in read names.
    """

    def __init__(self, length_tag):
        self.regex = re.compile(r"\b" + length_tag + r"[0-9]*\b")
        self.length_tag = length_tag

    def __repr__(self):
        return f"LengthTagModifier(length_tag='{self.length_tag}')"

    def __call__(self, read, info: ModificationInfo):
        read = read[:]
        if read.name.find(self.length_tag) >= 0:
            read.name = self.regex.sub(
                self.length_tag + str(len(read.sequence)), read.name
            )
        return read


class SuffixRemover(SingleEndModifier):
    """
    Remove a given suffix from read names.
    """

    def __init__(self, suffix):
        self.suffix = suffix

    def __repr__(self):
        return f"SuffixRemover('{self.suffix}')"

    def __call__(self, read, info: ModificationInfo):
        read = read[:]
        if read.name.endswith(self.suffix):
            read.name = read.name[: -len(self.suffix)]
        return read


class PrefixSuffixAdder(SingleEndModifier):
    """
    Add a suffix and a prefix to read names
    """

    def __init__(self, prefix, suffix):
        self.prefix = prefix
        self.suffix = suffix

    def __repr__(self):
        return f"PrefixSuffixAdder(prefix='{self.prefix}', suffix='{self.suffix}')"

    def __call__(self, read, info):
        read = read[:]
        adapter_name = info.matches[-1].adapter.name if info.matches else "no_adapter"
        read.name = (
            self.prefix.replace("{name}", adapter_name)
            + read.name
            + self.suffix.replace("{name}", adapter_name)
        )
        return read


class InvalidTemplate(Exception):
    pass


class Renamer(SingleEndModifier):
    """
    Rename reads using a template

    The template string can contain the following placeholders:

    - {header} -- full, unchanged header
    - {id} -- the part of the header before the first whitespace
    - {comment} -- the part of the header after the ID, excluding initial whitespace
    - {cut_prefix} -- prefix removed by UnconditionalCutter (with positive length argument)
    - {cut_suffix} -- suffix removed by UnconditionalCutter (with negative length argument)
    - {adapter_name} -- name of the *last* adapter match or no_adapter if there was none
    - {match_sequence} -- the sequence that matched the adapter (this includes possible errors)
          or an empty string if there was no match
    - {rc} -- the string 'rc' if the read was reverse complemented (with --revcomp) or '' otherwise
    """

    variables = {
        "header",
        "id",
        "comment",
        "cut_prefix",
        "cut_suffix",
        "adapter_name",
        "rc",
        "match_sequence",
    }

    def __init__(self, template: str):
        template = template.replace(r"\t", "\t")
        try:
            self._tokens = list(tokenize_braces(template))
        except TokenizeError as e:
            raise InvalidTemplate(f"Error in template '{template}': {e}")
        self.raise_if_invalid_variable(self._tokens, self.variables)
        self._template = template
        self._rename = self.compile_rename_function()

    def __repr__(self):
        return f"{self.__class__.__name__}('{self._template}')"

    def __reduce__(self):
        return Renamer, (self._template,)

    def compile_rename_function(self):
        """
        Create the function that computes a new name

        By creating the code dynamically, we can ensure that only those placeholder values are
        computed that are actually used in the template.
        """
        code = {
            "header": "read.name",
            "id": "id_",
            "comment": "comment",
            "cut_prefix": "info.cut_prefix if info.cut_prefix else ''",
            "cut_suffix": "info.cut_suffix if info.cut_suffix else ''",
            "adapter_name": "info.matches[-1].adapter.name if info.matches else 'no_adapter'",
            "rc": "'rc' if info.is_rc else ''",
            "match_sequence": "info.matches[-1].match_sequence() if info.matches else ''",
        }
        placeholders = set(
            token.value for token in self._tokens if isinstance(token, BraceToken)
        )
        lines = ["def rename(self, read, info):"]
        if "id" in placeholders or "header" in placeholders:
            lines.append("  id_, comment = self.parse_name(read.name)")
        lines.append("  return self._template.format(")
        for placeholder in placeholders:
            lines.append(f"    {placeholder}={code[placeholder]},")
        lines.append("  )")
        logger.debug("Generated code of rename function:\n%s", "\n".join(lines))
        namespace = dict()
        exec("\n".join(lines), namespace)
        return namespace["rename"]

    @staticmethod
    def raise_if_invalid_variable(tokens: List[Token], allowed: Set[str]) -> None:
        for token in tokens:
            if not isinstance(token, BraceToken):
                continue
            value = token.value
            if value not in allowed:
                raise InvalidTemplate(
                    f"Error in template: Variable '{value}' not recognized"
                )

    @staticmethod
    def parse_name(read_name: str) -> Tuple[str, str]:
        """Parse read header and return (id, comment) tuple"""
        fields = read_name.split(maxsplit=1)
        if len(fields) == 2:
            return (fields[0], fields[1])
        else:
            return (read_name, "")

    def __call__(self, read: SequenceRecord, info: ModificationInfo) -> SequenceRecord:
        read.name = self._rename(self, read, info)
        return read


class PairedEndRenamer(PairedEndModifier):
    """
    Rename paired-end reads using a template. The template is applied to both
    R1 and R2, and the same template variables as in the (single-end) renamer
    are allowed. However,
    these variables are evaluated separately for each read. For example, if `{comment}`
    is used, it gets replaced with the R1 comment in the R1 header, and with the R2
    comment in the R2 header.

    Additionally, all template variables except `id` can be used in the read-specific
    forms `{r1.variablename}` and `{r2.variablename}`. For example, `{r1.comment}`
    always gets replaced with the R1 comment, even in R2.
    """

    def __init__(self, template: str):
        try:
            self._tokens = list(tokenize_braces(template))
        except TokenizeError as e:
            raise InvalidTemplate(f"Error in template '{template}': {e}")
        Renamer.raise_if_invalid_variable(self._tokens, self._get_allowed_variables())
        self._template = template.replace(r"\t", "\t")

    @staticmethod
    def _get_allowed_variables() -> Set[str]:
        allowed = (Renamer.variables - {"rc"}) | {"rn"}
        for v in Renamer.variables - {"id", "rc"}:
            allowed.add("r1." + v)
            allowed.add("r2." + v)
        return allowed

    def __call__(
        self,
        read1: SequenceRecord,
        read2: SequenceRecord,
        info1: ModificationInfo,
        info2: ModificationInfo,
    ) -> Tuple[SequenceRecord, SequenceRecord]:
        if not record_names_match(read1.name, read2.name):
            id1 = Renamer.parse_name(read1.name)[0]
            id2 = Renamer.parse_name(read1.name)[1]
            raise ValueError(f"Input read IDs not identical: '{id1}' != '{id2}'")

        name1, name2 = self._rename(read1, read2, info1, info2)

        if not record_names_match(name1, name2):
            new_id1 = Renamer.parse_name(name1)[0]
            new_id2 = Renamer.parse_name(name2)[0]
            id1 = Renamer.parse_name(read1.name)[0]
            raise InvalidTemplate(
                "After renaming R1 and R2, their IDs are no longer identical: "
                f"'{new_id1}' != '{new_id2}'. Original read ID: '{id1}'. "
            )
        read1.name = name1
        read2.name = name2
        return read1, read2

    def _rename(
        self,
        read1: SequenceRecord,
        read2: SequenceRecord,
        info1: ModificationInfo,
        info2: ModificationInfo,
    ) -> Tuple[str, str]:
        id1, comment1 = Renamer.parse_name(read1.name)
        id2, comment2 = Renamer.parse_name(read2.name)
        header1 = read1.name
        header2 = read2.name

        d = []
        for id_, comment, header, info in (
            (id1, comment1, header1, info1),
            (id2, comment2, header2, info2),
        ):
            if info.matches:
                adapter_name = info.matches[-1].adapter.name
                match_sequence = info.matches[-1].match_sequence()
            else:
                adapter_name = "no_adapter"
                match_sequence = ""
            d.append(
                dict(
                    comment=comment,
                    header=header,
                    cut_prefix=info.cut_prefix if info.cut_prefix else "",
                    cut_suffix=info.cut_suffix if info.cut_suffix else "",
                    adapter_name=adapter_name,
                    match_sequence=match_sequence,
                )
            )
        name1 = self._template.format(
            id=id1,
            rn=1,
            **d[0],
            r1=SimpleNamespace(**d[0]),
            r2=SimpleNamespace(**d[1]),
        )
        name2 = self._template.format(
            id=id2,
            rn=2,
            **d[1],
            r1=SimpleNamespace(**d[0]),
            r2=SimpleNamespace(**d[1]),
        )
        return name1, name2


class ZeroCapper(SingleEndModifier):
    """
    Change negative quality values of a read to zero
    """

    def __init__(self, quality_base=33):
        self.quality_base = quality_base
        qb = quality_base
        self.zero_cap_trans = str.maketrans("".join(map(chr, range(qb))), chr(qb) * qb)

    def __repr__(self):
        return f"ZeroCapper(quality_base={self.quality_base})"

    def __call__(self, read, info: ModificationInfo):
        read = read[:]
        read.qualities = read.qualities.translate(self.zero_cap_trans)
        return read


class NextseqQualityTrimmer(SingleEndModifier):
    def __init__(self, cutoff: int, base: int = 33):
        self.cutoff = cutoff
        self.base = base
        self.trimmed_bases = 0

    def __repr__(self):
        return f"NextseqQualityTrimmer(cutoff={self.cutoff}, base={self.base})"

    def __call__(self, read, info: ModificationInfo):
        stop = nextseq_trim_index(read, self.cutoff, self.base)
        self.trimmed_bases += len(read) - stop
        return read[:stop]


class QualityTrimmer(SingleEndModifier):
    def __init__(self, cutoff_front: int, cutoff_back: int, base: int = 33):
        self.cutoff_front = cutoff_front
        self.cutoff_back = cutoff_back
        self.base = base
        self.trimmed_bases = 0

    def __repr__(self):
        return (
            f"QualityTrimmer(cutoff_front={self.cutoff_front}, "
            f"cutoff_back={self.cutoff_back}, base={self.base})"
        )

    def __call__(self, read, info: ModificationInfo):
        start, stop = quality_trim_index(
            read.qualities, self.cutoff_front, self.cutoff_back, self.base
        )
        self.trimmed_bases += len(read) - (stop - start)
        return read[start:stop]


class PolyATrimmer(SingleEndModifier):
    """Trim poly-A tails or poly-T heads"""

    def __init__(self, revcomp=False):
        self.trimmed_bases = defaultdict(int)
        self.revcomp = revcomp

    def __repr__(self):
        return "PolyATrimmer()"

    def __call__(self, record: SequenceRecord, info: ModificationInfo):
        if self.revcomp:
            index = poly_a_trim_index(record.sequence, revcomp=True)
            self.trimmed_bases[index] += 1
            return record[index:]
        else:
            index = poly_a_trim_index(record.sequence)
            self.trimmed_bases[len(record) - index] += 1
            return record[:index]


class Shortener(SingleEndModifier):
    """Unconditionally shorten a read to the given length

    If the length is positive, the bases are removed from the end of the read.
    If the length is negative, the bases are removed from the beginning of the read.
    """

    def __init__(self, length):
        self.length = length

    def __repr__(self):
        return f"Shortener(length={self.length})"

    def __call__(self, read, info: ModificationInfo):
        if self.length >= 0:
            return read[: self.length]
        else:
            return read[self.length :]


class NEndTrimmer(SingleEndModifier):
    """Trims Ns from the 3' and 5' end of reads"""

    def __init__(self):
        self.start_trim = re.compile(r"^N+")
        self.end_trim = re.compile(r"N+$")

    def __repr__(self):
        return "NEndTrimmer()"

    def __call__(self, read, info: ModificationInfo):
        sequence = read.sequence
        start_cut = self.start_trim.match(sequence)
        end_cut = self.end_trim.search(sequence)
        start_cut = start_cut.end() if start_cut else 0
        end_cut = end_cut.start() if end_cut else len(read)
        return read[start_cut:end_cut]
