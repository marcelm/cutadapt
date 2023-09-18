"""
Routines for printing a report.
"""
from dataclasses import dataclass
from io import StringIO
import textwrap
from collections import defaultdict, Counter
from typing import Any, Optional, List, Dict, Iterator, Tuple, Mapping
from .adapters import (
    EndStatistics,
    AdapterStatistics,
    FrontAdapter,
    BackAdapter,
    AnywhereAdapter,
    LinkedAdapter,
    SingleAdapter,
    LinkedAdapterStatistics,
    FrontAdapterStatistics,
    BackAdapterStatistics,
    AnywhereAdapterStatistics,
)
from .json import OneLine
from .modifiers import (
    QualityTrimmer,
    NextseqQualityTrimmer,
    AdapterCutter,
    PairedAdapterCutter,
    ReverseComplementer,
    PairedEndModifierWrapper,
    PolyATrimmer,
)
from .statistics import ReadLengthStatistics
from .steps import HasStatistics, HasFilterStatistics
from .utils import MICRO

FILTERS = {
    "too_short": "that were too short",
    "too_long": "that were too long",
    "too_many_n": "with too many N",
    "too_many_expected_errors": "with too many exp. errors",
    "casava_filtered": "failed CASAVA filter",
    "discard_trimmed": "discarded as trimmed",
    "discard_untrimmed": "discarded as untrimmed",
}


def safe_divide(numerator: Optional[int], denominator: int) -> float:
    if numerator is None or not denominator:
        return 0.0
    else:
        return numerator / denominator


def add_if_not_none(a: Optional[int], b: Optional[int]) -> Optional[int]:
    if a is None:
        return b
    if b is None:
        return a
    return a + b


class Statistics:
    def __init__(self) -> None:
        """ """
        self.paired: Optional[bool] = None
        # Map a filter name to the number of filtered reads/read pairs
        self.filtered: Dict[str, int] = defaultdict(int)
        self.reverse_complemented: Optional[int] = None
        self.n = 0
        self.total_bp = [0, 0]
        self.read_length_statistics = ReadLengthStatistics()
        self.with_adapters: List[Optional[int]] = [None, None]
        self.quality_trimmed_bp: List[Optional[int]] = [None, None]
        self.poly_a_trimmed_lengths: List[Optional[defaultdict[int, int]]] = [
            None,
            None,
        ]
        self.adapter_stats: List[List[AdapterStatistics]] = [[], []]
        self._collected: bool = False

    def __iadd__(self, other: Any):
        if not isinstance(other, Statistics):
            raise ValueError(f"Cannot add {other.__type__.__name__}")
        self.n += other.n
        self.read_length_statistics += other.read_length_statistics

        if self.paired is None:
            self.paired = other.paired
        elif self.paired != other.paired:
            raise ValueError("Incompatible Statistics: paired is not equal")

        self.reverse_complemented = add_if_not_none(
            self.reverse_complemented, other.reverse_complemented
        )

        for filter_name, count in other.filtered.items():
            self.filtered[filter_name] += count

        for i in (0, 1):
            self.total_bp[i] += other.total_bp[i]
            self.with_adapters[i] = add_if_not_none(
                self.with_adapters[i], other.with_adapters[i]
            )
            self.quality_trimmed_bp[i] = add_if_not_none(
                self.quality_trimmed_bp[i], other.quality_trimmed_bp[i]
            )
            if self.poly_a_trimmed_lengths[i] is None:
                self.poly_a_trimmed_lengths[i] = other.poly_a_trimmed_lengths[i]
            elif other.poly_a_trimmed_lengths[i] is not None:
                self.poly_a_trimmed_lengths[i] = defaultdict(
                    int,
                    Counter(self.poly_a_trimmed_lengths[i])
                    + Counter(other.poly_a_trimmed_lengths[i]),
                )

            if self.adapter_stats[i] and other.adapter_stats[i]:
                if len(self.adapter_stats[i]) != len(other.adapter_stats[i]):
                    raise ValueError(
                        "Incompatible Statistics objects (adapter_stats length)"
                    )
                for j in range(len(self.adapter_stats[i])):
                    self.adapter_stats[i][j] += other.adapter_stats[i][j]
            elif other.adapter_stats[i]:
                assert self.adapter_stats[i] == []
                self.adapter_stats[i] = other.adapter_stats[i]
        return self

    def collect(
        self, n: int, total_bp1: int, total_bp2: Optional[int], modifiers, steps
    ):
        """
        n -- total number of reads
        total_bp1 -- number of bases in first reads
        total_bp2 -- number of bases in second reads. None for single-end data.
        """
        if self._collected:
            raise ValueError("Cannot call Statistics.collect more than once")
        self.n = n
        self.total_bp[0] = total_bp1
        if total_bp2 is None:
            self.paired = False
        else:
            self.paired = True
            self.total_bp[1] = total_bp2

        for step in steps:
            self._collect_step(step)
        for modifier in modifiers:
            self._collect_modifier(modifier)
        self._collected = True

        # For chaining
        return self

    def _collect_step(self, step) -> None:
        if isinstance(step, HasStatistics):
            self.read_length_statistics += step.get_statistics()
        if isinstance(step, HasFilterStatistics):
            name = step.descriptive_identifier()
            self.filtered[name] = step.filtered()

    def _collect_modifier(self, m) -> None:
        if isinstance(m, PairedAdapterCutter):
            for i in 0, 1:
                self.with_adapters[i] = m.with_adapters
                self.adapter_stats[i] = list(m.adapter_statistics[i].values())
            return
        if isinstance(m, PairedEndModifierWrapper):
            modifiers_list = [(0, m._modifier1), (1, m._modifier2)]
        else:
            modifiers_list = [(0, m)]
        for i, modifier in modifiers_list:
            if isinstance(modifier, (QualityTrimmer, NextseqQualityTrimmer)):
                self.quality_trimmed_bp[i] = add_if_not_none(
                    self.quality_trimmed_bp[i], modifier.trimmed_bases
                )
            if isinstance(modifier, PolyATrimmer):
                self.poly_a_trimmed_lengths[i] = modifier.trimmed_bases
            elif isinstance(modifier, AdapterCutter):
                assert self.with_adapters[i] is None
                self.with_adapters[i] = modifier.with_adapters
                self.adapter_stats[i] = list(modifier.adapter_statistics.values())
            elif isinstance(modifier, ReverseComplementer):
                assert self.with_adapters[i] is None
                self.with_adapters[i] = modifier.adapter_cutter.with_adapters
                self.adapter_stats[i] = list(
                    modifier.adapter_cutter.adapter_statistics.values()
                )
                self.reverse_complemented = modifier.reverse_complemented

    def as_json(self, gc_content: float = 0.5, one_line: bool = False) -> Dict:
        """
        Return a dict representation suitable for dumping in JSON format

        To achieve a more compact representation, set one_line to True, which
        will wrap some items in a `cutadapt.json.OneLine` object, and use
        `cutadapt.json.dumps` instead of `json.dumps` to dump the dict.
        """
        filtered = {name: self.filtered.get(name) for name in FILTERS.keys()}
        filtered_total = sum(self.filtered.values())
        written_reads = self.read_length_statistics.written_reads()
        written_bp = self.read_length_statistics.written_bp()
        assert written_reads + filtered_total == self.n
        return {
            "read_counts": {  # pairs or reads
                "input": self.n,
                "filtered": filtered,
                "output": self.read_length_statistics.written_reads(),
                "reverse_complemented": self.reverse_complemented,
                "read1_with_adapter": self.with_adapters[0],
                "read2_with_adapter": self.with_adapters[1] if self.paired else None,
            },
            "basepair_counts": {
                "input": self.total,
                "input_read1": self.total_bp[0],
                "input_read2": self.total_bp[1] if self.paired else None,
                "quality_trimmed": self.quality_trimmed,
                "quality_trimmed_read1": self.quality_trimmed_bp[0],
                "quality_trimmed_read2": self.quality_trimmed_bp[1],
                "poly_a_trimmed": self.poly_a_trimmed,
                "poly_a_trimmed_read1": self.poly_a_trimmed_bp[0],
                "poly_a_trimmed_read2": self.poly_a_trimmed_bp[1],
                "output": self.total_written_bp,
                "output_read1": written_bp[0],
                "output_read2": written_bp[1] if self.paired else None,
            },
            "adapters_read1": [
                self._adapter_statistics_as_json(
                    astats, self.n, gc_content, one_line=one_line
                )
                for astats in self.adapter_stats[0]
            ],
            "adapters_read2": [
                self._adapter_statistics_as_json(
                    astats, self.n, gc_content, one_line=one_line
                )
                for astats in self.adapter_stats[1]
            ]
            if self.paired
            else None,
            "poly_a_trimmed_read1": self._poly_a_trimmed_as_json(
                self.poly_a_trimmed_lengths[0]
            ),
            "poly_a_trimmed_read2": self._poly_a_trimmed_as_json(
                self.poly_a_trimmed_lengths[1]
            ),
        }

    def _adapter_statistics_as_json(
        self,
        adapter_statistics: AdapterStatistics,
        n: int,
        gc_content: float,
        one_line: bool = False,
    ):
        adapter = adapter_statistics.adapter
        ends: List[Optional[Dict[str, Any]]] = []
        total_trimmed_reads = 0
        make_line = OneLine if one_line else list
        for end_statistics in adapter_statistics.end_statistics():
            if end_statistics is None:
                ends.append(None)
                continue
            total = sum(end_statistics.lengths.values())
            if end_statistics.allows_partial_matches:
                eranges = ErrorRanges(
                    length=end_statistics.effective_length,
                    error_rate=end_statistics.max_error_rate,
                ).lengths()
            else:
                eranges = None
            base_stats = AdjacentBaseStatistics(end_statistics.adjacent_bases)
            trimmed_lengths = [
                make_line(
                    {
                        "len": row.length,
                        "expect": round(row.expect, 1),
                        "counts": row.error_counts,
                    }
                )
                for row in histogram_rows(end_statistics, n, gc_content)
            ]
            ends.append(
                {
                    "type": end_statistics.adapter_type,
                    "sequence": end_statistics.sequence,
                    "error_rate": end_statistics.max_error_rate,
                    "indels": end_statistics.indels,
                    "error_lengths": make_line(eranges),
                    "matches": total,
                    "adjacent_bases": base_stats.as_json(),
                    "dominant_adjacent_base": base_stats.warnbase,
                    "trimmed_lengths": trimmed_lengths,
                }
            )
            total_trimmed_reads += total

        on_reverse_complement = (
            adapter_statistics.reverse_complemented
            if self.reverse_complemented
            else None
        )
        return {
            "name": adapter_statistics.name,
            "total_matches": total_trimmed_reads,
            "on_reverse_complement": on_reverse_complement,
            "linked": isinstance(adapter, LinkedAdapter),
            "five_prime_end": ends[0],
            "three_prime_end": ends[1],
        }

    @staticmethod
    def _poly_a_trimmed_as_json(poly_a):
        if poly_a is None:
            return None
        return [
            OneLine({"len": length, "count": poly_a[length]})
            for length in sorted(poly_a)
        ]

    @property
    def total(self) -> int:
        return sum(self.total_bp)

    @property
    def quality_trimmed(self) -> Optional[int]:
        return add_if_not_none(*self.quality_trimmed_bp)

    @property
    def poly_a_trimmed_bp(self) -> Tuple[Optional[int], Optional[int]]:
        def trimmed(i: int) -> Optional[int]:
            lengths = self.poly_a_trimmed_lengths[i]
            if lengths is None:
                return None
            return sum(length * count for length, count in lengths.items())

        return (trimmed(0), trimmed(1))

    @property
    def poly_a_trimmed(self) -> Optional[int]:
        return add_if_not_none(*self.poly_a_trimmed_bp)

    @property
    def total_written_bp(self) -> int:
        return sum(self.read_length_statistics.written_bp())

    @property
    def written(self) -> int:
        return self.read_length_statistics.written_reads()

    @property
    def written_fraction(self) -> float:
        return safe_divide(self.read_length_statistics.written_reads(), self.n)

    @property
    def with_adapters_fraction(self) -> List[float]:
        return [safe_divide(v, self.n) for v in self.with_adapters]

    @property
    def quality_trimmed_fraction(self) -> float:
        return safe_divide(self.quality_trimmed, self.total)

    @property
    def written_bp(self) -> Tuple[int, int]:
        return self.read_length_statistics.written_bp()

    @property
    def total_written_bp_fraction(self) -> float:
        return safe_divide(self.total_written_bp, self.total)

    @property
    def reverse_complemented_fraction(self) -> float:
        return safe_divide(self.reverse_complemented, self.n)

    def filtered_fraction(self, filter_name: str) -> float:
        return safe_divide(self.filtered.get(filter_name), self.n)

    @property
    def poly_a_trimmed_fraction(self) -> float:
        return safe_divide(self.poly_a_trimmed, self.total)


class ErrorRanges:
    """
    Representation of the lengths up to which a number of errors is allowed
    for partial adapter matches.

    >>> ErrorRanges(length=8, error_rate=0.1).lengths()
    [8]
    >>> ErrorRanges(length=19, error_rate=0.1).lengths()
    [9, 19]
    >>> ErrorRanges(length=20, error_rate=0.1).lengths()
    [9, 19, 20]
    >>> ErrorRanges(length=21, error_rate=0.1).lengths()
    [9, 19, 21]

    The entry at index i in the returned list is the length up to which
    i errors are allowed. For example, the list [9, 19, 23] describes that
    - 0 errors are allowed up to length 9
    - 1 error is allowed up to length 19
    - 2 errors are allowed up to length 23

    The last number in the list is always the length of the adapter sequence.
    """

    def __init__(self, length: int, error_rate: float):
        self.length = length
        self.error_rate = error_rate
        self._lengths = self._compute_lengths()

    def _compute_lengths(self) -> List[int]:
        lengths = [
            int(errors / self.error_rate) - 1
            for errors in range(1, int(self.error_rate * self.length) + 1)
        ]
        if not lengths or lengths[-1] < self.length:
            lengths.append(self.length)
        return lengths

    def __repr__(self):
        return (
            "ErrorRanges("
            f"length={self.length}, error_rate={self.error_rate}, _lengths={self._lengths})"
        )

    def __str__(self):
        """
        >>> str(ErrorRanges(length=8, error_rate=0.1))
        '1-8 bp: 0'
        >>> str(ErrorRanges(length=20, error_rate=0.1))
        '1-9 bp: 0; 10-19 bp: 1; 20 bp: 2'
        >>> str(ErrorRanges(length=23, error_rate=0.1))
        '1-9 bp: 0; 10-19 bp: 1; 20-23 bp: 2'
        """
        prev = 1
        s = ""
        for errors, r in enumerate(self._lengths[:-1]):
            s += f"{prev}-{r} bp: {errors}; "
            prev = r + 1
        if prev == self._lengths[-1]:
            s += f"{prev} bp: {len(self._lengths) - 1}"
        else:
            s += f"{prev}-{self._lengths[-1]} bp: {len(self._lengths) - 1}"
        return s

    def lengths(self):
        return self._lengths


def error_ranges(end_statistics: EndStatistics) -> str:
    length = end_statistics.effective_length
    error_rate = end_statistics.max_error_rate
    if end_statistics.allows_partial_matches:
        s = "\n" + str(ErrorRanges(length, error_rate))
    else:
        s = f" {int(error_rate * length)}"
    return "No. of allowed errors:" + s + "\n"


def histogram(end_statistics: EndStatistics, n: int, gc_content: float) -> str:
    """
    Return a formatted histogram. Include the no. of reads expected to be
    trimmed by chance (assuming a uniform distribution of nucleotides in the reads).

    adapter_statistics -- EndStatistics object
    adapter_length -- adapter length
    n -- total no. of reads.
    """
    sio = StringIO()

    print("length", "count", "expect", "max.err", "error counts", sep="\t", file=sio)
    for row in histogram_rows(end_statistics, n, gc_content):
        print(
            row.length,
            row.count,
            f"{row.expect:.1F}",
            row.max_err,
            " ".join(str(e) for e in row.error_counts),
            sep="\t",
            file=sio,
        )
    return sio.getvalue() + "\n"


@dataclass
class HistogramRow:
    """One row in the "trimmed lengths" histogram"""

    length: int
    count: int
    expect: float
    max_err: int
    error_counts: List[int]


def histogram_rows(
    end_statistics: EndStatistics,
    n: int,
    gc_content: float,
) -> Iterator[HistogramRow]:
    """
    Yield histogram rows

    Include the no. of reads expected to be
    trimmed by chance (assuming a uniform distribution of nucleotides in the reads).

    n -- total no. of reads.
    """
    d = end_statistics.lengths
    errors = end_statistics.errors

    match_probabilities = end_statistics.random_match_probabilities(
        gc_content=gc_content
    )
    for length in sorted(d):
        # when length surpasses adapter_length, the
        # probability does not increase anymore
        expect = n * match_probabilities[min(len(end_statistics.sequence), length)]
        count = d[length]
        max_errors = max(errors[length].keys())
        error_counts = [errors[length][e] for e in range(max_errors + 1)]
        row = HistogramRow(
            length=length,
            count=count,
            expect=expect,
            max_err=int(
                end_statistics.max_error_rate
                * min(length, end_statistics.effective_length)
            ),
            error_counts=error_counts,
        )
        yield row


class AdjacentBaseStatistics:
    def __init__(self, bases: Dict[str, int]):
        """ """
        self.bases: Dict[str, int] = bases
        self._warnbase: Optional[str] = None
        total = sum(self.bases.values())
        if total == 0:
            self._fractions = None
        else:
            self._fractions = []
            for base in ["A", "C", "G", "T", ""]:
                text = base if base != "" else "none/other"
                fraction = 1.0 * self.bases[base] / total
                self._fractions.append((text, 1.0 * self.bases[base] / total))
                if fraction > 0.8 and base != "":
                    self._warnbase = text
            if total < 20:
                self._warnbase = None

    def __repr__(self):
        return f"AdjacentBaseStatistics(bases={self.bases})"

    @property
    def should_warn(self) -> bool:
        return self._warnbase is not None

    @property
    def warnbase(self) -> Optional[str]:
        return self._warnbase

    def __str__(self) -> str:
        if not self._fractions:
            return ""
        sio = StringIO()
        print("Bases preceding removed adapters:", file=sio)
        for text, fraction in self._fractions:
            print(f"  {text}: {fraction:.1%}", file=sio)
        if self.should_warn:
            print("WARNING:", file=sio)
            print(
                f"    The adapter is preceded by '{self._warnbase}' extremely often.",
                file=sio,
            )
            print(
                "    The provided adapter sequence could be incomplete at its 5' end.",
                file=sio,
            )
            print("    Ignore this warning when trimming primers.", file=sio)
        return sio.getvalue()

    def as_json(self) -> Optional[Dict[str, int]]:
        if self._fractions:
            return {b: self.bases.get(b, 0) for b in ["A", "C", "G", "T", ""]}
        else:
            return None


def full_report(stats: Statistics, time: float, gc_content: float) -> str:  # noqa: C901
    """Print report to standard output."""
    if stats.n == 0:
        return "No reads processed!"
    if time == 0:
        time = 1e-6
    sio = StringIO()

    def print_s(*args, **kwargs):
        kwargs["file"] = sio
        print(*args, **kwargs)

    print_s(
        "Finished in {:.3F} s ({:.3F} {}s/read; {:.2F} M reads/minute).".format(
            time, 1e6 * time / stats.n, MICRO, stats.n / time * 60 / 1e6
        )
    )

    report = "\n=== Summary ===\n\n"
    if stats.paired:
        report += f"Total read pairs processed:      {stats.n:13,d}\n"
        for i in (0, 1):
            if stats.with_adapters[i] is not None:
                report += (
                    f"  Read {i+1} with adapter:           "
                    f"{stats.with_adapters[i]:13,d} ({stats.with_adapters_fraction[i]:.1%})\n"
                )
    else:
        report += f"Total reads processed:           {stats.n:13,d}\n"
        if stats.with_adapters[0] is not None:
            report += (
                f"Reads with adapters:             "
                f"{stats.with_adapters[0]:13,d} ({stats.with_adapters_fraction[0]:.1%})\n"
            )

    if stats.reverse_complemented is not None:
        report += (
            "Reverse-complemented:            "
            "{o.reverse_complemented:13,d} ({o.reverse_complemented_fraction:.1%})\n"
        )

    filter_report = format_filter_report(stats)
    if filter_report:
        report += "\n== Read fate breakdown ==\n"
        report += filter_report

    report += textwrap.dedent(
        """\
    {pairs_or_reads} written (passing filters): {o.written:13,d} ({o.written_fraction:.1%})

    Total basepairs processed: {o.total:13,d} bp
    """
    )
    if stats.paired:
        report += "  Read 1: {o.total_bp[0]:13,d} bp\n"
        report += "  Read 2: {o.total_bp[1]:13,d} bp\n"

    if stats.quality_trimmed is not None:
        report += (
            "Quality-trimmed:           "
            f"{stats.quality_trimmed:13,d} bp ({stats.quality_trimmed_fraction:.1%})\n"
        )
        if stats.paired:
            for i in (0, 1):
                if stats.quality_trimmed_bp[i] is not None:
                    report += f"  Read {i + 1}: {stats.quality_trimmed_bp[i]:13,d} bp\n"

    if stats.poly_a_trimmed is not None:
        report += (
            "Poly-A-trimmed:            "
            f"{stats.poly_a_trimmed:13,d} bp ({stats.poly_a_trimmed_fraction:.1%})\n"
        )
        if stats.paired:
            for i in (0, 1):
                if stats.poly_a_trimmed_bp[i] is not None:
                    report += f"  Read {i + 1}: {stats.poly_a_trimmed_bp[i]:13,d} bp\n"

    report += (
        "Total written (filtered):  "
        "{o.total_written_bp:13,d} bp ({o.total_written_bp_fraction:.1%})\n"
    )
    if stats.paired:
        report += "  Read 1: {o.written_bp[0]:13,d} bp\n"
        report += "  Read 2: {o.written_bp[1]:13,d} bp\n"
    pairs_or_reads = "Pairs" if stats.paired else "Reads"
    report = report.format(o=stats, pairs_or_reads=pairs_or_reads)
    print_s(report)

    warning = False
    for which_in_pair in (0, 1):
        for adapter_statistics in stats.adapter_stats[which_in_pair]:
            end_statistics = adapter_statistics.end_statistics()
            if end_statistics[0] is not None:
                total_front = sum(end_statistics[0].lengths.values())
            else:
                total_front = 0
            if end_statistics[1] is not None:
                total_back = sum(end_statistics[1].lengths.values())
            else:
                total_back = 0
            total = total_front + total_back
            reverse_complemented = adapter_statistics.reverse_complemented
            adapter = adapter_statistics.adapter
            if isinstance(adapter, BackAdapter):
                assert total_front == 0
            if isinstance(adapter, FrontAdapter):
                assert total_back == 0

            if stats.paired:
                extra = "First read: " if which_in_pair == 0 else "Second read: "
            else:
                extra = ""

            print_s("=" * 3, extra + "Adapter", adapter_statistics.name, "=" * 3)
            print_s()

            if isinstance(adapter_statistics, LinkedAdapterStatistics):
                print_s(
                    "Sequence: {}...{}; Type: linked; Length: {}+{}; "
                    "5' trimmed: {} times; 3' trimmed: {} times".format(
                        adapter_statistics.front.sequence,
                        adapter_statistics.back.sequence,
                        len(adapter_statistics.front.sequence),
                        len(adapter_statistics.back.sequence),
                        total_front,
                        total_back,
                    ),
                    end="",
                )
            else:
                assert isinstance(adapter, (SingleAdapter, AnywhereAdapter))
                print_s(
                    "Sequence: {}; Type: {}; Length: {}; Trimmed: {} times".format(
                        adapter.sequence,
                        adapter.description,
                        len(adapter.sequence),
                        total,
                    ),
                    end="",
                )
            if stats.reverse_complemented is not None:
                print_s(f"; Reverse-complemented: {reverse_complemented} times")
            else:
                print_s()
            if total == 0:
                print_s()
                continue
            if isinstance(adapter_statistics, AnywhereAdapterStatistics):
                assert isinstance(adapter, AnywhereAdapter)
                print_s(total_front, "times, it overlapped the 5' end of a read")
                print_s(
                    total_back, "times, it overlapped the 3' end or was within the read"
                )
                print_s()
                print_s("Minimum overlap:", adapter.min_overlap)
                print_s(error_ranges(adapter_statistics.front))
                print_s("Overview of removed sequences (5')")
                print_s(histogram(adapter_statistics.front, stats.n, gc_content))
                print_s()
                print_s("Overview of removed sequences (3' or within)")
                print_s(histogram(adapter_statistics.back, stats.n, gc_content))
            elif isinstance(adapter_statistics, LinkedAdapterStatistics):
                assert isinstance(adapter, LinkedAdapter)
                print_s()
                print_s(
                    f"Minimum overlap: "
                    f"{adapter.front_adapter.min_overlap}+{adapter.back_adapter.min_overlap}"
                )
                print_s(error_ranges(adapter_statistics.front))
                print_s(error_ranges(adapter_statistics.back))
                print_s("Overview of removed sequences at 5' end")
                print_s(histogram(adapter_statistics.front, stats.n, gc_content))
                print_s()
                print_s("Overview of removed sequences at 3' end")
                print_s(histogram(adapter_statistics.back, stats.n, gc_content))
            elif isinstance(adapter_statistics, FrontAdapterStatistics):
                assert isinstance(adapter, FrontAdapter)
                print_s()
                if adapter.allows_partial_matches:
                    print_s("Minimum overlap:", adapter.min_overlap)
                print_s(error_ranges(adapter_statistics.end))
                print_s("Overview of removed sequences")
                print_s(histogram(adapter_statistics.end, stats.n, gc_content))
            else:
                assert isinstance(adapter_statistics, BackAdapterStatistics)
                assert isinstance(adapter, BackAdapter)
                print_s()
                if adapter.allows_partial_matches:
                    print_s("Minimum overlap:", adapter.min_overlap)
                print_s(error_ranges(adapter_statistics.end))
                base_stats = AdjacentBaseStatistics(
                    adapter_statistics.end.adjacent_bases
                )
                warning = warning or base_stats.should_warn
                print_s(base_stats)
                print_s("Overview of removed sequences")
                print_s(histogram(adapter_statistics.end, stats.n, gc_content))

        poly_a = stats.poly_a_trimmed_lengths[which_in_pair]
        if poly_a is not None:
            print_s(poly_a_report(poly_a, which_in_pair if stats.paired else None))

    if warning:
        print_s("WARNING:")
        print_s("    One or more of your adapter sequences may be incomplete.")
        print_s("    Please see the detailed output above.")

    return sio.getvalue().rstrip()


def poly_a_report(poly_a: Mapping[int, int], which_in_pair: Optional[int]) -> str:
    sio = StringIO()
    if which_in_pair is None:
        title = "Poly-A"
    elif which_in_pair == 0:
        title = "R1 poly-A"
    else:
        assert which_in_pair == 1
        title = "R2 poly-A"

    print(f"=== {title} trimmed ===", file=sio)
    print(file=sio)
    print("length", "count", sep="\t", file=sio)
    for length in sorted(poly_a):
        count = poly_a[length]
        print(length, count, sep="\t", file=sio)

    return sio.getvalue() + "\n"


def format_filter_report(stats):
    report = ""
    for name, description in FILTERS.items():
        if name not in stats.filtered:
            continue
        value = stats.filtered[name]
        fraction = stats.filtered_fraction(name)
        line = (
            "{pairs_or_reads} "
            + (description + ":").ljust(27)
            + f"{value:13,d} ({fraction:.1%})\n"
        )
        report += line
    return report


def minimal_report(stats: Statistics, time: float, gc_content: float) -> str:
    """Create a minimal tabular report suitable for concatenation"""
    _ = time
    _ = gc_content

    fields = [
        "OK",
        stats.n,  # reads/pairs in
        stats.total,  # bases in
        stats.filtered.get("too_short", 0),  # reads/pairs
        stats.filtered.get("too_long", 0),  # reads/pairs
        stats.filtered.get("too_many_n", 0),  # reads/pairs
        stats.read_length_statistics.written_reads(),  # reads/pairs out
        stats.with_adapters[0] if stats.with_adapters[0] is not None else 0,  # reads
        stats.quality_trimmed_bp[0]
        if stats.quality_trimmed_bp[0] is not None
        else 0,  # bases
        stats.read_length_statistics.written_bp()[0],  # bases out
    ]
    if stats.paired:
        fields += [
            stats.with_adapters[1]
            if stats.with_adapters[1] is not None
            else 0,  # reads/pairs
            stats.quality_trimmed_bp[1]
            if stats.quality_trimmed_bp[1] is not None
            else 0,  # bases
            stats.read_length_statistics.written_bp()[1],  # bases
        ]

    warning = False
    for which_in_pair in (0, 1):
        for adapter_statistics in stats.adapter_stats[which_in_pair]:
            if isinstance(adapter_statistics, BackAdapterStatistics):
                if AdjacentBaseStatistics(
                    adapter_statistics.end.adjacent_bases
                ).should_warn:
                    warning = True
                    break
    if warning:
        fields[0] = "WARN"
    header = [
        "status",
        "in_reads",
        "in_bp",
        "too_short",
        "too_long",
        "too_many_n",
        "out_reads",
        "w/adapters",
        "qualtrim_bp",
        "out_bp",
    ]
    if stats.paired:
        header += ["w/adapters2", "qualtrim2_bp", "out2_bp"]
    return "\t".join(header) + "\n" + "\t".join(str(x) for x in fields)
