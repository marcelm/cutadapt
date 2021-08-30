"""
Routines for printing a report.
"""
import sys
from io import StringIO
import textwrap
from collections import Counter
from typing import Any, Optional, List, Dict, Tuple, Iterator
from .adapters import (
    EndStatistics, AdapterStatistics, FrontAdapter, NonInternalFrontAdapter, PrefixAdapter,
    BackAdapter, NonInternalBackAdapter, SuffixAdapter, AnywhereAdapter, LinkedAdapter,
)
from .modifiers import (QualityTrimmer, NextseqQualityTrimmer,
    AdapterCutter, PairedAdapterCutter, ReverseComplementer, PairedEndModifierWrapper)
from .steps import SingleEndFinalStep, PairedEndFinalStep


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


def sum_with_none(values):
    """Compute sum ignoring None values"""
    return sum(0 if v is None else v for v in values)


class Statistics:
    def __init__(self) -> None:
        """
        """
        self.paired: Optional[bool] = None
        self.did_quality_trimming: Optional[bool] = None
        self.too_short: Optional[int] = None
        self.too_long: Optional[int] = None
        self.too_many_n: Optional[int] = None
        self.too_many_expected_errors: Optional[int] = None
        self.casava_filtered: Optional[int] = None
        self.discard_trimmed: Optional[int] = None
        self.discard_untrimmed: Optional[int] = None
        self.reverse_complemented: Optional[int] = None
        self.n = 0
        self.written = 0
        self.total_bp = [0, 0]
        self.written_bp = [0, 0]
        self.written_lengths: List[Counter] = [Counter(), Counter()]
        self.with_adapters = [0, 0]
        self.quality_trimmed_bp = [0, 0]
        self.adapter_stats: List[List[AdapterStatistics]] = [[], []]

    def __iadd__(self, other: Any):
        self.n += other.n
        self.written += other.written

        if self.paired is None:
            self.paired = other.paired
        elif self.paired != other.paired:
            raise ValueError('Incompatible Statistics: paired is not equal')
        if self.did_quality_trimming is None:
            self.did_quality_trimming = other.did_quality_trimming
        elif self.did_quality_trimming != other.did_quality_trimming:
            raise ValueError('Incompatible Statistics: did_quality_trimming is not equal')

        self.reverse_complemented = add_if_not_none(
            self.reverse_complemented, other.reverse_complemented)
        self.too_short = add_if_not_none(self.too_short, other.too_short)
        self.too_long = add_if_not_none(self.too_long, other.too_long)
        self.too_many_n = add_if_not_none(self.too_many_n, other.too_many_n)
        self.too_many_expected_errors = add_if_not_none(
            self.too_many_expected_errors, other.too_many_expected_errors)
        self.discard_trimmed = add_if_not_none(self.discard_trimmed, other.discard_trimmed)
        self.discard_untrimmed = add_if_not_none(self.discard_untrimmed, other.discard_untrimmed)
        self.casava_filtered = add_if_not_none(self.casava_filtered, other.casava_filtered)
        for i in (0, 1):
            self.total_bp[i] += other.total_bp[i]
            self.written_bp[i] += other.written_bp[i]
            self.written_lengths[i] += other.written_lengths[i]
            self.with_adapters[i] += other.with_adapters[i]
            self.quality_trimmed_bp[i] += other.quality_trimmed_bp[i]
            if self.adapter_stats[i] and other.adapter_stats[i]:
                if len(self.adapter_stats[i]) != len(other.adapter_stats[i]):
                    raise ValueError('Incompatible Statistics objects (adapter_stats length)')
                for j in range(len(self.adapter_stats[i])):
                    self.adapter_stats[i][j] += other.adapter_stats[i][j]
            elif other.adapter_stats[i]:
                assert self.adapter_stats[i] == []
                self.adapter_stats[i] = other.adapter_stats[i]
        return self

    def collect(self, n: int, total_bp1: int, total_bp2: Optional[int], modifiers, writers):
        """
        n -- total number of reads
        total_bp1 -- number of bases in first reads
        total_bp2 -- number of bases in second reads. None for single-end data.
        """
        self.n = n
        self.total_bp[0] = total_bp1
        if total_bp2 is None:
            self.paired = False
        else:
            self.paired = True
            self.total_bp[1] = total_bp2

        for writer in writers:
            self._collect_writer(writer)
        assert self.written is not None
        for modifier in modifiers:
            self._collect_modifier(modifier)

        # For chaining
        return self

    def _collect_writer(self, w) -> None:
        if isinstance(w, (PairedEndFinalStep, SingleEndFinalStep)):
            self.written += w.statistics.written_reads()
            written_bp = w.statistics.written_bp()
            written_lengths = w.statistics.written_lengths()
            for i in 0, 1:
                self.written_bp[i] += written_bp[i]
                self.written_lengths[i] += written_lengths[i]
        if hasattr(w, "filter") and hasattr(w.filter, "name"):
            filter_name = w.filter.name
            if filter_name in {
                "too_short", "too_long", "too_many_n", "too_many_expected_errors",
                "casava_filtered", "discard_trimmed", "discard_untrimmed",
            }:
                setattr(self, filter_name, w.filtered)

    def _collect_modifier(self, m) -> None:
        if isinstance(m, PairedAdapterCutter):
            for i in 0, 1:
                self.with_adapters[i] += m.with_adapters
                self.adapter_stats[i] = list(m.adapter_statistics[i].values())
            return
        if isinstance(m, PairedEndModifierWrapper):
            modifiers_list = [(0, m._modifier1), (1, m._modifier2)]
        else:
            modifiers_list = [(0, m)]
        for i, modifier in modifiers_list:
            if isinstance(modifier, (QualityTrimmer, NextseqQualityTrimmer)):
                self.quality_trimmed_bp[i] = modifier.trimmed_bases
                self.did_quality_trimming = True
            elif isinstance(modifier, AdapterCutter):
                self.with_adapters[i] += modifier.with_adapters
                self.adapter_stats[i] = list(modifier.adapter_statistics.values())
            elif isinstance(modifier, ReverseComplementer):
                self.with_adapters[i] += modifier.adapter_cutter.with_adapters
                self.adapter_stats[i] = list(modifier.adapter_cutter.adapter_statistics.values())
                self.reverse_complemented = modifier.reverse_complemented

    def as_json(self, gc_content: float) -> Dict:
        """Return a dict representation suitable for dumping in JSON format"""
        filtered = {
            "too_short": self.too_short,
            "too_long": self.too_long,
            "too_many_n": self.too_many_n,
            "too_many_expected_errors": self.too_many_expected_errors,
            "casava_filtered": self.casava_filtered,
            "discard_trimmed": self.discard_trimmed,
            "discard_untrimmed": self.discard_untrimmed,
        }
        filtered_total = sum_with_none(filtered.values())
        assert self.written + filtered_total == self.n
        return {
            "read_counts": {  # pairs or reads
                "input": self.n,
                "filtered": filtered,
                "output": self.written,
                "reverse_complemented": self.reverse_complemented,
                "read1_with_adapter": self.with_adapters[0],
                "read2_with_adapter": self.with_adapters[1] if self.paired else None,
            },
            "basepair_counts": {
                "input": self.total,
                "input_read1": self.total_bp[0],
                "input_read2": self.total_bp[1] if self.paired else None,
                "quality_trimmed": self.quality_trimmed if self.did_quality_trimming else None,
                "quality_trimmed_read1": self.quality_trimmed_bp[0] if self.did_quality_trimming else None,
                "quality_trimmed_read2": self.quality_trimmed_bp[1] if self.paired else None,
                "output": self.total_written_bp,
                "output_read1": self.written_bp[0],
                "output_read2": self.written_bp[1] if self.paired else None,
            },
            "adapter_statistics_read1": [
                self._adapter_statistics_as_json(astats, self.n, gc_content)
                for astats in self.adapter_stats[0]
            ],
            "adapter_statistics_read2": [
                self._adapter_statistics_as_json(astats, self.n, gc_content)
                for astats in self.adapter_stats[1]
            ] if self.paired else None,
        }

    def _adapter_statistics_as_json(
        self, adapter_statistics: AdapterStatistics, n: int, gc_content: float
    ):
        adapter = adapter_statistics.adapter
        ends = []
        for which_end, end_statistics in (
            ("five_prime", adapter_statistics.front),
            ("three_prime", adapter_statistics.back),
        ):
            total = sum(end_statistics.lengths.values())
            if end_statistics.allows_partial_matches:
                eranges = ErrorRanges(
                    length=end_statistics.effective_length, error_rate=end_statistics.max_error_rate
                ).lengths()
            else:
                eranges = None
            base_stats = AdjacentBaseStatistics(end_statistics.adjacent_bases)
            ends.append({
                "which_end": which_end,
                "error_rate": end_statistics.max_error_rate,
                "error_lengths": eranges,
                "trimmed_reads": total,
                # "histogram": list(histogram_rows(end_statistics, n, gc_content)),
                "adjacent_bases": base_stats.as_json(),
                "dominant_adjacent_base": base_stats.warnbase,
            })

        on_reverse_complement = adapter_statistics.reverse_complemented if self.reverse_complemented else None
        return {
            "name": adapter_statistics.name,
            "type": adapter.description,
            "specification": adapter.spec(),
            "on_reverse_complement": on_reverse_complement,
            "total_trimmed_reads": ends[0]["trimmed_reads"] + ends[1]["trimmed_reads"],
            "ends": ends,
        }

    @property
    def total(self) -> int:
        return sum(self.total_bp)

    @property
    def quality_trimmed(self) -> int:
        return sum(self.quality_trimmed_bp)

    @property
    def total_written_bp(self) -> int:
        return sum(self.written_bp)

    @property
    def written_fraction(self) -> float:
        return safe_divide(self.written, self.n)

    @property
    def with_adapters_fraction(self) -> List[float]:
        return [safe_divide(v, self.n) for v in self.with_adapters]

    @property
    def quality_trimmed_fraction(self) -> float:
        return safe_divide(self.quality_trimmed, self.total)

    @property
    def total_written_bp_fraction(self) -> float:
        return safe_divide(self.total_written_bp, self.total)

    @property
    def reverse_complemented_fraction(self) -> float:
        return safe_divide(self.reverse_complemented, self.n)

    @property
    def too_short_fraction(self) -> float:
        return safe_divide(self.too_short, self.n)

    @property
    def too_long_fraction(self) -> float:
        return safe_divide(self.too_long, self.n)

    @property
    def too_many_n_fraction(self) -> float:
        return safe_divide(self.too_many_n, self.n)

    @property
    def too_many_expected_errors_fraction(self) -> float:
        return safe_divide(self.too_many_expected_errors, self.n)

    @property
    def casava_filtered_fraction(self) -> float:
        return safe_divide(self.casava_filtered, self.n)

    @property
    def discard_trimmed_fraction(self) -> float:
        return safe_divide(self.discard_trimmed, self.n)

    @property
    def discard_untrimmed_fraction(self) -> float:
        return safe_divide(self.discard_untrimmed, self.n)


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
        length, count, expect, max_err, error_counts = row
        print(
            length,
            count,
            "{:.1F}".format(expect),
            max_err,
            " ".join(str(e) for e in error_counts),
            sep="\t",
            file=sio,
        )
    return sio.getvalue() + "\n"


def histogram_rows(
    end_statistics: EndStatistics, n: int, gc_content: float,
) -> Iterator[Tuple[int, int, float, int, List[int]]]:
    """
    Yield tuples (length, count, expect, max_err, error_counts)

    Include the no. of reads expected to be
    trimmed by chance (assuming a uniform distribution of nucleotides in the reads).

    adapter_statistics -- EndStatistics object
    adapter_length -- adapter length
    n -- total no. of reads.
    """
    d = end_statistics.lengths
    errors = end_statistics.errors

    match_probabilities = end_statistics.random_match_probabilities(gc_content=gc_content)
    for length in sorted(d):
        # when length surpasses adapter_length, the
        # probability does not increase anymore
        expect = n * match_probabilities[min(len(end_statistics.sequence), length)]
        count = d[length]
        max_errors = max(errors[length].keys())
        error_counts = [errors[length][e] for e in range(max_errors + 1)]
        t = (
            length,
            count,
            expect,
            int(end_statistics.max_error_rate * min(length, len(end_statistics.sequence))),
            error_counts,
        )
        yield t


class AdjacentBaseStatistics:
    def __init__(self, bases: Dict[str, int]):
        """
        """
        self.bases: Dict[str, int] = bases
        self._warnbase: Optional[str] = None
        total = sum(self.bases.values())
        if total == 0:
            self._fractions = None
        else:
            self._fractions = []
            for base in ['A', 'C', 'G', 'T', '']:
                text = base if base != '' else 'none/other'
                fraction = 1.0 * self.bases[base] / total
                self._fractions.append((text, 1.0 * self.bases[base] / total))
                if fraction > 0.8 and base != '':
                    self._warnbase = text
            if total < 20:
                self._warnbase = None

    def __repr__(self):
        return 'AdjacentBaseStatistics(bases={})'.format(self.bases)

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
        print('Bases preceding removed adapters:', file=sio)
        for text, fraction in self._fractions:
            print('  {}: {:.1%}'.format(text, fraction), file=sio)
        if self.should_warn:
            print('WARNING:', file=sio)
            print('    The adapter is preceded by "{}" extremely often.'.format(self._warnbase), file=sio)
            print("    The provided adapter sequence could be incomplete at its 5' end.", file=sio)
            print("    Ignore this warning when trimming primers.", file=sio)
        return sio.getvalue()

    def as_json(self) -> Dict[str, int]:
        return {b: self.bases.get(b, 0) for b in ["A", "C", "G", "T", ""]}


def full_report(stats: Statistics, time: float, gc_content: float) -> str:  # noqa: C901
    """Print report to standard output."""
    if stats.n == 0:
        return "No reads processed!"
    if time == 0:
        time = 1E-6
    sio = StringIO()

    def print_s(*args, **kwargs):
        kwargs['file'] = sio
        print(*args, **kwargs)

    if sys.version_info[:2] <= (3, 6):
        micro = "u"
    else:
        micro = "µ"
    print_s("Finished in {:.2F} s ({:.0F} {}s/read; {:.2F} M reads/minute).".format(
        time, 1E6 * time / stats.n, micro, stats.n / time * 60 / 1E6))

    report = "\n=== Summary ===\n\n"
    if stats.paired:
        report += textwrap.dedent("""\
        Total read pairs processed:      {o.n:13,d}
          Read 1 with adapter:           {o.with_adapters[0]:13,d} ({o.with_adapters_fraction[0]:.1%})
          Read 2 with adapter:           {o.with_adapters[1]:13,d} ({o.with_adapters_fraction[1]:.1%})
        """)
    else:
        report += textwrap.dedent("""\
        Total reads processed:           {o.n:13,d}
        Reads with adapters:             {o.with_adapters[0]:13,d} ({o.with_adapters_fraction[0]:.1%})
        """)
    if stats.reverse_complemented is not None:
        report += "Reverse-complemented:            " \
                  "{o.reverse_complemented:13,d} ({o.reverse_complemented_fraction:.1%})\n"

    filter_report = format_filter_report(stats)
    if filter_report:
        report += "\n== Read fate breakdown ==\n"
        report += filter_report

    report += textwrap.dedent("""\
    {pairs_or_reads} written (passing filters): {o.written:13,d} ({o.written_fraction:.1%})

    Total basepairs processed: {o.total:13,d} bp
    """)
    if stats.paired:
        report += "  Read 1: {o.total_bp[0]:13,d} bp\n"
        report += "  Read 2: {o.total_bp[1]:13,d} bp\n"

    if stats.did_quality_trimming:
        report += "Quality-trimmed:           {o.quality_trimmed:13,d} bp ({o.quality_trimmed_fraction:.1%})\n"
        if stats.paired:
            report += "  Read 1: {o.quality_trimmed_bp[0]:13,d} bp\n"
            report += "  Read 2: {o.quality_trimmed_bp[1]:13,d} bp\n"

    report += "Total written (filtered):  {o.total_written_bp:13,d} bp ({o.total_written_bp_fraction:.1%})\n"
    if stats.paired:
        report += "  Read 1: {o.written_bp[0]:13,d} bp\n"
        report += "  Read 2: {o.written_bp[1]:13,d} bp\n"
    pairs_or_reads = "Pairs" if stats.paired else "Reads"
    report = report.format(o=stats, pairs_or_reads=pairs_or_reads)
    print_s(report)

    warning = False
    for which_in_pair in (0, 1):
        for adapter_statistics in stats.adapter_stats[which_in_pair]:
            total_front = sum(adapter_statistics.front.lengths.values())
            total_back = sum(adapter_statistics.back.lengths.values())
            total = total_front + total_back
            reverse_complemented = adapter_statistics.reverse_complemented
            adapter = adapter_statistics.adapter
            if isinstance(adapter, (BackAdapter, NonInternalBackAdapter, SuffixAdapter)):
                assert total_front == 0
            if isinstance(adapter, (FrontAdapter, NonInternalFrontAdapter, PrefixAdapter)):
                assert total_back == 0

            if stats.paired:
                extra = 'First read: ' if which_in_pair == 0 else 'Second read: '
            else:
                extra = ''

            print_s("=" * 3, extra + "Adapter", adapter_statistics.name, "=" * 3)
            print_s()

            if isinstance(adapter, LinkedAdapter):
                print_s("Sequence: {}...{}; Type: linked; Length: {}+{}; "
                    "5' trimmed: {} times; 3' trimmed: {} times".format(
                        adapter_statistics.front.sequence,
                        adapter_statistics.back.sequence,
                        len(adapter_statistics.front.sequence),
                        len(adapter_statistics.back.sequence),
                        total_front, total_back), end="")
            else:
                print_s("Sequence: {}; Type: {}; Length: {}; Trimmed: {} times".
                    format(adapter_statistics.front.sequence, adapter.description,
                        len(adapter_statistics.front.sequence), total), end="")
            if stats.reverse_complemented is not None:
                print_s("; Reverse-complemented: {} times".format(reverse_complemented))
            else:
                print_s()
            if total == 0:
                print_s()
                continue
            if isinstance(adapter, AnywhereAdapter):
                print_s(total_front, "times, it overlapped the 5' end of a read")
                print_s(total_back, "times, it overlapped the 3' end or was within the read")
                print_s()
                print_s("Minimum overlap:", adapter.min_overlap)
                print_s(error_ranges(adapter_statistics.front))
                print_s("Overview of removed sequences (5')")
                print_s(histogram(adapter_statistics.front, stats.n, gc_content))
                print_s()
                print_s("Overview of removed sequences (3' or within)")
                print_s(histogram(adapter_statistics.back, stats.n, gc_content))
            elif isinstance(adapter, LinkedAdapter):
                print_s()
                print_s(f"Minimum overlap: "
                        f"{adapter.front_adapter.min_overlap}+{adapter.back_adapter.min_overlap}")
                print_s(error_ranges(adapter_statistics.front))
                print_s(error_ranges(adapter_statistics.back))
                print_s("Overview of removed sequences at 5' end")
                print_s(histogram(adapter_statistics.front, stats.n, gc_content))
                print_s()
                print_s("Overview of removed sequences at 3' end")
                print_s(histogram(adapter_statistics.back, stats.n, gc_content))
            elif isinstance(adapter, (FrontAdapter, NonInternalFrontAdapter, PrefixAdapter)):
                print_s()
                if not isinstance(adapter, PrefixAdapter):
                    print_s("Minimum overlap:", adapter.min_overlap)
                print_s(error_ranges(adapter_statistics.front))
                print_s("Overview of removed sequences")
                print_s(histogram(adapter_statistics.front, stats.n, gc_content))
            else:
                assert isinstance(adapter, (BackAdapter, NonInternalBackAdapter, SuffixAdapter))
                print_s()
                if not isinstance(adapter, SuffixAdapter):
                    print_s("Minimum overlap:", adapter.min_overlap)
                print_s(error_ranges(adapter_statistics.back))
                base_stats = AdjacentBaseStatistics(adapter_statistics.back.adjacent_bases)
                warning = warning or base_stats.should_warn
                print_s(base_stats)
                print_s("Overview of removed sequences")
                print_s(histogram(adapter_statistics.back, stats.n, gc_content))

    if warning:
        print_s('WARNING:')
        print_s('    One or more of your adapter sequences may be incomplete.')
        print_s('    Please see the detailed output above.')

    return sio.getvalue().rstrip()


def format_filter_report(stats):
    report = ""
    if stats.too_short is not None:
        report += "{pairs_or_reads} that were too short:       {o.too_short:13,d} ({o.too_short_fraction:.1%})\n"
    if stats.too_long is not None:
        report += "{pairs_or_reads} that were too long:        {o.too_long:13,d} ({o.too_long_fraction:.1%})\n"
    if stats.too_many_n is not None:
        report += "{pairs_or_reads} with too many N:           {o.too_many_n:13,d} ({o.too_many_n_fraction:.1%})\n"
    if stats.too_many_expected_errors is not None:
        report += "{pairs_or_reads} with too many exp. errors: " \
                  "{o.too_many_expected_errors:13,d} ({o.too_many_expected_errors_fraction:.1%})\n"
    if stats.casava_filtered is not None:
        report += "{pairs_or_reads} failed CASAVA filter:      " \
                  "{o.casava_filtered:13,d} ({o.casava_filtered_fraction:.1%})\n"
    if stats.discard_trimmed is not None:
        report += "{pairs_or_reads} discarded as trimmed:      " \
                  "{o.discard_trimmed:13,d} ({o.discard_trimmed_fraction:.1%})\n"
    if stats.discard_untrimmed is not None:
        report += "{pairs_or_reads} discarded as untrimmed:    " \
                  "{o.discard_untrimmed:13,d} ({o.discard_untrimmed_fraction:.1%})\n"
    return report


def minimal_report(stats: Statistics, time: float, gc_content: float) -> str:
    """Create a minimal tabular report suitable for concatenation"""
    _ = time
    _ = gc_content

    def none(value):
        return 0 if value is None else value

    fields = [
        "OK",
        stats.n,  # reads/pairs in
        stats.total,  # bases in
        none(stats.too_short),  # reads/pairs
        none(stats.too_long),  # reads/pairs
        none(stats.too_many_n),  # reads/pairs
        stats.written,  # reads/pairs out
        stats.with_adapters[0],  # reads
        stats.quality_trimmed_bp[0],  # bases
        stats.written_bp[0],  # bases out
    ]
    if stats.paired:
        fields += [
            stats.with_adapters[1],  # reads/pairs
            stats.quality_trimmed_bp[1],  # bases
            stats.written_bp[1],  # bases
        ]

    warning = False
    for which_in_pair in (0, 1):
        for adapter_statistics in stats.adapter_stats[which_in_pair]:
            if isinstance(adapter_statistics.adapter, (BackAdapter, NonInternalBackAdapter, SuffixAdapter)):
                if AdjacentBaseStatistics(adapter_statistics.back.adjacent_bases).should_warn:
                    warning = True
                    break
    if warning:
        fields[0] = "WARN"
    header = [
        'status', 'in_reads', 'in_bp', 'too_short', 'too_long', 'too_many_n', 'out_reads',
        'w/adapters', 'qualtrim_bp', 'out_bp']
    if stats.paired:
        header += ['w/adapters2', 'qualtrim2_bp', 'out2_bp']
    return "\t".join(header) + "\n" + "\t".join(str(x) for x in fields)
