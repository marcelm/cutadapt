"""
Routines for printing a report.
"""
import sys
from io import StringIO
from contextlib import contextmanager
import textwrap
from .adapters import Where, EndStatistics
from .modifiers import QualityTrimmer, NextseqQualityTrimmer, AdapterCutter, PairedAdapterCutter
from .filters import (NoFilter, PairedNoFilter, TooShortReadFilter, TooLongReadFilter,
    PairedDemultiplexer, Demultiplexer, NContentFilter, InfoFileWriter, WildcardFileWriter,
    RestFileWriter)


def safe_divide(numerator, denominator):
    if numerator is None or not denominator:
        return 0.0
    else:
        return numerator / denominator


class Statistics:
    def __init__(self):
        """
        """
        self.paired = None
        self.too_short = None
        self.too_long = None
        self.too_many_n = None
        self.did_quality_trimming = None
        self.n = 0
        self.written = 0
        self.total_bp = [0, 0]
        self.written_bp = [0, 0]
        self.with_adapters = [0, 0]
        self.quality_trimmed_bp = [0, 0]
        self.adapter_stats = [[], []]

    def __iadd__(self, other):
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

        def add_if_not_none(a, b):
            if a is None:
                return b
            if b is None:
                return a
            return a + b
        self.too_short = add_if_not_none(self.too_short, other.too_short)
        self.too_long = add_if_not_none(self.too_long, other.too_long)
        self.too_many_n = add_if_not_none(self.too_many_n, other.too_many_n)
        for i in (0, 1):
            self.total_bp[i] += other.total_bp[i]
            self.written_bp[i] += other.written_bp[i]
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

    def collect(self, n, total_bp1, total_bp2, modifiers, writers):
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

        # Collect statistics from writers/filters
        for w in writers:
            if isinstance(w, (InfoFileWriter, RestFileWriter, WildcardFileWriter)):
                pass
            elif isinstance(w, (NoFilter, PairedNoFilter, PairedDemultiplexer, Demultiplexer)):
                self.written += w.written
                self.written_bp[0] += w.written_bp[0]
                self.written_bp[1] += w.written_bp[1]
            elif isinstance(w.filter, TooShortReadFilter):
                self.too_short = w.filtered
            elif isinstance(w.filter, TooLongReadFilter):
                self.too_long = w.filtered
            elif isinstance(w.filter, NContentFilter):
                self.too_many_n = w.filtered
        assert self.written is not None

        # Collect statistics from modifiers
        for m in modifiers:
            if isinstance(m, PairedAdapterCutter):
                for i in 0, 1:
                    self.with_adapters[i] += m.with_adapters
                    self.adapter_stats[i] = list(m.adapter_statistics[i].values())
                continue
            if getattr(m, 'paired', False):
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

        # For chaining
        return self

    @property
    def total(self):
        return sum(self.total_bp)

    @property
    def quality_trimmed(self):
        return sum(self.quality_trimmed_bp)

    @property
    def total_written_bp(self):
        return sum(self.written_bp)

    @property
    def written_fraction(self):
        return safe_divide(self.written, self.n)

    @property
    def with_adapters_fraction(self):
        return [safe_divide(v, self.n) for v in self.with_adapters]

    @property
    def quality_trimmed_fraction(self):
        return safe_divide(self.quality_trimmed, self.total)

    @property
    def total_written_bp_fraction(self):
        return safe_divide(self.total_written_bp, self.total)

    @property
    def too_short_fraction(self):
        return safe_divide(self.too_short, self.n)

    @property
    def too_long_fraction(self):
        return safe_divide(self.too_long, self.n)

    @property
    def too_many_n_fraction(self):
        return safe_divide(self.too_many_n, self.n)


ADAPTER_TYPES = {
    Where.BACK: "regular 3'",
    Where.BACK_NOT_INTERNAL: "non-internal 3'",
    Where.FRONT: "regular 5'",
    Where.FRONT_NOT_INTERNAL: "non-internal 5'",
    Where.PREFIX: "anchored 5'",
    Where.SUFFIX: "anchored 3'",
    Where.ANYWHERE: "variable 5'/3'",
    Where.LINKED: "linked",
}


def error_ranges(adapter_statistics: EndStatistics):
    length = adapter_statistics.effective_length
    error_rate = adapter_statistics.max_error_rate
    prev = 0
    s = ""
    for errors in range(1, int(error_rate * length) + 1):
        r = int(errors / error_rate)
        s += "{}-{} bp: {}; ".format(prev, r - 1, errors - 1)
        prev = r
    if prev == length:
        s += "{} bp: {}".format(length, int(error_rate * length))
    else:
        s += "{}-{} bp: {}".format(prev, length, int(error_rate * length))

    return "No. of allowed errors:\n" + s + "\n"


def histogram(end_statistics: EndStatistics, n: int, gc_content: float):
    """
    Return a formatted histogram. Include the no. of reads expected to be
    trimmed by chance (assuming a uniform distribution of nucleotides in the reads).

    adapter_statistics -- EndStatistics object
    adapter_length -- adapter length
    n -- total no. of reads.
    """
    sio = StringIO()
    d = end_statistics.lengths
    errors = end_statistics.errors

    match_probabilities = end_statistics.random_match_probabilities(gc_content=gc_content)
    print("length", "count", "expect", "max.err", "error counts", sep="\t", file=sio)
    for length in sorted(d):
        # when length surpasses adapter_length, the
        # probability does not increase anymore
        expect = n * match_probabilities[min(len(end_statistics.sequence), length)]
        count = d[length]
        max_errors = max(errors[length].keys())
        errs = ' '.join(str(errors[length][e]) for e in range(max_errors + 1))
        print(
            length,
            count,
            "{:.1F}".format(expect),
            int(end_statistics.max_error_rate * min(length, len(end_statistics.sequence))),
            errs,
            sep="\t",
            file=sio,
        )
    return sio.getvalue() + "\n"


class AdjacentBaseStatistics:
    def __init__(self, bases):
        """
        """
        self.bases = bases
        self._warnbase = None
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
    def should_warn(self):
        return self._warnbase is not None

    def __str__(self):
        if not self._fractions:
            return ""
        sio = StringIO()
        print('Bases preceding removed adapters:', file=sio)
        for text, fraction in self._fractions:
            print('  {}: {:.1%}'.format(text, fraction), file=sio)
        if self.should_warn:
            print('WARNING:', file=sio)
            print('    The adapter is preceded by "{}" extremely often.'.format(self._warnbase), file=sio)
            print("    The provided adapter sequence could be incomplete at its 3' end.", file=sio)
        return sio.getvalue()


def full_report(stats, time, gc_content) -> str:
    """Print report to standard output."""
    if stats.n == 0:
        return "No reads processed!"

    sio = StringIO()

    def print_s(*args, **kwargs):
        kwargs['file'] = sio
        print(*args, **kwargs)

    print_s("Finished in {:.2F} s ({:.0F} us/read; {:.2F} M reads/minute).".format(
        time, 1E6 * time / stats.n, stats.n / time * 60 / 1E6))

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
    if stats.too_short is not None:
        report += "{pairs_or_reads} that were too short:       {o.too_short:13,d} ({o.too_short_fraction:.1%})\n"
    if stats.too_long is not None:
        report += "{pairs_or_reads} that were too long:        {o.too_long:13,d} ({o.too_long_fraction:.1%})\n"
    if stats.too_many_n is not None:
        report += "{pairs_or_reads} with too many N:           {o.too_many_n:13,d} ({o.too_many_n_fraction:.1%})\n"

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
            where = adapter_statistics.where
            assert (where in (Where.ANYWHERE, Where.LINKED)
                or (where in (Where.BACK, Where.BACK_NOT_INTERNAL, Where.SUFFIX) and total_front == 0)
                or (where in (Where.FRONT, Where.FRONT_NOT_INTERNAL, Where.PREFIX) and total_back == 0)), (where, total_front)

            if stats.paired:
                extra = 'First read: ' if which_in_pair == 0 else 'Second read: '
            else:
                extra = ''

            print_s("=" * 3, extra + "Adapter", adapter_statistics.name, "=" * 3)
            print_s()

            if where is Where.LINKED:
                print_s("Sequence: {}...{}; Type: linked; Length: {}+{}; "
                    "5' trimmed: {} times; 3' trimmed: {} times".format(
                        adapter_statistics.front.sequence,
                        adapter_statistics.back.sequence,
                        len(adapter_statistics.front.sequence),
                        len(adapter_statistics.back.sequence),
                        total_front, total_back))
            else:
                print_s("Sequence: {}; Type: {}; Length: {}; Trimmed: {} times.".
                    format(adapter_statistics.front.sequence, ADAPTER_TYPES[adapter_statistics.where],
                        len(adapter_statistics.front.sequence), total))
            if total == 0:
                print_s()
                continue
            if where is Where.ANYWHERE:
                print_s(total_front, "times, it overlapped the 5' end of a read")
                print_s(total_back, "times, it overlapped the 3' end or was within the read")
                print_s()
                print_s(error_ranges(adapter_statistics.front))
                print_s("Overview of removed sequences (5')")
                print_s(histogram(adapter_statistics.front, stats.n, gc_content))
                print_s()
                print_s("Overview of removed sequences (3' or within)")
                print_s(histogram(adapter_statistics.back, stats.n, gc_content))
            elif where is Where.LINKED:
                print_s()
                print_s(error_ranges(adapter_statistics.front))
                print_s(error_ranges(adapter_statistics.back))
                print_s("Overview of removed sequences at 5' end")
                print_s(histogram(adapter_statistics.front, stats.n, gc_content))
                print_s()
                print_s("Overview of removed sequences at 3' end")
                print_s(histogram(adapter_statistics.back, stats.n, gc_content))
            elif where in (Where.FRONT, Where.PREFIX, Where.FRONT_NOT_INTERNAL):
                print_s()
                print_s(error_ranges(adapter_statistics.front))
                print_s("Overview of removed sequences")
                print_s(histogram(adapter_statistics.front, stats.n, gc_content))
            else:
                assert where in (Where.BACK, Where.SUFFIX, Where.BACK_NOT_INTERNAL)
                print_s()
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


def minimal_report(stats, time, gc_content) -> str:
    """Create a minimal tabular report suitable for concatenation"""

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
            if adapter_statistics.where in (Where.BACK, Where.SUFFIX, Where.BACK_NOT_INTERNAL):
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
