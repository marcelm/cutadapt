# coding: utf-8
"""
Routines for printing a report.
"""
from __future__ import print_function, division, absolute_import

import sys
from collections import namedtuple
from contextlib import contextmanager
import textwrap
from .adapters import BACK, FRONT, PREFIX, SUFFIX, ANYWHERE, LINKED
from .modifiers import ModType
from .filters import *

# TODO: Would be good to have a more generic way for filters
# and modifiers to report summary information

def collect_process_statistics(modifiers, filters):
    """
    Collect statistics from filters, and modifiers
    """
    stats = {}
    
    stats["too_short"] = None
    if FilterType.TOO_SHORT in filters:
        stats["too_short"] = filters[FilterType.TOO_SHORT].filtered
    
    stats["too_long"] = None
    if FilterType.TOO_LONG in filters:
        stats["too_long"] = filters[FilterType.TOO_LONG].filtered
    
    stats["too_many_n"] = None
    if FilterType.N_CONTENT in filters:
        stats["too_many_n"] = filters[FilterType.N_CONTENT].filtered
    
    stats["with_adapters"] = [0, 0]
    stats["quality_trimmed_bp"] = [0, 0]
    stats["did_quality_trimming"] = False
    for i in (0, 1):
        if ModType.TRIM_QUAL in modifiers[i]:
            stats["quality_trimmed_bp"][i] = modifiers[i][ModType.TRIM_QUAL].trimmed_bases
            stats["did_quality_trimming"] = True
            
        if ModType.ADAPTER in modifiers[i]:
            stats["with_adapters"][i] += modifiers[i][ModType.ADAPTER].with_adapters
    
    stats["quality_trimmed"] = sum(stats["quality_trimmed_bp"])
    
    return stats

def collect_adapter_statistics(adapters):
    """
    TODO: push this method into Adapter?
    """
    stats = {}
    for adapter in adapters:
        name = adapter.name
        total_front = sum(adapter.lengths_front.values())
        total_back = sum(adapter.lengths_back.values())
        where = adapter.where
        assert (where in (ANYWHERE, LINKED) or 
            (where in (BACK, SUFFIX) and total_front == 0) or 
            (where in (FRONT, PREFIX) and total_back == 0)
        )
        
        stats[name] = dict(
            name=name,
            where=where,
            total_front=total_front,
            total_back=total_back,
            total=total_front + total_back
        )
            
        if where == LINKED:
            stats[name]["front_sequence"] = adapter.front_adapter.sequence
            stats[name]["back_sequence"] = adapter.back_adapter.sequence
            stats[name]["front_max_error_rate"] = adapter.front_adapter.max_error_rate
            stats[name]["back_max_error_rate"] = adapter.back_adapter.max_error_rate
            stats[name]["front_lengths_front"] = adapter.front_adapter.lengths_front
            stats[name]["front_errors_front"] = adapter.front_adapter.errors_front
            stats[name]["front_lengths_back"] = adapter.front_adapter.lengths_back
            stats[name]["front_errors_back"] = adapter.front_adapter.errors_back
            stats[name]["back_lengths_front"] = adapter.back_adapter.lengths_front
            stats[name]["back_errors_front"] = adapter.back_adapter.errors_front
            stats[name]["back_lengths_back"] = adapter.back_adapter.lengths_back
            stats[name]["back_errors_back"] = adapter.back_adapter.errors_back
        else:
            stats[name]["sequence"] = adapter.sequence
            stats[name]["max_error_rate"] = adapter.max_error_rate
            if where in (ANYWHERE, FRONT, PREFIX):
                stats[name]["lengths_front"] = adapter.lengths_front
                stats[name]["errors_front"] = adapter.errors_front
            if where in (ANYWHERE, BACK, SUFFIX):
                stats[name]["lengths_back"] = adapter.lengths_back
                stats[name]["errors_back"] = adapter.errors_back
            if where in (BACK, SUFFIX):
                stats[name]["adjacent_bases"] = adapter.adjacent_bases
    
    return stats

def collect_writer_statistics(N, total_bp1, total_bp2, writers):
    """
    n -- total number of reads
    total_bp1 -- number of bases in first reads
    total_bp2 -- number of bases in second reads (set to None for single-end data)
    """
    written, written_bp = writers.summary()
    assert written is not None
    
    return dict(
        N=N,
        total_bp1=total_bp1,
        total_bp2=total_bp2,
        total_bp=total_bp1 + total_bp2,
        written=written,
        written_bp=written_bp,
        total_written_bp=sum(written_bp)
    )

class Summary(object):
    def __init__(self, writer_stats, process_stats={}, adapter_stats=[{},{}]):
        self.writer_stats = writer_stats
        self.process_stats = process_stats
        self.adapter_stats = adapter_stats
    
    def add_process_stats(self, stats):
        self._merge_stats(self.process_stats, stats)
    
    def add_adapter_stats(self, stats):
        """
        stats = [ {name : adapter_stats}, {name : adapter_stats} ]
        """
        for i in (0,1):
            for name in set(self.adapter_stats[i].keys()) + set(stats[i].keys()):
                if name in self.adapter_stats[i]:
                    self._merge_stats(self.adapter_stats[i][name], stats[i][name])
                else:
                    self.adapter_stats[i][name] = stats[i][name]
    
    def _merge_stats(self, dest, src):
        for k in set(dest.keys()) + set(src.keys()):
            v1 = dest.get(k, None)
            v2 = src.get(k, None)
            if v1 is None:
                dest[k] = v2
            elif v2 is not None:
                if isinstance(v1, list):
                    assert isinstance(v2, list)
                    dest[k] = [v1[0]+v2[0], v1[1]+v2[1]]
                else:
                    dest[k] = v1 + v2
    
    def finish(self):
        stats = self.writer_stats.copy()
        stats.update(self.process_stats)
        
        stats["written_fraction"] = 0
        stats["too_short_fraction"] = 0
        stats["too_long_fraction"] = 0
        stats["too_many_n_fraction"] = 0
        stats["with_adapters_fraction"] = [0, 0]
        stats["quality_trimmed_fraction"] = 0.0
        stats["total_written_bp_fraction"] = 0.0
        
        N = stats["N"]
        if N > 0:
            stats["written_fraction"] = stats["written"] / N if stats["written"] else 0
            stats["too_short_fraction"] = stats["too_short"] / N if stats["too_short"] else 0
            stats["too_long_fraction"] = stats["too_long"] / N if stats["too_long"] else 0
            stats["too_many_n_fraction"] = stats["too_many_n"] / N if stats["too_many_n"] else 0
            stats["with_adapters_fraction"] = [ (v / N) for v in stats["with_adapters"] ]
        
        if stats["total_bp"] > 0:
            stats["quality_trimmed_fraction"] = (stats["quality_trimmed"] / stats["total_bp"]) if stats["quality_trimmed"] else 0
            stats["total_written_bp_fraction"] = (stats["total_written_bp"] / stats["total_bp"]) if stats["total_written_bp"] else 0
        
        stats["adapters"] = [
            self.adapter_stats[0].values(),
            self.adapter_stats[1].values()
        ]
        
        return stats

def print_report(paired, options, time, stats):
    outfile = options.report_file
    close = False
    
    if outfile is not None:
        outfile = open(outfile, "w")
        close = True
    elif not options.quiet:
        outfile = sys.stderr if options.output is None else sys.stdout
    else:
        return
    
    stats["paired"] = paired
    stats["time"] = max(time, 0.01)
    stats['pairs_or_reads'] = "Pairs" if paired else "Reads"
    
    try:
        # TODO: rewrite generate_report so that this call 
        # receives a string result and writes it to outfile
        # outfile.write(generate_report(paired, time, stats))
        generate_report(stats, outfile)

    finally:
        if close and outfile is not None:
            outfile.close()
    
        def get_values(self):
            values = var(self)
            
            return values

def generate_report(stats, outfile):
    def _print(*args, **kwargs): print(*args, file=outfile, **kwargs)
    
    ADAPTER_TYPES = {
        BACK: "regular 3'",
        FRONT: "regular 5'",
        PREFIX: "anchored 5'",
        SUFFIX: "anchored 3'",
        ANYWHERE: "variable 5'/3'",
        LINKED: "linked",
    }

    # TODO: this should all be re-written to use a template
    # engine (either simple string formatting or e.g. mako) 
    # and generate a string result.

    def print_error_ranges(adapter_length, error_rate):
        _print("No. of allowed errors:")
        prev = 0
        for errors in range(1, int(error_rate * adapter_length) + 1):
            r = int(errors / error_rate)
            _print("{0}-{1} bp: {2};".format(prev, r - 1, errors - 1), end=' ')
            prev = r
        if prev == adapter_length:
            _print("{0} bp: {1}".format(adapter_length, int(error_rate * adapter_length)))
        else:
            _print("{0}-{1} bp: {2}".format(prev, adapter_length, int(error_rate * adapter_length)))
        _print()

    def print_histogram(d, adapter_length, n, error_rate, errors):
        """
        Print a histogram. Also, print the no. of reads expected to be
        trimmed by chance (assuming a uniform distribution of nucleotides in the reads).
        d -- a dictionary mapping lengths of trimmed sequences to their respective frequency
        adapter_length -- adapter length
        n -- total no. of reads.
        """
        h = []
        for length in sorted(d):
            # when length surpasses adapter_length, the
            # probability does not increase anymore
            estimated = n * 0.25 ** min(length, adapter_length)
            h.append( (length, d[length], estimated) )

        _print("length", "count", "expect", "max.err", "error counts", sep="\t")
        for length, count, estimate in h:
            max_errors = max(errors[length].keys())
            errs = ' '.join(str(errors[length][e]) for e in range(max_errors+1))
            _print(length, count, "{0:.1F}".format(estimate), int(error_rate*min(length, adapter_length)), errs, sep="\t")
        _print()

    def print_adjacent_bases(bases, sequence):
        """
        Print a summary of the bases preceding removed adapter sequences.
        Print a warning if one of the bases is overrepresented and there are
        at least 20 preceding bases available.

        Return whether a warning was printed.
        """
        total = sum(bases.values())
        if total == 0:
            return False
        _print('Bases preceding removed adapters:')
        warnbase = None
        for base in ['A', 'C', 'G', 'T', '']:
            b = base if base != '' else 'none/other'
            fraction = 1.0 * bases[base] / total
            _print('  {0}: {1:.1%}'.format(b, fraction))
            if fraction > 0.8 and base != '':
                warnbase = b
        if total >= 20 and warnbase is not None:
            _print('WARNING:')
            _print('    The adapter is preceded by "{0}" extremely often.'.format(warnbase))
            _print('    The provided adapter sequence may be incomplete.')
            _print('    To fix the problem, add "{0}" to the beginning of the adapter sequence.'.format(warnbase))
            _print()
            return True
        _print()
        return False
    
    """Print report to standard output."""
    if stats["N"] == 0:
        _print("No reads processed! Either your input file is empty or you used the wrong -f/--format parameter.")
        return
    _print("Finished in {0:.2F} s ({1:.0F} us/read; {2:.2F} M reads/minute).".format(
        stats["time"], 1E6 * stats["time"] / stats["N"], stats["N"] / stats["time"] * 60 / 1E6))

    report = "\n=== Summary ===\n\n"
    if stats["paired"]:
        report += textwrap.dedent("""\
        Total read pairs processed:      {N:13,d}
          Read 1 with adapter:           {with_adapters[0]:13,d} ({with_adapters_fraction[0]:.1%})
          Read 2 with adapter:           {with_adapters[1]:13,d} ({with_adapters_fraction[1]:.1%})
        """)
    else:
        report += textwrap.dedent("""\
        Total reads processed:           {N:13,d}
        Reads with adapters:             {with_adapters[0]:13,d} ({with_adapters_fraction[0]:.1%})
        """)
    if stats["too_short"] is not None:
        report += "{pairs_or_reads} that were too short:       {too_short:13,d} ({too_short_fraction:.1%})\n"
    if stats["too_long"] is not None:
        report += "{pairs_or_reads} that were too long:        {too_long:13,d} ({too_long_fraction:.1%})\n"
    if stats["too_many_n"] is not None:
        report += "{pairs_or_reads} with too many N:           {too_many_n:13,d} ({too_many_n_fraction:.1%})\n"

    report += textwrap.dedent("""\
    {pairs_or_reads} written (passing filters): {written:13,d} ({written_fraction:.1%})

    Total basepairs processed: {total_bp:13,d} bp
    """)
    if stats["paired"]:
        report += "  Read 1: {total_bp1:13,d} bp\n"
        report += "  Read 2: {total_bp2:13,d} bp\n"

    if stats["did_quality_trimming"]:
        report += "Quality-trimmed:           {quality_trimmed:13,d} bp ({quality_trimmed_fraction:.1%})\n"
        if stats["paired"]:
            report += "  Read 1: {quality_trimmed_bp[0]:13,d} bp\n"
            report += "  Read 2: {quality_trimmed_bp[1]:13,d} bp\n"

    report += "Total written (filtered):  {total_written_bp:13,d} bp ({total_written_bp_fraction:.1%})\n"
    if stats["paired"]:
        report += "  Read 1: {written_bp[0]:13,d} bp\n"
        report += "  Read 2: {written_bp[1]:13,d} bp\n"
    
    try:
        report = report.format(**stats)
    except ValueError:
        # Python 2.6 does not support the comma format specifier (PEP 378)
        report = report.replace(",d}", "d}").format(**stats)
    _print(report)

    warning = False
    for which_in_pair in (0, 1):
        for adapter in stats["adapters"][which_in_pair]:
            if stats["paired"]:
                extra = 'First read: ' if which_in_pair == 0 else 'Second read: '
            else:
                extra = ''

            _print("=" * 3, extra + "Adapter", adapter["name"], "=" * 3)
            _print()
            if adapter["where"] == LINKED:
                _print("Sequence: {0}...{1}; Type: linked; Length: {2}+{3}; Trimmed: {4} times; Half matches: {5}".
                    format(
                        adapter["front_sequence"],
                        adapter["back_sequence"],
                        len(adapter["front_sequence"]),
                        len(adapter["back_sequence"]),
                        adapter["total_front"], adapter["total_back"]
                    ))
            else:
                _print("Sequence: {0}; Type: {1}; Length: {2}; Trimmed: {3} times.".
                    format(adapter["sequence"], ADAPTER_TYPES[adapter["where"]],
                        len(adapter["sequence"]), adapter["total"]))
    
            if adapter["total"] == 0:
                _print()
                continue
            if adapter["where"] == ANYWHERE:
                _print(adapter["total_front"], "times, it overlapped the 5' end of a read")
                _print(adapter["total_back"], "times, it overlapped the 3' end or was within the read")
                _print()
                print_error_ranges(len(adapter["sequence"]), adapter["max_error_rate"])
                _print("Overview of removed sequences (5')")
                print_histogram(adapter["lengths_front"], len(adapter["sequence"]), stats["N"],
                    adapter["max_error_rate"], adapter["errors_front"])
                _print()
                _print("Overview of removed sequences (3' or within)")
                print_histogram(adapter["errors_back"], len(adapter["sequence"]), stats["N"],
                    adapter["max_error_rate"], adapter["errors_back"])
            elif adapter["where"] == LINKED:
                _print()
                print_error_ranges(len(adapter["front_sequence"]), adapter["front_max_error_rate"])
                print_error_ranges(len(adapter["back_sequence"]), adapter["back_max_error_rate"])
                _print("Overview of removed sequences at 5' end")
                print_histogram(adapter["front_lengths_front"], len(adapter["front_sequence"]), 
                    stats["N"], adapter["front_max_error_rate"], adapter["front_errors_front"])
                _print()
                _print("Overview of removed sequences at 3' end")
                print_histogram(adapter["back_lengths_back"], len(adapter["back_sequence"]),
                    stats["N"], adapter["back_max_error_rate"], adapter["back_errors_back"])
            elif adapter["where"] in (FRONT, PREFIX):
                _print()
                print_error_ranges(len(adapter["sequence"]), adapter["max_error_rate"])
                _print("Overview of removed sequences")
                print_histogram(adapter["lengths_front"], len(adapter["sequence"]),
                    stats["N"], adapter["max_error_rate"], adapter["errors_front"])
            else:
                assert adapter["where"] in (BACK, SUFFIX)
                _print()
                print_error_ranges(len(adapter["sequence"]), adapter["max_error_rate"])
                warning = warning or print_adjacent_bases(adapter["adjacent_bases"], adapter["sequence"])
                _print("Overview of removed sequences")
                print_histogram(adapter["lengths_back"], len(adapter["sequence"]),
                    stats["N"], adapter["max_error_rate"], adapter["errors_back"])
    
    if warning:
        _print('WARNING:')
        _print('    One or more of your adapter sequences may be incomplete.')
        _print('    Please see the detailed output above.')
