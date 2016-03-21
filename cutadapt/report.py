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
from .modifiers import QualityTrimmer, AdapterCutter
from .filters import (NoFilter, PairedNoFilter, TooShortReadFilter, TooLongReadFilter,
	DiscardTrimmedFilter, DiscardUntrimmedFilter, Demultiplexer, NContentFilter)


class Statistics:
	def __init__(self, n, total_bp1, total_bp2):
		"""
		n -- total number of reads
		total_bp1 -- number of bases in first reads
		total_bp2 -- number of bases in second reads (set to None for single-end data)
		"""
		self.n = n
		self.total_bp = total_bp1
		self.total_bp1 = total_bp1
		if total_bp2 is None:
			self.paired = False
		else:
			self.paired = True
			self.total_bp2 = total_bp2
			self.total_bp += total_bp2

	def collect(self, adapters_pair, time, modifiers, writers):
		self.time = max(time, 0.01)
		self.too_short = None
		self.too_long = None
		self.written = 0
		self.written_bp = [0, 0]
		self.too_many_n = None
		# Collect statistics from writers/filters
		for w in writers:
			if isinstance(w, (NoFilter, PairedNoFilter, Demultiplexer)) or isinstance(w.filter, (DiscardTrimmedFilter, DiscardUntrimmedFilter)):
				self.written += w.written
				if self.n > 0:
					self.written_fraction = self.written / self.n
				self.written_bp = self.written_bp[0] + w.written_bp[0], self.written_bp[1] + w.written_bp[1]
			elif isinstance(w.filter, TooShortReadFilter):
				self.too_short = w.filtered
			elif isinstance(w.filter, TooLongReadFilter):
				self.too_long = w.filtered
			elif isinstance(w.filter, NContentFilter):
				self.too_many_n = w.filtered
		assert self.written is not None

		# Collect statistics from modifiers
		self.with_adapters = [0, 0]
		self.quality_trimmed_bp = [0, 0]
		self.did_quality_trimming = False
		for i, mod_dict in [(0, modifiers.modifiers1), (1, modifiers.modifiers2)]:
			if QualityTrimmer in mod_dict:
				modifier = mod_dict[QualityTrimmer]
				self.quality_trimmed_bp[i] = modifier.trimmed_bases
				self.did_quality_trimming = True
			if AdapterCutter in mod_dict:
				self.with_adapters[i] += mod_dict[AdapterCutter].with_adapters
		self.with_adapters_fraction = [ (v / self.n if self.n > 0 else 0) for v in self.with_adapters ]
		self.quality_trimmed = sum(self.quality_trimmed_bp)
		self.quality_trimmed_fraction = self.quality_trimmed / self.total_bp if self.total_bp > 0 else 0.0

		self.total_written_bp = sum(self.written_bp)
		self.total_written_bp_fraction = self.total_written_bp / self.total_bp if self.total_bp > 0 else 0.0

		if self.n > 0:
			if self.too_short is not None:
				self.too_short_fraction = self.too_short / self.n
			if self.too_long is not None:
				self.too_long_fraction = self.too_long / self.n
			if self.too_many_n is not None:
				self.too_many_n_fraction = self.too_many_n / self.n


ADAPTER_TYPES = {
	BACK: "regular 3'",
	FRONT: "regular 5'",
	PREFIX: "anchored 5'",
	SUFFIX: "anchored 3'",
	ANYWHERE: "variable 5'/3'",
	LINKED: "linked",
}


def print_error_ranges(adapter_length, error_rate):
	print("No. of allowed errors:")
	prev = 0
	for errors in range(1, int(error_rate * adapter_length) + 1):
		r = int(errors / error_rate)
		print("{0}-{1} bp: {2};".format(prev, r - 1, errors - 1), end=' ')
		prev = r
	if prev == adapter_length:
		print("{0} bp: {1}".format(adapter_length, int(error_rate * adapter_length)))
	else:
		print("{0}-{1} bp: {2}".format(prev, adapter_length, int(error_rate * adapter_length)))
	print()


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

	print("length", "count", "expect", "max.err", "error counts", sep="\t")
	for length, count, estimate in h:
		max_errors = max(errors[length].keys())
		errs = ' '.join(str(errors[length][e]) for e in range(max_errors+1))
		print(length, count, "{0:.1F}".format(estimate), int(error_rate*min(length, adapter_length)), errs, sep="\t")
	print()


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
	print('Bases preceding removed adapters:')
	warnbase = None
	for base in ['A', 'C', 'G', 'T', '']:
		b = base if base != '' else 'none/other'
		fraction = 1.0 * bases[base] / total
		print('  {0}: {1:.1%}'.format(b, fraction))
		if fraction > 0.8 and base != '':
			warnbase = b
	if total >= 20 and warnbase is not None:
		print('WARNING:')
		print('    The adapter is preceded by "{0}" extremely often.'.format(warnbase))
		print('    The provided adapter sequence may be incomplete.')
		print('    To fix the problem, add "{0}" to the beginning of the adapter sequence.'.format(warnbase))
		print()
		return True
	print()
	return False


@contextmanager
def redirect_standard_output(file):
	if file is None:
		yield
		return
	old_stdout = sys.stdout
	sys.stdout = file
	yield
	sys.stdout = old_stdout


def print_report(stats, adapters_pair):
	"""Print report to standard output."""
	if stats.n == 0:
		print("No reads processed! Either your input file is empty or you used the wrong -f/--format parameter.")
		return
	print("Finished in {0:.2F} s ({1:.0F} us/read; {2:.2F} M reads/minute).".format(
		stats.time, 1E6 * stats.time / stats.n, stats.n / stats.time * 60 / 1E6))

	report = "\n=== Summary ===\n\n"
	if stats.paired:
		report += textwrap.dedent("""\
		Total read pairs processed:      {n:13,d}
		  Read 1 with adapter:           {with_adapters[0]:13,d} ({with_adapters_fraction[0]:.1%})
		  Read 2 with adapter:           {with_adapters[1]:13,d} ({with_adapters_fraction[1]:.1%})
		""")
	else:
		report += textwrap.dedent("""\
		Total reads processed:           {n:13,d}
		Reads with adapters:             {with_adapters[0]:13,d} ({with_adapters_fraction[0]:.1%})
		""")
	if stats.too_short is not None:
		report += "{pairs_or_reads} that were too short:       {too_short:13,d} ({too_short_fraction:.1%})\n"
	if stats.too_long is not None:
		report += "{pairs_or_reads} that were too long:        {too_long:13,d} ({too_long_fraction:.1%})\n"
	if stats.too_many_n is not None:
		report += "{pairs_or_reads} with too many N:           {too_many_n:13,d} ({too_many_n_fraction:.1%})\n"

	report += textwrap.dedent("""\
	{pairs_or_reads} written (passing filters): {written:13,d} ({written_fraction:.1%})

	Total basepairs processed: {total_bp:13,d} bp
	""")
	if stats.paired:
		report += "  Read 1: {total_bp1:13,d} bp\n"
		report += "  Read 2: {total_bp2:13,d} bp\n"

	if stats.did_quality_trimming:
		report += "Quality-trimmed:           {quality_trimmed:13,d} bp ({quality_trimmed_fraction:.1%})\n"
		if stats.paired:
			report += "  Read 1: {quality_trimmed_bp[0]:13,d} bp\n"
			report += "  Read 2: {quality_trimmed_bp[1]:13,d} bp\n"

	report += "Total written (filtered):  {total_written_bp:13,d} bp ({total_written_bp_fraction:.1%})\n"
	if stats.paired:
		report += "  Read 1: {written_bp[0]:13,d} bp\n"
		report += "  Read 2: {written_bp[1]:13,d} bp\n"
	v = vars(stats)
	v['pairs_or_reads'] = "Pairs" if stats.paired else "Reads"
	try:
		report = report.format(**v)
	except ValueError:
		# Python 2.6 does not support the comma format specifier (PEP 378)
		report = report.replace(",d}", "d}").format(**v)
	print(report)

	warning = False
	for which_in_pair in (0, 1):
		for adapter in adapters_pair[which_in_pair]:
			total_front = sum(adapter.lengths_front.values())
			total_back = sum(adapter.lengths_back.values())
			total = total_front + total_back
			where = adapter.where
			assert where in (ANYWHERE, LINKED) or (where in (BACK, SUFFIX) and total_front == 0) or (where in (FRONT, PREFIX) and total_back == 0)

			if stats.paired:
				extra = 'First read: ' if which_in_pair == 0 else 'Second read: '
			else:
				extra = ''

			print("=" * 3, extra + "Adapter", adapter.name, "=" * 3)
			print()
			if where == LINKED:
				print("Sequence: {0}...{1}; Type: linked; Length: {2}+{3}; Trimmed: {4} times; Half matches: {5}".
					format(adapter.front_adapter.sequence,
						adapter.back_adapter.sequence,
						len(adapter.front_adapter.sequence),
						len(adapter.back_adapter.sequence),
						total_front, total_back))
			else:
				print("Sequence: {0}; Type: {1}; Length: {2}; Trimmed: {3} times.".
					format(adapter.sequence, ADAPTER_TYPES[adapter.where],
						len(adapter.sequence), total))
			if total == 0:
				print()
				continue
			if where == ANYWHERE:
				print(total_front, "times, it overlapped the 5' end of a read")
				print(total_back, "times, it overlapped the 3' end or was within the read")
				print()
				print_error_ranges(len(adapter), adapter.max_error_rate)
				print("Overview of removed sequences (5')")
				print_histogram(adapter.lengths_front, len(adapter), stats.n, adapter.max_error_rate, adapter.errors_front)
				print()
				print("Overview of removed sequences (3' or within)")
				print_histogram(adapter.lengths_back, len(adapter), stats.n, adapter.max_error_rate, adapter.errors_back)
			elif where == LINKED:
				print()
				print_error_ranges(len(adapter.front_adapter), adapter.front_adapter.max_error_rate)
				print_error_ranges(len(adapter.back_adapter), adapter.back_adapter.max_error_rate)
				print("Overview of removed sequences at 5' end")
				print_histogram(adapter.front_adapter.lengths_front,
					len(adapter.front_adapter), stats.n,
					adapter.front_adapter.max_error_rate,
					adapter.front_adapter.errors_front)
				print()
				print("Overview of removed sequences at 3' end")
				print_histogram(adapter.back_adapter.lengths_back,
					len(adapter.back_adapter), stats.n,
					adapter.back_adapter.max_error_rate, adapter.back_adapter.errors_back)
			elif where in (FRONT, PREFIX):
				print()
				print_error_ranges(len(adapter), adapter.max_error_rate)
				print("Overview of removed sequences")
				print_histogram(adapter.lengths_front, len(adapter), stats.n, adapter.max_error_rate, adapter.errors_front)
			else:
				assert where in (BACK, SUFFIX)
				print()
				print_error_ranges(len(adapter), adapter.max_error_rate)
				warning = warning or print_adjacent_bases(adapter.adjacent_bases, adapter.sequence)
				print("Overview of removed sequences")
				print_histogram(adapter.lengths_back, len(adapter), stats.n, adapter.max_error_rate, adapter.errors_back)

	if warning:
		print('WARNING:')
		print('    One or more of your adapter sequences may be incomplete.')
		print('    Please see the detailed output above.')
