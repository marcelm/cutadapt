# coding: utf-8
"""
Routines for printing a report.
"""
from __future__ import print_function, division, absolute_import

import sys
from collections import namedtuple
import textwrap
from .adapters import BACK, FRONT, PREFIX, SUFFIX, ANYWHERE
from .modifiers import QualityTrimmer
from .writers import (TooShortReadFilter, TooLongReadFilter,
	ProcessedReadWriter, Demultiplexer, NContentTrimmer)

Statistics = namedtuple('Statistics', 'total_bp n')

ADAPTER_TYPES = {
	BACK: "regular 3'",
	FRONT: "regular 5'",
	PREFIX: "anchored 5'",
	SUFFIX: "anchored 3'",
	ANYWHERE: "variable 5'/3'"
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


def qtrimmed(modifiers):
	"""
	Look for a QualityTrimmer in the given list of modifiers and return its
	trimmed_bases attribute. If not found, return None.
	"""
	for m in modifiers:
		if isinstance(m, QualityTrimmer):
			return m.trimmed_bases
	return None


def print_statistics(adapters_pair, paired, time, stats,
		modifiers, modifiers2, writers, file=None):
	"""Print summary to file"""
	old_stdout = sys.stdout
	if file is not None:
		sys.stdout = file
	n = stats.n
	if stats.n == 0:
		print("No reads processed! Either your input file is empty or you used the wrong -f/--format parameter.")
		sys.stdout = old_stdout
		return
	time = max(time, 0.1)
	print("Finished in {0:.2F} s ({1:.0F} us/read; {2:.2F} M reads/minute).".format(
		time, 1E6 * time / n, n / time * 60 / 1E6))

	too_short = None
	too_long = None
	written = None
	written_bp = None
	too_many_n = None
	for w in writers:
		if isinstance(w, TooShortReadFilter):
			too_short = w.too_short
		elif isinstance(w, TooLongReadFilter):
			too_long = w.too_long
		elif isinstance(w, NContentTrimmer):
			too_many_n = w.too_many_n
		elif isinstance(w, (ProcessedReadWriter, Demultiplexer)):
			written = w.written
			written_fraction = written / n
			written_bp = w.written_bp
	assert written is not None

	with_adapters = [0, 0]
	for i in (0, 1):
		for adapter in adapters_pair[i]:
			with_adapters[i] += sum(adapter.lengths_front.values())
			with_adapters[i] += sum(adapter.lengths_back.values())
	with_adapters_fraction = [ v / n for v in with_adapters ]
	total_bp = sum(stats.total_bp)
	report = "\n=== Summary ===\n\n"
	if paired:
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
	pairs_or_reads = "Pairs" if paired else "Reads"
	if too_short is not None:
		too_short_fraction = too_short / n
		report += "{pairs_or_reads} that were too short:       {too_short:13,d} ({too_short_fraction:.1%})\n"
	if too_long is not None:
		too_long_fraction = too_long / n
		report += "{pairs_or_reads} that were too long:        {too_long:13,d} ({too_long_fraction:.1%})\n"
	if too_many_n is not None:
		too_many_n_fraction = too_many_n / n
		report += "{pairs_or_reads} with too many N:           {too_many_n:13,d} ({too_many_n_fraction:.1%})\n"

	report += textwrap.dedent("""\
	{pairs_or_reads} written (passing filters): {written:13,d} ({written_fraction:.1%})

	Total basepairs processed: {total_bp:13,d} bp
	""")
	if paired:
		report += "  Read 1: {stats.total_bp[0]:13,d} bp\n"
		report += "  Read 2: {stats.total_bp[1]:13,d} bp\n"

	quality_trimmed_bp = [qtrimmed(modifiers), qtrimmed(modifiers2)]
	if quality_trimmed_bp[0] is not None or quality_trimmed_bp[1] is not None:
		if quality_trimmed_bp[0] is None:
			quality_trimmed_bp[0] = 0
		if quality_trimmed_bp[1] is None:
			quality_trimmed_bp[1] = 0
		quality_trimmed = sum(quality_trimmed_bp)
		quality_trimmed_fraction = quality_trimmed / total_bp
		report += "Quality-trimmed:           {quality_trimmed:13,d} bp ({quality_trimmed_fraction:.1%})\n"
		if paired:
			report += "  Read 1: {quality_trimmed_bp[0]:13,d} bp\n"
			report += "  Read 2: {quality_trimmed_bp[1]:13,d} bp\n"

	total_written_bp = sum(written_bp)
	total_written_bp_fraction = total_written_bp / total_bp if total_bp > 0 else 0.0
	report += "Total written (filtered):  {total_written_bp:13,d} bp ({total_written_bp_fraction:.1%})\n"
	if paired:
		report += "  Read 1: {written_bp[0]:13,d} bp\n"
		report += "  Read 2: {written_bp[1]:13,d} bp\n"
	try:
		report = report.format(**vars())
	except ValueError:
		# Python 2.6 does not support the comma format specifier (PEP 378)
		report = report.replace(",d}", "d}").format(**vars())
	print(report)

	warning = False
	for which_in_pair in (0, 1):
		for index, adapter in enumerate(adapters_pair[which_in_pair]):
			total_front = sum(adapter.lengths_front.values())
			total_back = sum(adapter.lengths_back.values())
			total = total_front + total_back
			where = adapter.where
			assert where == ANYWHERE or (where in (BACK, SUFFIX) and total_front == 0) or (where in (FRONT, PREFIX) and total_back == 0)

			name = str(adapter.name)
			if not adapter.name_is_generated:
				name = "'{0}'".format(name)
			if paired:
				extra = 'First read: ' if which_in_pair == 0 else 'Second read: '
			else:
				extra = ''

			print("=" * 3, extra + "Adapter", name, "=" * 3)
			print()
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
				print_histogram(adapter.lengths_front, len(adapter), n, adapter.max_error_rate, adapter.errors_front)
				print()
				print("Overview of removed sequences (3' or within)")
				print_histogram(adapter.lengths_back, len(adapter), n, adapter.max_error_rate, adapter.errors_back)
			elif where in (FRONT, PREFIX):
				print()
				print_error_ranges(len(adapter), adapter.max_error_rate)
				print("Overview of removed sequences")
				print_histogram(adapter.lengths_front, len(adapter), n, adapter.max_error_rate, adapter.errors_front)
			else:
				assert where in (BACK, SUFFIX)
				print()
				print_error_ranges(len(adapter), adapter.max_error_rate)
				warning = warning or print_adjacent_bases(adapter.adjacent_bases, adapter.sequence)
				print("Overview of removed sequences")
				print_histogram(adapter.lengths_back, len(adapter), n, adapter.max_error_rate, adapter.errors_back)

	if warning:
		print('WARNING:')
		print('    One or more of your adapter sequences may be incomplete.')
		print('    Please see the detailed output above.')

	sys.stdout = old_stdout
