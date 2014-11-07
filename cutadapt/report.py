"""
Routines for printing a report.
"""
from __future__ import print_function, division, absolute_import

import sys
import platform
from collections import namedtuple
from . import __version__
from .adapters import BACK, FRONT, PREFIX, SUFFIX, ANYWHERE

Statistics = namedtuple('Statistics', 'total_bp n quality_trimmed_bases')

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
		print('    The adapter is preceded by "{}" extremely often.'.format(warnbase))
		print('    The provided adapter sequence may be incomplete.')
		print('    To fix the problem, add "{}" to the beginning of the adapter sequence.'.format(warnbase))
		print()
		return True
	print()
	return False


def print_statistics(adapters, time, stats, action, reads_matched,
		error_rate, too_short, too_long, args, file=None):
	"""Print summary to file"""
	old_stdout = sys.stdout
	if file is not None:
		sys.stdout = file
	n = stats.n
	total_bp = stats.total_bp
	quality_trimmed = stats.quality_trimmed_bases
	print("You are running cutadapt", __version__, "with Python", platform.python_version())
	print("Command line parameters:", " ".join(args))
	print("Maximum error rate: {0:.2%}".format(error_rate))
	print("   No. of adapters:", len(adapters))
	print("   Processed reads: {0:12}".format(n))
	print("   Processed bases: {0:12} bp ({1:.1F} Mbp)".format(total_bp, total_bp / 1E6))
	trimmed_bp = 0
	for adapter in adapters:
		for d in (adapter.lengths_front, adapter.lengths_back):
			trimmed_bp += sum(seqlen * count for (seqlen, count) in d.items())

	if n > 0:
		operation = "Trimmed" if action == 'trim' else "Matched"
		print("     {0} reads: {1:12} ({2:.1%})".format(operation, reads_matched, reads_matched / n))
		t = [ ("Quality-trimmed", quality_trimmed), ("  Trimmed bases", trimmed_bp)]
		if quality_trimmed < 0:
			del t[0]
		for what, bp in t:
			s = " ({0:.2%} of total)".format(float(bp)/total_bp) if total_bp > 0 else ''
			print("   {0}: {1:12} bp ({2:.1F} Mbp){3}".format(what, bp, bp/1E6, s))
		print("   Too short reads: {0:12} ({1:.1%} of processed reads)".format(too_short, too_short / n))
		print("    Too long reads: {0:12} ({1:.1%} of processed reads)".format(too_long, too_long / n))
	print("        Total time: {0:9.2F} s".format(time))
	if n > 0:
		print("     Time per read: {0:10.3F} ms".format(1000. * time / n))
	print()

	warning = False
	for index, adapter in enumerate(adapters):
		total_front = sum(adapter.lengths_front.values())
		total_back = sum(adapter.lengths_back.values())
		total = total_front + total_back
		where = adapter.where
		assert where == ANYWHERE or (where in (BACK, SUFFIX) and total_front == 0) or (where in (FRONT, PREFIX) and total_back == 0)

		name = str(adapter.name)
		if not adapter.name_is_generated:
			name = "'{0}'".format(name)
		print("=" * 3, "Adapter", name, "=" * 3)
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
	if n == 0:
		print("No reads were read! Either your input file is empty or you used the wrong -f/--format parameter.")
	sys.stdout = old_stdout
