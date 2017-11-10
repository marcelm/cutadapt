# coding: utf-8
"""
Routines for printing a report.
"""
from __future__ import print_function, division, absolute_import

import sys
from contextlib import contextmanager
import textwrap
from .adapters import BACK, FRONT, PREFIX, SUFFIX, ANYWHERE, LINKED
from .modifiers import QualityTrimmer, AdapterCutter
from .filters import (NoFilter, PairedNoFilter, TooShortReadFilter, TooLongReadFilter,
	DiscardTrimmedFilter, DiscardUntrimmedFilter, PairedEndDemultiplexer, Demultiplexer,
	NContentFilter, InfoFileWriter, WildcardFileWriter, RestFileWriter)


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

	def collect(self, n, total_bp1, total_bp2, modifiers, modifiers2, writers):
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
			elif isinstance(w, (NoFilter, PairedNoFilter, PairedEndDemultiplexer, Demultiplexer)) or \
					isinstance(w.filter, (DiscardTrimmedFilter, DiscardUntrimmedFilter)):
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
		for i, modifiers_list in [(0, modifiers), (1, modifiers2)]:
			for modifier in modifiers_list:
				if isinstance(modifier, QualityTrimmer):
					self.quality_trimmed_bp[i] = modifier.trimmed_bases
					self.did_quality_trimming = True
				elif isinstance(modifier, AdapterCutter):
					self.with_adapters[i] += modifier.with_adapters
					self.adapter_stats[i] = list(modifier.adapter_statistics.values())

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


def print_histogram(end_statistics, n, gc_content):
	"""
	Print a histogram. Also, print the no. of reads expected to be
	trimmed by chance (assuming a uniform distribution of nucleotides in the reads).

	adapter_statistics -- EndStatistics object
	adapter_length -- adapter length
	n -- total no. of reads.
	"""
	d = end_statistics.lengths
	errors = end_statistics.errors

	match_probabilities = end_statistics.random_match_probabilities(gc_content=gc_content)
	print("length", "count", "expect", "max.err", "error counts", sep="\t")
	for length in sorted(d):
		# when length surpasses adapter_length, the
		# probability does not increase anymore
		expect = n * match_probabilities[min(len(end_statistics.sequence), length)]
		count = d[length]
		max_errors = max(errors[length].keys())
		errs = ' '.join(str(errors[length][e]) for e in range(max_errors+1))
		print(
			length,
			count,
			"{0:.1F}".format(expect),
			int(end_statistics.max_error_rate*min(length, len(end_statistics.sequence))),
			errs,
			sep="\t")
	print()


def print_adjacent_bases(bases):
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


def print_report(stats, time, gc_content):
	"""Print report to standard output."""
	if stats.n == 0:
		print("No reads processed! Either your input file is empty or you used the wrong -f/--format parameter.")
		return
	print("Finished in {0:.2F} s ({1:.0F} us/read; {2:.2F} M reads/minute).".format(
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
	print(report)

	warning = False
	for which_in_pair in (0, 1):
		for adapter_statistics in stats.adapter_stats[which_in_pair]:
			total_front = sum(adapter_statistics.front.lengths.values())
			total_back = sum(adapter_statistics.back.lengths.values())
			total = total_front + total_back
			where = adapter_statistics.where
			assert where in (ANYWHERE, LINKED) or (where in (BACK, SUFFIX) and total_front == 0) or (where in (FRONT, PREFIX) and total_back == 0)

			if stats.paired:
				extra = 'First read: ' if which_in_pair == 0 else 'Second read: '
			else:
				extra = ''

			print("=" * 3, extra + "Adapter", adapter_statistics.name, "=" * 3)
			print()

			if where == LINKED:
				print("Sequence: {0}...{1}; Type: linked; Length: {2}+{3}; "
					"5' trimmed: {4} times; 3' trimmed: {5} times".format(
						adapter_statistics.front.sequence,
						adapter_statistics.back.sequence,
						len(adapter_statistics.front.sequence),
						len(adapter_statistics.back.sequence),
						total_front, total_back))
			else:
				print("Sequence: {0}; Type: {1}; Length: {2}; Trimmed: {3} times.".
					format(adapter_statistics.front.sequence, ADAPTER_TYPES[adapter_statistics.where],
						len(adapter_statistics.front.sequence), total))
			if total == 0:
				print()
				continue
			if where == ANYWHERE:
				print(total_front, "times, it overlapped the 5' end of a read")
				print(total_back, "times, it overlapped the 3' end or was within the read")
				print()
				print_error_ranges(len(adapter_statistics.front.sequence), adapter_statistics.front.max_error_rate)
				print("Overview of removed sequences (5')")
				print_histogram(adapter_statistics.front, stats.n, gc_content)
				print()
				print("Overview of removed sequences (3' or within)")
				print_histogram(adapter_statistics.back, stats.n, gc_content)
			elif where == LINKED:
				print()
				print_error_ranges(len(adapter_statistics.front.sequence), adapter_statistics.front.max_error_rate)
				print_error_ranges(len(adapter_statistics.back.sequence), adapter_statistics.back.max_error_rate)
				print("Overview of removed sequences at 5' end")
				print_histogram(adapter_statistics.front, stats.n, gc_content)
				print()
				print("Overview of removed sequences at 3' end")
				print_histogram(adapter_statistics.back, stats.n, gc_content)
			elif where in (FRONT, PREFIX):
				print()
				print_error_ranges(len(adapter_statistics.front.sequence), adapter_statistics.front.max_error_rate)
				print("Overview of removed sequences")
				print_histogram(adapter_statistics.front, stats.n, gc_content)
			else:
				assert where in (BACK, SUFFIX)
				print()
				print_error_ranges(len(adapter_statistics.back.sequence), adapter_statistics.back.max_error_rate)
				warning = warning or print_adjacent_bases(adapter_statistics.back.adjacent_bases)
				print("Overview of removed sequences")
				print_histogram(adapter_statistics.back, stats.n, gc_content)

	if warning:
		print('WARNING:')
		print('    One or more of your adapter sequences may be incomplete.')
		print('    Please see the detailed output above.')
