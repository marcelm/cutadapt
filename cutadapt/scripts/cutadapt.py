#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kate: word-wrap off; remove-trailing-space on; replace-trailing-space-save on;
#
# Copyright (c) 2010-2012 Marcel Martin <marcel.martin@tu-dortmund.de>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""%prog [options] <FASTA/FASTQ FILE> [<QUALITY FILE>]

Reads a FASTA or FASTQ file, finds and removes adapters,
and writes the changed sequence to standard output.
When finished, statistics are printed to standard error.

Use a dash "-" as file name to read from standard input
(FASTA/FASTQ is autodetected).

If two file names are given, the first must be a .fasta or .csfasta
file and the second must be a .qual file. This is the file format
used by some 454 software and by the SOLiD sequencer.
If you have color space data, you still need to provide the -c option
to correctly deal with color space!

If the name of any input or output file ends with '.gz', it is
assumed to be gzip-compressed.

If you want to search for the reverse complement of an adapter, you must
provide an additional adapter sequence using another -a, -b or -g parameter.

If the input sequences are in color space, the adapter
can be given in either color space (as a string of digits 0, 1, 2, 3) or in
nucleotide space.

EXAMPLE

Assuming your sequencing data is available as a FASTQ file, use this
command line:
$ cutadapt -e ERROR-RATE -a ADAPTER-SEQUENCE input.fastq > output.fastq

See the README file for more help and examples."""

from __future__ import print_function, division

import sys
import re
import time
import errno
if sys.version_info[0] < 3:
	from string import maketrans
else:
	maketrans = bytes.maketrans
	xrange = range
from optparse import OptionParser, OptionGroup

from .. import seqio, __version__
from ..xopen import xopen
from ..qualtrim import quality_trim_index
from ..adapters import Adapter, ColorspaceAdapter, BACK, FRONT, PREFIX, ANYWHERE


class HelpfulOptionParser(OptionParser):
	"""An OptionParser that prints full help on errors."""
	def error(self, msg):
		self.print_help(sys.stderr)
		self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


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


def print_histogram(d, adapter_length, n, error_rate):
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

	print("length", "count", "expected", "max. errors", sep="\t")
	for length, count, estimate in h:
		print(length, count, "{0:.1F}".format(estimate), int(error_rate*min(length, adapter_length)), sep="\t")
	print()


class Statistics(object):
	"""Store statistics about reads and adapters"""

	def __init__(self, adapters):
		self._start_time = time.clock()
		self.time = None
		self.adapters = adapters

	def stop_clock(self):
		"""Stop the timer that was automatically started when the class was instantiated."""
		self.time = time.clock() - self._start_time

	def print_statistics(self, n, total_bp, quality_trimmed, reads_changed,
			error_rate, too_short, too_long, file=None):
		"""Print summary to file"""
		if self.time is None:
			self.stop_clock()
		old_stdout = sys.stdout
		if file is not None:
			sys.stdout = file
		print("cutadapt version", __version__)
		print("Command line parameters:", " ".join(sys.argv[1:]))
		print("Maximum error rate: {0:.2%}".format(error_rate))
		print("   No. of adapters:", len(self.adapters))
		print("   Processed reads: {0:12}".format(n))
		print("   Processed bases: {0:12} bp ({1:.1F} Mbp)".format(total_bp, total_bp/1E6))
		trimmed_bp = 0
		for adapter in self.adapters:
			for d in (adapter.lengths_front, adapter.lengths_back):
				trimmed_bp += sum( seqlen*count for (seqlen, count) in d.items() )

		if n > 0:
			print("     Trimmed reads: {0:12} ({1:.1%})".format(reads_changed, reads_changed / n))
			t = [ ("Quality-trimmed", quality_trimmed), ("  Trimmed bases", trimmed_bp)]
			if quality_trimmed < 0:
				del t[0]
			for what, bp in t:
				s = " ({0:.2%} of total)".format(float(bp)/total_bp) if total_bp > 0 else ''
				print("   {0}: {1:12} bp ({2:.1F} Mbp){3}".format(what, bp, trimmed_bp/1E6, s))
			print("   Too short reads: {0:12} ({1:.1%} of processed reads)".format(too_short, too_short / n))
			print("    Too long reads: {0:12} ({1:.1%} of processed reads)".format(too_long, too_long / n))
		print("        Total time: {0:9.2F} s".format(self.time))
		if n > 0:
			print("     Time per read: {0:9.2F} ms".format(1000. * self.time / n))
		print()
		for index, adapter in enumerate(self.adapters):
			total_front = sum(adapter.lengths_front.values())
			total_back = sum(adapter.lengths_back.values())
			total = total_front + total_back
			where = adapter.where
			assert where == ANYWHERE or (where == BACK and total_front == 0) or (where in (FRONT, PREFIX) and total_back == 0)

			print("=" * 3, "Adapter", index+1, "=" * 3)
			print()
			if adapter.name:
				name = "'{0}' ({1})".format(adapter.name, adapter.sequence)
			else:
				name = "'{0}'".format(adapter.sequence)
			print("Adapter {0}, length {1}, was trimmed {2} times.".format(name, len(adapter.sequence), total))
			if where == ANYWHERE:
				print(total_front, "times, it overlapped the 5' end of a read")
				print(total_back, "times, it overlapped the 3' end or was within the read")
				print()
				print_error_ranges(len(adapter), adapter.max_error_rate)
				print("Lengths of removed sequences (5')")
				print_histogram(adapter.lengths_front, len(adapter), n, adapter.max_error_rate)
				print()
				print("Lengths of removed sequences (3' or within)")
				print_histogram(adapter.lengths_back, len(adapter), n, adapter.max_error_rate)
			elif where in (FRONT, PREFIX):
				print()
				print_error_ranges(len(adapter), adapter.max_error_rate)
				print("Lengths of removed sequences")
				print_histogram(adapter.lengths_front, len(adapter), n, adapter.max_error_rate)
			else:
				assert where == BACK
				print()
				print_error_ranges(len(adapter), adapter.max_error_rate)
				print("Lengths of removed sequences")
				print_histogram(adapter.lengths_back, len(adapter), n, adapter.max_error_rate)

		if n == 0:
			print("No reads were read! Either your input file is empty or you used the wrong -f/--format parameter.")
		sys.stdout = old_stdout


# TODO make this a class and add a trim_primer parameter
def write_read(read, outfile, twoheaders=False):
	"""
	Write read in either FASTA or FASTQ format
	(depending on whether qualities is None or not) to outfile

	If twoheaders is True and the output is FASTQ, then the sequence name
	(description) is also written after the "+" character in the third line.
	"""
	initial = getattr(read, 'primer', '')
	if read.qualities is None:
		# FASTA
		print('>%s\n%s%s' % (read.name, initial, read.sequence), file=outfile)
	else:
		# FASTQ
		tmp = read.name if twoheaders else ''
		print('@%s\n%s%s\n+%s\n%s' % (read.name, initial, read.sequence, tmp, read.qualities), file=outfile)


def read_sequences(seqfilename, qualityfilename, colorspace, fileformat):
	"""
	Read sequences and (if available) quality information from either:
	* seqfilename in FASTA format (qualityfilename must be None)
	* seqfilename in FASTQ format (qualityfilename must be None)
	* seqfilename in .csfasta format and qualityfilename in .qual format
	  (SOLiD color space)

	Return a generator over tuples (description, sequence, qualities).
	qualities is None if no qualities are available.
	qualities are ASCII-encoded (chr(quality) + 33).
	"""
	#if ftype == 'FASTQ' and qualityfilename is not None:
		#raise ValueError("If a FASTQ file is given, no quality file can be provided.")

	if qualityfilename is not None:
		if colorspace:
		# read from .(CS)FASTA/.QUAL
			return seqio.ColorspaceFastaQualReader(seqfilename, qualityfilename)
		else:
			return seqio.FastaQualReader(seqfilename, qualityfilename)
	else:
		# read from FASTA or FASTQ
		return seqio.SequenceReader(seqfilename, colorspace, fileformat)


class ReadFilter(object):
	"""Filter reads according to length and according to whether any adapter matches."""

	def __init__(self, minimum_length, maximum_length, too_short_outfile, discard_trimmed, discard_untrimmed, trim_primer):
		self.minimum_length = minimum_length
		self.maximum_length = maximum_length
		self.too_short_outfile = too_short_outfile
		self.discard_trimmed = discard_trimmed
		self.discard_untrimmed = discard_untrimmed
		self.trim_primer = trim_primer
		self.too_long = 0
		self.too_short = 0

	def keep(self, read, trimmed):
		"""
		Return whether to keep the given read.
		"""
		if self.discard_trimmed and trimmed:
			return False
		if self.discard_untrimmed and not trimmed:
			return False
		if len(read.sequence) < self.minimum_length:
			self.too_short += 1
			if self.too_short_outfile is not None:
				if self.trim_primer: # TODO refactor
					read = read[1:]
					read.primer = ''
				write_read(read, self.too_short_outfile)
			return False
		if len(read.sequence) > self.maximum_length:
			self.too_long += 1
			return False
		return True


class LengthTagModifier:
	"""
	Replace "length=..." strings in read names.
	"""
	def __init__(self, length_tag):
		self.regex = re.compile(r"\b" + length_tag + r"[0-9]*\b")
		self.length_tag = length_tag

	def apply(self, read):
		read = read[:]
		if read.name.find(self.length_tag) >= 0:
			read.name = self.regex.sub(self.length_tag + str(len(read.sequence)), read.name)
		return read


class SuffixRemover:
	"""
	Remove a given suffix from read names.
	"""
	def __init__(self, suffix):
		self.suffix = suffix

	def apply(self, read):
		read = read[:]
		if read.name.endswith(self.suffix):
			read.name = read.name[:-len(self.suffix)]
		return read


class PrefixSuffixAdder:
	"""
	Add a suffix and a prefix to read names
	"""
	def __init__(self, prefix, suffix):
		self.prefix = prefix
		self.suffix = suffix

	def apply(self, read):
		read = read[:]
		read.name = self.prefix + read.name + self.suffix
		return read


class DoubleEncoder:
	"""
	Double-encode colorspace reads, using characters ACGTN to represent colors.
	"""
	def __init__(self):
		self.DOUBLE_ENCODE_TRANS = maketrans(b'0123.', b'ACGTN')

	def apply(self, read):
		read = read[:]
		read.sequence = read.sequence.translate(self.DOUBLE_ENCODE_TRANS)
		return read


class ZeroCapper:
	"""
	Change negative quality values of a read to zero
	"""
	def __init__(self, quality_base=33):
		qb = quality_base
		if sys.version_info[0] < 3:
			self.ZERO_CAP_TRANS = maketrans(''.join(map(chr, range(qb))), chr(qb) * qb)
		else:
			self.ZERO_CAP_TRANS = maketrans(bytes(range(qb)), bytes([qb] * qb))

	def apply(self, read):
		read = read[:]
		read.qualities = read.qualities.translate(self.ZERO_CAP_TRANS)
		return read


class RestFileWriter(object):
	def __init__(self, file):
		self.file = file

	def write(self, match):
		rest = match.rest()
		if len(rest) > 0:
			print(rest, match.read.name, file=self.file)


class RepeatedAdapterMatcher(object):
	"""
	Repeatedly find one of multiple adapters in reads.
	The number of times the search is repeated is specified by the
	times parameter.
	"""

	def __init__(self, adapters, times=1, wildcard_file=None, info_file=None):
		"""
		adapters -- list of Adapter objects
		"""
		self.adapters = adapters
		self.times = times
		self.info_file = info_file
		self.wildcard_file = wildcard_file
		self.reads_changed = 0

	def _best_match(self, read):
		"""
		Find the best matching adapter.

		read -- The read to which each adapter will be aligned

		Return an AdapterMatch instance or None if there are no matches.
		"""
		best = None
		for adapter in self.adapters:
			match = adapter.match(read)
			if match is None:
				continue

			# the no. of matches determines which adapter fits best
			if best is None or match.matches > best.matches:
				best = match
		return best

	def _write_info(self, match):
		"""write one line to the info file"""
		# TODO move to separate class
		if not self.info_file:
			return
		seq = match.read.sequence
		if match is None:
			print(match.read.name, -1, seq, sep='\t', file=self.info_file)
		else:
			print(match.read.name, match.errors, match.rstart, match.rstop, seq[0:match.rstart], seq[match.rstart:match.rstop], seq[match.rstop:], sep='\t', file=self.info_file)

	def find_match(self, read):
		"""
		Determine the adapter that best matches the given read.
		Since the best adapter is searched repeatedly, a list
		of AdapterMatch instances is returned, which
		need to be applied consecutively to the read.
		The list is empty if there are no adapter matches.

		The read will be converted to uppercase
		before it is compared to the adapter sequences.
		"""
		matches = []

		# try at most self.times times to remove an adapter
		for t in xrange(self.times):
			match = self._best_match(read)
			if match is None:
				# nothing found
				break
			self._write_info(match) # FIXME move to cut() or somewhere else
			assert match.length > 0
			assert match.errors / match.length <= match.adapter.max_error_rate
			assert match.length - match.errors > 0

			if self.wildcard_file: # FIXME move to cut() or somewhere else
				print(match.wildcards(), read.name, file=self.wildcard_file)

			matches.append(match)
			if t != self.times - 1:
				read = match.adapter.trimmed(match)
		return matches


	def cut(self, matches):
		"""
		Cut found adapters from a single read.

		matches -- a list of AdapterMatch instances
		"""
		# TODO move these lines out of here
		read = matches[0].read

		if __debug__:
			old_length = len(read.sequence)
		assert matches

		# The last match contains a copy of the read it was matched to.
		# No iteration is necessary.
		read = matches[-1].adapter.trimmed(matches[-1])

		# if an adapter was found, then the read should now be shorter
		assert len(read.sequence) < old_length
		self.reads_changed += 1 # TODO move to filter class

		return read


class QualityTrimmer:
	def __init__(self, cutoff, base):
		self.cutoff = cutoff
		self.base = base
		self.trimmed_bases = 0  # statistics

	def trimmed(self, read):
		index = quality_trim_index(read.qualities, self.cutoff, self.base)
		self.trimmed_bases += len(read.qualities) - index
		return read[:index]


def process_reads(reader, adapter_matcher, quality_trimmer, modifiers,
		readfilter, trimmed_outfile, untrimmed_outfile, rest_writer, trim_primer):
	"""
	Loop over reads, find adapters, trim reads, apply modifiers and
	output modified reads.
	Return a tuple (number_of_processed_reads, number_of_processed_basepairs)
	"""
	n = 0 # no. of processed reads
	total_bp = 0
	twoheaders = None
	for read in reader:
		n += 1
		total_bp += len(read.sequence)
		if quality_trimmer:
			read = quality_trimmer.trimmed(read)
		matches = adapter_matcher.find_match(read)
		if len(matches) > 0:
			read = adapter_matcher.cut(matches)
			trimmed = True
		else:
			trimmed = False
		if rest_writer:
			rest_writer.write(matches[-1])
		for modifier in modifiers:
			read = modifier.apply(read)
		if twoheaders is None:
			try:
				twoheaders = reader.twoheaders
			except AttributeError:
				twoheaders = False
		if not readfilter.keep(read, trimmed):
			continue
		if trim_primer:
			read = read[1:]
			read.primer = ''
		write_read(read, trimmed_outfile if trimmed else untrimmed_outfile, twoheaders)
	return (n, total_bp)


def main(cmdlineargs=None, trimmed_outfile=sys.stdout):
	"""
	Main function that evaluates command-line parameters and iterates
	over all reads.

	trimmed_outfile is the default output file to which trimmed reads
	are sent. It can be overriden by using the '-o' parameter.
	"""
	parser = HelpfulOptionParser(usage=__doc__, version=__version__)

	parser.add_option("-f", "--format", default=None,
		help="Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. "
		"Ignored when reading csfasta/qual files (default: auto-detect "
		"from file name extension).")

	group = OptionGroup(parser, "Options that influence how the adapters are found",
		description="Each of the following three parameters (-a, -b, -g) can be used "
			"multiple times and in any combination to search for an entire set of "
			"adapters of possibly different types. All of the "
			"given adapters will be searched for in each read, but only the best "
			"matching one will be trimmed (but see the --times option).")
	group.add_option("-a", "--adapter", action="append", metavar="ADAPTER", dest="adapters", default=[],
		help="Sequence of an adapter that was ligated to the 3' end. The adapter itself and anything that follows is trimmed.")
	group.add_option("-b", "--anywhere", action="append", metavar="ADAPTER", default=[],
		help="Sequence of an adapter that was ligated to the 5' or 3' end. If the adapter is found within the read or overlapping the 3' end of the read, the behavior is the same as for the -a option. If the adapter overlaps the 5' end (beginning of the read), the initial portion of the read matching the adapter is trimmed, but anything that follows is kept.")
	group.add_option("-g", "--front", action="append", metavar="ADAPTER", default=[],
		help="Sequence of an adapter that was ligated to the 5' end. If the "
		"adapter sequence starts with the character '^', the adapter is "
		"'anchored'. An anchored adapter must appear in its entirety at the "
		"5' end of the read (it is a prefix of the read). A non-anchored adapter may "
		"appear partially at the 5' end, or it may occur within the read. If it is "
		"found within a read, the sequence preceding the adapter is also trimmed. "
		"In all cases, the adapter itself is trimmed.")
	group.add_option("-e", "--error-rate", type=float, default=0.1,
		help="Maximum allowed error rate (no. of errors divided by the length of the matching region) (default: %default)")
	group.add_option("-n", "--times", type=int, metavar="COUNT", default=1,
		help="Try to remove adapters at most COUNT times. Useful when an adapter gets appended multiple times (default: %default).")
	group.add_option("-O", "--overlap", type=int, metavar="LENGTH", default=3,
		help="Minimum overlap length. If the overlap between the read and the adapter is shorter than LENGTH, the read is not modified."
			"This reduces the no. of bases trimmed purely due to short random adapter matches (default: %default).")
	group.add_option("--match-read-wildcards", action="store_true", default=False,
		help="Allow 'N's in the read as matches to the adapter (default: %default).")
	group.add_option("-N", "--no-match-adapter-wildcards", action="store_false",
		default=True, dest='match_adapter_wildcards',
		help="Do not treat 'N' in the adapter sequence as wildcards. This is needed when you want to search for literal 'N' characters.")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Options for filtering of processed reads")
	group.add_option("--discard-trimmed", "--discard", action='store_true', default=False,
		help="Discard reads that contain the adapter instead of trimming them. Also use -O in order to avoid throwing away too many randomly matching reads!")
	group.add_option("--discard-untrimmed", "--trimmed-only", action='store_true', default=False,
		help="Discard reads that do not contain the adapter.")
	group.add_option("-m", "--minimum-length", type=int, default=0, metavar="LENGTH",
		help="Discard trimmed reads that are shorter than LENGTH. Reads that are too short even before adapter removal are also discarded. In colorspace, an initial primer is not counted (default: 0).")
	group.add_option("-M", "--maximum-length", type=int, default=sys.maxsize, metavar="LENGTH",
		help="Discard trimmed reads that are longer than LENGTH. "
			"Reads that are too long even before adapter removal "
			"are also discarded. In colorspace, an initial primer "
			"is not counted (default: no limit).")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Options that influence what gets output to where")
	group.add_option("-o", "--output", default=None, metavar="FILE",
		help="Write the modified sequences to this file instead of standard output and send the summary report to standard output. "
		     "The format is FASTQ if qualities are available, FASTA otherwise. (default: standard output)")
	group.add_option("--info-file", metavar="FILE",
		help="Write information about each read and its adapter matches into FILE. "
			"Currently experimental: Expect the file format to change!")
	group.add_option("-r", "--rest-file", default=None, metavar="FILE",
		help="When the adapter matches in the middle of a read, write the rest (after the adapter) into a file. Use - for standard output.")
	group.add_option("--wildcard-file", default=None, metavar="FILE",
		help="When the adapter has wildcard bases ('N's) write adapter bases matching wildcard "
		     "positions to FILE. Use - for standard output. "
		     "When there are indels in the alignment, this may occasionally "
		     "not be quite accurate.")
	group.add_option("--too-short-output", default=None, metavar="FILE",
		help="Write reads that are too short (according to length specified by -m) to FILE. (default: discard reads)")
	group.add_option("--untrimmed-output", default=None, metavar="FILE",
		help="Write reads that do not contain the adapter to FILE, instead "
			"of writing them to the regular output file. (default: output "
			"to same file as trimmed)")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Additional modifications to the reads")
	group.add_option("-q", "--quality-cutoff", type=int, default=0, metavar="CUTOFF",
		help="Trim low-quality ends from reads before adapter removal. "
			"The algorithm is the same as the one used by BWA "
			"(Subtract CUTOFF from all qualities; "
			"compute partial sums from all indices to the end of the "
			"sequence; cut sequence at the index at which the sum "
			"is minimal) (default: %default)")
	group.add_option("--quality-base", type=int, default=33,
		help="Assume that quality values are encoded as ascii(quality + QUALITY_BASE). The default (33) is usually correct, "
		     "except for reads produced by some versions of the Illumina pipeline, where this should be set to 64. (default: %default)")
	group.add_option("-x", "--prefix", default='',
		help="Add this prefix to read names")
	group.add_option("-y", "--suffix", default='',
		help="Add this suffix to read names")
	group.add_option("--strip-suffix", action='append', default=[],
		help="Remove this suffix from read names if present. Can be given multiple times.")
	group.add_option("-c", "--colorspace", action='store_true', default=False,
		help="Colorspace mode: Also trim the color that is adjacent to the found adapter.")
	group.add_option("-d", "--double-encode", action='store_true', default=False,
		help="When in color space, double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).")
	group.add_option("-t", "--trim-primer", action='store_true', default=False,
		help="When in color space, trim primer base and the first color "
			"(which is the transition to the first nucleotide)")
	group.add_option("--strip-f3", action='store_true', default=False,
		help="For color space: Strip the _F3 suffix of read names")
	group.add_option("--maq", "--bwa", action='store_true', default=False,
		help="MAQ- and BWA-compatible color space output. This enables -c, -d, -t, --strip-f3, -y '/1' and -z.")
	group.add_option("--length-tag", default=None, metavar="TAG",
		help="Search for TAG followed by a decimal number in the name of the read "
			"(description/comment field of the FASTA or FASTQ file). Replace the "
			"decimal number with the correct length of the trimmed read. "
			"For example, use --length-tag 'length=' to correct fields "
			"like 'length=123'.")
	group.add_option("--zero-cap", "-z", action='store_true', default=False,
		help="Change negative quality values to zero (workaround to avoid segmentation faults in BWA)")
	parser.add_option_group(group)

	options, args = parser.parse_args(args=cmdlineargs)

	if len(args) == 0:
		parser.error("At least one parameter needed: name of a FASTA or FASTQ file.")
	elif len(args) > 2:
		parser.error("Too many parameters.")

	input_filename = args[0]
	quality_filename = None
	if len(args) == 2:
		quality_filename = args[1]
	if input_filename.endswith('.qual') and quality_filename.endswith('fasta'):
		parser.error("FASTA and QUAL file given, but the FASTA file must be first.")

	if options.format is not None and options.format.lower() not in ['fasta', 'fastq', 'sra-fastq']:
		parser.error("The input file format must be either 'fasta', 'fastq' or 'sra-fastq' (not '{0}').".format(options.format))

	# TODO should this really be an error?
	if options.format is not None and quality_filename is not None:
		parser.error("If a pair of .fasta and .qual files is given, the -f/--format parameter cannot be used.")

	# default output files (overwritten below)
	too_short_outfile = None # too short reads go here
	#too_long_outfile = None # too long reads go here

	if options.output is not None:
		trimmed_outfile = xopen(options.output, 'w')
	untrimmed_outfile = trimmed_outfile # reads without adapters go here
	if options.untrimmed_output is not None:
		untrimmed_outfile = xopen(options.untrimmed_output, 'w')
	if options.too_short_output is not None:
		too_short_outfile = xopen(options.too_short_output, 'w')
	#if options.too_long_output is not None:
		#too_long_outfile = xopen(options.too_long_output, 'w')

	if options.maq:
		options.colorspace = True
		options.double_encode = True
		options.trim_primer = True
		options.strip_suffix.append('_F3')
		options.suffix = "/1"
		options.zero_cap = True
	if options.trim_primer and not options.colorspace:
		parser.error("Trimming the primer makes only sense in color space.")
	if options.double_encode and not options.colorspace:
		parser.error("Double-encoding makes only sense in color space.")
	if options.anywhere and options.colorspace:
		parser.error("Using --anywhere with color space reads is currently not supported  (if you think this may be useful, contact the author).")
	if not (0 <= options.error_rate <= 1.):
		parser.error("The maximum error rate must be between 0 and 1.")
	if options.overlap < 1:
		parser.error("The overlap must be at least 1.")

	if options.rest_file is not None:
		options.rest_file = xopen(options.rest_file, 'w')
		rest_writer = RestFileWriter(options.rest_file)
	else:
		rest_writer = None
	if options.info_file is not None:
		options.info_file = xopen(options.info_file, 'w')
	if options.wildcard_file is not None:
		options.wildcard_file = xopen(options.wildcard_file, 'w')

	adapters = []

	ADAPTER_CLASS = ColorspaceAdapter if options.colorspace else Adapter
	def append_adapters(adapter_list, where):
		for seq in adapter_list:
			fields = seq.split('=', 1)
			if len(fields) > 1:
				name, seq = fields
				name = name.strip()
			else:
				name = None
			seq = seq.strip()
			w = where
			if w == FRONT and seq.startswith('^'):
				seq = seq[1:]
				w = PREFIX
			if len(seq) == 0:
				parser.error("The adapter sequence is empty")
			adapter = ADAPTER_CLASS(seq, w, options.error_rate,
				options.overlap, options.match_read_wildcards,
				options.match_adapter_wildcards, name=name)
			adapters.append(adapter)

	append_adapters(options.adapters, BACK)
	append_adapters(options.anywhere, ANYWHERE)
	append_adapters(options.front, FRONT)


	# make sure these aren't used by accident
	del options.adapters
	del options.anywhere
	del options.front

	if not adapters and options.quality_cutoff == 0:
		print("You need to provide at least one adapter sequence.", file=sys.stderr)
		return 1

	modifiers = []
	if options.length_tag:
		modifiers.append(LengthTagModifier(options.length_tag))
	if options.strip_f3:
		options.strip_suffix.append('_F3')
	for suffix in options.strip_suffix:
		modifiers.append(SuffixRemover(suffix))
	if options.prefix or options.suffix:
		modifiers.append(PrefixSuffixAdder(options.prefix, options.suffix))
	if options.double_encode:
		modifiers.append(DoubleEncoder())
	if options.zero_cap:
		modifiers.append(ZeroCapper(quality_base=options.quality_base))
	if options.quality_cutoff > 0:
		quality_trimmer = QualityTrimmer(options.quality_cutoff, options.quality_base)
	else:
		quality_trimmer = None

	adapter_matcher = RepeatedAdapterMatcher(adapters, options.times,
				options.wildcard_file, options.info_file)
	readfilter = ReadFilter(options.minimum_length, options.maximum_length,
		too_short_outfile, options.discard_trimmed, options.discard_untrimmed,
		options.trim_primer)
	stats = Statistics(adapters)
	try:
		reader = read_sequences(input_filename, quality_filename, colorspace=options.colorspace, fileformat=options.format)
		(n, total_bp) = process_reads(reader, adapter_matcher, quality_trimmer, modifiers, readfilter, trimmed_outfile, untrimmed_outfile, rest_writer, options.trim_primer)
	except IOError as e:
		if e.errno == errno.EPIPE:
			return 1
		raise
	except seqio.FormatError as e:
		print("Error:", e, file=sys.stderr)
		return 1
	if options.rest_file is not None:
		options.rest_file.close()
	if options.wildcard_file is not None:
		options.wildcard_file.close()
	if options.info_file is not None:
		options.info_file.close()
	# send statistics to stderr if result was sent to stdout
	stat_file = sys.stderr if options.output is None else None

	total_quality_trimmed = quality_trimmer.trimmed_bases if quality_trimmer else -1
	stats.print_statistics(
		n, total_bp, total_quality_trimmed, adapter_matcher.reads_changed,
		options.error_rate, readfilter.too_short, readfilter.too_long, file=stat_file)

	return 0


if __name__ == '__main__':
	if len(sys.argv) > 1 and sys.argv[1] == '--profile':
		del sys.argv[1]
		import cProfile as profile
		profile.run('main()', 'cutadapt.prof')
	else:
		sys.exit(main())
