#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kate: word-wrap off; remove-trailing-spaces all;
#
# Copyright (c) 2010-2015 Marcel Martin <marcel.martin@scilifelab.se>
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

"""
cutadapt version %version
Copyright (C) 2010-2015 Marcel Martin <marcel.martin@scilifelab.se>

cutadapt removes adapter sequences from high-throughput sequencing reads.

Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. The reverse complement is *not* automatically
searched. All reads from input.fastq will be written to output.fastq with the
adapter sequence removed. Adapter matching is error-tolerant. Multiple adapter
sequences can be given (use further -a options), but only the best-matching
adapter will be removed.

Input may also be in FASTA format. Compressed input and output is supported and
auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for
standard input/output. Without the -o option, output is sent to standard output.

Some other available features are:
  * Various other adapter types (5' adapters, "mixed" 5'/3' adapters etc.)
  * Trimming a fixed number of bases
  * Quality trimming
  * Trimming colorspace reads
  * Filtering reads by various criteria

Use "cutadapt --help" to see all command-line options.
See http://cutadapt.readthedocs.org/ for full documentation.
"""

from __future__ import print_function, division, absolute_import

# Print a helpful error message if the extension modules cannot be imported.
from cutadapt import check_importability
check_importability()

import sys
import time
import errno
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
import logging
import platform
import textwrap

from cutadapt import seqio, __version__
from cutadapt.xopen import xopen
from cutadapt.adapters import (Adapter, ColorspaceAdapter, gather_adapters,
	BACK, FRONT, PREFIX, SUFFIX, ANYWHERE)
from cutadapt.modifiers import (LengthTagModifier, SuffixRemover, PrefixSuffixAdder,
	DoubleEncoder, ZeroCapper, PrimerTrimmer, QualityTrimmer, UnconditionalCutter,
	NEndTrimmer)
from cutadapt.filters import (TooShortReadFilter, TooLongReadFilter,
	Demultiplexer, NContentFilter, DiscardUntrimmedFilter, DiscardTrimmedFilter)
from cutadapt.report import Statistics, print_report, redirect_standard_output
from cutadapt.compat import next

logger = logging.getLogger(__name__)

class CutadaptOptionParser(OptionParser):
	def get_usage(self):
		return self.usage.lstrip().replace('%version', __version__)


class RestFileWriter(object):
	def __init__(self, file):
		self.file = file

	def write(self, match):
		rest = match.rest()
		if len(rest) > 0:
			print(rest, match.read.name, file=self.file)


class AdapterCutter(object):
	"""
	Repeatedly find one of multiple adapters in reads.
	The number of times the search is repeated is specified by the
	times parameter.
	"""

	def __init__(self, adapters, times=1, wildcard_file=None, info_file=None,
			rest_writer=None, action='trim'):
		"""
		adapters -- list of Adapter objects

		action -- What to do with a found adapter: None, 'trim', or 'mask'
		"""
		self.adapters = adapters
		self.times = times
		self.wildcard_file = wildcard_file
		self.info_file = info_file
		self.rest_writer = rest_writer
		self.action = action
		self.reads_matched = 0

	def _best_match(self, read):
		"""
		Find the best matching adapter in the given read.

		Return either an AdapterMatch instance or None if there are no matches.
		"""
		best = None
		for adapter in self.adapters:
			match = adapter.match_to(read)
			if match is None:
				continue

			# the no. of matches determines which adapter fits best
			if best is None or match.matches > best.matches:
				best = match
		return best

	def _write_info(self, read, matches):
		"""
		Write to the info, wildcard and rest files.
		# TODO
		# This design with a read having a .match attribute and
		# a match having a .read attribute is really confusing.
		"""
		match = read.match
		if self.rest_writer and match:
			self.rest_writer.write(match)

		if self.wildcard_file and match:
			print(match.wildcards(), read.name, file=self.wildcard_file)

		if self.info_file:
			if match:
				for m in matches:
					seq = m.read.sequence
					qualities = m.read.qualities
					if qualities is None:
						qualities = ''
					print(
						m.read.name,
						m.errors,
						m.rstart,
						m.rstop,
						seq[0:m.rstart],
						seq[m.rstart:m.rstop],
						seq[m.rstop:],
						m.adapter.name,
						qualities[0:m.rstart],
						qualities[m.rstart:m.rstop],
						qualities[m.rstop:],
						sep='\t', file=self.info_file
					)
			else:
				seq = read.sequence
				qualities = read.qualities if read.qualities is not None else ''
				print(read.name, -1, seq, qualities, sep='\t', file=self.info_file)

	def __call__(self, read):
		"""
		Determine the adapter that best matches the given read.
		Since the best adapter is searched repeatedly, a list
		of AdapterMatch instances is returned, which
		need to be applied consecutively to the read.
		The list is empty if there are no adapter matches.

		The read is converted to uppercase before it is compared to the adapter
		sequences.

		Cut found adapters from a single read. Return modified read.
		"""
		matches = []

		# try at most self.times times to remove an adapter
		trimmed_read = read
		for t in range(self.times):
			match = self._best_match(trimmed_read)
			if match is None:
				# nothing found
				break
			assert match.length > 0
			assert match.errors / match.length <= match.adapter.max_error_rate
			assert match.length - match.errors > 0
			matches.append(match)
			trimmed_read = match.adapter.trimmed(match)

		trimmed_read.match = matches[-1] if matches else None
		self._write_info(trimmed_read, matches)

		if not matches:
			return trimmed_read

		if __debug__:
			assert len(trimmed_read) < len(read), "Trimmed read isn't shorter than original"

		if self.action == 'trim':
			# read is already trimmed, nothing to do
			pass
		elif self.action == 'mask':
			# add N from last modification
			masked_sequence = trimmed_read.sequence
			for match in sorted(matches, reverse=True, key=lambda m: m.astart):
				ns = 'N' * (len(match.read.sequence) -
							len(match.adapter.trimmed(match).sequence))
				# add N depending on match position
				if match.front:
					masked_sequence = ns + masked_sequence
				else:
					masked_sequence += ns
			# set masked sequence as sequence with original quality
			trimmed_read.sequence = masked_sequence
			trimmed_read.qualities = matches[0].read.qualities

			assert len(trimmed_read.sequence) == len(read)
		elif self.action is None:
			trimmed_read = read
			trimmed_read.match = matches[-1]

		self.reads_matched += 1  # TODO move to filter class

		return trimmed_read


def process_single_reads(reader, modifiers, filters):
	"""
	Loop over reads, find adapters, trim reads, apply modifiers and
	output modified reads.

	Return a Statistics object.
	"""
	n = 0  # no. of processed reads
	total_bp = 0
	for read in reader:
		n += 1
		total_bp += len(read.sequence)
		for modifier in modifiers:
			read = modifier(read)
		for filter in filters:
			if filter(read):
				break

	return Statistics(n=n, total_bp1=total_bp, total_bp2=None)


def process_paired_reads(paired_reader, modifiers1, modifiers2, filters):
	"""
	Loop over reads, find adapters, trim reads, apply modifiers and
	output modified reads.

	Return a Statistics object.
	"""
	n = 0  # no. of processed reads
	total1_bp = 0
	total2_bp = 0
	for read1, read2 in paired_reader:
		n += 1
		total1_bp += len(read1.sequence)
		total2_bp += len(read2.sequence)
		for modifier in modifiers1:
			read1 = modifier(read1)
		for modifier in modifiers2:
			read2 = modifier(read2)
		for filter in filters:
			# Stop writing as soon as one of the filters was successful.
			if filter(read1, read2):
				break
	return Statistics(n=n, total_bp1=total1_bp, total_bp2=total2_bp)


def trimmed_and_untrimmed_writers(
		default_output,
		output_path, output_paired_path,
		untrimmed_path, untrimmed_paired_path,
		discard_trimmed,
		discard_untrimmed,
		fileformat,
		colorspace,
		):
	"""
	Figure out (from command-line parameters) where trimmed and untrimmed reads
	should be written.

	Return a pair (trimmed_writer, untrimmed_writer). The values are either
	Fasta/FastqWriters (depending on fileformat) or None, in which case no
	output should be produced.

	Later parameters have higher precedence.

	default_output -- If nothing else is specified below, this file-like object
		is returned for both trimmed and untrimmed output.
	output_path -- Path to output file for both trimmed and untrimmed output.
	untrimmed_path -- Path to an output file for untrimmed reads.
	discard_trimmed -- bool, overrides earlier options.
	discard_untrimmed -- bool, overrides earlier options.
	"""
	def open_writer(path_or_file, path_or_file2=None):
		return seqio.open(path_or_file, path_or_file2,
			mode='w', colorspace=colorspace, fileformat=fileformat)

	if discard_trimmed:
		if discard_untrimmed:
			untrimmed = None
		elif untrimmed_path is not None:
			untrimmed = open_writer(untrimmed_path, untrimmed_paired_path)
		elif output_path is not None:
			untrimmed = open_writer(output_path, output_paired_path)
		else:
			# should only happen for single-end data
			untrimmed = open_writer(default_output)
		return (None, untrimmed)

	if discard_untrimmed:
		trimmed = open_writer(default_output)
		if output_path is not None:
			trimmed = open_writer(output_path, output_paired_path)
		return (trimmed, None)

	trimmed = None
	untrimmed = None
	if output_path is not None:
		trimmed = untrimmed = open_writer(output_path, output_paired_path)
	if untrimmed_path is not None:
		untrimmed = open_writer(untrimmed_path, untrimmed_paired_path)
	if trimmed is None or untrimmed is None:
		default_writer = open_writer(default_output)
		if trimmed is None:
			trimmed = default_writer
		if untrimmed is None:
			untrimmed = default_writer
	return (trimmed, untrimmed)


def get_option_parser():
	parser = CutadaptOptionParser(usage=__doc__, version=__version__)

	parser.add_option("--debug", action='store_true', default=False,
		help="Print debugging information.")
	parser.add_option("-f", "--format",
		help="Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. "
			"Ignored when reading csfasta/qual files (default: auto-detect "
			"from file name extension).")

	group = OptionGroup(parser, "Options that influence how the adapters are found",
		description="Each of the three parameters -a, -b, -g can be used "
			"multiple times and in any combination to search for an entire set of "
			"adapters of possibly different types. Only the best matching "
			"adapter is trimmed from each read (but see the --times option). "
			"Instead of giving an adapter directly, you can also write "
			"file:FILE and the adapter sequences will be read from the given "
			"FASTA FILE.")
	group.add_option("-a", "--adapter", action="append", default=[], metavar="ADAPTER",
		dest="adapters",
		help="Sequence of an adapter that was ligated to the 3' end. The "
			"adapter itself and anything that follows is trimmed. If the "
			"adapter sequence ends with the '$' character, the adapter is "
			"anchored to the end of the read and only found if it is a "
			"suffix of the read.")
	group.add_option("-g", "--front", action="append", default=[], metavar="ADAPTER",
		help="Sequence of an adapter that was ligated to the 5' end. If the "
		"adapter sequence starts with the character '^', the adapter is "
		"'anchored'. An anchored adapter must appear in its entirety at the "
		"5' end of the read (it is a prefix of the read). A non-anchored adapter may "
		"appear partially at the 5' end, or it may occur within the read. If it is "
		"found within a read, the sequence preceding the adapter is also trimmed. "
		"In all cases, the adapter itself is trimmed.")
	group.add_option("-b", "--anywhere", action="append", default=[], metavar="ADAPTER",
		help="Sequence of an adapter that was ligated to the 5' or 3' end. If "
			"the adapter is found within the read or overlapping the 3' end of "
			"the read, the behavior is the same as for the -a option. If the "
			"adapter overlaps the 5' end (beginning of the read), the initial "
			"portion of the read matching the adapter is trimmed, but anything "
			"that follows is kept.")
	group.add_option("-e", "--error-rate", type=float, default=0.1,
		help="Maximum allowed error rate (no. of errors divided by the length "
			"of the matching region) (default: %default)")
	group.add_option("--no-indels", action='store_false', dest='indels', default=True,
		help="Do not allow indels in the alignments (allow only mismatches). "
			"(default: allow both mismatches and indels)")
	group.add_option("-n", "--times", type=int, metavar="COUNT", default=1,
		help="Remove up to COUNT adapters from each read (default: %default)")
	group.add_option("-O", "--overlap", type=int, metavar="LENGTH", default=3,
		help="Minimum overlap length. If the overlap between the read and the "
			"adapter is shorter than LENGTH, the read is not modified. "
			"This reduces the no. of bases trimmed purely due to short random "
			"adapter matches (default: %default).")
	group.add_option("--match-read-wildcards", action="store_true", default=False,
		help="Allow IUPAC wildcards in reads (default: %default).")
	group.add_option("-N", "--no-match-adapter-wildcards", action="store_false",
		default=True, dest='match_adapter_wildcards',
		help="Do not interpret IUPAC wildcards in adapters.")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Options for filtering of processed reads")
	group.add_option("--discard-trimmed", "--discard", action='store_true', default=False,
		help="Discard reads that contain an adapter. Also use -O to avoid "
			"discarding too many randomly matching reads!")
	group.add_option("--discard-untrimmed", "--trimmed-only", action='store_true', default=False,
		help="Discard reads that do not contain the adapter.")
	group.add_option("-m", "--minimum-length", type=int, default=0, metavar="LENGTH",
		help="Discard trimmed reads that are shorter than LENGTH. Reads that "
			"are too short even before adapter removal are also discarded. In "
			"colorspace, an initial primer is not counted (default: 0).")
	group.add_option("-M", "--maximum-length", type=int, default=sys.maxsize, metavar="LENGTH",
		help="Discard trimmed reads that are longer than LENGTH. "
			"Reads that are too long even before adapter removal "
			"are also discarded. In colorspace, an initial primer "
			"is not counted (default: no limit).")
	group.add_option("--no-trim", dest='action', action='store_const', const=None,
		help="Match and redirect reads to output/untrimmed-output as usual, "
			"but do not remove adapters.")
	group.add_option("--max-n", type=float, default=-1.0, metavar="LENGTH",
		help="The max proportion of N's allowed in a read. A number < 1 will be treated as a proportion while"
			 " a number > 1 will be treated as the maximum number of N's contained.")
	group.add_option("--mask-adapter", dest='action', action='store_const', const='mask',
		help="Mask adapters with 'N' characters instead of trimming them.")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Options that influence what gets output to where")
	group.add_option("--quiet", default=False, action='store_true',
		help="Do not print a report at the end.")
	group.add_option("-o", "--output", metavar="FILE",
		help="Write modified reads to FILE. FASTQ or FASTA format is chosen "
			"depending on input. The summary report is sent to standard output. "
			"Use '{name}' in FILE to demultiplex reads into multiple "
			"files. (default: trimmed reads are written to standard output)")
	group.add_option("--info-file", metavar="FILE",
		help="Write information about each read and its adapter matches into FILE. "
			"See the documentation for the file format.")
	group.add_option("-r", "--rest-file", metavar="FILE",
		help="When the adapter matches in the middle of a read, write the "
			"rest (after the adapter) into FILE.")
	group.add_option("--wildcard-file", metavar="FILE",
		help="When the adapter has N bases (wildcards), write adapter bases "
			"matching wildcard positions to FILE. When there are indels in the "
			"alignment, this will often not be accurate.")
	group.add_option("--too-short-output", metavar="FILE",
		help="Write reads that are too short (according to length specified by "
		"-m) to FILE. (default: discard reads)")
	group.add_option("--too-long-output", metavar="FILE",
		help="Write reads that are too long (according to length specified by "
		"-M) to FILE. (default: discard reads)")
	group.add_option("--untrimmed-output", default=None, metavar="FILE",
		help="Write reads that do not contain the adapter to FILE. (default: "
			"output to same file as trimmed reads)")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Additional read modifications")
	group.add_option("-u", "--cut", action='append', default=[], type=int, metavar="LENGTH",
		help="Remove LENGTH bases from the beginning or end of each read. "
			"If LENGTH is positive, bases are removed from the beginning of each read. "
			"If LENGTH is negative, bases are removed from the end of each read. "
			"This option can be specified twice if the LENGTHs have different signs.")
	group.add_option("-q", "--quality-cutoff", default=None, metavar="[5'CUTOFF,]3'CUTOFF",
		help="Trim low-quality bases from 5' and/or 3' ends of reads before "
			"adapter removal. If one value is given, only the 3' end is trimmed. "
			"If two comma-separated cutoffs are given, the 5' end is trimmed with "
			"the first cutoff, the 3' end with the second. See documentation for "
			"the algorithm. (default: no trimming)")
	group.add_option("--quality-base", type=int, default=33,
		help="Assume that quality values in FASTQ are encoded as ascii(quality "
			"+ QUALITY_BASE). This needs to be set to 64 for some old Illumina "
			"FASTQ files. Default: %default")
	group.add_option("--trim-n", action='store_true', default=False,
		help="Trim N's on ends of reads.")
	group.add_option("-x", "--prefix", default='',
		help="Add this prefix to read names. Use {name} to insert the name of the matching adapter.")
	group.add_option("-y", "--suffix", default='',
		help="Add this suffix to read names; can also include {name}")
	group.add_option("--strip-suffix", action='append', default=[],
		help="Remove this suffix from read names if present. Can be given multiple times.")
	group.add_option("--length-tag", metavar="TAG",
		help="Search for TAG followed by a decimal number in the description "
			"field of the read. Replace the decimal number with the correct "
			"length of the trimmed read. For example, use --length-tag 'length=' "
			"to correct fields like 'length=123'.")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Colorspace options")
	group.add_option("-c", "--colorspace", action='store_true', default=False,
		help="Enable colorspace mode: Also trim the color that is adjacent to the found adapter.")
	group.add_option("-d", "--double-encode", action='store_true', default=False,
		help="Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).")
	group.add_option("-t", "--trim-primer", action='store_true', default=False,
		help="Trim primer base and the first color (which is the transition "
			"to the first nucleotide)")
	group.add_option("--strip-f3", action='store_true', default=False,
		help="Strip the _F3 suffix of read names")
	group.add_option("--maq", "--bwa", action='store_true', default=False,
		help="MAQ- and BWA-compatible colorspace output. This enables -c, -d, "
			"-t, --strip-f3 and -y '/1'.")
	group.add_option("--no-zero-cap", dest='zero_cap', action='store_false',
		help="Do not change negative quality values to zero in colorspace "
			"data. By default, they are changed to zero since many tools have "
			"problems with negative qualities.")
	group.add_option("--zero-cap", "-z", action='store_true',
		help="Change negative quality values to zero. This is enabled "
		"by default when -c/--colorspace is also enabled. Use the above option "
		"to disable it.")
	parser.set_defaults(zero_cap=None, action='trim')
	parser.add_option_group(group)

	group = OptionGroup(parser, "Paired-end options", description="The "
		"-A/-G/-B/-U options work like their -a/-b/-g/-u counterparts.")
	group.add_option("-A", dest='adapters2', action='append', default=[], metavar='ADAPTER',
		help="3' adapter to be removed from second read in a pair.")
	group.add_option("-G", dest='front2', action='append', default=[], metavar='ADAPTER',
		help="5' adapter to be removed from second read in a pair.")
	group.add_option("-B", dest='anywhere2', action='append', default=[], metavar='ADAPTER',
		help="5'/3 adapter to be removed from second read in a pair.")
	group.add_option("-U", dest='cut2', action='append', default=[], type=int, metavar="LENGTH",
		help="Remove LENGTH bases from the beginning or end of each second read (see --cut).")
	group.add_option("-p", "--paired-output", metavar="FILE",
		help="Write second read in a pair to FILE.")
	group.add_option("--interleaved", action='store_true', default=False,
		help="Read and write interleaved paired-end reads.")
	group.add_option("--untrimmed-paired-output", metavar="FILE",
		help="Write the second read in a pair to this FILE when no adapter "
			"was found in the first read. Use this option together with "
			"--untrimmed-output when trimming paired-end reads. (Default: output "
			"to same file as trimmed reads.)")
	parser.add_option_group(group)

	return parser


def main(cmdlineargs=None, default_outfile=sys.stdout):
	"""
	Main function that evaluates command-line parameters and iterates
	over all reads.

	default_outfile is the file to which trimmed reads are sent if the ``-o``
	parameter is not used.
	"""
	logging.basicConfig(level=logging.INFO, format='%(message)s')  #  %(levelname)s
	parser = get_option_parser()
	if cmdlineargs is None:
		cmdlineargs = sys.argv[1:]
	options, args = parser.parse_args(args=cmdlineargs)

	if len(args) == 0:
		parser.error("At least one parameter needed: name of a FASTA or FASTQ file.")
	elif len(args) > 2:
		parser.error("Too many parameters.")
	input_filename = args[0]

	# Find out which 'mode' we need to use.
	# Default: single-read trimming (neither -p nor -A/-G/-B/-U given)
	paired = False
	if options.paired_output:
		# Modify first read only, keep second in sync (-p given, but not -A/-G/-B/-U).
		# This exists for backwards compatibility ('legacy mode').
		paired = 'first'
	if options.adapters2 or options.front2 or options.anywhere2 or options.cut2 or options.interleaved:
		# Full paired-end trimming when both -p and -A/-G/-B/-U given
		# Also the read modifications (such as quality trimming) are applied
		# to second read.
		paired = 'both'

	if paired and len(args) == 1 and not options.interleaved:
		parser.error("When paired-end trimming is enabled via -A/-G/-B/-U or -p, "
			"two input files are required.")
	if options.interleaved and len(args) != 1:
		parser.error("When reading interleaved files, only one input file may "
			"be given.")
	if paired:
		if options.interleaved:
			input_paired_filename = None
		else:
			input_paired_filename = args[1]
		quality_filename = None
	else:
		input_paired_filename = None
		if len(args) == 2:
			if args[0].endswith('.qual'):
				parser.error("The QUAL file must be the second argument.")
			quality_filename = args[1]
		else:
			quality_filename = None

	if paired:
		if not options.interleaved:
			if not options.paired_output:
				parser.error("When paired-end trimming is enabled via -A/-G/-B/-U, "
					"a second output file needs to be specified via -p (--paired-output).")
			if bool(options.untrimmed_output) != bool(options.untrimmed_paired_output):
				parser.error("When trimming paired-end reads, you must use either none "
					"or both of the --untrimmed-output/--untrimmed-paired-output options.")
	else:
		if options.untrimmed_paired_output:
			parser.error("Option --untrimmed-paired-output can only be used when "
				"trimming paired-end reads (with option -p).")
		if input_filename.endswith('.qual'):
			parser.error("Need a FASTA file in addition to the QUAL file.")
		if options.format is not None and quality_filename is not None:
			parser.error("If a pair of .fasta and .qual files is given, the -f/--format parameter cannot be used.")

	if options.format is not None and options.format.lower() not in ['fasta', 'fastq', 'sra-fastq']:
		parser.error("The input file format must be either 'fasta', 'fastq' or "
			"'sra-fastq' (not '{0}').".format(options.format))

	# Open input file(s)
	try:
		reader = seqio.open(input_filename, file2=input_paired_filename,
				qualfile=quality_filename, colorspace=options.colorspace,
				fileformat=options.format, interleaved=options.interleaved)
	except (seqio.UnknownFileType, IOError) as e:
		parser.error(e)
	fileformat = 'fastq' if reader.delivers_qualities else 'fasta'

	if options.quality_cutoff is not None:
		cutoffs = options.quality_cutoff.split(',')
		if len(cutoffs) == 1:
			try:
				cutoffs = [0, int(cutoffs[0])]
			except ValueError as e:
				parser.error("Quality cutoff value not recognized: {0}".format(e))
		elif len(cutoffs) == 2:
			try:
				cutoffs = [int(cutoffs[0]), int(cutoffs[1])]
			except ValueError as e:
				parser.error("Quality cutoff value not recognized: {0}".format(e))
		else:
			parser.error("Expected one value or two values separated by comma for the quality cutoff")
	else:
		cutoffs = None
	filters = []
	too_short_writer = None  # too short reads go here
	too_short_filter = None
	# TODO pass file name to TooShortReadFilter, add a .close() method?
	if options.minimum_length > 0:
		if options.too_short_output:
			too_short_writer = seqio.open(options.too_short_output, mode='w',
				fileformat=fileformat, colorspace=options.colorspace)
		else:
			too_short_writer = None
		too_short_filter = TooShortReadFilter(options.minimum_length,
			too_short_writer, paired=='both')
		filters.append(too_short_filter)
	too_long_writer = None  # too long reads go here
	too_long_filter = None
	if options.maximum_length < sys.maxsize:
		if options.too_long_output is not None:
			too_long_writer = seqio.open(options.too_long_output,
				mode='w', fileformat=fileformat, colorspace=options.colorspace)
		else:
			too_long_writer = None
		too_long_filter = TooLongReadFilter(options.maximum_length,
			too_long_writer, check_second=paired=='both')
		filters.append(too_long_filter)

	if options.max_n != -1:
		filters.append(NContentFilter(options.max_n, check_second=paired=='both'))

	demultiplexer = None
	if options.output is not None and '{name}' in options.output:
		if options.discard_trimmed:
			parser.error("Do not use --discard-trimmed when demultiplexing.")
		if paired:
			parser.error("Demultiplexing not supported for paired-end files, yet.")
		untrimmed = options.output.replace('{name}', 'unknown')
		if options.untrimmed_output:
			untrimmed = options.untrimmed_output
		if options.discard_untrimmed:
			untrimmed = None
		demultiplexer = Demultiplexer(options.output, untrimmed,
			fileformat=fileformat, colorspace=colorspace)
		filters.append(demultiplexer)
		trimmed_writer = None
		untrimmed_writer = None
	else:
		trimmed_writer, untrimmed_writer = trimmed_and_untrimmed_writers(
			default_outfile,
			options.output, options.paired_output,
			options.untrimmed_output, options.untrimmed_paired_output,
			options.discard_trimmed,
			options.discard_untrimmed,
			fileformat=fileformat,
			colorspace=options.colorspace)

		filters.append(DiscardUntrimmedFilter(
			untrimmed_writer,
			check_second=paired=='both'))
		filters.append(DiscardTrimmedFilter(
			trimmed_writer,
			check_second=paired=='both'))

	if options.maq:
		options.colorspace = True
		options.double_encode = True
		options.trim_primer = True
		options.strip_suffix.append('_F3')
		options.suffix = "/1"
	if options.zero_cap is None:
		options.zero_cap = options.colorspace
	if options.trim_primer and not options.colorspace:
		parser.error("Trimming the primer makes only sense in colorspace.")
	if options.double_encode and not options.colorspace:
		parser.error("Double-encoding makes only sense in colorspace.")
	if options.anywhere and options.colorspace:
		parser.error("Using --anywhere with colorspace reads is currently not supported (if you think this may be useful, contact the author).")
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

	if options.colorspace:
		if options.match_read_wildcards:
			parser.error('IUPAC wildcards not supported in colorspace')
		options.match_adapter_wildcards = False

	ADAPTER_CLASS = ColorspaceAdapter if options.colorspace else Adapter
	try:
		# TODO refactor this a bit
		def collect(back, anywhere, front):
			adapters = []
			for name, seq, where in gather_adapters(back, anywhere, front):
				if not seq:
					parser.error("The adapter sequence is empty.")
				adapter = ADAPTER_CLASS(seq, where, options.error_rate,
					options.overlap, options.match_read_wildcards,
					options.match_adapter_wildcards, name=name, indels=options.indels)
				if options.debug:
					adapter.enable_debug()
				adapters.append(adapter)
			return adapters

		adapters = collect(options.adapters, options.anywhere, options.front)
		adapters2 = collect(options.adapters2, options.anywhere2, options.front2)
	except IOError as e:
		if e.errno == errno.ENOENT:
			parser.error(e)
		raise

	if not adapters and not adapters2 and not cutoffs and \
			options.cut == [] and options.cut2 == [] and \
			options.minimum_length == 0 and \
			options.maximum_length == sys.maxsize and \
			quality_filename is None and \
			options.max_n == -1:
		parser.error("You need to provide at least one adapter sequence.")

	# Create the single-end processing pipeline (a list of "modifiers")
	modifiers = []
	if options.cut:
		if len(options.cut) > 2:
			parser.error("You cannot remove bases from more than two ends.")
		if len(options.cut) == 2 and options.cut[0] * options.cut[1] > 0:
			parser.error("You cannot remove bases from the same end twice.")
		for cut in options.cut:
			if cut != 0:
				modifiers.append(UnconditionalCutter(cut))

	if cutoffs:
		modifiers.append(QualityTrimmer(cutoffs[0], cutoffs[1], options.quality_base))
	if adapters:
		adapter_cutter = AdapterCutter(adapters, options.times,
				options.wildcard_file, options.info_file,
				rest_writer, options.action)
		modifiers.append(adapter_cutter)
	else:
		adapter_cutter = None

	# Modifiers that apply to both reads of paired-end reads unless in legacy mode
	modifiers_both = []
	if options.trim_n:
		modifiers_both.append(NEndTrimmer())
	if options.length_tag:
		modifiers_both.append(LengthTagModifier(options.length_tag))
	if options.strip_f3:
		options.strip_suffix.append('_F3')
	for suffix in options.strip_suffix:
		modifiers_both.append(SuffixRemover(suffix))
	if options.prefix or options.suffix:
		modifiers_both.append(PrefixSuffixAdder(options.prefix, options.suffix))
	if options.double_encode:
		modifiers_both.append(DoubleEncoder())
	if options.zero_cap and reader.delivers_qualities:
		modifiers_both.append(ZeroCapper(quality_base=options.quality_base))
	if options.trim_primer:
		modifiers_both.append(PrimerTrimmer)
	modifiers.extend(modifiers_both)

	# For paired-end data, create a second processing pipeline.
	# However, if no second-read adapters were given (via -A/-G/-B/-U), we need to
	# be backwards compatible and *no modifications* are done to the second read.
	modifiers2 = []
	if paired == 'both':
		if options.cut2:
			if len(options.cut2) > 2:
				parser.error("You cannot remove bases from more than two ends.")
			if len(options.cut2) == 2 and options.cut2[0] * options.cut2[1] > 0:
				parser.error("You cannot remove bases from the same end twice.")
			for cut in options.cut2:
				if cut != 0:
					modifiers2.append(UnconditionalCutter(cut))

		if cutoffs:
			modifiers2.append(QualityTrimmer(cutoffs[0], cutoffs[1], options.quality_base))
		if adapters2:
			adapter_cutter2 = AdapterCutter(adapters2, options.times,
					None, None, None, options.action)
			modifiers2.append(adapter_cutter2)
		else:
			adapter_cutter2 = None
		modifiers2.extend(modifiers_both)

	# Due to backwards compatibility, from here on logging output needs to be
	# sent to standard output instead of standard error if the -o option is used.
	if options.output:
		logger.root.handlers = []
		logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stdout)
	logger.info("This is cutadapt %s with Python %s", __version__, platform.python_version())
	logger.info("Command line parameters: %s", " ".join(cmdlineargs))
	logger.info("Trimming %s adapter%s with at most %.1f%% errors in %s mode ...",
		len(adapters) + len(adapters2), 's' if len(adapters) + len(adapters2) != 1 else '',
		options.error_rate * 100,
		{ False: 'single-end', 'first': 'paired-end legacy', 'both': 'paired-end' }[paired])

	if paired == 'first' and (modifiers_both or cutoffs):
		logger.warn('\n'.join(textwrap.wrap('WARNING: Requested read '
			'modifications are applied only to the first '
			'read since backwards compatibility mode is enabled. '
			'To modify both reads, also use any of the -A/-B/-G/-U options. '
			'Use a dummy adapter sequence when necessary: -A XXX')))

	start_time = time.clock()
	try:
		if paired:
			stats = process_paired_reads(reader, modifiers, modifiers2, filters)
		else:
			stats = process_single_reads(reader, modifiers, filters)
	except KeyboardInterrupt as e:
		print("Interrupted", file=sys.stderr)
		sys.exit(130)
	except IOError as e:
		if e.errno == errno.EPIPE:
			sys.exit(1)
		raise
	except (seqio.FormatError, EOFError) as e:
		sys.exit("cutadapt: error: {0}".format(e))

	# close open files
	for f in [trimmed_writer, untrimmed_writer,
			options.rest_file, options.wildcard_file,
			options.info_file, too_short_writer, too_long_writer,
			options.info_file, demultiplexer]:
		if f is not None and f is not sys.stdin and f is not sys.stdout:
			f.close()

	elapsed_time = time.clock() - start_time
	if not options.quiet:
		stats.collect((adapters, adapters2), elapsed_time,
			modifiers, modifiers2, filters)
		# send statistics to stderr if result was sent to stdout
		stat_file = sys.stderr if options.output is None else None
		with redirect_standard_output(stat_file):
			print_report(stats, (adapters, adapters2))


if __name__ == '__main__':
	main()
