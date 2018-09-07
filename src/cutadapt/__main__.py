#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kate: word-wrap off; remove-trailing-spaces all;
#
# Copyright (c) 2010-2018 Marcel Martin <marcel.martin@scilifelab.se>
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

Copyright (C) 2010-2018 Marcel Martin <marcel.martin@scilifelab.se>

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

Citation:

Marcel Martin. Cutadapt removes adapter sequences from high-throughput
sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011.
http://dx.doi.org/10.14806/ej.17.1.200

Run "cutadapt --help" to see all command-line options.
See https://cutadapt.readthedocs.io/ for full documentation.
"""

from __future__ import print_function, division, absolute_import

import sys
import errno
import time
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
import logging
import platform
import textwrap
from xopen import xopen

from cutadapt import seqio, __version__
from cutadapt.adapters import AdapterParser
from cutadapt.modifiers import (LengthTagModifier, SuffixRemover, PrefixSuffixAdder,
	DoubleEncoder, ZeroCapper, PrimerTrimmer, QualityTrimmer, UnconditionalCutter,
	NEndTrimmer, AdapterCutter, NextseqQualityTrimmer, Shortener)
from cutadapt.report import print_report, print_minimal_report, redirect_standard_output
from cutadapt.pipeline import SingleEndPipeline, PairedEndPipeline, OutputFiles, ParallelPipelineRunner
from cutadapt.utils import available_cpu_count
from cutadapt.compat import PY3

logger = logging.getLogger()


class CutadaptOptionParser(OptionParser):
	def get_usage(self):
		return self.usage.lstrip().replace('%version', __version__)

	def error(self, msg):
		print('Run "cutadapt --help" to see command-line options.', file=sys.stderr)
		print('See https://cutadapt.readthedocs.io/ for full documentation.', file=sys.stderr)
		self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


class CommandLineError(Exception):
	pass


class NiceFormatter(logging.Formatter):
	"""
	Do not prefix "INFO:" to info-level log messages (but do it for all other
	levels).

	Based on http://stackoverflow.com/a/9218261/715090 .
	"""
	def format(self, record):
		if record.levelno != logging.INFO:
			record.msg = '{}: {}'.format(record.levelname, record.msg)
		return super(NiceFormatter, self).format(record)


def setup_logging(stdout=False, quiet=False):
	"""
	Attach handler to the global logger object
	"""
	# Due to backwards compatibility, logging output is sent to standard output
	# instead of standard error if the -o option is used.
	stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
	stream_handler.setFormatter(NiceFormatter())
	stream_handler.setLevel(logging.ERROR if quiet else logging.INFO)
	logger.setLevel(logging.INFO)
	logger.addHandler(stream_handler)


def get_option_parser():
	parser = CutadaptOptionParser(usage=__doc__, version=__version__)

	parser.add_option("--debug", action='store_true', default=False,
		help="Print debugging information.")
	parser.add_option("-f", "--format",
		help="Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. "
			"Ignored when reading csfasta/qual files. Default: auto-detect "
			"from file name extension.")
	parser.add_option('-j', '--cores', type=int, default=1,
		help='Number of CPU cores to use. Use 0 to auto-detect. Default: %default')

	# Hidden options
	parser.add_option("--gc-content", type=float, default=50,  # it's a percentage
		help=SUPPRESS_HELP)
	parser.add_option("--buffer-size", type=int, default=4000000,
		help=SUPPRESS_HELP)  # buffer size for the reader process when running in parallel

	group = OptionGroup(parser, "Finding adapters",
		description="Parameters -a, -g, -b specify adapters to be removed from "
			"each read (or from the first read in a pair if data is paired). "
			"If specified multiple times, only the best matching adapter is "
			"trimmed (but see the --times option). When the special notation "
			"'file:FILE' is used, adapter sequences are read from the given "
			"FASTA file.")
	group.add_option("-a", "--adapter", action="append", default=[], metavar="ADAPTER",
		dest="adapters",
		help="Sequence of an adapter ligated to the 3' end (paired data: of the "
			"first read). The adapter and subsequent bases are trimmed. If a "
			"'$' character is appended ('anchoring'), the adapter is only "
			"found if it is a suffix of the read.")
	group.add_option("-g", "--front", action="append", default=[], metavar="ADAPTER",
		help="Sequence of an adapter ligated to the 5' end (paired data: of the "
			"first read). The adapter and any preceding bases are trimmed. "
			"Partial matches at the 5' end are allowed. If a '^' character is "
			"prepended ('anchoring'), the adapter is only found if it is a "
			"prefix of the read.")
	group.add_option("-b", "--anywhere", action="append", default=[], metavar="ADAPTER",
		help="Sequence of an adapter that may be ligated to the 5' or 3' end "
			"(paired data: of the first read). Both types of matches as "
			"described under -a und -g are allowed. If the first base of the "
			"read is part of the match, the behavior is as with -g, otherwise "
			"as with -a. This option is mostly for rescuing failed library "
			"preparations - do not use if you know which end your adapter was "
			"ligated to!")
	group.add_option("-e", "--error-rate", type=float, default=0.1, metavar="RATE",
		help="Maximum allowed error rate as value between 0 and 1 (no. of "
			"errors divided by length of matching region). Default: %default (=10%)")
	group.add_option("--no-indels", action='store_false', dest='indels', default=True,
		help="Allow only mismatches in alignments. "
			"Default: allow both mismatches and indels")
	group.add_option("-n", "--times", type=int, metavar="COUNT", default=1,
		help="Remove up to COUNT adapters from each read. Default: %default")
	group.add_option("-O", "--overlap", type=int, metavar="MINLENGTH", default=3,
		help="Require MINLENGTH overlap between read and adapter for an adapter "
			"to be found. Default: %default")
	group.add_option("--match-read-wildcards", action="store_true", default=False,
		help="Interpret IUPAC wildcards in reads. Default: %default")
	group.add_option("-N", "--no-match-adapter-wildcards", action="store_false",
		default=True, dest='match_adapter_wildcards',
		help="Do not interpret IUPAC wildcards in adapters.")
	group.add_option("--action", choices=('mask', 'trim', 'none'),
		help="What to do with found adapters. trim: remove; "
			"mask: replace with 'N' characters; "
			"none: leave unchanged (useful with "
			"--discard-untrimmed). Default: trim")
	group.add_option("--no-trim", dest='action', action='store_const', const='none',
		help="Deprecated synonym for --action=none")
	group.add_option("--mask-adapter", dest='action', action='store_const', const='mask',
		help="Deprecated synonym for --action=mask")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Additional read modifications")
	group.add_option("-u", "--cut", action='append', default=[], type=int, metavar="LENGTH",
		help="Remove bases from each read (first read only if paired). "
			"If LENGTH is positive, remove bases from the beginning. "
			"If LENGTH is negative, remove bases from the end. "
			"Can be used twice if LENGTHs have different signs. "
			"This is applied *before* adapter trimming.")
	group.add_option("--nextseq-trim", type=int, default=None, metavar="3'CUTOFF",
		help="NextSeq-specific quality trimming (each read). Trims also dark "
			"cycles appearing as high-quality G bases.")
	group.add_option("-q", "--quality-cutoff", default=None, metavar="[5'CUTOFF,]3'CUTOFF",
		help="Trim low-quality bases from 5' and/or 3' ends of each read before "
			"adapter removal. Applied to both reads if data is paired. If one "
			"value is given, only the 3' end is trimmed. If two "
			"comma-separated cutoffs are given, the 5' end is trimmed with "
			"the first cutoff, the 3' end with the second.")
	group.add_option("--quality-base", type=int, default=33, metavar='N',
		help="Assume that quality values in FASTQ are encoded as ascii(quality "
			"+ N). This needs to be set to 64 for some old Illumina "
			"FASTQ files. Default: %default")
	group.add_option("--length", "-l", type=int, default=None, metavar="LENGTH",
			help="Shorten reads to LENGTH. Positive values remove bases at the end "
			"while negative ones remove bases at the beginning. This and the following modifications "
			"are applied after adapter trimming.")
	group.add_option("--trim-n", action='store_true', default=False,
		help="Trim N's on ends of reads.")
	group.add_option("--length-tag", metavar="TAG",
		help="Search for TAG followed by a decimal number in the description "
			"field of the read. Replace the decimal number with the correct "
			"length of the trimmed read. For example, use --length-tag 'length=' "
			"to correct fields like 'length=123'.")
	group.add_option("--strip-suffix", action='append', default=[],
		help="Remove this suffix from read names if present. Can be given multiple times.")
	group.add_option("-x", "--prefix", default='',
		help="Add this prefix to read names. Use {name} to insert the name of the matching adapter.")
	group.add_option("-y", "--suffix", default='',
		help="Add this suffix to read names; can also include {name}")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Filtering of processed reads",
		description="Filters are applied after above read modifications. "
			"Paired-end reads are always discarded pairwise (see also "
			"--pair-filter).")
	group.add_option("-m", "--minimum-length", default=None, metavar="LEN[:LEN2]",
		help="Discard reads shorter than LEN. Default: 0")
	group.add_option("-M", "--maximum-length", default=None, metavar="LEN[:LEN2]",
		help="Discard reads longer than LEN. Default: no limit")
	group.add_option("--max-n", type=float, default=None, metavar="COUNT",
		help="Discard reads with more than COUNT 'N' bases. If COUNT is a number "
			"between 0 and 1, it is interpreted as a fraction of the read length.")
	group.add_option("--discard-trimmed", "--discard", action='store_true', default=False,
		help="Discard reads that contain an adapter. Also use -O to avoid "
			"discarding too many randomly matching reads!")
	group.add_option("--discard-untrimmed", "--trimmed-only", action='store_true', default=False,
		help="Discard reads that do not contain an adapter.")
	group.add_option("--discard-casava", action='store_true', default=False,
		help="Discard reads that did not pass CASAVA filtering (header has :Y:).")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Output")
	group.add_option("--quiet", default=False, action='store_true',
		help="Print only error messages.")
	group.add_option("--report", choices=('full', 'minimal'), default=None,
		help="Which type of report to print. Default: full")
	group.add_option("-o", "--output", metavar="FILE",
		help="Write trimmed reads to FILE. FASTQ or FASTA format is chosen "
			"depending on input. The summary report is sent to standard output. "
			"Use '{name}' in FILE to demultiplex reads into multiple "
			"files. Default: write to standard output")
	group.add_option("--info-file", metavar="FILE",
		help="Write information about each read and its adapter matches into FILE. "
			"See the documentation for the file format.")
	group.add_option("-r", "--rest-file", metavar="FILE",
		help="When the adapter matches in the middle of a read, write the "
			"rest (after the adapter) to FILE.")
	group.add_option("--wildcard-file", metavar="FILE",
		help="When the adapter has N wildcard bases, write adapter bases "
			"matching wildcard positions to FILE. (Inaccurate with indels.)")
	group.add_option("--too-short-output", metavar="FILE",
		help="Write reads that are too short (according to length specified by "
		"-m) to FILE. Default: discard reads")
	group.add_option("--too-long-output", metavar="FILE",
		help="Write reads that are too long (according to length specified by "
		"-M) to FILE. Default: discard reads")
	group.add_option("--untrimmed-output", default=None, metavar="FILE",
		help="Write reads that do not contain any adapter to FILE. Default: "
			"output to same file as trimmed reads")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Colorspace options")
	group.add_option("-c", "--colorspace", action='store_true', default=False,
		help="Enable colorspace mode")
	group.add_option("-d", "--double-encode", action='store_true', default=False,
		help="Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).")
	group.add_option("-t", "--trim-primer", action='store_true', default=False,
		help="Trim primer base and the first color")
	group.add_option("--strip-f3", action='store_true', default=False,
		help="Strip the _F3 suffix of read names")
	group.add_option("--maq", "--bwa", action='store_true', default=False,
		help="MAQ- and BWA-compatible colorspace output. This enables -c, -d, "
			"-t, --strip-f3 and -y '/1'.")
	group.add_option("--zero-cap", "-z", action='store_true',
		help="Change negative quality values to zero. Enabled by default "
			"in colorspace mode since many tools have problems with "
			"negative qualities")
	group.add_option("--no-zero-cap", dest='zero_cap', action='store_false',
		help="Disable zero capping")
	parser.set_defaults(zero_cap=None, action='trim')
	parser.add_option_group(group)

	group = OptionGroup(parser, "Paired-end options", description="The "
		"-A/-G/-B/-U options work like their -a/-b/-g/-u counterparts, but "
		"are applied to the second read in each pair.")
	group.add_option("-A", dest='adapters2', action='append', default=[], metavar='ADAPTER',
		help="3' adapter to be removed from second read in a pair.")
	group.add_option("-G", dest='front2', action='append', default=[], metavar='ADAPTER',
		help="5' adapter to be removed from second read in a pair.")
	group.add_option("-B", dest='anywhere2', action='append', default=[], metavar='ADAPTER',
		help="5'/3 adapter to be removed from second read in a pair.")
	group.add_option("-U", dest='cut2', action='append', default=[], type=int, metavar="LENGTH",
		help="Remove LENGTH bases from second read in a pair.")
	group.add_option("-p", "--paired-output", metavar="FILE",
		help="Write second read in a pair to FILE.")
	# Setting the default for pair_filter to None allows us to find out whether
	# the option was used at all.
	group.add_option("--pair-filter", metavar='(any|both|first)', default=None,
		choices=("any", "both", "first"),
		help="Which of the reads in a paired-end read have to match the "
			"filtering criterion in order for the pair to be filtered. "
			"Default: any")
	group.add_option("--interleaved", action='store_true', default=False,
		help="Read and write interleaved paired-end reads.")
	group.add_option("--untrimmed-paired-output", metavar="FILE",
		help="Write second read in a pair to this FILE when no adapter "
			"was found. Use with --untrimmed-output. Default: output "
			"to same file as trimmed reads")
	group.add_option("--too-short-paired-output", metavar="FILE", default=None,
		help="Write second read in a pair to this file if pair is too short. "
			"Use also --too-short-output.")
	group.add_option("--too-long-paired-output", metavar="FILE", default=None,
		help="Write second read in a pair to this file if pair is too long. "
			"Use also --too-long-output.")
	parser.add_option_group(group)

	return parser


def parse_cutoffs(s):
	"""Parse a string INT[,INT] into a two-element list of integers"""
	cutoffs = s.split(',')
	if len(cutoffs) == 1:
		try:
			cutoffs = [0, int(cutoffs[0])]
		except ValueError as e:
			raise CommandLineError("Quality cutoff value not recognized: {0}".format(e))
	elif len(cutoffs) == 2:
		try:
			cutoffs = [int(cutoffs[0]), int(cutoffs[1])]
		except ValueError as e:
			raise CommandLineError("Quality cutoff value not recognized: {0}".format(e))
	else:
		raise CommandLineError("Expected one value or two values separated by comma for "
			"the quality cutoff")
	return cutoffs


def parse_lengths(s):
	"""Parse [INT][:[INT]] into a pair of integers. If a value is omitted, use None

	>>> parse_lengths('25')
	(25,)
	>>> parse_lengths('17:25')
	(17, 25)
	>>> parse_lengths('25:')
	(25, None)
	>>> parse_lengths(':25')
	(None, 25)
	"""
	fields = s.split(':')
	if len(fields) not in (1, 2):
		raise CommandLineError("Only at most one colon is allowed")
	try:
		values = tuple(int(f) if f != '' else None for f in fields)
	except ValueError as e:
		raise CommandLineError("Value not recognized: {0}".format(e))
	if len(values) == 2 and values[0] is None and values[1] is None:
		raise CommandLineError("Cannot parse {!r}: At least one length needs to be given".format(s))
	return tuple(values)


def open_output_files(options, default_outfile, interleaved):
	"""
	Return an OutputFiles instance. If demultiplex is True, the untrimmed, untrimmed2, out and out2
	attributes are not opened files, but paths (out and out2 with the '{name}' template).
	"""
	rest_file = info_file = wildcard = None
	if options.rest_file is not None:
		rest_file = xopen(options.rest_file, 'w')
	if options.info_file is not None:
		info_file = xopen(options.info_file, 'w')
	if options.wildcard_file is not None:
		wildcard = xopen(options.wildcard_file, 'w')

	def open2(path1, path2):
		file1 = file2 = None
		if path1 is not None:
			file1 = xopen(path1, 'w')
			if path2 is not None:
				file2 = xopen(path2, 'w')
		return file1, file2

	too_short = too_short2 = None
	if options.minimum_length is not None:
		too_short, too_short2 = open2(options.too_short_output, options.too_short_paired_output)

	too_long = too_long2 = None
	if options.maximum_length is not None:
		too_long, too_long2 = open2(options.too_long_output, options.too_long_paired_output)

	if int(options.discard_trimmed) + int(options.discard_untrimmed) + int(
					options.untrimmed_output is not None) > 1:
		raise CommandLineError("Only one of the --discard-trimmed, --discard-untrimmed "
			"and --untrimmed-output options can be used at the same time.")

	demultiplex = options.output is not None and '{name}' in options.output
	if options.paired_output is not None and (demultiplex != ('{name}' in options.paired_output)):
		raise CommandLineError('When demultiplexing paired-end data, "{name}" must appear in '
			'both output file names (-o and -p)')

	if demultiplex:
		if options.discard_trimmed:
			raise CommandLineError("Do not use --discard-trimmed when demultiplexing.")

		out = options.output
		untrimmed = options.output.replace('{name}', 'unknown')
		if options.untrimmed_output:
			untrimmed = options.untrimmed_output
		if options.discard_untrimmed:
			untrimmed = None

		if options.paired_output is not None:
			out2 = options.paired_output
			untrimmed2 = options.paired_output.replace('{name}', 'unknown')
			if options.untrimmed_paired_output:
				untrimmed2 = options.untrimmed_paired_output
			if options.discard_untrimmed:
				untrimmed2 = None

		else:
			untrimmed2 = out2 = None
	else:
		untrimmed, untrimmed2 = open2(options.untrimmed_output, options.untrimmed_paired_output)
		out, out2 = open2(options.output, options.paired_output)
		if out is None:
			out = default_outfile

	if demultiplex:
		assert out is not None and '{name}' in out and (out2 is None or '{name}' in out2)
	return OutputFiles(
		rest=rest_file,
		info=info_file,
		wildcard=wildcard,
		too_short=too_short,
		too_short2=too_short2,
		too_long=too_long,
		too_long2=too_long2,
		untrimmed=untrimmed,
		untrimmed2=untrimmed2,
		out=out,
		out2=out2,
		demultiplex=demultiplex,
		interleaved=interleaved,
	)


def determine_paired_mode(options):
	"""
	Determine the paired-end mode: single-end, paired-end or legacy paired-end.

	Return False, 'first' or 'both'.

	False -- single-end
	'first' -- Backwards-compatible "legacy" mode in which read modifications apply only to read 1
	'both' -- normal paired-end mode in which read modifications apply to read 1 and 2

	Legacy mode is deactivated as soon as any option is used that exists only in cutadapt 1.8 or
	later, such as -A/-G/-B/-U/--interleaved/--nextseq-trim.
	"""
	paired = False
	if options.paired_output:
		paired = 'first'

	# Switch off legacy mode if certain options given
	if paired and options.nextseq_trim:
		paired = 'both'
	if (options.adapters2 or options.front2 or options.anywhere2 or
			options.cut2 or options.interleaved or options.pair_filter or
			options.too_short_paired_output or options.too_long_paired_output):
		paired = 'both'
	return paired


def determine_interleaved(options, args):
	is_interleaved_input = False
	is_interleaved_output = False
	if options.interleaved:
		is_interleaved_input = len(args) == 1
		is_interleaved_output = not options.paired_output
		if not is_interleaved_input and not is_interleaved_output:
			raise CommandLineError("When --interleaved is used, you cannot provide both two "
				"input files and two output files")
	return is_interleaved_input, is_interleaved_output


def input_files_from_parsed_args(args, paired, interleaved):
	"""
	Return tuple (input_filename, input_paired_filename, quality_filename)
	"""
	if len(args) == 0:
		raise CommandLineError("Please give me something to do!")
	elif len(args) > 2:
		raise CommandLineError("Too many parameters.")
	input_filename = args[0]
	if input_filename.endswith('.qual'):
		raise CommandLineError("If a .qual file is given, it must be the second argument.")
	if paired and len(args) == 1 and not interleaved:
		raise CommandLineError("When paired-end trimming is enabled via -A/-G/-B/-U/"
			"--interleaved or -p, two input files are required.")

	input_paired_filename = None
	quality_filename = None
	if paired:
		if not interleaved:
			input_paired_filename = args[1]
	elif len(args) == 2:
		quality_filename = args[1]

	return input_filename, input_paired_filename, quality_filename


def pipeline_from_parsed_args(options, paired, pair_filter_mode, quality_filename, is_interleaved_output):
	"""
	Setup a processing pipeline from parsed command-line options.

	If there are any problems parsing the arguments, a CommandLineError is thrown.

	Return an instance of Pipeline (SingleEndPipeline or PairedEndPipeline)
	"""

	if not paired:
		if options.untrimmed_paired_output:
			raise CommandLineError("Option --untrimmed-paired-output can only be used when "
				"trimming paired-end reads (with option -p).")

	if paired:
		if not is_interleaved_output:
			if not options.paired_output:
				raise CommandLineError("When paired-end trimming is enabled via -A/-G/-B/-U, "
					"a second output file needs to be specified via -p (--paired-output).")
			if not options.output:
				raise CommandLineError("When you use -p or --paired-output, you must also "
					"use the -o option.")

		if bool(options.untrimmed_output) != bool(options.untrimmed_paired_output):
			raise CommandLineError("When trimming paired-end reads, you must use either none "
				"or both of the --untrimmed-output/--untrimmed-paired-output options.")
		if options.too_short_output and not options.too_short_paired_output:
			raise CommandLineError("When using --too-short-output with paired-end "
				"reads, you also need to use --too-short-paired-output")
		if options.too_long_output and not options.too_long_paired_output:
			raise CommandLineError("When using --too-long-output with paired-end "
				"reads, you also need to use --too-long-paired-output")
	elif quality_filename is not None:
		if options.format is not None:
			raise CommandLineError('If a pair of .fasta and .qual files is given, the -f/--format '
				'parameter cannot be used.')

	if options.format is not None and options.format.lower() not in ['fasta', 'fastq', 'sra-fastq']:
		raise CommandLineError("The input file format must be either 'fasta', 'fastq' or "
			"'sra-fastq' (not '{0}').".format(options.format))

	if options.maq:
		options.colorspace = True
		options.double_encode = True
		options.trim_primer = True
		options.strip_suffix.append('_F3')
		options.suffix = "/1"
	if options.zero_cap is None:
		options.zero_cap = options.colorspace
	if options.trim_primer and not options.colorspace:
		raise CommandLineError("Trimming the primer makes only sense in colorspace.")
	if options.double_encode and not options.colorspace:
		raise CommandLineError("Double-encoding makes only sense in colorspace.")
	if options.anywhere and options.colorspace:
		raise CommandLineError("Using --anywhere with colorspace reads is currently not supported "
			"(if you think this may be useful, contact the author).")
	if not (0 <= options.error_rate <= 1.):
		raise CommandLineError("The maximum error rate must be between 0 and 1.")
	if options.overlap < 1:
		raise CommandLineError("The overlap must be at least 1.")
	if not (0 <= options.gc_content <= 100):
		raise CommandLineError("GC content must be given as percentage between 0 and 100")
	if options.action == 'none':
		options.action = None

	if options.colorspace:
		if options.match_read_wildcards:
			raise CommandLineError('IUPAC wildcards not supported in colorspace')
		options.match_adapter_wildcards = False

	adapter_parser = AdapterParser(
		colorspace=options.colorspace,
		max_error_rate=options.error_rate,
		min_overlap=options.overlap,
		read_wildcards=options.match_read_wildcards,
		adapter_wildcards=options.match_adapter_wildcards,
		indels=options.indels)

	try:
		adapters = adapter_parser.parse_multi(options.adapters, options.anywhere, options.front)
		adapters2 = adapter_parser.parse_multi(options.adapters2, options.anywhere2, options.front2)
	except IOError as e:
		if e.errno == errno.ENOENT:
			raise CommandLineError(e)
		raise
	except ValueError as e:
		raise CommandLineError(e)
	if options.debug:
		for adapter in adapters + adapters2:
			adapter.enable_debug()

	# Create the processing pipeline.
	# If no second-read adapters were given (via -A/-G/-B/-U), we need to
	# be backwards compatible and *no modifications* are done to the second read.
	if paired:
		pipeline = PairedEndPipeline(pair_filter_mode, modify_first_read_only=paired == 'first')
	else:
		pipeline = SingleEndPipeline()

	if options.cut:
		if len(options.cut) > 2:
			raise CommandLineError("You cannot remove bases from more than two ends.")
		if len(options.cut) == 2 and options.cut[0] * options.cut[1] > 0:
			raise CommandLineError("You cannot remove bases from the same end twice.")
		for cut in options.cut:
			if cut != 0:
				pipeline.add1(UnconditionalCutter(cut))

	if options.cut2:
		if len(options.cut2) > 2:
			raise CommandLineError("You cannot remove bases from more than two ends.")
		if len(options.cut2) == 2 and options.cut2[0] * options.cut2[1] > 0:
			raise CommandLineError("You cannot remove bases from the same end twice.")
		for cut in options.cut2:
			if cut != 0:
				pipeline.add2(UnconditionalCutter(cut))

	if options.nextseq_trim is not None:
		pipeline.add(NextseqQualityTrimmer(options.nextseq_trim, options.quality_base))
	if options.quality_cutoff is not None:
		cutoffs = parse_cutoffs(options.quality_cutoff)
		pipeline.add(QualityTrimmer(cutoffs[0], cutoffs[1], options.quality_base))

	if adapters:
		adapter_cutter = AdapterCutter(adapters, options.times, options.action)
		pipeline.add1(adapter_cutter)
	if adapters2:
		adapter_cutter2 = AdapterCutter(adapters2, options.times, options.action)
		pipeline.add2(adapter_cutter2)

	# Modifiers that apply to both reads of paired-end reads unless in legacy mode
	if options.length is not None:
		pipeline.add(Shortener(options.length))
	if options.trim_n:
		pipeline.add(NEndTrimmer())
	if options.length_tag:
		pipeline.add(LengthTagModifier(options.length_tag))
	if options.strip_f3:
		options.strip_suffix.append('_F3')
	for suffix in options.strip_suffix:
		pipeline.add(SuffixRemover(suffix))
	if options.prefix or options.suffix:
		pipeline.add(PrefixSuffixAdder(options.prefix, options.suffix))
	if options.double_encode:
		pipeline.add(DoubleEncoder())
	if options.zero_cap:
		pipeline.add(ZeroCapper(quality_base=options.quality_base))
	if options.trim_primer:
		pipeline.add(PrimerTrimmer())

	# Set filtering parameters
	# Minimum/maximum length
	for attr in 'minimum_length', 'maximum_length':
		param = getattr(options, attr)
		if param is not None:
			lengths = parse_lengths(param)
			if not paired and len(lengths) == 2:
				raise CommandLineError('Two minimum or maximum lengths given for single-end data')
			if paired and len(lengths) == 1:
				lengths = (lengths[0], lengths[0])
			setattr(pipeline, attr, lengths)
	pipeline.max_n = options.max_n
	pipeline.discard_casava = options.discard_casava
	pipeline.discard_trimmed = options.discard_trimmed
	pipeline.discard_untrimmed = options.discard_untrimmed

	return pipeline


def main(cmdlineargs=None, default_outfile=sys.stdout):
	"""
	Main function that sets up a processing pipeline and runs it.

	default_outfile is the file to which trimmed reads are sent if the ``-o``
	parameter is not used.
	"""
	start_time = time.time()
	parser = get_option_parser()
	if cmdlineargs is None:
		cmdlineargs = sys.argv[1:]
	options, args = parser.parse_args(args=cmdlineargs)
	# Setup logging only if there are not already any handlers (can happen when
	# this function is being called externally such as from unit tests)
	if not logging.root.handlers:
		setup_logging(stdout=bool(options.output), quiet=options.quiet or options.report == 'minimal')
	if options.quiet and options.report:
		parser.error("Options --quiet and --report cannot be used at the same time")

	paired = determine_paired_mode(options)
	assert paired in (False, 'first', 'both')

	if paired == 'first':
		# legacy mode
		assert options.pair_filter is None
		pair_filter_mode = 'first'
	elif options.pair_filter is None:
		# default
		pair_filter_mode = 'any'
	else:
		# user-provided behavior
		pair_filter_mode = options.pair_filter

	try:
		is_interleaved_input, is_interleaved_output = determine_interleaved(options, args)
		input_filename, input_paired_filename, quality_filename = input_files_from_parsed_args(args,
			paired, is_interleaved_input)
		pipeline = pipeline_from_parsed_args(options, paired, pair_filter_mode, quality_filename, is_interleaved_output)
		outfiles = open_output_files(options, default_outfile, is_interleaved_output)
	except CommandLineError as e:
		parser.error(e)
		return  # avoid IDE warnings below

	if options.cores < 0:
		parser.error('Value for --cores cannot be negative')
	cores = available_cpu_count() if options.cores == 0 else options.cores
	if cores > 1:
		if (
			PY3
			and ParallelPipelineRunner.can_output_to(outfiles)
			and quality_filename is None
			and not options.colorspace
			and options.format is None
		):
			runner = ParallelPipelineRunner(pipeline, cores, options.buffer_size)
		else:
			if not PY3:
				logger.error('Running in parallel is not supported on Python 2')
			else:
				logger.error('Running in parallel is currently not supported for '
					'the given combination of command-line parameters.\nThese '
					'options are not supported: --info-file, --rest-file, '
					'--wildcard-file, --untrimmed-output, '
					'--untrimmed-paired-output, --too-short-output, '
					'--too-short-paired-output, --too-long-output, '
					'--too-long-paired-output, --format, --colorspace')
			sys.exit(1)
	else:
		runner = pipeline
	try:
		runner.set_input(input_filename, file2=input_paired_filename,
			qualfile=quality_filename, colorspace=options.colorspace,
			fileformat=options.format, interleaved=is_interleaved_input)
		runner.set_output(outfiles)
	except (seqio.UnknownFileType, IOError) as e:
		parser.error(e)

	implementation = platform.python_implementation()
	opt = ' (' + implementation + ')' if implementation != 'CPython' else ''
	logger.info("This is cutadapt %s with Python %s%s", __version__,
		platform.python_version(), opt)
	logger.info("Command line parameters: %s", " ".join(cmdlineargs))
	logger.info("Processing reads on %d core%s in %s mode ...",
		cores, 's' if cores > 1 else '',
		{False: 'single-end', 'first': 'paired-end legacy', 'both': 'paired-end'}[pipeline.paired])

	if pipeline.should_warn_legacy:
		logger.warning('\n'.join(textwrap.wrap('Legacy mode is '
			'enabled. Read modification and filtering options *ignore* '
			'the second read. To switch to regular paired-end mode, '
			'provide the --pair-filter=any option or use any of the '
			'-A/-B/-G/-U/--interleaved options.')))

	try:
		stats = runner.run()
		# cProfile.runctx('stats=runner.run()', globals(), locals(), 'profile_main.prof')
		runner.close()
	except KeyboardInterrupt:
		print("Interrupted", file=sys.stderr)
		sys.exit(130)
	except IOError as e:
		if e.errno == errno.EPIPE:
			sys.exit(1)
		raise
	except (seqio.FormatError, seqio.UnknownFileType, EOFError) as e:
		sys.exit("cutadapt: error: {0}".format(e))

	elapsed = time.time() - start_time
	if not options.quiet:
		# send statistics to stderr if result was sent to stdout
		stat_file = sys.stderr if options.output is None else None
		with redirect_standard_output(stat_file):
			if options.report == 'minimal':
				print_minimal_report(stats, elapsed, options.gc_content / 100)
			else:
				print_report(stats, elapsed, options.gc_content / 100)


if __name__ == '__main__':
	main()
