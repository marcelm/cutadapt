#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kate: word-wrap off; remove-trailing-spaces all;
#
# Copyright (c) 2010-2016 Marcel Martin <marcel.martin@scilifelab.se>
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
Copyright (C) 2010-2016 Marcel Martin <marcel.martin@scilifelab.se>

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

Use "cutadapt --help" to see all command-line options.
See http://cutadapt.readthedocs.org/ for full documentation.
"""

# This is a fork of Cutadapt that can use multiple threads to process
# reads in parallel.
# This program incorporates code/ideas from the following sources:
# TrimGalore (http://www.bioinformatics.babraham.ac.uk/projects/download.html#trim_galore)
# miRge (https://github.com/BarasLab/miRge/blob/master/trim_file.py)

# TODO: support reading from a pair of files but writing to interleaved file (including stdout)
# TODO: finish adding defaults from TrimGalore
# TODO: implement an autodetect option for adapters similar to TrimGalore? - read first
# 1M reads, search for common adapters, and pick the one that appears most frequently.

from __future__ import print_function, division, absolute_import

# Print a helpful error message if the extension modules cannot be imported.
from cutadapt import *
check_importability()

from cutadapt import __version__
from cutadapt.scripts import cutadapt as cutadapt_script
from cutadapt import seqio
from cutadapt.xopen import xopen
from cutadapt.adapters import AdapterParser
from cutadapt.modifiers import ModType, create_modifier
from cutadapt.filters import FilterType, FilterFactory
from cutadapt.report import *
from cutadapt.compat import next

from collections import OrderedDict
from heapq import nsmallest
import io
import logging
from multiprocessing import Process, JoinableQueue, Value, cpu_count
from optparse import OptionParser, OptionGroup
import os
import platform
from queue import Empty, Full
import sys
import time

MAGNITUDE = dict(
	G=(1E9, "000000000"),
	M=(1E6, "000000"),
	K=(1E3, "000")
)

#from memory_profiler import profile
#fp=open('memory_profiler.log','w+')
#@profile(stream=fp)
def main(cmdlineargs=None, default_outfile="-"):
	"""
	Main function that evaluates command-line parameters and iterates
	over all reads.
	
	default_outfile is the file to which trimmed reads are sent if the ``-o``
	parameter is not used.
	"""
	# For CPU profiling
	#import yappi
	#yappi.start()
	
	parser = get_option_parser()
	
	if cmdlineargs is None:
		cmdlineargs = sys.argv[1:]
	
	options, args = parser.parse_args(args=cmdlineargs)
	
	# Setup logging only if there are not already any handlers (can happen when
	# this function is being called externally such as from unit tests)
	if not logging.root.handlers:
		setup_logging(stdout=bool(options.output), quiet=options.quiet)
	
	paired, multiplexed = validate_options(options, args, parser)
	reader, qualities, has_qual_file = create_reader(args, options, parser)
	writers = Writers(options, multiplexed, qualities, default_outfile)
	modifiers, num_adapters = create_modifiers(options, paired, qualities, has_qual_file, parser)
	min_affected = 2 if options.pair_filter == 'both' else 1
	filters = create_filters(options, paired, min_affected)
	
	logger = logging.getLogger()
	logger.info("This is cutadapt %s with Python %s", __version__, platform.python_version())
	logger.info("Command line parameters: %s", " ".join(cmdlineargs))
	logger.info("Trimming %s adapter%s with at most %.1f%% errors in %s mode ...",
		num_adapters, 's' if num_adapters > 1 else '', options.error_rate * 100,
		{ False: 'single-end', 'first': 'paired-end legacy', 'both': 'paired-end' }[paired])
	if paired == 'first' and (len(modifiers.mod2) > 0 or options.quality_cutoff):
		import textwrap
		logger.warning('\n'.join(textwrap.wrap('WARNING: Requested read '
			'modifications are applied only to the first '
			'read since backwards compatibility mode is enabled. '
			'To modify both reads, also use any of the -A/-B/-G/-U options. '
			'Use a dummy adapter sequence when necessary: -A XXX')))
	
	start_time = time.clock()
	
	if options.threads is None:
		# Run cutadapt normally
		summary = run_cutadapt_serial(reader, writers, modifiers, filters)
	else:
		# Run multiprocessing version
		summary = run_cutadapt_parallel(reader, writers, modifiers, filters,
			options.threads, options.thread_timeout, options.preserve_order, 
			options.read_queue_size, options.result_queue_size)

	report = print_report(paired, options, time.clock() - start_time, summary)
	
	#yappi.stop()
	#yappi.get_func_stats().print_all()
	#yappi.get_thread_stats().print_all()
	
def get_option_parser():
	# TODO: upgrade to argparse
	class CutadaptOptionParser(OptionParser):
		def get_usage(self):
			return self.usage.lstrip().replace('%version', __version__)
	
	parser = CutadaptOptionParser(usage=__doc__, version=__version__)
	
	parser.add_option("--debug", action='store_true', default=False,
		help="Print debugging information.")
	parser.add_option("-f", "--format",
		help="Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. "
			"Ignored when reading csfasta/qual files. Default: auto-detect "
			"from file name extension.")
	parser.add_option("--max-reads", default=None,
		help="Maximum number of reads/pairs to process")
	parser.add_option("--no-progress", action="store_true", default=False,
		help="Don't show a progress bar")
	
	group = OptionGroup(parser, "Finding adapters:",
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
	group.add_option("-e", "--error-rate", type=float, default=None,
		help="Maximum allowed error rate (no. of errors divided by the length "
			"of the matching region). Default: %default")
	group.add_option("--no-indels", action='store_false', dest='indels', default=True,
		help="Allow only mismatches in alignments. "
			"Default: allow both mismatches and indels")
	group.add_option("-n", "--times", type=int, metavar="COUNT", default=1,
		help="Remove up to COUNT adapters from each read. Default: %default")
	group.add_option("-O", "--overlap", type=int, metavar="MINLENGTH", default=3,
		help="If the overlap between the read and the adapter is shorter than "
			"MINLENGTH, the read is not modified. Reduces the no. of bases "
			"trimmed due to random adapter matches. Default: %default")
	group.add_option("--match-read-wildcards", action="store_true", default=False,
		help="Interpret IUPAC wildcards in reads. Default: %default")
	group.add_option("-N", "--no-match-adapter-wildcards", action="store_false",
		default=True, dest='match_adapter_wildcards',
		help="Do not interpret IUPAC wildcards in adapters.")
	group.add_option("--no-trim", dest='action', action='store_const', const=None,
		help="Match and redirect reads to output/untrimmed-output as usual, "
			"but do not remove adapters.")
	group.add_option("--mask-adapter", dest='action', action='store_const', const='mask',
		help="Mask adapters with 'N' characters instead of trimming them.")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Additional read modifications")
	group.add_option("-u", "--cut", action='append', default=[], type=int, metavar="LENGTH",
		help="Remove bases from each read (first read only if paired). "
			"If LENGTH is positive, remove bases from the beginning. "
			"If LENGTH is negative, remove bases from the end. "
			"Can be used twice if LENGTHs have different signs.")
	group.add_option("-q", "--quality-cutoff", default=None, metavar="[5'CUTOFF,]3'CUTOFF",
		help="Trim low-quality bases from 5' and/or 3' ends of each read before "
			"adapter removal. Applied to both reads if data is paired. If one "
			"value is given, only the 3' end is trimmed. If two "
			"comma-separated cutoffs are given, the 5' end is trimmed with "
			"the first cutoff, the 3' end with the second.")
	group.add_option("--nextseq-trim", type=int, default=None, metavar="3'CUTOFF",
		help="NextSeq-specific quality trimming (each read). Trims also dark "
			"cycles appearing as high-quality G bases (EXPERIMENTAL).")
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
	
	group = OptionGroup(parser, "Filtering of processed reads")
	group.add_option("--discard-trimmed", "--discard", action='store_true', default=False,
		help="Discard reads that contain an adapter. Also use -O to avoid "
			"discarding too many randomly matching reads!")
	group.add_option("--discard-untrimmed", "--trimmed-only", action='store_true', default=False,
		help="Discard reads that do not contain the adapter.")
	group.add_option("-m", "--minimum-length", type=int, default=None, metavar="LENGTH",
		help="Discard trimmed reads that are shorter than LENGTH. Reads that "
			"are too short even before adapter removal are also discarded. In "
			"colorspace, an initial primer is not counted. Default: 0")
	group.add_option("-M", "--maximum-length", type=int, default=sys.maxsize, metavar="LENGTH",
		help="Discard trimmed reads that are longer than LENGTH. "
			"Reads that are too long even before adapter removal "
			"are also discarded. In colorspace, an initial primer "
			"is not counted. Default: no limit")
	group.add_option("--max-n", type=float, default=-1.0, metavar="COUNT",
		help="Discard reads with too many N bases. If COUNT is an integer, it "
			"is treated as the absolute number of N bases. If it is between 0 "
			"and 1, it is treated as the proportion of N's allowed in a read.")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Output")
	group.add_option("--quiet", default=False, action='store_true',
		help="Print only error messages.")
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
			"rest (after the adapter) into FILE.")
	group.add_option("--wildcard-file", metavar="FILE",
		help="When the adapter has N bases (wildcards), write adapter bases "
			"matching wildcard positions to FILE. When there are indels in the "
			"alignment, this will often not be accurate.")
	group.add_option("--too-short-output", metavar="FILE",
		help="Write reads that are too short (according to length specified by "
		"-m) to FILE. Default: discard reads")
	group.add_option("--too-long-output", metavar="FILE",
		help="Write reads that are too long (according to length specified by "
		"-M) to FILE. Default: discard reads")
	group.add_option("--untrimmed-output", default=None, metavar="FILE",
		help="Write reads that do not contain the adapter to FILE. Default: "
			"output to same file as trimmed reads")
	group.add_option("--report-file", default=None, metavar="FILE",
		help="Write report to file rather than stdout/stderr.")
	group.add_option("--file-buffer-size", type=int, default=1024*64, metavar="SIZE",
		help="Buffer size to use for file I/O; a larger buffer can potentially improve "
			"speed by reducing the number of I/O operations.")
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
			"data. By default, they are since many tools have problems with "
			"negative qualities.")
	group.add_option("--zero-cap", "-z", action='store_true',
		help="Change negative quality values to zero. This is enabled "
		"by default when -c/--colorspace is also enabled. Use the above option "
		"to disable it.")
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
		help="Remove LENGTH bases from second read in a pair (see --cut).")
	group.add_option("-p", "--paired-output", metavar="FILE",
		help="Write second read in a pair to FILE.")
	# Setting the default for pair_filter to None allows us to find out whether
	# the option was used at all.
	group.add_option("--pair-filter", metavar='(any|both)', default=None,
		choices=("any", "both"),
		help="Which of the reads in a paired-end read have to match the "
			"filtering criterion in order for it to be filtered. "
			"Default: any")
	group.add_option("--interleaved", action='store_true', default=False,
		help="Read and write interleaved paired-end reads.")
	group.add_option("--untrimmed-paired-output", metavar="FILE",
		help="Write second read in a pair to this FILE when no adapter "
			"was found in the first read. Use this option together with "
			"--untrimmed-output when trimming paired-end reads. Default: output "
			"to same file as trimmed reads")
	group.add_option("--too-short-paired-output", metavar="FILE", default=None,
		help="Write second read in a pair to this file if pair is too short. "
			"Use together with --too-short-output.")
	group.add_option("--too-long-paired-output", metavar="FILE", default=None,
		help="Write second read in a pair to this file if pair is too long. "
			"Use together with --too-long-output.")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Parallel options")
	group.add_option("--threads", type=int, default=None, metavar="THREADS",
		help="Number of threads to use for read trimming. Set to 0 to use max available threads.")
	group.add_option("--batch-size", type=int, default=1000, metavar="SIZE",
		help="Number of records to process in each batch.")
	group.add_option("--preserve-order", action="store_true", default=False,
		help="Preserve order of reads in input files.")
	group.add_option("--thread-timeout", type=int, default=60, metavar="SECONDS",
		help="Number of seconds threads should wait before exiting if task queue is empty.")
	group.add_option("--read-queue-size", type=int, default=None, metavar="SIZE",
		help="Size of queue for batches of reads to be processed.")
	group.add_option("--result-queue-size", type=int, default=None, metavar="SIZE",
		help="Size of queue for batches of results to be written.")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Method-specific options")
	group.add_option("--bisulfite", action="store_true", default=False,
		help="Set default option values for bisulfite-treated data.")
	group.add_option("--non-directional", action="store_true", default=False,
		help="Bisulfite sequencing libraries are non-directional.")
	group.add_option("--mirna", action="store_true", default=False,
		help="Set default option values for miRNA data.")
	parser.add_option_group(group)
	
	return parser

def setup_logging(stdout=False, quiet=False):
	"""
	Attach handler to the global logger object
	"""
	# Due to backwards compatibility, logging output is sent to standard output
	# instead of standard error if the -o option is used.
	stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
	stream_handler.setFormatter(logging.Formatter('%(message)s'))
	stream_handler.setLevel(logging.ERROR if quiet else logging.INFO)
	logging.getLogger().setLevel(logging.INFO)
	logging.getLogger().addHandler(stream_handler)

def validate_options(options, args, parser):
	if options.mirna:
		if options.adapter is None and options.front is None and options.anywhere is None:
			options.adapter = ['TGGAATTCTCGG'] # illumina small RNA adapter
		if options.quality_cutoff is None:
			options.quality_cutoff = "20"
		if options.minimum_length is None:
			options.minimum_length = 16
		if options.error_rate is None:
			options.error_rate = 0.12
	elif options.bisulfite and options.quality_cutoff is None:
		options.quality_cutoff = str((options.quality_base - 33) + 20)
	
	if options.debug and options.threads is not None:
		parser.error("Cannot use debug mode with multiple threads")
		
	if len(args) == 0:
		parser.error("At least one parameter needed: name of a FASTA or FASTQ file.")
	elif len(args) > 2:
		parser.error("Too many parameters.")
	if args[0].endswith('.qual'):
		parser.error("If a .qual file is given, it must be the second argument.")
	
	# Find out which 'mode' we need to use.
	# Default: single-read trimming (neither -p nor -A/-G/-B/-U/--interleaved given)
	paired = False
	if options.paired_output:
		# Modify first read only, keep second in sync (-p given, but not -A/-G/-B/-U).
		# This exists for backwards compatibility ('legacy mode').
		paired = 'first'
	# Any of these options switch off legacy mode
	if (options.adapters2 or options.front2 or options.anywhere2 or
		options.cut2 or options.interleaved or options.pair_filter or
		options.too_short_paired_output or options.too_long_paired_output):
		# Full paired-end trimming when both -p and -A/-G/-B/-U given
		# Read modifications (such as quality trimming) are applied also to second read.
		paired = 'both'

	if paired and len(args) == 1 and not options.interleaved:
		parser.error("When paired-end trimming is enabled via -A/-G/-B/-U or -p, "
			"two input files are required.")
	if options.interleaved and len(args) != 1:
		parser.error("When reading interleaved files, only one input file may "
			"be given.")
	if not paired:
		if options.untrimmed_paired_output:
			parser.error("Option --untrimmed-paired-output can only be used when "
				"trimming paired-end reads (with option -p).")

	if paired:
		if not options.interleaved:
			if not options.paired_output:
				parser.error("When paired-end trimming is enabled via -A/-G/-B/-U, "
					"a second output file needs to be specified via -p (--paired-output).")
			if not options.output:
				parser.error("When you use -p or --paired-output, you must also "
					"use the -o option.")
			if bool(options.untrimmed_output) != bool(options.untrimmed_paired_output):
				parser.error("When trimming paired-end reads, you must use either none "
					"or both of the --untrimmed-output/--untrimmed-paired-output options.")
			if options.too_short_output and not options.too_short_paired_output:
				parser.error("When using --too-short-output with paired-end "
					"reads, you also need to use --too-short-paired-output")
			if options.too_long_output and not options.too_long_paired_output:
				parser.error("When using --too-long-output with paired-end "
					"reads, you also need to use --too-long-paired-output")
	elif len(args) == 2:
		if options.format is not None:
			parser.error("If a pair of .fasta and .qual files is given, the -f/--format "
				"parameter cannot be used.")

	if options.format is not None and options.format.lower() not in ['fasta', 'fastq', 'sra-fastq']:
		parser.error("The input file format must be either 'fasta', 'fastq' or "
			"'sra-fastq' (not '{0}').".format(options.format))

	if options.quality_cutoff is not None:
		cutoffs = options.quality_cutoff.split(',')
		if len(cutoffs) == 1:
			try:
				options.quality_cutoff = [0, int(cutoffs[0])]
			except ValueError as e:
				parser.error("Quality cutoff value not recognized: {0}".format(e))
		elif len(cutoffs) == 2:
			try:
				options.quality_cutoff = [int(cutoffs[0]), int(cutoffs[1])]
			except ValueError as e:
				parser.error("Quality cutoff value not recognized: {0}".format(e))
		else:
			parser.error("Expected one value or two values separated by comma for the quality cutoff")

	if options.pair_filter is None:
		options.pair_filter = 'any'

	if (options.discard_trimmed or options.discard_untrimmed) and (options.untrimmed_output is not None):
		parser.error("Only one of the --discard-trimmed, --discard-untrimmed "
			"and --untrimmed-output options can be used at the same time.")
	
	multiplexed = False
	if options.output is not None and '{name}' in options.output:
		multiplexed = True
		if options.discard_trimmed:
			parser.error("Do not use --discard-trimmed when demultiplexing.")
		if paired:
			parser.error("Demultiplexing not supported for paired-end files, yet.")
	
	if options.maq:
		options.colorspace = True
		options.double_encode = True
		options.trim_primer = True
		options.suffix = "/1"
	
	if options.strip_f3 or options.maq:
		options.strip_suffix.append('_F3')
	
	if options.zero_cap is None:
		options.zero_cap = options.colorspace
	
	if options.trim_primer and not options.colorspace:
		parser.error("Trimming the primer makes only sense in colorspace.")
	
	if options.double_encode and not options.colorspace:
		parser.error("Double-encoding makes only sense in colorspace.")
	
	if options.anywhere and options.colorspace:
		parser.error("Using --anywhere with colorspace reads is currently not supported (if you "
			"think this may be useful, contact the author).")
	
	if options.error_rate is None:
		options.error_rate = 0.1
	elif not (0 <= options.error_rate <= 1.):
		parser.error("The maximum error rate must be between 0 and 1.")
	
	if options.overlap < 1:
		parser.error("The overlap must be at least 1.")

	if options.cut:
		if len(options.cut) > 2:
			parser.error("You cannot remove bases from more than two ends.")
		if len(options.cut) == 2 and options.cut[0] * options.cut[1] > 0:
			parser.error("You cannot remove bases from the same end twice.")

	if paired == 'both' and options.cut2:
		if len(options.cut2) > 2:
			parser.error("You cannot remove bases from more than two ends.")
		if len(options.cut2) == 2 and options.cut2[0] * options.cut2[1] > 0:
			parser.error("You cannot remove bases from the same end twice.")

	if options.colorspace:
		if options.match_read_wildcards:
			parser.error('IUPAC wildcards not supported in colorspace')
		options.match_adapter_wildcards = False
	
	if options.quiet or options.output is None or options.output == "-":
		options.no_progress = True
	elif not options.no_progress:
		try:
			import progressbar
		except:
			options.no_progress = True
	
	# TODO: once we switch to argparse, int_or_str can be passed
	# as the argument type
	def int_or_str(x):
		if x is None or isinstance(x, int):
			return x
		elif isinstance(x, str):
			x = x.upper()
			for a, mag in MAGNITUDE.items():
				x = x.replace(a, mag[1])
			return int(x)
		else:
			raise Exception("Unsupported type {}".format(x))
	
	options.max_reads = int_or_str(options.max_reads)
	
	if options.threads is not None:
		if options.threads <= 0:
			options.threads = cpu_count() - 2
			if options.threads <= 0:
				# TODO: warn user
				options.threads = 1
	
		options.batch_size = int_or_str(options.batch_size)
		
		if options.thread_timeout is None:
			# set timeout based on batch size
			options.thread_timeout = max(60, int(options.batch_size * 0.05))
		
		if options.file_buffer_size < io.DEFAULT_BUFFER_SIZE:
			parser.error("File buffer size must be at least {}".format(io.DEFAULT_BUFFER_SIZE))
		
		default_queue_size = options.threads * 10
		
		if options.read_queue_size is None:
			options.read_queue_size = default_queue_size
		elif options.read_queue_size > 0:
			assert options.read_queue_size >= options.threads
	
		if options.result_queue_size is None:
			options.result_queue_size = default_queue_size
		elif options.result_queue_size > 0:
			assert options.result_queue_size > options.threads
	
	return (paired, multiplexed)

def create_reader(input_files, options, parser, counter_magnitude="M"):
	input_filename = input_files[0]
	input_paired_filename = None
	quality_filename = None
	if len(input_files) > 1:
		if paired:
			input_paired_filename = input_files[1]
		else:
			quality_filename = input_files[1]
	
	try:
		reader = seqio.open(input_filename, file2=input_paired_filename,
			qualfile=quality_filename, colorspace=options.colorspace, 
			fileformat=options.format, interleaved=options.interleaved)
	except (seqio.UnknownFileType, IOError) as e:
		parser.error(e)
	
	qualities = reader.delivers_qualities
	
	# Wrap reader in batch iterator
	batch_size = options.batch_size or 1000
	reader = BatchIterator(reader, batch_size, options.max_reads)
	
	# Wrap iterator in progress bar, unless requested not to
	if not options.no_progress:
		import progressbar
		import progressbar.widgets
		import math

		class BatchProgressBar(progressbar.ProgressBar):
			def __init__(self, iterable, widgets, max_value=None):
				super(BatchProgressBar, self).__init__(widgets=widgets, max_value=max_value)
				self._iterable = iterable
			
			def __next__(self):
				try:
					value = next(self._iterable)
					if self.start_time is None:
						self.start()
					self.update(self.value + value[0])
					return value
				except StopIteration:
					self.finish()
					raise
			
			def close(self):
				self._iterable.close()
		
		class MagCounter(progressbar.widgets.WidgetBase):
			def __init__(self, magnitude):
				suffix = ""
				if magnitude is None:
					div = 1.0
				else:
					div = float(MAGNITUDE[magnitude][0])
					suffix = magnitude
			
				self._format = lambda val: "{:.1f} {}".format(val / div, suffix)
		
			def __call__(self, progress, data):
				return self._format(data["value"])
			
		if options.max_reads:
			reader = BatchProgressBar(reader, [
				MagCounter(counter_magnitude), " Reads (", progressbar.Percentage(), ") ", 
				progressbar.Timer(), " ", progressbar.Bar(), progressbar.AdaptiveETA()
			], options.max_reads)
		else:
			reader = BatchProgressBar(reader, [
				MagCounter(counter_magnitude), " Reads", progressbar.Timer(), 
				progressbar.AnimatedMarker()
			])
	
	return (reader, qualities, quality_filename is not None)

class BatchIterator(object):
	def __init__(self, reader, size, max_reads=None):
		self.reader = enumerate(reader, 1)
		self.size = size
		self.max_reads = max_reads
		self.done = False
		self._empty_batch = [None] * size
	
	def __iter__(self):
		return self
	
	def __next__(self):
		if self.done:
			raise StopIteration()
		
		try:
			read_index, record = next(self.reader)
		except StopIteration:
			self.done = True
			raise
		
		batch = self._empty_batch.copy()
		batch[0] = record
		batch_index = 1
		max_size = self.size
		if self.max_reads:
			max_size = min(max_size, self.max_reads - read_index)
		
		while batch_index < max_size:
			try:
				read_index, record = next(self.reader)
				batch[batch_index] = record
				batch_index += 1
			except StopIteration:
				self.done = True
				break
		
		if self.max_reads and read_index >= self.max_reads:
			self.done = True
		
		if batch_index == self.size:
			return (batch_index, batch)
		else:
			return (batch_index, batch[0:batch_index])
	
	# py2x alias
	next = __next__

class Filters(object):
	def __init__(self, filter_factory):
		self.filters = OrderedDict()
		self.filter_factory = filter_factory
	
	def add_filter(self, filter_type, *args, **kwargs):
		self.filters[filter_type] = self.filter_factory(filter_type, *args, **kwargs)
	
	def filter(self, read1, read2=None):
		dest = FilterType.NONE
		for filter_type, f in self.filters.items():
			if f(read1, read2):
				dest = filter_type
				# Stop writing as soon as one of the filters was successful.
				break
		return dest
	
	def __contains__(self, filter_type):
		return filter_type in self.filters
	
	def __getitem__(self, filter_type):
		return self.filters[filter_type]

def create_filters(options, paired, min_affected):
	filters = Filters(FilterFactory(paired, min_affected))
	
	if options.minimum_length is not None and options.minimum_length > 0:
		filters.add_filter(FilterType.TOO_SHORT, options.minimum_length)
	
	if options.maximum_length < sys.maxsize:
		filters.add_filter(FilterType.TOO_LONG, options.maximum_length)
	
	if options.max_n >= 0:
		filters.add_filter(FilterType.N_CONTENT, options.max_n)
	
	if options.discard_trimmed:
		filters.add_filter(FilterType.TRIMMED)
	
	if options.discard_untrimmed or options.untrimmed_output:
		filters.add_filter(FilterType.UNTRIMMED)
	
	return filters

class Modifiers(object):
	def __init__(self, paired):
		self.mod1 = OrderedDict()
		self.mod2 = OrderedDict()
		self.paired = paired
	
	def add_modifier_pair(self, mod_type, *args, **kwargs):
		args1 = [a[0] for a in args]
		if not any(a is None for a in args1):
			self.mod1[mod_type] = create_modifier(mod_type, *args1, **kwargs)
		
		if self.paired == "both":
			args2 = [a[1] for a in args]
			if not any(a is None for a in args2):
				self.mod2[mod_type] = create_modifier(mod_type, *args2, **kwargs)
	
	def add_modifier(self, mod_type, *args, _first_only=False, **kwargs):
		self.mod1[mod_type] = create_modifier(mod_type, *args, **kwargs)
		if self.paired == "both" and not _first_only:
			self.mod2[mod_type] = create_modifier(mod_type, *args, **kwargs)
	
	def get_modifier_pair(self, mod_type):
		return [self.mod1.get(mod_type, None), self.mod2.get(mod_type, None)]
	
	def modify(self, record):
		bp = [0,0]
		if self.paired:
			read1, read2 = record
			bp[0] = len(read1.sequence)
			bp[1] = len(read2.sequence)
			for mod in self.mod1.values():
				read1 = mod(read1)
			for mod in self.mod2.values():
				read2 = mod(read2)
			reads = [read1, read2]
		else:
			read = record
			bp[0] = len(read.sequence)
			for mod in self.mod1.values():
				read = mod(read)
			reads = [read]
		return (reads, bp)
	
	def __getitem__(self, i):
		if i == 0:
			return self.mod1
		else:
			return self.mod2

def create_modifiers(options, paired, qualities, has_qual_file, parser):
	adapter_parser = AdapterParser(
		colorspace=options.colorspace,
		max_error_rate=options.error_rate,
		min_overlap=options.overlap,
		read_wildcards=options.match_read_wildcards,
		adapter_wildcards=options.match_adapter_wildcards,
		indels=options.indels)

	try:
		adapters1 = adapter_parser.parse_multi(options.adapters, options.anywhere, options.front)
		adapters2 = adapter_parser.parse_multi(options.adapters2, options.anywhere2, options.front2)
	except IOError as e:
		if e.errno == errno.ENOENT:
			parser.error(e)
		raise
	except ValueError as e:
		parser.error(e)
	
	if not adapters1 and not adapters2 and not options.quality_cutoff and \
			options.nextseq_trim is None and \
			options.cut == [] and options.cut2 == [] and \
			(options.minimum_length is None or options.minimum_length <= 0) and \
			options.maximum_length == sys.maxsize and \
			not has_qual_file and \
			options.max_n == -1 and not options.trim_n:
		parser.error("You need to provide at least one adapter sequence.")
	
	if options.debug:
		for adapter in adapters1 + adapters2:
			adapter.enable_debug()
	
	modifiers = Modifiers(paired)
	
	if options.cut or options.cut2:
		modifiers.add_modifier_pair(ModType.CUT, (options.cut, options.cut2))
	
	if options.nextseq_trim is not None:
		modifiers.add_modifier(ModType.TRIM_NEXTSEQ_QUAL, options.nextseq_trim, 
			options.quality_base, _first_only=True)
	
	if options.quality_cutoff:
		modifiers.add_modifier(ModType.TRIM_QUAL, options.quality_cutoff[0], 
			options.quality_cutoff[1], options.quality_base)
	
	if adapters1 or adapters2:
		modifiers.add_modifier_pair(ModType.ADAPTER, (adapters1, adapters2), 
			times=options.times, action=options.action)
	
	if options.bisulfite:
		# For non-directional RRBS/WGBS libraries (which implies that they were digested
		# using MspI), sequences starting with either 'CAA' or 'CGA' will have 2 bp 
		# trimmed off either end to remove potential methylation-biased bases from the 
		# end-repair reaction. 
		# For all RRBS/WGBS libraries, sequences that are adapter trimmed and are either
		# directional or do not start with CAA/CGA, 2 bp are are removed from the 3'
		# end to remove potential methylation-biased bases from the end-repair reaction.
		# TODO: log messge?
		modifiers.add_modifier(ModType.BISULFITE, non_directional=options.non_directional)
	
	if options.trim_n:
		modifiers.add_modifier(ModType.TRIM_END_N)
	
	if options.length_tag:
		modifiers.add_modifier(ModType.LENGTH_TAG, options.length_tag)
	
	if options.strip_suffix:
		modifiers.add_modifier(ModType.REMOVE_SUFFIX, options.strip_suffix)
	
	if options.prefix or options.suffix:
		modifiers.add_modifier(ModType.ADD_PREFIX_SUFFIX, options.prefix, options.suffix)
	
	if options.double_encode:
		modifiers.add_modifier(ModType.CS_DOUBLE_ENCODE)
	
	if options.zero_cap and qualities:
		modifiers.add_modifier(ModType.ZERO_CAP, quality_base=options.quality_base)
	
	if options.trim_primer:
		modifiers.add_modifier(ModType.CS_TRIM_PRIMER)
	
	return (modifiers, len(adapters1) + len(adapters2))

class SingleEndWriter(object):
	def __init__(self, handle):
		self.handle = handle
		self.written = 0
		self.read1_bp = 0
		self.read2_bp = 0
	
	def write(self, read1, read2=None):
		assert read2 is None
		self.handle.write(read1)
		self.written += 1
		self.read1_bp += len(read1)
	
	@property
	def written_bp(self):
		return (self.read1_bp, self.read2_bp)
	
	def close(self):
		if self.handle is not sys.stdout and self.handle is not sys.stderr:
			self.handle.close()

class PairedEndWriter(SingleEndWriter):
	def write(self, read1, read2):
		self.handle.write(read1, read2)
		self.written += 1
		self.read1_bp += len(read1)
		self.read2_bp += len(read2)

class Writers(object):
	def __init__(self, options, multiplexed, qualities, default_outfile):
		self.multiplexed = multiplexed
		self.output = options.output
		self.seqio_open_args = dict(
			qualities=qualities,
			colorspace=options.colorspace, 
			interleaved=options.interleaved,
			buffer_size=options.file_buffer_size
		)
		self.seqfile_paths = {}
		self.force_create = {}
		self.rest_path = options.rest_file
		self.info_path = options.info_file
		self.wildcard_path = options.wildcard_file
		self.writers = {}
		self._rest_writer = None
		self.discarded = 0
		self.buffer_size = options.file_buffer_size
		
		if (options.minimum_length is not None 
				and options.minimum_length > 0 
				and options.too_short_output):
			self.add_seq_writer(FilterType.TOO_SHORT,
				options.too_short_output, options.too_short_paired_output)
	
		if options.maximum_length < sys.maxsize and options.too_long_output is not None:
			self.add_seq_writer(FilterType.TOO_LONG,
				options.too_long_output, options.too_long_paired_output)
		
		if not self.multiplexed:
			if options.output is not None:
				self.add_seq_writer(FilterType.NONE, options.output, 
					options.paired_output, force_create=True)
			elif not (options.discard_trimmed and options.untrimmed_output):
				self.add_seq_writer(FilterType.NONE, default_outfile, 
					force_create=True)
		
		if not options.discard_untrimmed:
			if self.multiplexed:
				untrimmed = options.untrimmed_output or options.output.format(name='unknown')
				self.add_seq_writer(FilterType.UNTRIMMED, untrimmed)
				self.add_seq_writer(FilterType.NONE, untrimmed)
			elif options.untrimmed_output:
				self.add_seq_writer(FilterType.UNTRIMMED,
					options.untrimmed_output, options.untrimmed_paired_output)
	
	def add_seq_writer(self, filter_type, file1, file2=None, force_create=False):
		self.seqfile_paths[filter_type] = (file1, file2)
		if force_create and isinstance(file1, str) and file1 != "-":
			self.force_create[file1] = False
			if file2 is not None:
				self.force_create[file2] = False
	
	def has_seq_writer(self, filter_type):
		return filter_type in self.seqfile_paths
	
	def get_seq_writer(self, filter_type):
		paths = self.seqfile_paths[filter_type]
		if paths not in self.writers:
			self.writers[paths] = self._create_seq_writer(*paths)
		return self.writers[paths]
	
	def get_mux_writer(self, name):
		assert self.multiplexed
		path = self.output.format(name=name)
		if path not in self.writers:
			self.writers[path] = self._create_seq_writer(path)
		return self.writers[path]
		
	def _create_seq_writer(self, file1, file2=None):
		seqfile = seqio.open(file1, file2, mode='w', **self.seqio_open_args)
		if file1 in self.force_create:
			self.force_create[file1] = True
		if file2 is not None and file2 in self.force_create:
			self.force_create[file2] = True
		if isinstance(seqfile, seqio.SingleRecordWriter):
			return SingleEndWriter(seqfile)
		elif isinstance(seqfile, seqio.PairRecordWriter):
			return PairedEndWriter(seqfile)
		else:
			raise Exception("Unrecognized type of writer {}".format(writer.__class__))
	
	def get_seq_writers(self):
		seq_writers = set()
		for paths in self.seqfile_paths.values():
			if paths in self.writers:
				seq_writers.add(self.writers[paths])
		return seq_writers
	
	def summary(self):
		seq_writers = self.get_seq_writers()
		written = sum(w.written for w in seq_writers)
		written_bp = (
			sum(w.read1_bp for w in seq_writers),
			sum(w.read2_bp for w in seq_writers),
		)
		return (written, written_bp)
		
	def write(self, dest, read1, read2=None):
		writer = None
		
		if dest == FilterType.NONE and self.multiplexed and read1.match:
			name = read1.match.adapter.name
			writer = self.get_mux_writer(name)
		
		if writer is None and self.has_seq_writer(dest):
			writer = self.get_seq_writer(dest)
		
		if writer is not None:
			writer.write(read1, read2)
		else:
			self.discarded += 1
		
		self._write_info(read1)
		if read2:
			self._write_info(read2)
	
	def _write_info(self, read):
		match = read.match
		
		if self.rest_path is not None and match:
			self.rest_writer.write(match)
		
		if self.wildcard_path is not None and match:
			print(match.wildcards(), read.name, file=self.wildcard_writer)
		
		if self.info_path is not None:
			if match:
				for m in read.matches:
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
						sep='\t', file=self.info_writer
					)
			else:
				seq = read.sequence
				qualities = read.qualities if read.qualities is not None else ''
				print(read.name, -1, seq, qualities, sep='\t', file=self.info_writer)
	
	@property
	def rest_writer(self):
		if self._rest_writer is None:
			class RestFileWriter(object):
				def __init__(self, file):
					self.file = file
	
				def write(self, match):
					rest = match.rest()
					if len(rest) > 0:
						print(rest, match.read.name, file=self.file)
			
			self._rest_writer = RestFileWriter(
				self._get_file_writer(self.rest_path))
		
		return self._rest_writer
		
	@property
	def info_writer(self):
		return self._get_file_writer(self.info_path)
	
	@property
	def wildcard_writer(self):
		return self._get_file_writer(self.wildcard_path)	

	def _get_file_writer(self, path):
		assert path is not None
		if path not in self.writers:
			writer = xopen(path, 'w', buffer_size=self.buffer_size)
			self.writers[path] = writer
		return self.writers[path]
	
	def close(self):
		# touch any files in force_create that haven't been created
		for path, created in self.force_create.items():
			if not created:
				with open(path, 'w'): os.utime(path, None)
		# close any open writers
		for writer in self.writers.values():
			if writer is not None:
				writer.close()

def summarize_adapters(modifiers):
	mods = modifiers.get_modifier_pair(ModType.ADAPTER)
	summary = [{}, {}]
	if mods[0] is not None:
		summary[0] = collect_adapter_statistics(mods[0].adapters)
	if mods[1] is not None:
		summary[1] = collect_adapter_statistics(mods[1].adapters)
	return summary

def run_cutadapt_serial(reader, writers, modifiers, filters):
	total_bp1 = 0
	total_bp2 = 0
	
	try:
		n = 0
		for batch_size, batch in reader:
			n += batch_size
			for record in batch:
				reads, bp = modifiers.modify(record)
				total_bp1 += bp[0]
				total_bp2 += bp[1]
				dest = filters.filter(*reads)
				writers.write(dest, *reads)
		
		return Summary(
			collect_writer_statistics(n, total_bp1, total_bp2, writers),
			collect_process_statistics(modifiers, filters),
			summarize_adapters(modifiers)
		).finish()
	
	except KeyboardInterrupt as e:
		print("Interrupted", file=sys.stderr)
		sys.exit(130)
	except IOError as e:
		if e.errno == errno.EPIPE:
			sys.exit(1)
		raise
	except (seqio.FormatError, EOFError) as e:
		sys.exit("cutadapt: error: {0}".format(e))
	finally:
		reader.close()
		writers.close()

def run_cutadapt_parallel(reader, writers, modifiers, filters, threads, timeout=0, 
						  preserve_order=False, read_queue_size=0, result_queue_size=0):
	
	# TODO: should probably log the run_cutadapt_parallel arguments
	# TODO: would it help to use watermarking queues - i.e. queues that block
	# at some max level X and don't unblock until some min level Y?
	
	# timeout <= 0 means block forever
	if timeout <= 0: timeout = None
	
	# Limit queue sizes to prevent filling up memory
	read_queue = JoinableQueue(read_queue_size) 
	result_queue = JoinableQueue(result_queue_size)
	# Queue for threads to send summary information back to main process
	worker_summary_queue = JoinableQueue()
	writer_summary_queue = JoinableQueue()
	# Shared variable for telling threads to stop
	control = Value('l', 0)
	
	worker_threads = set()
	for i in range(threads):
		worker = WorkerThread(i, modifiers, filters, read_queue, 
			result_queue, worker_summary_queue, control, timeout)
		worker_threads.add(worker)
		worker.start()
	
	if preserve_order:
		writer_thread = OrderPreservingWriterThread(writers, 
			result_queue, writer_summary_queue, control, timeout)
	else:
		writer_thread = WriterThread(writers, 
			result_queue, writer_summary_queue, control, timeout)
	writer_thread.start()
	
	try:
		# Main loop - add batches of records to the queue
		for batch in enumerate(reader, 1):
			read_queue.put(batch, block=True)
		
		# Tell the writer thread the max number of batches to expect
		control.value = num_batches
		
		# Done writing to the input queue
		read_queue.close()
		read_queue.join()
		
		# Wait for the writer summary
		writer_summary = writer_summary_queue.get(block=True)
		if writer_summary is None:
			raise Excpetion("Writer thread died unexpectedly")
		summary = Summary(writer_summary)
		writer_summary_queue.task_done()
		writer_summary_queue.close()
		
		# Wait for the writer thread to end
		writer_thread.join()
		
		# Wait for summary information from worker threads
		seen_summaries = set()
		while len(seen_summaries) < threads:
			worker_summary = worker_summary_queue.get(block=True)
			if worker_summary is None:
				raise Exception("Worker thread died unexpectedly")
			seen_summaries.add(worker_summary[0])
			summary.add_process_stats(worker_summary[1])
			summary.add_adapter_stats(worker_summary[2])
			worker_summary_queue.task_done()
		assert worker_summary_queue.empty()
		worker_summary_queue.close()
		
		# wait for worker threads to complete and get
		for worker in worker_threads:
			worker.join()
		
		return summary.finish()
	
	except KeyboardInterrupt as e:
		print("Interrupted", file=sys.stderr)
		sys.exit(130)
	
	except IOError as e:
		if e.errno == errno.EPIPE:
			sys.exit(1)
		else:
			raise
	
	except (seqio.FormatError, EOFError) as e:
		sys.exit("cutadapt: error: {0}".format(e))
	
	finally:
		reader.close()
		
		# notify all threads that they should stop
		with control.get_lock():
			control.value = -1
		
		# Wait for all threads to finish
		def kill(t):
			if t.is_alive():
				# first try to be nice by waiting
				t.join(60)
			if t.is_alive():
				# if it takes too long, use harsher measures
				t.terminate()
				
		kill(writer_thread)
		for worker in worker_threads:
			kill(worker)

class PendingQueue(object):
	def _init(self):
		self.queue = {}
		self.min = None

	def push(self, priority, value):
		self.queue[priority] = value
		if self.min is None or priority < self.min:
			self.min = priority

	def pop(self):
		size = len(self.queue)
		if size == 0:
			raise Exception("PendingQueue is empty")
		value = self.queue.pop(self.min)
		if size == 1:
			self.min = None
		else:
			self.min = min(self.queue.keys())
		return value

	@property
	def empty(self):
		return len(self.queue) == 0

class WorkerThread(Process):
	"""
	Thread that processes batches of reads.
	input_queue: queue with batches of records to process
	output_queue: queue where results are written
	summary_queue: queue where summary information is written
	control: shared value with value 0 unless the thread
	should die (< 0) or there are no more batches coming (> 0).
	"""
	def __init__(self, index, modifiers, filters, input_queue, output_queue, 
				summary_queue, control, timeout=60):
		super(WorkerThread, self).__init__()
		self.index = index
		self.modifiers = modifiers
		self.filters = filters
		self.input_queue = input_queue
		self.output_queue = output_queue
		self.summary_queue = summary_queue
		self.control = control
		self.timeout = timeout
		self.processed_reads = 0
		self.processed_batches = 0
	
	def run(self):
		def process_next_batch():
			waiting = None
			while self.control.value >= 0:
				try:
					batch_num, (batch_size, records) = self.input_queue.get(block=True, timeout=1)
					result = [self._modify_and_filter(record) for record in records]
					self.processed_reads += batch_size
					self.input_queue.task_done()
					return (batch_num, result)
				except Empty:
					# TODO: log
					#print("Worker {} waiting for batch".format(self.index))
					if self.control.value == 0 and not waiting:
						waiting = time.time()
					elif (self.control.value > 0) or (time.time() - waiting >= self.timeout):
						break
			return None
		
		def send_result(result):
			if result is not None:
				waiting = None
				while self.control.value >= 0:
					try:
						self.output_queue.put(result, block=True, timeout=1)
						return True
					except Full:
						if not waiting:
							waiting = time.time()
						elif time.time() - waiting >= self.timeout:
							break
			return False
		
		while send_result(process_next_batch()):
			self.processed_batches += 1
		
		self.output_queue.close()
		self.output_queue.join()
		
		if self.control.value == 0:
			# TODO: log error
			print("Warning: worker thread {} exiting early "
				"due to no batches available".format(self.index))
			self.summary_queue.put(None)
		elif self.control.value > 0:
			process_stats = collect_process_statistics(self.modifiers, self.filters)
			adapter_stats = summarize_adapters(self.modifiers)
			self.summary_queue.put((self.index, process_stats, adapter_stats))
		
		self.summary_queue.close()
	
	def _modify_and_filter(self, record):
		reads, bp = self.modifiers.modify(record)
		dest = self.filters.filter(*reads)
		return (reads, bp, dest)

class WriterThread(Process):
	"""
	Thread that accepts results from the worker threads 
	and writes the results to output file(s). Not guaranteed
	to preserve the original order of sequence records.
	
	queue: input queue
	writers: Writers object used to write results to output files
	control: a shared value that is 0 unless the thread should
	die (< 0) or the main process has finished reading records
	and communicates the total number of batches (> 0).
	summary: a shared array to hold n (number of records processed),
	total_bp1 and total_bp2 (total read1 and read2 bp written).
	"""
	def __init__(self, writers, input_queue, summary_queue, control, timeout=60):
		super(WriterThread, self).__init__()
		self.writers = writers
		self.queue = input_queue
		self.summary_queue = summary_queue
		self.control = control
		self.n = 0	# no. of processed reads
		self.total_bp1 = 0
		self.total_bp2 = 0
		self.seen_batches = set()
		self.timeout = timeout
	
	def run(self):
		try:
			self._init()
			
			waiting = None
			while self.control.value == 0 or self.control.value > len(self.seen_batches):
				try:
					batch_num, records = self.queue.get(block=True, timeout=1)
					waiting = None
					self.seen_batches.add(batch_num)
					self.n += len(records)
					self._process_batch(batch_num, records)
					self.queue.task_done()
				except Empty:
					# TODO: log
					#print("Writer waiting for results")
					if not waiting:
						waiting = time.time()
					elif time.time() - waiting >= self.timeout:
						break
			
			if self.control.value == 0:
				# TODO: log error
				print("Warning: writer thread exiting early due to no results available")
				self.summary_queue.put(None)
			elif self.control.value > 0:
				self._no_more_batches()
				self.summary_queue.put(collect_writer_statistics(
					self.n, self.total_bp1, self.total_bp2, self.writers))
		finally:
			self.writers.close()
		
		self.summary_queue.close()
		self.summary_queue.join()
	
	def _init(self):
		pass

	def _process_batch(self, batch_num, records):
		self._write_batch(records)
		
	def _no_more_batches(self):
		pass
	
	def _write_batch(self, records):
		for reads, bp, dest in records:
			self.writers.write(dest, *reads)
			self.total_bp1 += bp[0]
			self.total_bp2 += bp[1]

class OrderPreservingWriterThread(WriterThread):
	"""
	Writer thread that is less time/memory efficient, but is
	guaranteed to preserve the original order of records.
	"""
	def _init(self):
		self.pending = PendingQueue()
		self.cur_batch = 1
	
	def _process_batch(self, batch_num, records):
		if batch_num == self.cur_batch:
			self._write_batch(batch_num, records)
			self.cur_batch += 1
			self._consume_pending()
		else:
			self.pending.push(batch_num, records)
	
	def _no_more_batches(self):
		self._consume_pending()
		if not self.pending.empty:
			raise Exception("Did not receive all batches")
	
	def _consume_pending(self):
		while (not self.pending.empty) and (self.cur_batch == pending.min):
			self.write_batch(pending.pop())
			self.cur_batch += 1

if __name__ == '__main__':
	main()
