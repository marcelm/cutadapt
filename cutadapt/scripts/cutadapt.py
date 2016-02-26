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
from cutadapt.report import Statistics, print_report, redirect_standard_output
from cutadapt.compat import next

from collections import OrderedDict
from heapq import nsmallest
import logging
from multiprocessing import Process, Queue, Value, cpu_count
from queue import Empty
from optparse import OptionParser, OptionGroup
import os
import platform
import sys
import time

def main(cmdlineargs=None, default_outfile="-"):
    """
    Main function that evaluates command-line parameters and iterates
    over all reads.

    default_outfile is the file to which trimmed reads are sent if the ``-o``
    parameter is not used.
    """
    parser = get_option_parser()
    
    if cmdlineargs is None:
        cmdlineargs = sys.argv[1:]
    
    options, args = parser.parse_args(args=cmdlineargs)
    
    # Setup logging only if there are not already any handlers (can happen when
    # this function is being called externally such as from unit tests)
    if not logging.root.handlers:
        setup_logging(stdout=bool(options.output), quiet=options.quiet)
    
    if options.rrbs:
        set_rrbs_defaults(options)
    elif options.wgbs:
        set_wgbs_defaults(options)
    elif options.mirna:
        set_mirna_defaults(options)

    paired, multiplexed = validate_options(options, args, parser)
    input_filename = args[0]
    input_paired_filename = None
    quality_filename = None
    if len(args) > 1:
        if paired:
            input_paired_filename = args[1]
        else:
            quality_filename = args[1]

    reader = writers = None
    try:
        try:
            # Open input file(s)
            reader = seqio.open(input_filename, file2=input_paired_filename,
                qualfile=quality_filename, colorspace=options.colorspace, 
                fileformat=options.format, interleaved=options.interleaved)
        except (seqio.UnknownFileType, IOError) as e:
            parser.error(e)
        
        qualities = reader.delivers_qualities
        has_qual_file = quality_filename is not None
        writers = Writers(options, multiplexed, qualities, default_outfile)
        modifiers, adapters = create_modifiers(options, paired, qualities, has_qual_file, parser)
        min_affected = 2 if options.pair_filter == 'both' else 1
        filters = create_filters(options, paired, min_affected)
        
        logger = logging.getLogger()
        logger.info("This is cutadapt %s with Python %s", __version__, platform.python_version())
        logger.info("Command line parameters: %s", " ".join(cmdlineargs))
        num_adapters = len(adapters[0]) + len(adapters[1])
        logger.info("Trimming %s adapter%s with at most %.1f%% errors in %s mode ...",
            num_adapters, 's' if num_adapters > 1 else '', options.error_rate * 100,
            { False: 'single-end', 'first': 'paired-end legacy', 'both': 'paired-end' }[paired])
        if paired == 'first' and (len(modifiers[1]) > 0 or options.quality_cutoff):
            import textwrap
            logger.warning('\n'.join(textwrap.wrap('WARNING: Requested read '
                'modifications are applied only to the first '
                'read since backwards compatibility mode is enabled. '
                'To modify both reads, also use any of the -A/-B/-G/-U options. '
                'Use a dummy adapter sequence when necessary: -A XXX')))
        
        start_time = time.clock()
        
        if options.threads == 1:
            # Run cutadapt normally
            stats = run_cutadapt_serial(paired, reader, writers, modifiers, filters)
    
        else:
            # Run multiprocessing version
            stats = run_cutadapt_parallel(paired, reader, writers, modifiers, filters, 
                options.threads, options.batch_size, options.preserve_order)
        
        elapsed_time = time.clock() - start_time
        
        if not options.quiet:
            stats.collect(adapters, elapsed_time, modifiers, filters, writers)
            # send statistics to stderr if result was sent to stdout
            stat_file = sys.stderr if options.output is None else None
            with redirect_standard_output(stat_file):
                print_report(stats, adapters)
    
    finally:
        # close open files
        if reader:
            reader.close()
        if writers:
            writers.close()

def get_option_parser():
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
    group.add_option("-e", "--error-rate", type=float, default=0.1,
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
    group.add_option("-m", "--minimum-length", type=int, default=0, metavar="LENGTH",
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
    group.add_option("--threads", type=int, default=1, metavar="THREADS",
        help="Number of threads to use for read trimming. Set to 0 to use max available threads.")
    group.add_option("--batch-size", type=int, default=1000, metavar="SIZE",
        help="Number of records to process in each batch")
    group.add_option("--preserve-order", action="store_true", default=False,
        help="Preserve order of reads in input files")
    parser.add_option_group(group)
    
    group = OptionGroup(parser, "Method-specific options")
    group.add_option("--rrbs", action="store_true", default=False,
        help="Set default option values for RRBS data")
    group.add_option("--wgbs", action="store_true", default=False,
        help="Set default option values for WGBS data")
    group.add_option("--non-directional", action="store_true", default=False,
        help="Bisulfite sequencing libraries are non-directional")
    group.add_option("--mirna", action="store_true", default=False,
        help="Set default option values for miRNA data")
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
    logger.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

def set_rrbs_defaults(options):
    pass

def set_wgbs_defaults(options):
    pass

def set_mirna_defaults(options):
    if options.minimum_length == 0 or options.minimum_length > 16:
        options.minimum_length = 16
    if options.error_rate < 0.12:
        options.error_rate = 0.12

def validate_options(options, args, parser):
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
    if not (0 <= options.error_rate <= 1.):
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

    return (paired, multiplexed)

def create_filters(options, paired, min_affected):
    filters = OrderedDict()
    filter_factory = FilterFactory(paired, min_affected)
    
    def add_filter(filter_type, *args, **kwargs):
        f = filter_factory(filter_type, *args, **kwargs)
        filters[filter_type] = f
    
    if options.minimum_length > 0:
        add_filter(FilterType.TOO_SHORT, options.minimum_length)
    
    if options.maximum_length < sys.maxsize:
        add_filter(FilterType.TOO_LONG, options.maximum_length)
    
    if options.max_n >= 0:
        add_filter(FilterType.N_CONTENT, options.max_n)
    
    if options.discard_trimmed:
        add_filter(FilterType.TRIMMED)
    
    if options.discard_untrimmed or options.untrimmed_output:
        add_filter(FilterType.UNTRIMMED)
    
    return filters

def create_modifiers(options, paired, qualities, has_qual_file, parser):
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
            parser.error(e)
        raise
    except ValueError as e:
        parser.error(e)
    if options.debug:
        for adapter in adapters + adapters2:
            adapter.enable_debug()
    
    if not adapters and not adapters2 and not options.quality_cutoff and \
            options.nextseq_trim is None and \
            options.cut == [] and options.cut2 == [] and \
            options.minimum_length == 0 and \
            options.maximum_length == sys.maxsize and \
            not has_qual_file and \
            options.max_n == -1 and not options.trim_n:
        parser.error("You need to provide at least one adapter sequence.")
    
    # Create the single-end processing pipeline (a list of "modifiers")
    modifiers = [OrderedDict(), OrderedDict()]
    
    def add_modifier(mod_type, *args, _dest=None, **kwargs):
        if _dest is None:
            mod = create_modifier(mod_type, *args, **kwargs)
            for d in modifiers:
                d[mod_type] = mod
        else:
            modifiers[_dest-1][mod_type] = create_modifier(mod_type, *args, **kwargs)
    
    if options.cut:
        add_modifier(ModType.CUT, options.cut, _dest=1)
    
    if options.nextseq_trim is not None:
        add_modifier(ModType.TRIM_NEXTSEQ_QUAL, options.nextseq_trim, options.quality_base, _dest=1)
    
    if options.quality_cutoff:
        add_modifier(ModType.TRIM_QUAL, options.quality_cutoff[0], 
            options.quality_cutoff[1], options.quality_base, _dest=1)
    
    if adapters:
        add_modifier(ModType.ADAPTER, adapters, times=options.times, action=options.action, _dest=1)
    
    # For paired-end data, create a second processing pipeline.
    # However, if no second-read adapters were given (via -A/-G/-B/-U), we need to
    # be backwards compatible and *no modifications* are done to the second read.
    if paired == 'both':
        if options.cut2:
            add_modifier(ModType.CUT, options.cut2, _dest=2)
    
        if options.quality_cutoff:
            add_modifier(ModType.TRIM_QUAL, options.quality_cutoff[0], 
                options.quality_cutoff[1], options.quality_base, _dest=2),
        
        if adapters2:
            add_modifier(ModType.ADAPTER, adapters2, times=options.times, action=options.action, _dest=2)

    # Modifiers that apply to both reads of paired-end reads unless in legacy mode
    
    if options.rrbs or options.wgbs:
        # trim two additional bases from adapter-trimmed sequences
        add_modifier(ModType.CLIP_ADAPTER_TRIMMED, modifiers, -2)
    
    if options.trim_n:
        add_modifier(ModType.TRIM_END_N)
    
    if options.length_tag:
        add_modifier(ModType.LENGTH_TAG, options.length_tag)
    
    if options.strip_suffix:
        add_modifier(ModType.REMOVE_SUFFIX, options.strip_suffix)
    
    if options.prefix or options.suffix:
        add_modifier(ModType.ADD_PREFIX_SUFFIX, options.prefix, options.suffix)
    
    if options.double_encode:
        add_modifier(ModType.CS_DOUBLE_ENCODE)
    
    if options.zero_cap and qualities:
        add_modifier(ModType.ZERO_CAP, quality_base=options.quality_base)
    
    if options.trim_primer:
        add_modifier(ModType.CS_TRIM_PRIMER)
    
    return (modifiers, (adapters, adapters2))

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
        self.open_args = dict(
            qualities=qualities,
            colorspace=options.colorspace, 
            interleaved=options.interleaved
        )
        self.seqfile_paths = dict()
        self.force_create = set()
        self.rest_path = options.rest_file
        self.info_path = options.info_file
        self.wildcard_path = options.wildcard_file
        self.writers = dict()
        self._rest_writer = None
        self.discarded = 0
        
        if options.minimum_length > 0 and options.too_short_output:
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
            self.force_create.add(file1)
            if file2 is not None:
                self.force_create.add(file2)
    
    def has_seq_writer(self, filter_type):
        return filter_type in self.seqfile_paths
    
    def get_seq_writer(self, filter_type):
        paths = self.seqfile_paths[filter_type]
        if paths not in self.writers:
            self.writers[paths] = self._create_seq_writer(*paths)
        return self.writers[paths]
    
    def get_mux_writer(self, name):
        path = self.output.format(name=name)
        if path not in self.writers:
            self.writers[path] = self._create_seq_writer(path)
        return self.writers[path]
        
    def _create_seq_writer(self, file1, file2=None):
        seqfile = seqio.open(file1, file2, mode='w', **self.open_args)
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
            writer = xopen(path, 'w')
            self.writers[path] = writer
        return self.writers[path]
    
    def close(self):
        # touch any files in force_create that haven't been created
        for f in self.force_create:
            if not os.path.exists(f):
                with open(f, 'a'): os.utime(f, None)
        # close any open writers
        for writer in self.writers.values():
            if writer is not None:
                writer.close()

def modify_and_filter(record, paired, modifiers, filters):
    dest = FilterType.NONE
    bp = [0,0]
    if paired:
        read1, read2 = record
        bp[0] = len(read1.sequence)
        bp[1] = len(read2.sequence)
        for mod in modifiers[0].values():
            read1 = mod(read1)
        for mod in modifiers[1].values():
            read2 = mod(read2)
        reads = [read1, read2]
    else:
        read = record
        bp[0] = len(read.sequence)
        for mod in modifiers[0].values():
            read = mod(read)
        reads = [read]
    
    for filter_type, f in filters.items():
        if f(*reads):
            dest = filter_type
            # Stop writing as soon as one of the filters was successful.
            break
    
    return (dest, reads, bp)

def run_cutadapt_serial(paired, reader, writers, modifiers, filters):
    try:
        n = 0  # no. of processed reads
        total1_bp = 0
        total2_bp = 0 if paired else None
        for record in reader:
            n += 1
            dest, reads, bp = modify_and_filter(record, paired, modifiers, filters)
            writers.write(dest, *reads)
            total1_bp += bp[0]
            if paired:
                total2_bp += bp[1]
        
        return Statistics(n=n, total_bp1=total1_bp, total_bp2=total2_bp)

    except KeyboardInterrupt as e:
        print("Interrupted", file=sys.stderr)
        sys.exit(130)
    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(1)
        raise
    except (seqio.FormatError, EOFError) as e:
        sys.exit("cutadapt: error: {0}".format(e))

def run_cutadapt_parallel(paired, reader, writers, modifiers, filters, threads=None, 
                          batch_size=1000, preserve_order=False, max_wait=60):
    if threads is None or threads <= 0:
        threads = min(cpu_count() - 1, 1)
    
    read_queue = Queue()
    result_queue = Queue()
    control = Value('l', 0)
    
    # start WorkerThreads
    workers = set()
    for i in range(threads):
        worker = Worker(read_queue, result_queue, control, paired, modifiers, filters)
        workers.add(worker)
        worker.start()
    
    # start WriterThread
    if preserve_order:
        writer = OrderPreservingWriterThread(result_queue, writers, control)
    else:
        writer = WriterThread(result_queue, writers, control)
    writer.start()
    
    try:
        # Main loop - add batches of records to the queue
        emtpy_batch = [None] * batch_size
        batch = empty_batch.copy()
        num_batches = 0
        batch_index = 0
        for record in reader:
            batch[batch_index] = record
            batch_index += 1
            if batch_index == batch_size:
                read_queue.put(batch)
                num_batches += 1
                batch = empty_batch.copy()
                batch_index = 0
        
        if len(batch) > 0:
            read_queue.put(batch)
            num_batches += 1
        
        # Tell the writer thread the max number of
        # batches to expect
        with control.get_lock():
            control.value = num_batches
        
        # Wait for the writer to complete
        writer_thread.join()
        
        return Statistics(
            n=writer_thread.n, 
            total_bp1=writer_thread.total1_bp, 
            total_bp2=writer_thread.total2_bp
        )

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
        # notify all threads that they should stop
        with control.get_lock():
            control.value = -1
        
        # Wait for all threads to finish
        print("Waiting for threads to die...")
        
        def kill(t):
            if t.is_alive():
                # first try to be nice by waiting
                t.join(max_wait)
            if t.is_alive():
                # if it takes too long, use harsher measures
                t.terminate()
                
        kill(writer_thread)
        for worker in workers:
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
    """
    def __init__(self, queue, writers, control):
        super(WriterThread, self).__init__()
        self.queue = queue
        self.writers = writers
        self.control = control
        self.n = 0  # no. of processed reads
        self.total1_bp = 0
        self.total2_bp = 0
        self.seen_batches = set()
    
    def run(self):
        self._init()
        while self.control.value >= 0:
            done = False
            if 0 < self.control.value <= len(self.seen_batches):
                done = True
            else:
                try:
                    batch = self.queue.get(timeout=1)
                except Empty:
                    continue
                
                if batch is None:
                    done = True
                else:
                    batch_num, records = batch
                    self.n += len(records)
                    self.seen_batches.add(batch_num)
                    self._process_batch(batch_num, records)
            
            if done:
                self._no_more_batches()
                break
    
    def _init(self):
        pass

    def _process_batch(self, batch_num, records):
        self._write_batch(records)
        
    def _no_more_batches(self):
        pass
    
    def _write_batch(records):
        for dest, reads, bp in records:
            self.writers.write(dest, *reads)
            self.total1_bp += bp[0]
            self.total2_bp += bp[1]

class OrderPreservingWriterThread(WriterThread):
    """
    Writer thread that is less time/memory efficient, but is
    guaranteed to preserve the original order of records.
    """
    def _init(self):
        self.pending = PendingQueue()
        self.cur_batch = 0
    
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

class WorkerThread(Process):
    """
    Thread that processes batches of reads.
    input_queue: queue with batches of records to process
    output_queue: queue where results are written
    control: shared value with value 0 unless the thread
    should die (< 0) or there are no more batches coming (> 0).
    """
    def __init__(self, input_queue, output_queue, control, paired, modifiers, filters):
        self.input_queue = input_queue
        self.output_queue = output_queue
        self.control = control
        self.paired = paired
        self.modifiers = modifiers
        self.filters = filters
    
    def run(self):
        while self.control.value >= 0:
            try:
                batch = self.input_queue.get(timeout=1)
                result = modify_and_filter(record, self.paired, self.modifiers, self.filters)
                self.output_queue.put(result)
            except Empty:
                if self.control.value > 0:
                    break

if __name__ == '__main__':
    main()
