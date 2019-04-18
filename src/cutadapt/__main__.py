#!/usr/bin/env python
#
# Copyright (c) 2010-2019 Marcel Martin <marcel.martin@scilifelab.se>
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
cutadapt version {version}

Copyright (C) 2010-2019 Marcel Martin <marcel.martin@scilifelab.se>

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

import sys
import errno
import time
from argparse import ArgumentParser, SUPPRESS, HelpFormatter
import logging
import platform
from xopen import xopen
import dnaio

from cutadapt import __version__
from cutadapt.adapters import AdapterParser
from cutadapt.modifiers import (LengthTagModifier, SuffixRemover, PrefixSuffixAdder,
    ZeroCapper, QualityTrimmer, UnconditionalCutter, NEndTrimmer, AdapterCutter,
    PairedAdapterCutterError, PairedAdapterCutter, NextseqQualityTrimmer, Shortener)
from cutadapt.report import full_report, minimal_report
from cutadapt.pipeline import (SingleEndPipeline, PairedEndPipeline, InputFiles, OutputFiles,
    SerialPipelineRunner, ParallelPipelineRunner)
from cutadapt.utils import available_cpu_count

logger = logging.getLogger()


class CutadaptArgumentParser(ArgumentParser):
    """
    This ArgumentParser customizes two things:
    - The usage message is not prefixed with 'usage:'
    - A brief message is shown on errors, not full usage
    """
    class CustomUsageHelpFormatter(HelpFormatter):
        def add_usage(self, usage, actions, groups, prefix=None):
            if usage is not SUPPRESS:
                args = usage, actions, groups, ''
                self._add_item(self._format_usage, args)

    def __init__(self, *args, **kwargs):
        kwargs['formatter_class'] = self.CustomUsageHelpFormatter
        kwargs['usage'] = kwargs['usage'].replace("{version}", __version__)
        super().__init__(*args, **kwargs)

    def error(self, message):
        """
        If you override this in a subclass, it should not return -- it
        should either exit or raise an exception.
        """
        print('Run "cutadapt --help" to see command-line options.', file=sys.stderr)
        print('See https://cutadapt.readthedocs.io/ for full documentation.', file=sys.stderr)
        self.exit(2, "\n{prog}: error: {message}\n".format(prog=self.prog, message=message))


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
        return super().format(record)


def setup_logging(stdout=False, quiet=False, debug=False):
    """
    Attach handler to the global logger object
    """
    # Due to backwards compatibility, logging output is sent to standard output
    # instead of standard error if the -o option is used.
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setFormatter(NiceFormatter())
    # debug overrides quiet
    if debug:
        level = logging.DEBUG
    elif quiet:
        level = logging.ERROR
    else:
        level = logging.INFO
    stream_handler.setLevel(level)
    logger.setLevel(level)
    logger.addHandler(stream_handler)


def get_argument_parser():
    parser = CutadaptArgumentParser(usage=__doc__, add_help=False)
    group = parser.add_argument_group("Options")
    group.add_argument("-h", "--help", action="help", help="Show this help message and exit")
    group.add_argument("--version", action="version", help="Show version number and exit",
        version=__version__)
    group.add_argument("--debug", action='store_true', default=False,
        help="Print debugging information.")
    group.add_argument("--profile", action="store_true", default=False, help=SUPPRESS)
    group.add_argument('-j', '--cores', type=int, default=1,
        help='Number of CPU cores to use. Use 0 to auto-detect. Default: %(default)s')

    # Hidden options
    # GC content as a percentage
    group.add_argument("--gc-content", type=float, default=50,
        help=SUPPRESS)
    # Buffer size for the reader process when running in parallel
    group.add_argument("--buffer-size", type=int, default=4000000,
        help=SUPPRESS)
    # Deprecated: The input format is always auto-detected
    group.add_argument("-f", "--format", help=SUPPRESS)

    group = parser.add_argument_group("Finding adapters",
        description="Parameters -a, -g, -b specify adapters to be removed from "
            "each read (or from the first read in a pair if data is paired). "
            "If specified multiple times, only the best matching adapter is "
            "trimmed (but see the --times option). When the special notation "
            "'file:FILE' is used, adapter sequences are read from the given "
            "FASTA file.")
    group.add_argument("-a", "--adapter", action="append", default=[], metavar="ADAPTER",
        dest="adapters",
        help="Sequence of an adapter ligated to the 3' end (paired data: of the "
            "first read). The adapter and subsequent bases are trimmed. If a "
            "'$' character is appended ('anchoring'), the adapter is only "
            "found if it is a suffix of the read.")
    group.add_argument("-g", "--front", action="append", default=[], metavar="ADAPTER",
        help="Sequence of an adapter ligated to the 5' end (paired data: of the "
            "first read). The adapter and any preceding bases are trimmed. "
            "Partial matches at the 5' end are allowed. If a '^' character is "
            "prepended ('anchoring'), the adapter is only found if it is a "
            "prefix of the read.")
    group.add_argument("-b", "--anywhere", action="append", default=[], metavar="ADAPTER",
        help="Sequence of an adapter that may be ligated to the 5' or 3' end "
            "(paired data: of the first read). Both types of matches as "
            "described under -a und -g are allowed. If the first base of the "
            "read is part of the match, the behavior is as with -g, otherwise "
            "as with -a. This option is mostly for rescuing failed library "
            "preparations - do not use if you know which end your adapter was "
            "ligated to!")
    group.add_argument("-e", "--error-rate", type=float, default=0.1, metavar="RATE",
        help="Maximum allowed error rate as value between 0 and 1 (no. of "
            "errors divided by length of matching region). Default: %(default)s (=10%%)")
    group.add_argument("--no-indels", action='store_false', dest='indels', default=True,
        help="Allow only mismatches in alignments. "
            "Default: allow both mismatches and indels")
    group.add_argument("-n", "--times", type=int, metavar="COUNT", default=1,
        help="Remove up to COUNT adapters from each read. Default: %(default)s")
    group.add_argument("-O", "--overlap", type=int, metavar="MINLENGTH", default=3,
        help="Require MINLENGTH overlap between read and adapter for an adapter "
            "to be found. Default: %(default)s")
    group.add_argument("--match-read-wildcards", action="store_true", default=False,
        help="Interpret IUPAC wildcards in reads. Default: %(default)s")
    group.add_argument("-N", "--no-match-adapter-wildcards", action="store_false",
        default=True, dest='match_adapter_wildcards',
        help="Do not interpret IUPAC wildcards in adapters.")
    group.add_argument("--action", choices=('trim', 'mask', 'lowercase', 'none'), default='trim',
        help="What to do with found adapters. "
            "mask: replace with 'N' characters; "
            "lowercase: convert to lowercase; "
            "none: leave unchanged (useful with "
            "--discard-untrimmed). Default: trim")
    group.add_argument("--no-trim", dest='action', action='store_const', const='none',
        help=SUPPRESS)  # Deprecated, use --action=none
    group.add_argument("--mask-adapter", dest='action', action='store_const', const='mask',
        help=SUPPRESS)  # Deprecated, use --action=mask

    group = parser.add_argument_group("Additional read modifications")
    group.add_argument("-u", "--cut", action='append', default=[], type=int, metavar="LENGTH",
        help="Remove bases from each read (first read only if paired). "
            "If LENGTH is positive, remove bases from the beginning. "
            "If LENGTH is negative, remove bases from the end. "
            "Can be used twice if LENGTHs have different signs. "
            "This is applied *before* adapter trimming.")
    group.add_argument("--nextseq-trim", type=int, default=None, metavar="3'CUTOFF",
        help="NextSeq-specific quality trimming (each read). Trims also dark "
            "cycles appearing as high-quality G bases.")
    group.add_argument("-q", "--quality-cutoff", default=None, metavar="[5'CUTOFF,]3'CUTOFF",
        help="Trim low-quality bases from 5' and/or 3' ends of each read before "
            "adapter removal. Applied to both reads if data is paired. If one "
            "value is given, only the 3' end is trimmed. If two "
            "comma-separated cutoffs are given, the 5' end is trimmed with "
            "the first cutoff, the 3' end with the second.")
    group.add_argument("--quality-base", type=int, default=33, metavar='N',
        help="Assume that quality values in FASTQ are encoded as ascii(quality "
            "+ N). This needs to be set to 64 for some old Illumina "
            "FASTQ files. Default: %(default)s")
    group.add_argument("--length", "-l", type=int, default=None, metavar="LENGTH",
            help="Shorten reads to LENGTH. Positive values remove bases at the end "
            "while negative ones remove bases at the beginning. This and the "
            "following modifications are applied after adapter trimming.")
    group.add_argument("--trim-n", action='store_true', default=False,
        help="Trim N's on ends of reads.")
    group.add_argument("--length-tag", metavar="TAG",
        help="Search for TAG followed by a decimal number in the description "
            "field of the read. Replace the decimal number with the correct "
            "length of the trimmed read. For example, use --length-tag 'length=' "
            "to correct fields like 'length=123'.")
    group.add_argument("--strip-suffix", action='append', default=[],
        help="Remove this suffix from read names if present. Can be given multiple times.")
    group.add_argument("-x", "--prefix", default='',
        help="Add this prefix to read names. Use {name} to insert the name of the matching "
            "adapter.")
    group.add_argument("-y", "--suffix", default='',
        help="Add this suffix to read names; can also include {name}")
    group.add_argument("--zero-cap", "-z", action='store_true', default=False,
        help="Change negative quality values to zero.")

    group = parser.add_argument_group("Filtering of processed reads",
        description="Filters are applied after above read modifications. "
            "Paired-end reads are always discarded pairwise (see also "
            "--pair-filter).")
    group.add_argument("-m", "--minimum-length", default=None, metavar="LEN[:LEN2]",
        help="Discard reads shorter than LEN. Default: 0")
    group.add_argument("-M", "--maximum-length", default=None, metavar="LEN[:LEN2]",
        help="Discard reads longer than LEN. Default: no limit")
    group.add_argument("--max-n", type=float, default=None, metavar="COUNT",
        help="Discard reads with more than COUNT 'N' bases. If COUNT is a number "
             "between 0 and 1, it is interpreted as a fraction of the read length.")
    group.add_argument("--discard-trimmed", "--discard", action='store_true', default=False,
        help="Discard reads that contain an adapter. Use also -O to avoid "
            "discarding too many randomly matching reads.")
    group.add_argument("--discard-untrimmed", "--trimmed-only", action='store_true', default=False,
        help="Discard reads that do not contain an adapter.")
    group.add_argument("--discard-casava", action='store_true', default=False,
        help="Discard reads that did not pass CASAVA filtering (header has :Y:).")

    group = parser.add_argument_group("Output")
    group.add_argument("--quiet", default=False, action='store_true',
        help="Print only error messages.")
    group.add_argument("--report", choices=('full', 'minimal'), default=None,
        help="Which type of report to print: 'full' or 'minimal'. Default: full")
    group.add_argument("-o", "--output", metavar="FILE",
        help="Write trimmed reads to FILE. FASTQ or FASTA format is chosen "
            "depending on input. The summary report is sent to standard output. "
            "Use '{name}' in FILE to demultiplex reads into multiple "
            "files. Default: write to standard output")
    group.add_argument("--info-file", metavar="FILE",
        help="Write information about each read and its adapter matches into FILE. "
            "See the documentation for the file format.")
    group.add_argument("-r", "--rest-file", metavar="FILE",
        help="When the adapter matches in the middle of a read, write the "
            "rest (after the adapter) to FILE.")
    group.add_argument("--wildcard-file", metavar="FILE",
        help="When the adapter has N wildcard bases, write adapter bases "
            "matching wildcard positions to FILE. (Inaccurate with indels.)")
    group.add_argument("--too-short-output", metavar="FILE",
        help="Write reads that are too short (according to length specified by "
        "-m) to FILE. Default: discard reads")
    group.add_argument("--too-long-output", metavar="FILE",
        help="Write reads that are too long (according to length specified by "
        "-M) to FILE. Default: discard reads")
    group.add_argument("--untrimmed-output", default=None, metavar="FILE",
        help="Write reads that do not contain any adapter to FILE. Default: "
            "output to same file as trimmed reads")

    group = parser.add_argument_group("Paired-end options", description="The "
        "-A/-G/-B/-U options work like their -a/-b/-g/-u counterparts, but "
        "are applied to the second read in each pair.")
    group.add_argument("-A", dest='adapters2', action='append', default=[], metavar='ADAPTER',
        help="3' adapter to be removed from second read in a pair.")
    group.add_argument("-G", dest='front2', action='append', default=[], metavar='ADAPTER',
        help="5' adapter to be removed from second read in a pair.")
    group.add_argument("-B", dest='anywhere2', action='append', default=[], metavar='ADAPTER',
        help="5'/3 adapter to be removed from second read in a pair.")
    group.add_argument("-U", dest='cut2', action='append', default=[], type=int, metavar="LENGTH",
        help="Remove LENGTH bases from second read in a pair.")
    group.add_argument("-p", "--paired-output", metavar="FILE",
        help="Write second read in a pair to FILE.")
    group.add_argument("--pair-adapters", action="store_true",
        help="Treat adapters given with -a/-A etc. as pairs. Either both "
             "or none are removed from each read pair.")
    # Setting the default for pair_filter to None allows us to find out whether
    # the option was used at all.
    group.add_argument("--pair-filter", metavar='(any|both|first)', default=None,
        choices=("any", "both", "first"),
        help="Which of the reads in a paired-end read have to match the "
            "filtering criterion in order for the pair to be filtered. "
            "Default: any")
    group.add_argument("--interleaved", action='store_true', default=False,
        help="Read and write interleaved paired-end reads.")
    group.add_argument("--untrimmed-paired-output", metavar="FILE",
        help="Write second read in a pair to this FILE when no adapter "
            "was found. Use with --untrimmed-output. Default: output "
            "to same file as trimmed reads")
    group.add_argument("--too-short-paired-output", metavar="FILE", default=None,
        help="Write second read in a pair to this file if pair is too short. "
            "Use also --too-short-output.")
    group.add_argument("--too-long-paired-output", metavar="FILE", default=None,
        help="Write second read in a pair to this file if pair is too long. "
            "Use also --too-long-output.")

    for arg in ("--colorspace", "-c", "-d", "--double-encode", "-t", "--trim-primer",
            "--strip-f3", "--maq", "--bwa", "--no-zero-cap"):
        group.add_argument(arg, dest='colorspace', action='store_true', default=False,
        help=SUPPRESS)
    parser.set_defaults(colorspace=False)

    # We could have two positional arguments here, with the second one optional, but
    # we want custom, more helpful error messages.
    parser.add_argument("inputs", nargs='*', help=SUPPRESS)

    return parser


def parse_cutoffs(s):
    """Parse a string INT[,INT] into a two-element list of integers"""
    cutoffs = s.split(',')
    if len(cutoffs) == 1:
        try:
            cutoffs = [0, int(cutoffs[0])]
        except ValueError as e:
            raise CommandLineError("Quality cutoff value not recognized: {}".format(e))
    elif len(cutoffs) == 2:
        try:
            cutoffs = [int(cutoffs[0]), int(cutoffs[1])]
        except ValueError as e:
            raise CommandLineError("Quality cutoff value not recognized: {}".format(e))
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
        raise CommandLineError("Value not recognized: {}".format(e))
    if len(values) == 2 and values[0] is None and values[1] is None:
        raise CommandLineError("Cannot parse {!r}: At least one length needs to be given".format(s))
    return tuple(values)


def open_output_files(args, default_outfile, interleaved):
    """
    Return an OutputFiles instance. If demultiplex is True, the untrimmed, untrimmed2, out and out2
    attributes are not opened files, but paths (out and out2 with the '{name}' template).
    """
    rest_file = info_file = wildcard = None
    if args.rest_file is not None:
        rest_file = xopen(args.rest_file, 'w')
    if args.info_file is not None:
        info_file = xopen(args.info_file, 'w')
    if args.wildcard_file is not None:
        wildcard = xopen(args.wildcard_file, 'w')

    def open2(path1, path2):
        file1 = file2 = None
        if path1 is not None:
            file1 = xopen(path1, 'wb')
            if path2 is not None:
                file2 = xopen(path2, 'wb')
        return file1, file2

    too_short = too_short2 = None
    if args.minimum_length is not None:
        too_short, too_short2 = open2(args.too_short_output, args.too_short_paired_output)

    too_long = too_long2 = None
    if args.maximum_length is not None:
        too_long, too_long2 = open2(args.too_long_output, args.too_long_paired_output)

    if int(args.discard_trimmed) + int(args.discard_untrimmed) + int(
        args.untrimmed_output is not None) > 1:
        raise CommandLineError("Only one of the --discard-trimmed, --discard-untrimmed "
            "and --untrimmed-output options can be used at the same time.")

    demultiplex = args.output is not None and '{name}' in args.output
    if args.paired_output is not None and (demultiplex != ('{name}' in args.paired_output)):
        raise CommandLineError('When demultiplexing paired-end data, "{name}" must appear in '
            'both output file names (-o and -p)')

    if demultiplex:
        if args.discard_trimmed:
            raise CommandLineError("Do not use --discard-trimmed when demultiplexing.")

        out = args.output
        untrimmed = args.output.replace('{name}', 'unknown')
        if args.untrimmed_output:
            untrimmed = args.untrimmed_output
        if args.discard_untrimmed:
            untrimmed = None

        if args.paired_output is not None:
            out2 = args.paired_output
            untrimmed2 = args.paired_output.replace('{name}', 'unknown')
            if args.untrimmed_paired_output:
                untrimmed2 = args.untrimmed_paired_output
            if args.discard_untrimmed:
                untrimmed2 = None

        else:
            untrimmed2 = out2 = None
    else:
        untrimmed, untrimmed2 = open2(args.untrimmed_output, args.untrimmed_paired_output)
        out, out2 = open2(args.output, args.paired_output)
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


def determine_paired_mode(args):
    """
    Determine whether we should work in paired-end mode.
    """
    # Usage of any of these options enables paired-end mode
    return bool(
        args.paired_output
        or args.interleaved
        or args.adapters2
        or args.front2
        or args.anywhere2
        or args.cut2
        or args.pair_filter
        or args.too_short_paired_output
        or args.too_long_paired_output)


def determine_interleaved(args):
    is_interleaved_input = False
    is_interleaved_output = False
    if args.interleaved:
        is_interleaved_input = len(args.inputs) == 1
        is_interleaved_output = not args.paired_output
        if not is_interleaved_input and not is_interleaved_output:
            raise CommandLineError("When --interleaved is used, you cannot provide both two "
                "input files and two output files")
    return is_interleaved_input, is_interleaved_output


def input_files_from_parsed_args(inputs, paired, interleaved):
    """
    Return tuple (input_filename, input_paired_filename)
    """
    if len(inputs) == 0:
        raise CommandLineError(
            "You did not provide any input file names. Please give me something to do!")
    elif len(inputs) > 2:
        raise CommandLineError(
            "You provided {} input file names, but either one or two are expected. ".format(
                len(inputs))
            + "The file names were:\n - "
            + "\n - ".join("{!r}".format(p) for p in inputs)
            + "\nHint: If your path contains spaces, you need to enclose it in quotes")
    input_filename = inputs[0]
    if paired and not interleaved:
        # Two file names required
        if len(inputs) == 1:
            raise CommandLineError(
                "You used an option that enabled paired-end mode (such as -p, -A, -G, -B, -U), "
                "but then you also need to provide two input files (you provided one) or "
                "use --interleaved.")
        else:
            input_paired_filename = inputs[1]
    else:
        if len(inputs) == 2:
            raise CommandLineError(
                "It appears you want to trim paired-end data because you provided two input files, "
                "but then you also need to provide two output files (with -o and -p) or use the "
                "--interleaved option.")
        input_paired_filename = None

    return input_filename, input_paired_filename


def pipeline_from_parsed_args(args, paired, is_interleaved_output):
    """
    Setup a processing pipeline from parsed command-line arguments.

    If there are any problems parsing the arguments, a CommandLineError is thrown.

    Return an instance of Pipeline (SingleEndPipeline or PairedEndPipeline)
    """

    if not paired:
        if args.untrimmed_paired_output:
            raise CommandLineError("Option --untrimmed-paired-output can only be used when "
                "trimming paired-end reads (with option -p).")

    if paired:
        if not is_interleaved_output:
            if not args.paired_output:
                raise CommandLineError("When a paired-end trimming option such as -A/-G/-B/-U, "
                    "is used, a second output file needs to be specified via -p (--paired-output).")
            if not args.output:
                raise CommandLineError("When you use -p or --paired-output, you must also "
                    "use the -o option.")

        if bool(args.untrimmed_output) != bool(args.untrimmed_paired_output):
            raise CommandLineError("When trimming paired-end reads, you must use either none "
                "or both of the --untrimmed-output/--untrimmed-paired-output options.")
        if args.too_short_output and not args.too_short_paired_output:
            raise CommandLineError("When using --too-short-output with paired-end "
                "reads, you also need to use --too-short-paired-output")
        if args.too_long_output and not args.too_long_paired_output:
            raise CommandLineError("When using --too-long-output with paired-end "
                "reads, you also need to use --too-long-paired-output")

    if args.format is not None:
        logger.warning("Option --format is deprecated and ignored because the input file format is "
            "always auto-detected")

    if not (0 <= args.error_rate < 1.):
        raise CommandLineError("The maximum error rate must be at least 0 and less than 1.")
    if args.overlap < 1:
        raise CommandLineError("The overlap must be at least 1.")
    if not (0 <= args.gc_content <= 100):
        raise CommandLineError("GC content must be given as percentage between 0 and 100")
    if args.action == 'none':
        args.action = None

    adapter_parser = AdapterParser(
        max_error_rate=args.error_rate,
        min_overlap=args.overlap,
        read_wildcards=args.match_read_wildcards,
        adapter_wildcards=args.match_adapter_wildcards,
        indels=args.indels,
    )
    try:
        adapters = adapter_parser.parse_multi(args.adapters, args.anywhere, args.front)
        adapters2 = adapter_parser.parse_multi(args.adapters2, args.anywhere2, args.front2)
    except IOError as e:
        if e.errno == errno.ENOENT:
            raise CommandLineError(e)
        raise
    except ValueError as e:
        raise CommandLineError(e)
    if args.debug:
        for adapter in adapters + adapters2:
            adapter.enable_debug()

    # Create the processing pipeline
    if paired:
        pair_filter_mode = 'any' if args.pair_filter is None else args.pair_filter
        pipeline = PairedEndPipeline(pair_filter_mode)
    else:
        pipeline = SingleEndPipeline()

    # When adapters are being trimmed only in R1 or R2, override the pair filter mode
    # as using the default of 'any' would regard all read pairs as untrimmed.
    if paired and (not adapters2 or not adapters) and (
            args.discard_untrimmed or args.untrimmed_output or args.untrimmed_paired_output):
        pipeline.override_untrimmed_pair_filter = True

    for i, cut_arg in enumerate([args.cut, args.cut2]):
        # cut_arg is a list
        if not cut_arg:
            continue
        if len(cut_arg) > 2:
            raise CommandLineError("You cannot remove bases from more than two ends.")
        if len(cut_arg) == 2 and cut_arg[0] * cut_arg[1] > 0:
            raise CommandLineError("You cannot remove bases from the same end twice.")
        for c in cut_arg:
            if c == 0:
                continue
            if i == 0:  # R1
                if paired:
                    pipeline.add(UnconditionalCutter(c), None)
                else:
                    pipeline.add(UnconditionalCutter(c))
            else:
                # R2
                assert isinstance(pipeline, PairedEndPipeline)
                pipeline.add(None, UnconditionalCutter(c))

    pipeline_add = pipeline.add_both if paired else pipeline.add

    if args.nextseq_trim is not None:
        pipeline_add(NextseqQualityTrimmer(args.nextseq_trim, args.quality_base))
    if args.quality_cutoff is not None:
        cutoffs = parse_cutoffs(args.quality_cutoff)
        pipeline_add(QualityTrimmer(cutoffs[0], cutoffs[1], args.quality_base))

    if args.pair_adapters:
        if not paired:
            raise CommandLineError("Option --pair-adapters can only be used when trimming "
                "paired-end reads")
        if args.times != 1:
            raise CommandLineError("--pair-adapters cannot be used with --times")
        try:
            cutter = PairedAdapterCutter(adapters, adapters2, args.action)
        except PairedAdapterCutterError as e:
            raise CommandLineError("--pair-adapters: " + str(e))
        pipeline.add_paired_modifier(cutter)
    else:
        adapter_cutter, adapter_cutter2 = None, None
        if adapters:
            adapter_cutter = AdapterCutter(adapters, args.times, args.action)
        if adapters2:
            adapter_cutter2 = AdapterCutter(adapters2, args.times, args.action)
        if paired:
            if adapter_cutter or adapter_cutter2:
                pipeline.add(adapter_cutter, adapter_cutter2)
        else:
            if adapter_cutter:
                pipeline.add(adapter_cutter)

    # Remaining modifiers that apply to both reads of paired-end reads
    if args.length is not None:
        pipeline_add(Shortener(args.length))
    if args.trim_n:
        pipeline_add(NEndTrimmer())
    if args.length_tag:
        pipeline_add(LengthTagModifier(args.length_tag))
    for suffix in args.strip_suffix:
        pipeline_add(SuffixRemover(suffix))
    if args.prefix or args.suffix:
        pipeline_add(PrefixSuffixAdder(args.prefix, args.suffix))
    if args.zero_cap:
        pipeline_add(ZeroCapper(quality_base=args.quality_base))

    # Set filtering parameters
    # Minimum/maximum length
    for attr in 'minimum_length', 'maximum_length':
        param = getattr(args, attr)
        if param is not None:
            lengths = parse_lengths(param)
            if not paired and len(lengths) == 2:
                raise CommandLineError('Two minimum or maximum lengths given for single-end data')
            if paired and len(lengths) == 1:
                lengths = (lengths[0], lengths[0])
            setattr(pipeline, attr, lengths)
    pipeline.max_n = args.max_n
    pipeline.discard_casava = args.discard_casava
    pipeline.discard_trimmed = args.discard_trimmed
    pipeline.discard_untrimmed = args.discard_untrimmed

    return pipeline


def log_header(cmdlineargs):
    """Print the "This is cutadapt ..." header"""

    implementation = platform.python_implementation()
    opt = ' (' + implementation + ')' if implementation != 'CPython' else ''
    logger.info("This is cutadapt %s with Python %s%s", __version__,
        platform.python_version(), opt)
    logger.info("Command line parameters: %s", " ".join(cmdlineargs))


def main(cmdlineargs=None, default_outfile=sys.stdout.buffer):
    """
    Main function that sets up a processing pipeline and runs it.

    default_outfile is the file to which trimmed reads are sent if the ``-o``
    parameter is not used.
    """
    start_time = time.time()
    parser = get_argument_parser()
    if cmdlineargs is None:
        cmdlineargs = sys.argv[1:]
    args = parser.parse_args(args=cmdlineargs)
    # log to stderr if results are to be sent to stdout
    log_to_stdout = args.output is not None and args.output != "-" and args.paired_output != "-"
    # Setup logging only if there are not already any handlers (can happen when
    # this function is being called externally such as from unit tests)
    if not logging.root.handlers:
        setup_logging(stdout=log_to_stdout,
            quiet=args.quiet or args.report == 'minimal', debug=args.debug)
    if args.profile:
        import cProfile
        profiler = cProfile.Profile()
        profiler.enable()

    if args.quiet and args.report:
        parser.error("Options --quiet and --report cannot be used at the same time")

    if args.colorspace:
        parser.error(
            "These colorspace-specific options are no longer supported: "
            "--colorspace, -c, -d, --double-encode, -t, --trim-primer, "
            "--strip-f3, --maq, --bwa, --no-zero-cap. "
            "Use Cutadapt 1.18 or earlier to work with colorspace data.")
    paired = determine_paired_mode(args)
    assert paired in (False, True)

    # Print the header now because some of the functions below create logging output
    log_header(cmdlineargs)
    try:
        is_interleaved_input, is_interleaved_output = determine_interleaved(args)
        input_filename, input_paired_filename = input_files_from_parsed_args(args.inputs,
            paired, is_interleaved_input)
        pipeline = pipeline_from_parsed_args(args, paired, is_interleaved_output)
        outfiles = open_output_files(args, default_outfile, is_interleaved_output)
    except CommandLineError as e:
        parser.error(e)
        return  # avoid IDE warnings below

    if args.cores < 0:
        parser.error('Value for --cores cannot be negative')
    cores = available_cpu_count() if args.cores == 0 else args.cores
    if cores > 1:
        if ParallelPipelineRunner.can_output_to(outfiles):
            runner_class = ParallelPipelineRunner
            runner_kwargs = dict(n_workers=cores, buffer_size=args.buffer_size)
        else:
            logger.error('Running in parallel is currently not supported for '
                'the given combination of command-line parameters.\nThese '
                'options are not supported: --info-file, --rest-file, '
                '--wildcard-file, --untrimmed-output, '
                '--untrimmed-paired-output, --too-short-output, '
                '--too-short-paired-output, --too-long-output, '
                '--too-long-paired-output, --format\n'
                'Also, demultiplexing is not supported.\n'
                'Omit --cores/-j to continue.')
            sys.exit(1)
    else:
        runner_class = SerialPipelineRunner
        runner_kwargs = dict()
    infiles = InputFiles(input_filename, file2=input_paired_filename,
            interleaved=is_interleaved_input)
    try:
        runner = runner_class(pipeline, infiles, outfiles, **runner_kwargs)
    except (dnaio.UnknownFileFormat, IOError) as e:
        parser.error(e)

    logger.info("Processing reads on %d core%s in %s mode ...",
        cores, 's' if cores > 1 else '',
        {False: 'single-end', True: 'paired-end'}[pipeline.paired])
    try:
        stats = runner.run()
        runner.close()
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        sys.exit(130)
    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(1)
        raise
    except (dnaio.FileFormatError, dnaio.UnknownFileFormat, EOFError) as e:
        sys.exit("cutadapt: error: {}".format(e))

    elapsed = time.time() - start_time
    if args.report == 'minimal':
        report = minimal_report
    else:
        report = full_report
    logger.info('%s', report(stats, elapsed, args.gc_content / 100))
    if args.profile:
        import pstats
        profiler.disable()
        pstats.Stats(profiler).sort_stats('time').print_stats(20)


if __name__ == '__main__':
    main()
