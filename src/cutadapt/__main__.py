#!/usr/bin/env python
#
# Copyright (c) 2010-2020 Marcel Martin <marcel.martin@scilifelab.se>
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

Copyright (C) 2010-2020 Marcel Martin <marcel.martin@scilifelab.se>

cutadapt removes adapter sequences from high-throughput sequencing reads.

Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. All reads from input.fastq will be written to
output.fastq with the adapter sequence removed. Adapter matching is
error-tolerant. Multiple adapter sequences can be given (use further -a
options), but only the best-matching adapter will be removed.

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
import time
import logging
import platform
from typing import Tuple, Optional, Sequence, List, Any, Iterator, Union, Type
from argparse import ArgumentParser, SUPPRESS, HelpFormatter

import dnaio

from cutadapt import __version__
from cutadapt.adapters import warn_duplicate_adapters, Adapter
from cutadapt.parser import AdapterParser
from cutadapt.modifiers import (SingleEndModifier, LengthTagModifier, SuffixRemover, PrefixSuffixAdder,
    ZeroCapper, QualityTrimmer, UnconditionalCutter, NEndTrimmer, AdapterCutter,
    PairedAdapterCutterError, PairedAdapterCutter, NextseqQualityTrimmer, Shortener,
    ReverseComplementer)
from cutadapt.report import full_report, minimal_report
from cutadapt.pipeline import (Pipeline, SingleEndPipeline, PairedEndPipeline, InputFiles,
    OutputFiles, PipelineRunner, SerialPipelineRunner, ParallelPipelineRunner)
from cutadapt.utils import available_cpu_count, Progress, DummyProgress, FileOpener
from cutadapt.log import setup_logging, REPORT

logger = logging.getLogger()


class CutadaptArgumentParser(ArgumentParser):
    """
    This ArgumentParser customizes two things:
    - The usage message is not prefixed with 'usage:'
    - A brief message is shown on errors, not full usage
    """
    class CustomUsageHelpFormatter(HelpFormatter):
        def __init__(self, *args, **kwargs):
            try:
                import shutil
                kwargs['width'] = min(24 + 80, shutil.get_terminal_size().columns)
            except ImportError:
                pass
            super().__init__(*args, **kwargs)

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


def get_argument_parser() -> ArgumentParser:
    # noqa: E131
    parser = CutadaptArgumentParser(usage=__doc__, add_help=False)
    group = parser.add_argument_group("Options")
    group.add_argument("-h", "--help", action="help", help="Show this help message and exit")
    group.add_argument("--version", action="version", help="Show version number and exit",
        version=__version__)
    group.add_argument("--debug", nargs="?", const=True, default=False,
        choices=("trace",),
        help="Print debug log. 'trace' prints also DP matrices")
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
    # Compression level for gzipped output files. Not exposed since we have -Z
    group.add_argument("--compression-level", type=int, default=6,
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
    group.add_argument("-a", "--adapter", type=lambda x: ("back", x), action="append",
        default=[], metavar="ADAPTER", dest="adapters",
        help="Sequence of an adapter ligated to the 3' end (paired data: of the "
            "first read). The adapter and subsequent bases are trimmed. If a "
            "'$' character is appended ('anchoring'), the adapter is only "
            "found if it is a suffix of the read.")
    group.add_argument("-g", "--front", type=lambda x: ("front", x), action="append",
        default=[], metavar="ADAPTER", dest="adapters",
        help="Sequence of an adapter ligated to the 5' end (paired data: of the "
            "first read). The adapter and any preceding bases are trimmed. "
            "Partial matches at the 5' end are allowed. If a '^' character is "
            "prepended ('anchoring'), the adapter is only found if it is a "
            "prefix of the read.")
    group.add_argument("-b", "--anywhere", type=lambda x: ("anywhere", x), action="append",
        default=[], metavar="ADAPTER", dest="adapters",
        help="Sequence of an adapter that may be ligated to the 5' or 3' end "
            "(paired data: of the first read). Both types of matches as "
            "described under -a and -g are allowed. If the first base of the "
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
            "--discard-untrimmed). Default: %(default)s")
    group.add_argument("--rc", "--revcomp", dest="reverse_complement", default=False,
        action="store_true",
        help="Check both the read and its reverse complement for adapter matches. If "
            "match is on reverse-complemented version, output that one. "
            "Default: check only read")
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
    group.add_argument("--max-expected-errors", "--max-ee", type=float, default=None,
        metavar="ERRORS",
        help="Discard reads whose expected number of errors (computed "
            "from quality values) exceeds ERRORS.")
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
            "depending on input. Summary report is sent to standard output. "
            "Use '{name}' for demultiplexing (see docs). "
            "Default: write to standard output")
    group.add_argument("--fasta", default=False, action='store_true',
        help="Output FASTA to standard output even on FASTQ input.")
    group.add_argument("-Z", action="store_const", const=1, dest="compression_level",
        help="Use compression level 1 for gzipped output files (faster, but uses more space)")
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
    group.add_argument("-A", type=lambda x: ("back", x), dest='adapters2',
        action='append', default=[], metavar='ADAPTER',
        help="3' adapter to be removed from second read in a pair.")
    group.add_argument("-G", type=lambda x: ("front", x), dest='adapters2',
        action='append', default=[], metavar='ADAPTER',
        help="5' adapter to be removed from second read in a pair.")
    group.add_argument("-B", type=lambda x: ("anywhere", x), dest='adapters2',
        action='append', default=[], metavar='ADAPTER',
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
        help="Read and/or write interleaved paired-end reads.")
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


def parse_cutoffs(s: str) -> Tuple[int, int]:
    """Parse a string INT[,INT] into a pair of integers

    >>> parse_cutoffs("5")
    (0, 5)
    >>> parse_cutoffs("6,7")
    (6, 7)
    """
    try:
        cutoffs = [int(value) for value in s.split(",")]
    except ValueError as e:
        raise CommandLineError("Quality cutoff value not recognized: {}".format(e))

    if len(cutoffs) == 1:
        cutoffs = [0, cutoffs[0]]
    elif len(cutoffs) != 2:
        raise CommandLineError("Expected one value or two values separated by comma for "
            "the quality cutoff")

    return (cutoffs[0], cutoffs[1])


def parse_lengths(s: str) -> Tuple[Optional[int], ...]:
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


def open_output_files(
    args, default_outfile, file_opener: FileOpener
) -> OutputFiles:
    """
    Return an OutputFiles instance. If demultiplex is True, the untrimmed, untrimmed2, out and out2
    attributes are not opened files, but paths (out and out2 with the '{name}' template).
    """

    rest_file = file_opener.xopen_or_none(args.rest_file, "wb")
    info_file = file_opener.xopen_or_none(args.info_file, "wb")
    wildcard = file_opener.xopen_or_none(args.wildcard_file, "wb")

    too_short = too_short2 = None
    if args.minimum_length is not None:
        too_short, too_short2 = file_opener.xopen_pair(
            args.too_short_output, args.too_short_paired_output, "wb")

    too_long = too_long2 = None
    if args.maximum_length is not None:
        too_long, too_long2 = file_opener.xopen_pair(
            args.too_long_output, args.too_long_paired_output, "wb")

    if int(args.discard_trimmed) + int(args.discard_untrimmed) + int(
            args.untrimmed_output is not None) > 1:
        raise CommandLineError("Only one of the --discard-trimmed, --discard-untrimmed "
            "and --untrimmed-output options can be used at the same time.")

    demultiplex_mode = determine_demultiplex_mode(args)
    if demultiplex_mode and args.discard_trimmed:
        raise CommandLineError("Do not use --discard-trimmed when demultiplexing.")

    if demultiplex_mode == "normal":
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

        assert out is not None and '{name}' in out and (out2 is None or '{name}' in out2)
    elif demultiplex_mode == "combinatorial":
        out = args.output
        out2 = args.paired_output
        if args.untrimmed_output or args.untrimmed_paired_output:
            raise CommandLineError("Combinatorial demultiplexing (with {name1} and {name2})"
                " cannot be combined with --untrimmed-output or --untrimmed-paired-output")
        if args.discard_untrimmed:
            untrimmed = untrimmed2 = None
        else:
            untrimmed = untrimmed2 = 'unknown'
    else:
        untrimmed, untrimmed2 = file_opener.xopen_pair(
            args.untrimmed_output, args.untrimmed_paired_output, "wb")
        out, out2 = file_opener.xopen_pair(args.output, args.paired_output, "wb")
        if out is None:
            out = default_outfile

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
        demultiplex=bool(demultiplex_mode),
        force_fasta=args.fasta,
    )


def determine_demultiplex_mode(args) -> Union[str, bool]:
    """Return one of "normal", "combinatorial" or False"""

    demultiplex = args.output is not None and '{name}' in args.output

    if args.paired_output is not None and (demultiplex != ('{name}' in args.paired_output)):
        raise CommandLineError('When demultiplexing paired-end data, "{name}" must appear in '
                               'both output file names (-o and -p)')

    demultiplex_combinatorial = (
        args.output is not None
        and args.paired_output is not None
        and '{name1}' in args.output
        and '{name2}' in args.output
        and '{name1}' in args.paired_output
        and '{name2}' in args.paired_output
    )
    if demultiplex and demultiplex_combinatorial:
        raise CommandLineError("You cannot combine {name} with {name1} and {name2}")

    if demultiplex:
        return "normal"
    elif demultiplex_combinatorial:
        return "combinatorial"
    else:
        return False


def determine_paired(args) -> bool:
    """
    Determine whether we should work in paired-end mode.
    """
    # Usage of any of these options enables paired-end mode
    return bool(
        args.paired_output
        or args.interleaved
        or args.adapters2
        or args.cut2
        or args.pair_filter
        or args.untrimmed_paired_output
        or args.too_short_paired_output
        or args.too_long_paired_output)


def setup_input_files(
    inputs: Sequence[str], paired: bool, interleaved: bool
) -> Tuple[str, Optional[str]]:
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
            + "\n - ".join("'{}'".format(p) for p in inputs)
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
            input_paired_filename = inputs[1]  # type: Optional[str]
    else:
        if len(inputs) == 2:
            raise CommandLineError(
                "It appears you want to trim paired-end data because you provided two input files, "
                "but then you also need to provide two output files (with -o and -p) or use the "
                "--interleaved option.")
        input_paired_filename = None

    return input_filename, input_paired_filename


def check_arguments(args, paired: bool) -> None:
    if not paired:
        if args.untrimmed_paired_output:
            raise CommandLineError("Option --untrimmed-paired-output can only be used when "
                "trimming paired-end reads.")

        if args.pair_adapters:
            raise CommandLineError("Option --pair-adapters can only be used when trimming "
                "paired-end reads")

    if paired and not args.interleaved:
        if not args.paired_output:
            raise CommandLineError("When a paired-end trimming option such as -A/-G/-B/-U, "
                "is used, a second output file needs to be specified via -p (--paired-output).")
        if not args.output:
            raise CommandLineError("When you use -p or --paired-output, you must also "
                "use the -o option.")
        for out, paired_out, argname in [
            (args.untrimmed_output, args.untrimmed_paired_output, "untrimmed"),
            (args.too_short_output, args.too_short_paired_output, "too-short"),
            (args.too_long_output, args.too_long_paired_output, "too-long"),
        ]:
            if bool(out) != bool(paired_out):
                raise CommandLineError(
                    "When trimming paired-end data, you must use either none or both of the"
                    " --{name}-output/--{name}-paired-output options.".format(name=argname)
                )

    if args.format is not None:
        logger.warning("Option --format is deprecated and ignored because the input file format is "
            "always auto-detected")

    if not (0 <= args.error_rate < 1.):
        raise CommandLineError("The maximum error rate must be at least 0 and less than 1.")
    if args.overlap < 1:
        raise CommandLineError("The overlap must be at least 1.")
    if not (0 <= args.gc_content <= 100):
        raise CommandLineError("GC content must be given as percentage between 0 and 100")

    if args.pair_adapters and args.times != 1:
        raise CommandLineError("--pair-adapters cannot be used with --times")


def pipeline_from_parsed_args(args, paired, file_opener) -> Pipeline:
    """
    Setup a processing pipeline from parsed command-line arguments.

    If there are any problems parsing the arguments, a CommandLineError is raised.

    Return an instance of Pipeline (SingleEndPipeline or PairedEndPipeline)
    """
    if args.action == 'none':
        args.action = None

    adapters, adapters2 = adapters_from_args(args)

    # Create the processing pipeline
    if paired:
        pair_filter_mode = 'any' if args.pair_filter is None else args.pair_filter
        pipeline = PairedEndPipeline(
            pair_filter_mode, file_opener
        )  # type: Any
    else:
        pipeline = SingleEndPipeline(file_opener)

    # When adapters are being trimmed only in R1 or R2, override the pair filter mode
    # as using the default of 'any' would regard all read pairs as untrimmed.
    if isinstance(pipeline, PairedEndPipeline) and (not adapters2 or not adapters) and (
            args.discard_untrimmed or args.untrimmed_output or args.untrimmed_paired_output):
        pipeline.override_untrimmed_pair_filter = True

    add_unconditional_cutters(pipeline, args.cut, args.cut2, paired)

    pipeline_add = pipeline.add_both if paired else pipeline.add

    if args.nextseq_trim is not None:
        pipeline_add(NextseqQualityTrimmer(args.nextseq_trim, args.quality_base))
    if args.quality_cutoff is not None:
        cutoffs = parse_cutoffs(args.quality_cutoff)
        pipeline_add(QualityTrimmer(cutoffs[0], cutoffs[1], args.quality_base))

    add_adapter_cutter(
        pipeline,
        adapters,
        adapters2,
        paired,
        args.pair_adapters,
        args.action,
        args.times,
        args.reverse_complement,
    )

    for modifier in modifiers_applying_to_both_ends_if_paired(args):
        pipeline_add(modifier)

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
    pipeline.max_expected_errors = args.max_expected_errors
    pipeline.discard_casava = args.discard_casava
    pipeline.discard_trimmed = args.discard_trimmed
    pipeline.discard_untrimmed = args.discard_untrimmed

    return pipeline


def adapters_from_args(args) -> Tuple[List[Adapter], List[Adapter]]:
    adapter_parser = AdapterParser(
        max_error_rate=args.error_rate,
        min_overlap=args.overlap,
        read_wildcards=args.match_read_wildcards,
        adapter_wildcards=args.match_adapter_wildcards,
        indels=args.indels,
    )
    try:
        adapters = adapter_parser.parse_multi(args.adapters)
        adapters2 = adapter_parser.parse_multi(args.adapters2)
    except (FileNotFoundError, ValueError) as e:
        raise CommandLineError(e)
    warn_duplicate_adapters(adapters)
    warn_duplicate_adapters(adapters2)
    if args.debug == "trace":
        for adapter in adapters + adapters2:
            adapter.enable_debug()
    return adapters, adapters2


def add_unconditional_cutters(pipeline: Pipeline, cut1: List[int], cut2: List[int], paired: bool):
    for i, cut_arg in enumerate([cut1, cut2]):
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
                    assert isinstance(pipeline, PairedEndPipeline)
                    pipeline.add(UnconditionalCutter(c), None)
                else:
                    assert isinstance(pipeline, SingleEndPipeline)
                    pipeline.add(UnconditionalCutter(c))
            else:
                # R2
                assert isinstance(pipeline, PairedEndPipeline)
                pipeline.add(None, UnconditionalCutter(c))


def add_adapter_cutter(
    pipeline,
    adapters,
    adapters2,
    paired: bool,
    pair_adapters: bool,
    action: str,
    times: int,
    reverse_complement: bool,
):
    if pair_adapters:
        if reverse_complement:
            raise CommandLineError("Cannot use --revcomp with --pair-adapters")
        try:
            cutter = PairedAdapterCutter(adapters, adapters2, action)
        except PairedAdapterCutterError as e:
            raise CommandLineError("--pair-adapters: " + str(e))
        pipeline.add_paired_modifier(cutter)
    else:
        adapter_cutter, adapter_cutter2 = None, None
        if adapters:
            adapter_cutter = AdapterCutter(adapters, times, action)
        if adapters2:
            adapter_cutter2 = AdapterCutter(adapters2, times, action)
        if paired:
            if reverse_complement:
                raise CommandLineError("--revcomp not implemented for paired-end reads")
            if adapter_cutter or adapter_cutter2:
                pipeline.add(adapter_cutter, adapter_cutter2)
        elif adapter_cutter:
            if reverse_complement:
                modifier = ReverseComplementer(
                    adapter_cutter
                )  # type: Union[AdapterCutter,ReverseComplementer]
            else:
                modifier = adapter_cutter
            pipeline.add(modifier)


def modifiers_applying_to_both_ends_if_paired(args) -> Iterator[SingleEndModifier]:
    if args.length is not None:
        yield Shortener(args.length)
    if args.trim_n:
        yield NEndTrimmer()
    if args.length_tag:
        yield LengthTagModifier(args.length_tag)
    for suffix in args.strip_suffix:
        yield SuffixRemover(suffix)
    if args.prefix or args.suffix:
        yield PrefixSuffixAdder(args.prefix, args.suffix)
    if args.zero_cap:
        yield ZeroCapper(quality_base=args.quality_base)


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
    args, leftover_args = parser.parse_known_args(args=cmdlineargs)
    # log to stderr if results are to be sent to stdout
    log_to_stdout = args.output is not None and args.output != "-" and args.paired_output != "-"
    # Setup logging only if there are not already any handlers (can happen when
    # this function is being called externally such as from unit tests)
    if not logging.root.handlers:
        setup_logging(logger, stdout=log_to_stdout,
            quiet=args.quiet, minimal=args.report == 'minimal', debug=args.debug)
    profiler = setup_profiler_if_requested(args.profile)

    if args.quiet and args.report:
        parser.error("Options --quiet and --report cannot be used at the same time")

    if args.colorspace:
        parser.error(
            "These colorspace-specific options are no longer supported: "
            "--colorspace, -c, -d, --double-encode, -t, --trim-primer, "
            "--strip-f3, --maq, --bwa, --no-zero-cap. "
            "Use Cutadapt 1.18 or earlier to work with colorspace data.")

    paired = determine_paired(args)
    assert paired in (False, True)

    # Print the header now because some of the functions below create logging output
    log_header(cmdlineargs)
    if leftover_args:
        warn_if_en_dashes(cmdlineargs)
        parser.error("unrecognized arguments: " + " ".join(leftover_args))

    if args.cores < 0:
        parser.error('Value for --cores cannot be negative')

    cores = available_cpu_count() if args.cores == 0 else args.cores
    file_opener = FileOpener(
        compression_level=args.compression_level, threads=0 if cores == 1 else 1)
    if sys.stderr.isatty() and not args.quiet:
        progress = Progress()
    else:
        progress = DummyProgress()

    try:
        is_interleaved_input = args.interleaved and len(args.inputs) == 1
        input_filename, input_paired_filename = setup_input_files(args.inputs,
            paired, is_interleaved_input)
        check_arguments(args, paired)
        pipeline = pipeline_from_parsed_args(args, paired, file_opener)
        outfiles = open_output_files(args, default_outfile, file_opener)
        infiles = InputFiles(input_filename, file2=input_paired_filename,
            interleaved=is_interleaved_input)
        runner = setup_runner(pipeline, infiles, outfiles, progress, cores, args.buffer_size)
    except CommandLineError as e:
        parser.error(str(e))
        return  # avoid IDE warnings below

    logger.info("Processing reads on %d core%s in %s mode ...",
        cores, 's' if cores > 1 else '',
        {False: 'single-end', True: 'paired-end'}[pipeline.paired])
    try:
        with runner as r:
            stats = r.run()
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        sys.exit(130)
    except BrokenPipeError:
        sys.exit(1)
    except (dnaio.FileFormatError, dnaio.UnknownFileFormat, EOFError) as e:
        sys.exit("cutadapt: error: {}".format(e))

    elapsed = time.time() - start_time
    if args.report == 'minimal':
        report = minimal_report
    else:
        report = full_report
    logger.log(REPORT, '%s', report(stats, elapsed, args.gc_content / 100))
    if profiler is not None:
        import pstats
        profiler.disable()
        pstats.Stats(profiler).sort_stats('time').print_stats(20)


def setup_runner(pipeline: Pipeline, infiles, outfiles, progress, cores, buffer_size):
    if cores > 1:
        if ParallelPipelineRunner.can_output_to(outfiles):
            runner_class = ParallelPipelineRunner  # type: Type[PipelineRunner]
            runner_kwargs = dict(n_workers=cores, buffer_size=buffer_size)
        else:
            raise CommandLineError("Running in parallel is currently not supported "
                 "when using --format or when demultiplexing.\n"
                 "Omit --cores/-j to continue.")
            # return  # avoid IDE warnings below
    else:
        runner_class = SerialPipelineRunner
        runner_kwargs = dict()
    try:
        runner = runner_class(pipeline, infiles, outfiles, progress, **runner_kwargs)
    except (dnaio.UnknownFileFormat, dnaio.FileFormatError, OSError) as e:
        raise CommandLineError(e)

    return runner


def setup_profiler_if_requested(requested):
    if requested:
        import cProfile
        profiler = cProfile.Profile()
        profiler.enable()
    else:
        profiler = None
    return profiler


def warn_if_en_dashes(args):
    for arg in args:
        if arg.startswith("–"):
            logger.warning(
                "The first character in argument '%s' is '–' (an en-dash, Unicode U+2013)"
                " and will therefore be interpreted as a file name. If you wanted to"
                " provide an option, use a regular hyphen '-'.", arg
            )


if __name__ == '__main__':
    main()
