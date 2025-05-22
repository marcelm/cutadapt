#!/usr/bin/env python
#
# Copyright (c) 2010 Marcel Martin <marcel.martin@scilifelab.se> and contributors
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

Copyright (C) 2010 Marcel Martin <marcel.martin@scilifelab.se> and contributors

Cutadapt removes adapter sequences from high-throughput sequencing reads.

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
import copy
import sys
import time
import shutil
import logging
import platform
import itertools
import multiprocessing
from pathlib import Path
from typing import Tuple, Optional, Sequence, List, Iterator, Union, Dict
from argparse import ArgumentParser, SUPPRESS, HelpFormatter

import dnaio
import xopen as xopen

from cutadapt import __version__
from cutadapt.adapters import warn_duplicate_adapters, Adapter, InvalidCharacter
from cutadapt.json import OneLine, dumps as json_dumps
from cutadapt.parser import make_adapters_from_specifications
from cutadapt.modifiers import (
    SingleEndModifier,
    LengthTagModifier,
    SuffixRemover,
    PrefixSuffixAdder,
    ZeroCapper,
    QualityTrimmer,
    UnconditionalCutter,
    NEndTrimmer,
    AdapterCutter,
    PairedAdapterCutterError,
    PairedAdapterCutter,
    NextseqQualityTrimmer,
    Shortener,
    ReverseComplementer,
    PairedEndRenamer,
    Renamer,
    InvalidTemplate,
    PolyATrimmer,
    PairedReverseComplementer,
)
from cutadapt.predicates import (
    TooShort,
    TooLong,
    TooManyN,
    TooManyExpectedErrors,
    TooHighAverageErrorRate,
    CasavaFiltered,
    IsTrimmed,
    IsUntrimmed,
)
from cutadapt.report import full_report, minimal_report, Statistics
from cutadapt.pipeline import SingleEndPipeline, PairedEndPipeline
from cutadapt.runners import make_runner
from cutadapt.files import InputPaths, OutputFiles, FileOpener
from cutadapt.steps import (
    InfoFileWriter,
    PairedInfoFileWriter,
    PairedSingleEndStep,
    RestFileWriter,
    WildcardFileWriter,
    SingleEndFilter,
    PairedEndFilter,
    Demultiplexer,
    CombinatorialDemultiplexer,
    PairedDemultiplexer,
    PairedEndSink,
    SingleEndSink,
)
from cutadapt.utils import available_cpu_count, Progress, DummyProgress
from cutadapt.log import setup_logging, REPORT
from cutadapt.qualtrim import HasNoQualities

logger = logging.getLogger()


class CutadaptArgumentParser(ArgumentParser):
    """
    This ArgumentParser customizes two things:
    - The usage message is not prefixed with 'usage:'
    - A brief message is shown on errors, not full usage
    """

    class CustomUsageHelpFormatter(HelpFormatter):
        def __init__(self, *args, **kwargs):
            kwargs["width"] = min(24 + 80, shutil.get_terminal_size().columns)
            super().__init__(*args, **kwargs)

        def add_usage(self, usage, actions, groups, prefix=None):
            if usage is not SUPPRESS:  # pragma: no cover
                args = usage, actions, groups, ""
                self._add_item(self._format_usage, args)

    def __init__(self, *args, **kwargs):
        kwargs["formatter_class"] = self.CustomUsageHelpFormatter
        kwargs["usage"] = kwargs["usage"].replace("{version}", __version__)
        super().__init__(*args, **kwargs)

    def error(self, message):
        """
        If you override this in a subclass, it should not return -- it
        should either exit or raise an exception.
        """
        print('Run "cutadapt --help" to see command-line options.', file=sys.stderr)
        print(
            "See https://cutadapt.readthedocs.io/ for full documentation.",
            file=sys.stderr,
        )
        self.exit(2, f"\n{self.prog}: error: {message}\n")


class CommandLineError(Exception):
    pass


# fmt: off
def get_argument_parser() -> ArgumentParser:
    # noqa: E131
    parser = CutadaptArgumentParser(usage=__doc__, add_help=False)
    group = parser.add_argument_group("Options")
    group.add_argument("-h", "--help", action="help", help="Show this help message and exit")
    group.add_argument("--version", action="version", help="Show version number and exit",
        version=__version__)
    group.add_argument("--debug", action="count", default=0,
        help="Print debug log. Use twice to also print DP matrices")
    group.add_argument("--profile", action="store_true", default=False, help=SUPPRESS)
    group.add_argument("-j", "--cores", type=int, default=1,
        help='Number of CPU cores to use. Use 0 to auto-detect. Default: %(default)s')

    # Hidden options
    # GC content as a percentage
    group.add_argument("--gc-content", type=float, default=50,
        help=SUPPRESS)
    # Buffer size for the reader process when running in parallel
    group.add_argument("--buffer-size", type=int, default=4000000,
        help=SUPPRESS)
    # Disable adapter index creation
    group.add_argument("--no-index", dest="index", default=True, action="store_false", help=SUPPRESS)

    group = parser.add_argument_group("Finding adapters",
        description="Parameters -a, -g, -b specify adapters to be removed from "
            "each read (or from R1 if data is paired-end. "
            "If specified multiple times, only the best matching adapter is "
            "trimmed (but see the --times option). Use notation "
            "'file:FILE' to read adapter sequences from a FASTA file.")
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
    group.add_argument("-e", "--error-rate", "--errors",
        type=float, metavar="E", default=0.1,
        help="Maximum allowed error rate (if 0 <= E < 1), or absolute number of errors "
            "for full-length adapter match (if E is an integer >= 1). Error rate = "
            "no. of errors divided by length of matching region. Default: %(default)s (10%%)")
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
        default=True, dest="match_adapter_wildcards",
        help="Do not interpret IUPAC wildcards in adapters.")
    group.add_argument("--action", choices=("trim", "retain", "mask", "lowercase", "crop", "none"),
        default="trim",
        help="What to do if a match was found. "
            "trim: trim adapter and up- or downstream sequence; "
            "retain: trim, but retain adapter; "
            "mask: replace with 'N' characters; "
            "lowercase: convert to lowercase; "
            "crop: trim up and downstream sequence; "
            "none: leave unchanged. Default: %(default)s")
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
    group.add_argument("-u", "--cut", action='append', default=[], type=int, metavar="LEN",
        help="Remove LEN bases from each read (or R1 if paired; use -U option for R2). "
            "If LEN is positive, remove bases from the beginning. "
            "If LEN is negative, remove bases from the end. "
            "Can be used twice if LENs have different signs. "
            "Applied *before* adapter trimming.")
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
    group.add_argument("--poly-a", action="store_true", default=False,
        help="Trim poly-A tails")
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
    group.add_argument("--rename", metavar="TEMPLATE",
        help="Rename reads using TEMPLATE containing variables such as {id}, {adapter_name} "
            "etc. (see documentation)")
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
    group.add_argument("--max-average-error-rate", "--max-aer", type=float, default=None,
        metavar="ERROR_RATE",
        help="as --max-expected-errors (see above), but divided by length to "
             "account for reads of varying length.")
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
    group.add_argument("--json", metavar="FILE",
        help="Dump report in JSON format to FILE")
    group.add_argument("-o", "--output", metavar="FILE",
        help="Write trimmed reads to FILE. FASTQ or FASTA format is chosen "
            "depending on input. Summary report is sent to standard output. "
            "Use '{name}' for demultiplexing (see docs). "
            "Default: write to standard output")
    group.add_argument("--fasta", default=False, action='store_true',
        help="Output FASTA to standard output even on FASTQ input.")
    group.add_argument("--compression-level", type=int, default=1, metavar="N",
        help="Compression level for compressed output files. Default: %(default)s")
    group.add_argument("-Z", action="store_const", const=1, dest="compression_level",
        help=SUPPRESS)  # deprecated because compression level 1 is now the default
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
        "-A/-G/-B/-U/-Q options work like their lowercase counterparts, but "
        "are applied to R2 (second read in pair)")
    group.add_argument("-A", type=lambda x: ("back", x), dest='adapters2',
        action='append', default=[], metavar='ADAPTER',
        help="3' adapter to be removed from R2")
    group.add_argument("-G", type=lambda x: ("front", x), dest='adapters2',
        action='append', default=[], metavar='ADAPTER',
        help="5' adapter to be removed from R2")
    group.add_argument("-B", type=lambda x: ("anywhere", x), dest='adapters2',
        action='append', default=[], metavar='ADAPTER',
        help="5'/3 adapter to be removed from R2")
    group.add_argument("-U", dest='cut2', action='append', default=[], type=int, metavar="LENGTH",
        help="Remove LENGTH bases from R2")
    group.add_argument("-Q", dest="quality_cutoff2", default=None, metavar="[5'CUTOFF,]3'CUTOFF",
        help="Quality-trimming cutoff for R2. Default: same as for R1")
    group.add_argument("-L", dest="length2", type=int, default=None, metavar="LENGTH",
        help="Shorten R2 to LENGTH. Default: same as for R1")
    group.add_argument("-p", "--paired-output", metavar="FILE",
        help="Write R2 to FILE.")
    group.add_argument("--info-file-paired", dest="info_file2", metavar="FILE",
        help="Write info about R2 to FILE (see --info-file)")
    group.add_argument("--pair-adapters", action="store_true",
        help="Treat adapters given with -a/-A etc. as pairs. Either both "
             "or none are removed from each read pair.")
    # Setting the default for pair_filter to None allows us to find out whether
    # the option was used at all.
    group.add_argument("--pair-filter", default=None,
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
        help="Write second read in a pair to this file if pair is too short.")
    group.add_argument("--too-long-paired-output", metavar="FILE", default=None,
        help="Write second read in a pair to this file if pair is too long.")

    # We could have two positional arguments here, with the second one optional, but
    # we want custom, more helpful error messages.
    parser.add_argument("inputs", nargs='*', help=SUPPRESS)

    return parser
# fmt: on


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
        raise CommandLineError(f"Quality cutoff value not recognized: {e}")

    if len(cutoffs) == 1:
        cutoffs = [0, cutoffs[0]]
    elif len(cutoffs) != 2:
        raise CommandLineError(
            "Expected one value or two values separated by comma for "
            "the quality cutoff"
        )

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
    fields = s.split(":")
    if len(fields) not in (1, 2):
        raise CommandLineError("Only at most one colon is allowed")
    try:
        values = tuple(int(f) if f != "" else None for f in fields)
    except ValueError as e:
        raise CommandLineError(f"Value not recognized: {e}")
    if len(values) == 2 and values[0] is None and values[1] is None:
        raise CommandLineError(
            f"Cannot parse '{s}': At least one length needs to be given"
        )
    return tuple(values)


def complain_about_duplicate_paths(paths: List[str]):
    if sys.platform == "win32" and sys.version_info < (3, 8):
        # Bug in handling of NUL
        return
    seen = set()
    for path in paths:
        if path is None:
            continue
        p = Path(path)
        if p.exists() and not p.is_file():
            # assumed to be FIFO, /dev/null etc.
            continue
        if path in seen:
            raise CommandLineError(
                f"Path {path} specified more than once as an output file. "
                f"This is not supported at the moment."
            )
        seen.add(path)


def determine_demultiplex_mode(
    output: Optional[str], paired_output: Optional[str]
) -> Union[str, bool]:
    """Return one of "normal", "combinatorial" or False"""

    demultiplex = output is not None and "{name}" in output

    if paired_output is not None and (demultiplex != ("{name}" in paired_output)):
        raise CommandLineError(
            'When demultiplexing paired-end data, "{name}" must appear in '
            "both output file names (-o and -p)"
        )

    demultiplex_combinatorial = (
        output is not None
        and paired_output is not None
        and "{name1}" in output
        and "{name2}" in output
        and "{name1}" in paired_output
        and "{name2}" in paired_output
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
    # Any of these options enable paired-end mode
    return bool(
        args.paired_output
        or args.interleaved
        or args.adapters2
        or args.cut2
        or args.length2
        or args.pair_filter
        or args.untrimmed_paired_output
        or args.too_short_paired_output
        or args.too_long_paired_output
        or args.quality_cutoff2
        or args.info_file2
    )


def make_input_paths(
    inputs: Sequence[str], paired: bool, interleaved: bool
) -> InputPaths:
    """
    Do some other error checking of the input file names and return InputPaths.
    """
    if len(inputs) == 0:
        raise CommandLineError(
            "You did not provide any input file names. Please give me something to do!"
        )
    elif len(inputs) > 2:
        raise CommandLineError(
            f"You provided {len(inputs)} input file names, but either one or two are expected. "
            + "The file names were:\n - "
            + "\n - ".join(f"'{p}'" for p in inputs)
            + "\nHint: If your path contains spaces, you need to enclose it in quotes"
        )
    input_filename = inputs[0]
    if paired and not interleaved:
        # Two file names required
        if len(inputs) == 1:
            raise CommandLineError(
                "You used an option that enables paired-end mode (such as -p, -A, -G, -B, -U), but "
                "only provided one input file. Please either provide two input files or use "
                "use --interleaved as appropriate."
            )
        else:
            input_paired_filename = inputs[1]  # type: Optional[str]
    else:
        if len(inputs) == 2:
            raise CommandLineError(
                "It appears you want to trim paired-end data because you provided two input files, "
                "but then you also need to provide two output files (with -o and -p) or use the "
                "--interleaved option."
            )
        input_paired_filename = None

    if input_paired_filename:
        return InputPaths(
            input_filename, input_paired_filename, interleaved=interleaved
        )
    else:
        return InputPaths(input_filename, interleaved=interleaved)


def check_arguments(args, paired: bool) -> None:
    if not paired:
        if args.untrimmed_paired_output:
            raise CommandLineError(
                "Option --untrimmed-paired-output can only be used when "
                "trimming paired-end reads."
            )

        if args.pair_adapters:
            raise CommandLineError(
                "Option --pair-adapters can only be used when trimming "
                "paired-end reads"
            )

    if paired and not args.interleaved:
        if not args.paired_output:
            raise CommandLineError(
                "When a paired-end trimming option such as -A/-G/-B/-U, "
                "is used, a second output file needs to be specified via -p (--paired-output)."
            )
        if not args.output:
            raise CommandLineError(
                "When you use -p or --paired-output, you must also "
                "use the -o option."
            )
        for out, paired_out, argname in [
            (args.untrimmed_output, args.untrimmed_paired_output, "untrimmed"),
            (args.too_short_output, args.too_short_paired_output, "too-short"),
            (args.too_long_output, args.too_long_paired_output, "too-long"),
        ]:
            if bool(out) != bool(paired_out):
                raise CommandLineError(
                    "When trimming paired-end data, you must use either none or both of the"
                    " --{name}-output/--{name}-paired-output options.".format(
                        name=argname
                    )
                )

    if args.overlap < 1:
        raise CommandLineError("The overlap must be at least 1.")
    if not (0 <= args.gc_content <= 100):
        raise CommandLineError(
            "GC content must be given as percentage between 0 and 100"
        )

    if args.pair_adapters and args.times != 1:
        raise CommandLineError("--pair-adapters cannot be used with --times")


def make_pipeline_from_args(  # noqa: C901
    args, input_file_format, outfiles, paired, adapters, adapters2
):
    """
    Set up a processing pipeline from parsed command-line arguments.

    If there are any problems parsing the arguments, raise a CommandLineError.

    Return an instance of Pipeline (SingleEndPipeline or PairedEndPipeline)
    """
    action = None if args.action == "none" else args.action
    pair_filter_mode = None
    if paired:
        pair_filter_mode = "any" if args.pair_filter is None else args.pair_filter

    def make_filter(
        predicate1, predicate2, path1, path2, pair_filter_mode=pair_filter_mode
    ):
        record_writer = None
        if path1 or path2:
            paths = [path1, path2] if paired else [path1]
            if paired and path2 is None:
                interleaved = True
                paths = paths[:1]
            else:
                interleaved = False
            record_writer = outfiles.open_record_writer(*paths, interleaved=interleaved)
        if paired:
            step = PairedEndFilter(
                predicate1, predicate2, record_writer, pair_filter_mode=pair_filter_mode
            )
        else:
            step = SingleEndFilter(predicate1, record_writer)
        return step

    adapter_names: List[Optional[str]] = [a.name for a in adapters]
    adapter_names2: List[Optional[str]] = [a.name for a in adapters2]

    steps = []

    if args.rest_file is not None:
        step = RestFileWriter(outfiles.open_text(args.rest_file))
        if paired:
            step = PairedSingleEndStep(step)
        steps.append(step)

    if args.info_file is not None:
        if paired and args.info_file2 is not None:
            step = PairedInfoFileWriter(
                outfiles.open_text(args.info_file), outfiles.open_text(args.info_file2)
            )
        else:
            step = InfoFileWriter(outfiles.open_text(args.info_file))
            if paired:
                step = PairedSingleEndStep(step)
        steps.append(step)

    if args.wildcard_file is not None:
        step = WildcardFileWriter(outfiles.open_text(args.wildcard_file))
        if paired:
            step = PairedSingleEndStep(step)
        steps.append(step)

    # Add filtering steps

    for length, path1, path2, predicate_class in [
        (
            args.minimum_length,
            args.too_short_output,
            args.too_short_paired_output,
            TooShort,
        ),
        (
            args.maximum_length,
            args.too_long_output,
            args.too_long_paired_output,
            TooLong,
        ),
    ]:
        if length is None:
            if path1 or path2:
                if predicate_class is TooShort:
                    raise CommandLineError(
                        "When --too-short-output or --too-short-paired-output are used, "
                        "a minimum length must be provided with -m/--minimum-length"
                    )
                if predicate_class is TooLong:
                    raise CommandLineError(
                        "When --too-long-output or --too-long-paired-output are used, "
                        "a maximum length must be provided with -M/--maximum-length"
                    )
            continue
        if not paired and path2:
            raise CommandLineError(
                "--too-short/long-paired-output cannot be used with single-end data"
            )
        lengths = parse_lengths(length)
        if not paired and len(lengths) == 2:
            raise CommandLineError(
                "Two minimum or maximum lengths given for single-end data"
            )
        if paired and len(lengths) == 1:
            lengths = (lengths[0], lengths[0])
        predicate1 = predicate_class(lengths[0]) if lengths[0] is not None else None
        if len(lengths) == 2 and lengths[1] is not None:
            predicate2 = predicate_class(lengths[1])
        else:
            predicate2 = None

        steps.append(make_filter(predicate1, predicate2, path1, path2))

    if args.max_n is not None:
        predicate = TooManyN(args.max_n)
        if paired:
            step = PairedEndFilter(
                predicate, predicate, pair_filter_mode=pair_filter_mode
            )
        else:
            step = SingleEndFilter(predicate)
        steps.append(step)

    if args.max_expected_errors is not None:
        if not input_file_format.has_qualities():
            logger.warning(
                "Ignoring option --max-ee because input does not provide quality values"
            )
        else:
            predicate = TooManyExpectedErrors(args.max_expected_errors)
            if paired:
                step = PairedEndFilter(
                    predicate, predicate, pair_filter_mode=pair_filter_mode
                )
            else:
                step = SingleEndFilter(predicate)
            steps.append(step)

    if args.max_average_error_rate is not None:
        if not input_file_format.has_qualities():
            logger.warning(
                "Ignoring option --max-er because input does not contain quality values"
            )
        else:
            predicate = TooHighAverageErrorRate(args.max_average_error_rate)
            if paired:
                step = PairedEndFilter(
                    predicate, predicate, pair_filter_mode=pair_filter_mode
                )
            else:
                step = SingleEndFilter(predicate)
            steps.append(step)

    if args.discard_casava:
        predicate = CasavaFiltered()
        if paired:
            step = PairedEndFilter(
                predicate, predicate, pair_filter_mode=pair_filter_mode
            )
        else:
            step = SingleEndFilter(predicate)
        steps.append(step)

    # Add the last step that writes the records that made it through the pipeline
    # to the final output(s)

    if (
        int(args.discard_trimmed)
        + int(args.discard_untrimmed)
        + int(
            args.untrimmed_output is not None
            or args.untrimmed_paired_output is not None
        )
        > 1
    ):
        raise CommandLineError(
            "Only one of the --discard-trimmed, --discard-untrimmed "
            "and --untrimmed-output options can be used at the same time."
        )

    demultiplex_mode = determine_demultiplex_mode(args.output, args.paired_output)
    if demultiplex_mode and args.discard_trimmed:
        raise CommandLineError("Do not use --discard-trimmed when demultiplexing.")
    if demultiplex_mode == "combinatorial" and args.pair_adapters:
        raise CommandLineError(
            "With --pair-adapters, you can only use {name} in your output file name template, "
            "not {name1} and {name2} (no combinatorial demultiplexing)."
        )
    if demultiplex_mode == "normal":
        if paired:
            step = PairedDemultiplexer(
                adapter_names,
                template1=args.output,
                template2=args.paired_output,
                untrimmed_output=args.untrimmed_output,
                untrimmed_paired_output=args.untrimmed_paired_output,
                discard_untrimmed=args.discard_untrimmed,
                outfiles=outfiles,
            )
        else:
            step = Demultiplexer(
                adapter_names,
                template=args.output,
                untrimmed_output=args.untrimmed_output,
                discard_untrimmed=args.discard_untrimmed,
                outfiles=outfiles,
            )
        steps.append(step)
    elif demultiplex_mode == "combinatorial":
        assert "{name1}" in args.output and "{name2}" in args.output
        assert "{name1}" in args.paired_output and "{name2}" in args.paired_output
        if args.untrimmed_output or args.untrimmed_paired_output:
            raise CommandLineError(
                "Combinatorial demultiplexing (with {name1} and {name2})"
                " cannot be combined with --untrimmed-output or --untrimmed-paired-output"
            )
        step = CombinatorialDemultiplexer(
            adapter_names,
            adapter_names2,
            template1=args.output,
            template2=args.paired_output,
            discard_untrimmed=args.discard_untrimmed,
            outfiles=outfiles,
        )
        steps.append(step)
    else:
        # When adapters are being trimmed only in R1 or R2, override the pair filter mode
        # as using the default of 'any' would regard all read pairs as untrimmed.
        override_pair_filter_mode = (
            paired
            and (not adapters2 or not adapters)
            and (
                args.discard_untrimmed
                or args.untrimmed_output
                or args.untrimmed_paired_output
            )
        )

        # Set up the remaining filters to deal with --discard-trimmed,
        # --discard-untrimmed and --untrimmed-output. These options
        # are mutually exclusive to help prevent brain damage.
        if args.discard_trimmed:
            predicate = IsTrimmed()
            if paired:
                step = PairedEndFilter(
                    predicate, predicate, pair_filter_mode=pair_filter_mode
                )
            else:
                step = SingleEndFilter(predicate)

            steps.append(step)
        elif args.discard_untrimmed:
            predicate = IsUntrimmed()
            if paired:
                step = PairedEndFilter(
                    predicate,
                    predicate,
                    pair_filter_mode=(
                        "both" if override_pair_filter_mode else pair_filter_mode
                    ),
                )
            else:
                step = SingleEndFilter(predicate)
            steps.append(step)
        elif args.untrimmed_output or args.untrimmed_paired_output:
            predicate1 = IsUntrimmed()
            predicate2 = IsUntrimmed()
            steps.append(
                make_filter(
                    predicate1,
                    predicate2 if paired else None,
                    args.untrimmed_output,
                    args.untrimmed_paired_output,
                    pair_filter_mode=(
                        "both" if override_pair_filter_mode else pair_filter_mode
                    ),
                )
            )

        if paired:
            paths = [args.output, args.paired_output]
            if args.paired_output is None:
                interleaved = True
                paths = paths[:1]
            else:
                interleaved = False
            steps.append(
                PairedEndSink(
                    outfiles.open_record_writer(*paths, interleaved=interleaved)
                )
            )
        else:
            if args.output is None:
                out = outfiles.open_stdout_record_writer(
                    interleaved=paired and args.interleaved,
                    force_fasta=args.fasta,
                )
            else:
                out = outfiles.open_record_writer(args.output, force_fasta=args.fasta)
            steps.append(SingleEndSink(out))

    logger.debug("Pipeline steps:")
    for step in steps:
        logger.debug("- %s", step)
    modifiers = []
    modifiers.extend(make_unconditional_cutters(args.cut, args.cut2, paired))

    if args.nextseq_trim is not None:
        trimmer = NextseqQualityTrimmer(args.nextseq_trim, args.quality_base)
        if paired:
            modifiers.append((trimmer, copy.copy(trimmer)))
        else:
            modifiers.append(trimmer)

    modifiers.extend(
        make_quality_trimmers(
            args.quality_cutoff,
            args.quality_cutoff2,
            args.quality_base,
            paired,
        )
    )
    modifiers.extend(
        make_adapter_cutter(
            adapters,
            adapters2,
            paired,
            args.pair_adapters,
            action,
            args.times,
            args.reverse_complement,
            not args.rename,  # no "rc" suffix if --rename is used
            args.index,
        )
    )

    if args.poly_a:
        if paired:
            modifiers.append((PolyATrimmer(), PolyATrimmer(revcomp=True)))
        else:
            modifiers.append(PolyATrimmer())

    modifiers.extend(make_shortener(args.length, args.length2, paired))
    for modifier in modifiers_applying_to_both_ends_if_paired(args):
        if paired:
            modifiers.append((modifier, copy.copy(modifier)))
        else:
            modifiers.append(modifier)

    if args.rename and (args.prefix or args.suffix):
        raise CommandLineError(
            "Option --rename cannot be combined with --prefix (-x) or --suffix (-y)"
        )
    if args.rename and args.rename != "{header}":
        try:
            renamer = PairedEndRenamer(args.rename) if paired else Renamer(args.rename)
            modifiers.append(renamer)
        except InvalidTemplate as e:
            raise CommandLineError(e)

    # Create the processing pipeline
    if paired:
        pipeline = PairedEndPipeline(modifiers, steps)  # type: Any
    else:
        pipeline = SingleEndPipeline(modifiers, steps)

    return pipeline


def adapters_from_args(args) -> Tuple[List[Adapter], List[Adapter]]:
    search_parameters = dict(
        max_errors=args.error_rate,
        min_overlap=args.overlap,
        read_wildcards=args.match_read_wildcards,
        adapter_wildcards=args.match_adapter_wildcards,
        indels=args.indels,
    )
    try:
        adapters = make_adapters_from_specifications(args.adapters, search_parameters)
        adapters2 = make_adapters_from_specifications(args.adapters2, search_parameters)
    except (
        KeyError,
        ValueError,
        InvalidCharacter,
    ) as e:
        raise CommandLineError(e.args[0])
    warn_duplicate_adapters(adapters)
    warn_duplicate_adapters(adapters2)
    if args.debug > 1:
        for adapter in adapters + adapters2:
            adapter.enable_debug()
    return adapters, adapters2


def make_unconditional_cutters(cut1: List[int], cut2: List[int], paired: bool):
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
                    yield (UnconditionalCutter(c), None)
                else:
                    yield UnconditionalCutter(c)
            else:
                # R2
                assert paired
                yield (None, UnconditionalCutter(c))


def make_quality_trimmers(
    cutoff1: Optional[str],
    cutoff2: Optional[str],
    quality_base: int,
    paired: bool,
):
    qtrimmers = [
        (
            QualityTrimmer(*parse_cutoffs(cutoff), quality_base)
            if cutoff is not None and cutoff != "0"
            else None
        )
        for cutoff in (cutoff1, cutoff2)
    ]
    if paired:
        if cutoff1 is not None and cutoff2 is None:
            qtrimmers[1] = copy.copy(qtrimmers[0])
        if qtrimmers[0] is not None or qtrimmers[1] is not None:
            yield tuple(qtrimmers)
    elif qtrimmers[0] is not None:
        assert not paired
        yield qtrimmers[0]


def make_adapter_cutter(
    adapters,
    adapters2,
    paired: bool,
    pair_adapters: bool,
    action: Optional[str],
    times: int,
    reverse_complement: bool,
    add_rc_suffix: bool,
    allow_index: bool,
):
    if pair_adapters:
        if reverse_complement:
            raise CommandLineError("Cannot use --revcomp with --pair-adapters")
        try:
            cutter = PairedAdapterCutter(adapters, adapters2, action)
        except PairedAdapterCutterError as e:
            raise CommandLineError("--pair-adapters: " + str(e))
        yield cutter
    else:
        adapter_cutter, adapter_cutter2 = None, None
        try:
            if adapters:
                adapter_cutter = AdapterCutter(adapters, times, action, allow_index)
            if adapters2:
                adapter_cutter2 = AdapterCutter(adapters2, times, action, allow_index)
        except ValueError as e:
            raise CommandLineError(e)
        if paired:
            if adapter_cutter or adapter_cutter2:
                if reverse_complement:
                    yield PairedReverseComplementer(
                        adapter_cutter,
                        adapter_cutter2,
                        rc_suffix=" rc" if add_rc_suffix else None,
                    )
                else:
                    yield (adapter_cutter, adapter_cutter2)
        elif adapter_cutter:
            if reverse_complement:
                yield ReverseComplementer(
                    adapter_cutter,
                    rc_suffix=" rc" if add_rc_suffix else None,
                )
            else:
                yield adapter_cutter


def make_shortener(length1: Optional[int], length2: Optional[int], paired: bool):
    if paired:
        if length1 is not None and length2 is not None:
            yield Shortener(length1), Shortener(length2)
        elif length1 is not None and length2 is None:
            # If -L not given, use same setting for both
            yield Shortener(length1), Shortener(length1)
        elif length1 is None and length2 is not None:
            yield None, Shortener(length2)
    else:
        if length1 is not None:
            yield Shortener(length1)


def modifiers_applying_to_both_ends_if_paired(args) -> Iterator[SingleEndModifier]:
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
    opt = " (" + implementation + ")" if implementation != "CPython" else ""
    logger.info(
        "This is cutadapt %s with Python %s%s",
        __version__,
        platform.python_version(),
        opt,
    )
    logger.info("Command line parameters: %s", " ".join(cmdlineargs))


def main_cli():  # pragma: no cover
    """Entry point for command-line script"""
    multiprocessing.freeze_support()
    main(sys.argv[1:])
    return 0


def main(cmdlineargs) -> Statistics:
    """
    Set up a processing pipeline from the command-line arguments, run it and return
    a Statistics object.
    """
    start_time = time.time()
    parser = get_argument_parser()
    if not cmdlineargs:
        parser.print_usage()
        sys.exit(2)

    cmdlineargs = [str(arg) if isinstance(arg, Path) else arg for arg in cmdlineargs]
    args, leftover_args = parser.parse_known_args(args=cmdlineargs)
    # Setup logging only if there are not already any handlers (can happen when
    # this function is being called externally such as from unit tests)
    if not logging.root.handlers:
        setup_logging(
            logger,
            log_to_stderr=is_any_output_stdout(args),
            quiet=args.quiet,
            minimal=args.report == "minimal",
            debug=args.debug,
        )
    log_header(cmdlineargs)
    profiler = setup_profiler_if_requested(args.profile)

    log_system_info()
    if args.quiet and args.report:
        parser.error("Options --quiet and --report cannot be used at the same time")

    if leftover_args:
        warn_if_en_dashes(cmdlineargs)
        parser.error("unrecognized arguments: " + " ".join(leftover_args))

    if args.cores < 0:
        parser.error("Value for --cores cannot be negative")

    cores = available_cpu_count() if args.cores == 0 else args.cores
    file_opener = FileOpener(
        compression_level=args.compression_level,
        threads=estimate_compression_threads(cores),
    )
    if sys.stderr.isatty() and not args.quiet and not args.debug:
        progress = Progress()
    else:
        progress = DummyProgress()
    paired = determine_paired(args)

    try:
        is_interleaved_input = args.interleaved and len(args.inputs) == 1
        input_paths = make_input_paths(args.inputs, paired, is_interleaved_input)
        check_arguments(args, paired)
        adapters, adapters2 = adapters_from_args(args)
        log_adapters(adapters, adapters2 if paired else None)
        complain_about_duplicate_paths(
            [
                args.rest_file,
                args.info_file,
                args.wildcard_file,
                args.too_short_output,
                args.too_short_paired_output,
                args.too_long_output,
                args.too_long_paired_output,
                args.untrimmed_output,
                args.untrimmed_paired_output,
                args.output,
                args.paired_output,
            ]
        )

        with make_runner(input_paths, cores, args.buffer_size) as runner:
            outfiles = OutputFiles(
                proxied=cores > 1,
                qualities=runner.input_file_format().has_qualities(),
                file_opener=file_opener,
                interleaved=args.interleaved,
            )
            pipeline = make_pipeline_from_args(
                args,
                runner.input_file_format(),
                outfiles,
                paired,
                adapters,
                adapters2,
            )
            logger.info(
                "Processing %s reads on %d core%s ...",
                {False: "single-end", True: "paired-end"}[pipeline.paired],
                cores,
                "s" if cores > 1 else "",
            )
            stats = runner.run(pipeline, progress, outfiles)
    except KeyboardInterrupt:
        if args.debug:
            raise
        else:
            print("Interrupted", file=sys.stderr)
            sys.exit(130)
    except BrokenPipeError:
        sys.exit(1)
    except (
        OSError,
        EOFError,
        HasNoQualities,
        dnaio.UnknownFileFormat,
        dnaio.FileFormatError,
        CommandLineError,
    ) as e:
        logger.debug("Command line error. Traceback:", exc_info=True)
        logger.error("%s", e)
        exit_code = 2 if isinstance(e, CommandLineError) else 1
        sys.exit(exit_code)
    finally:
        # TODO ...
        try:
            outfiles.close()
        except UnboundLocalError:
            pass

    elapsed = time.time() - start_time
    if args.report == "minimal":
        report = minimal_report
    else:
        report = full_report
    logger.log(REPORT, "%s", report(stats, elapsed, args.gc_content / 100.0))
    if args.json is not None:
        with open(args.json, "w") as f:
            json_dict = json_report(
                stats=stats,
                cmdlineargs=cmdlineargs,
                path1=input_paths.paths[0],
                path2=input_paths.paths[1] if len(input_paths.paths) > 1 else None,
                cores=cores,
                paired=paired,
                gc_content=args.gc_content / 100.0,
            )
            f.write(json_dumps(json_dict))
            f.write("\n")
    if profiler is not None:
        import pstats

        profiler.disable()
        pstats.Stats(profiler).sort_stats("time").print_stats(20)
    return stats


def log_system_info():
    logger.debug("Python executable: %s", sys.executable)
    logger.debug("dnaio version: %s", dnaio.__version__)
    logger.debug("xopen version: %s", xopen.__version__)


def log_adapters(adapters, adapters2):
    paired = adapters2 is not None
    logger.debug("R1 adapters (%d):" if paired else "Adapters (%d):", len(adapters))
    for a in itertools.islice(adapters, 20):
        logger.debug("- %s", a)
    if len(adapters) > 20:
        logger.debug("- (%d more)", len(adapters) - 20)
    if paired:
        logger.debug("R2 adapters (%d):", len(adapters2))
        for a in itertools.islice(adapters2, 20):
            logger.debug("- %s", a)
        if len(adapters2) > 20:
            logger.debug("- (%d more)", len(adapters2) - 20)


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
                " provide an option, use a regular hyphen '-'.",
                arg,
            )


def estimate_compression_threads(cores: int) -> Optional[int]:
    return max(0, min(cores - 1, 4))


def is_any_output_stdout(args):
    return any(
        [
            args.output is None,
            args.output == "-",
            args.paired_output == "-",
            args.untrimmed_output == "-",
            args.untrimmed_paired_output == "-",
            args.too_short_output == "-",
            args.too_short_paired_output == "-",
            args.too_long_output == "-",
            args.too_long_paired_output == "-",
            args.rest_file == "-",
            args.info_file == "-",
            args.wildcard_file == "-",
        ]
    )


def json_report(
    stats: Statistics,
    cmdlineargs: List[str],
    path1: str,
    path2: Optional[str],
    cores: int,
    paired: bool,
    gc_content: float,
) -> Dict:
    d = {
        "tag": "Cutadapt report",
        "schema_version": OneLine([0, 3]),
        "cutadapt_version": __version__,
        "python_version": platform.python_version(),
        "command_line_arguments": cmdlineargs,
        "cores": cores,
        "input": {
            "path1": path1,
            "path2": path2,
            "paired": paired,
        },
    }
    d.update(stats.as_json(gc_content, one_line=True))
    return d


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main_cli())
