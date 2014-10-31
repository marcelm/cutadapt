#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kate: word-wrap off; remove-trailing-spaces all;
#
# Copyright (c) 2010-2014 Marcel Martin <marcel.martin@scilifelab.se>
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
Copyright (C) 2010-2014 Marcel Martin <marcel.martin@scilifelab.se>

cutadapt removes adapter sequences from high-throughput sequencing reads.

Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. The reverse complement is *not* automatically
searched. All reads from input.fastq will be written to output.fastq with the
adapter sequence removed. Adapter matching is error-tolerant. Multiple adapter
sequences can be given (use further -a options), but only the best-matching
adapter will be removed.

Input may also be in FASTA format. Compressed input and output is supported and
auto-detected from the file name (.gz or .bz2). Use the file name '-' for
standard input/output. Without the -o option, output is sent to standard output.

Some other available features are:
  * Various other adapter types (5' adapters, "mixed" 5'/3' adapters etc.)
  * Trimming a fixed number of bases
  * Quality trimming
  * Trimming paired-end reads
  * Trimming colorspace reads
  * Filtering reads by various criteria

Use "cutadapt --help" to see all command-line options.
See http://cutadapt.readthedocs.org/ for the full documentation.
"""

from __future__ import print_function, division, absolute_import

import sys
import time
import errno
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP

from cutadapt import check_importability

# Print a helpful error message if the extension modules cannot be imported.
check_importability()

from cutadapt import seqio, __version__
from cutadapt.xopen import xopen
from cutadapt.adapters import Adapter, ColorspaceAdapter, BACK, FRONT, PREFIX, SUFFIX, ANYWHERE
from cutadapt.modifiers import (LengthTagModifier, SuffixRemover, PrefixSuffixAdder,
	DoubleEncoder, ZeroCapper, PrimerTrimmer, QualityTrimmer, UnconditionalCutter)
from cutadapt.report import Statistics, print_statistics
from cutadapt.compat import next


class CutadaptOptionParser(OptionParser):
	def get_usage(self):
		return self.usage.lstrip().replace('%version', __version__)


def read_sequences(seqfilename, qualityfilename, colorspace, fileformat):
	"""
	Read sequences and (if available) quality information from either:
	* seqfilename in FASTA format (qualityfilename must be None)
	* seqfilename in FASTQ format (qualityfilename must be None)
	* seqfilename in .csfasta format and qualityfilename in .qual format
	  (SOLiD colorspace)

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


class ProcessedReadWriter(object):
	"""
	Write reads that have been processed (adapter trimming, quality trimming and
	so on) to the proper output file(s).

	Filter reads according to length and according to whether any adapter matches.
	"""

	def __init__(self,
			minimum_length,
			maximum_length,
			too_short_outfile,
			too_long_outfile,
			trimmed_outfile,
			trimmed_paired_outfile,
			untrimmed_outfile,
			untrimmed_paired_outfile):
		
		self.minimum_length = minimum_length
		self.maximum_length = maximum_length
		self.too_short_outfile = too_short_outfile
		self.too_long_outfile = too_long_outfile
		self.trimmed_outfile = trimmed_outfile
		self.untrimmed_outfile = untrimmed_outfile
		self.trimmed_paired_outfile = trimmed_paired_outfile
		self.untrimmed_paired_outfile = untrimmed_paired_outfile
		self.too_long = 0
		self.too_short = 0

	def write(self, read1, read2=None):
		"""
		Write this read to the proper file.

		If read2 is not None, this is a paired-end read.
		"""
		if len(read1.sequence) < self.minimum_length:
			self.too_short += 1
			if self.too_short_outfile is not None:
				read1.write(self.too_short_outfile)
			# TODO read2 is silently discarded
			return
		if len(read1.sequence) > self.maximum_length:
			self.too_long += 1
			if self.too_long_outfile is not None:
				read1.write(self.too_long_outfile)
			# TODO read2 is silently discarded
			return

		if read2 is None:
			# single end
			if read1.match is not None and self.trimmed_outfile:
				read1.write(self.trimmed_outfile)
			if read1.match is None and self.untrimmed_outfile:
				read1.write(self.untrimmed_outfile)
		else:
			# paired end
			if read1.match is not None:
				if self.trimmed_outfile:
					read1.write(self.trimmed_outfile)
				if self.trimmed_paired_outfile:
					read2.write(self.trimmed_paired_outfile)
			else:
				if self.untrimmed_outfile:
					read1.write(self.untrimmed_outfile)
				if self.untrimmed_paired_outfile:
					read2.write(self.untrimmed_paired_outfile)


class RestFileWriter(object):
	def __init__(self, file):
		self.file = file

	def write(self, match):
		rest = match.rest()
		if len(rest) > 0:
			print(rest, match.read.name, file=self.file)


class RepeatedAdapterCutter(object):
	"""
	Repeatedly find one of multiple adapters in reads.
	The number of times the search is repeated is specified by the
	times parameter.
	"""

	def __init__(self, adapters, times=1, wildcard_file=None, info_file=None, trim=True, rest_writer=None, mask_adapter=False):
		"""
		adapters -- list of Adapter objects

		trim -- whether to remove a found adapter from the read
		"""
		self.adapters = adapters
		self.times = times
		self.info_file = info_file
		self.wildcard_file = wildcard_file
		self.trim = trim
		self.reads_matched = 0
		self.rest_writer = rest_writer
		self.mask_adapter = mask_adapter

	def _best_match(self, read):
		"""
		Find the best matching adapter.

		read -- The read to which each adapter will be aligned

		Return an AdapterMatch instance or None if there are no matches.
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

	def _write_info(self, match):
		"""write one line to the info file"""
		# TODO move to separate class
		if not self.info_file:
			return
		seq = match.read.sequence
		if match is None:
			print(match.read.name, -1, seq, sep='\t', file=self.info_file)
		else:
			print(
				match.read.name,
				match.errors,
				match.rstart,
				match.rstop,
				seq[0:match.rstart],
				seq[match.rstart:match.rstop],
				seq[match.rstop:],
				match.adapter.name,
				sep='\t', file=self.info_file
			)

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
		for t in range(self.times):
			match = self._best_match(read)
			if match is None:
				# nothing found
				break
			self._write_info(match)  # FIXME move to cut() or somewhere else
			assert match.length > 0
			assert match.errors / match.length <= match.adapter.max_error_rate
			assert match.length - match.errors > 0

			if self.wildcard_file:  # FIXME move to cut() or somewhere else
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
		if self.trim:
			# The last match contains a copy of the read it was matched to.
			# No iteration is necessary.
			read = matches[-1].adapter.trimmed(matches[-1])

			# if an adapter was found, then the read should now be shorter
			assert len(read.sequence) < old_length

			if self.mask_adapter:
				# add N from last modification
				masked_sequence = read.sequence
				for match in sorted(matches, reverse=True, key=lambda m: m.astart):
					ns = 'N' * (len(match.read.sequence) -
								len(match.adapter.trimmed(match).sequence))
					# add N depending on match position
					if match.front:
						masked_sequence = ns + masked_sequence
					else:
						masked_sequence += ns
				# set masked sequence as sequence with original quality
				read.sequence = masked_sequence
				read.qualities = matches[0].read.qualities
				read.match = matches[-1]

				assert len(read.sequence) == old_length

		self.reads_matched += 1  # TODO move to filter class

		if self.rest_writer:
			self.rest_writer.write(matches[-1])

		return read

	def __call__(self, read):
		matches = self.find_match(read)
		if len(matches) > 0:
			read = self.cut(matches)
			read.match = matches[-1]
		else:
			read.match = None
		return read


def qtrimmed(modifiers):
	"""
	Look for a QualityTrimmer in the given list of modifiers and return its
	trimmed_bases attribute. If not found, return -1.
	"""
	for m in modifiers:
		if isinstance(m, QualityTrimmer):
			return m.trimmed_bases
	return -1


def process_single_reads(reader, modifiers, read_writer):
	"""
	Loop over reads, find adapters, trim reads, apply modifiers and
	output modified reads.

	Return a Statistics object.
	"""
	n = 0  # no. of processed reads
	total_bp = 0
	twoheaders = None
	for read in reader:
		n += 1
		total_bp += len(read.sequence)
		for modifier in modifiers:
			read = modifier(read)
		read_writer.write(read)

	return Statistics(total_bp=total_bp, n=n, quality_trimmed_bases=qtrimmed(modifiers))


def process_paired_reads(paired_reader, modifiers, read_writer):
	"""
	Loop over reads, find adapters, trim reads, apply modifiers and
	output modified reads.
	
	Note that all trimming is only done on the first read of each pair!

	Return a Statistics object.
	"""
	n = 0  # no. of processed reads
	total1_bp = 0
	total2_bp = 0
	twoheaders1, twoheaders2 = None, None
	for read1, read2 in paired_reader:
		n += 1
		total1_bp += len(read1.sequence)
		total2_bp += len(read2.sequence)
		for modifier in modifiers:
			read1 = modifier(read1)
		read_writer.write(read1, read2)

	return Statistics(total_bp=total1_bp, n=n, quality_trimmed_bases=qtrimmed(modifiers))


def parse_adapter_name(seq):
	"""
	Parse an adapter given as 'name=adapt' into 'name' and 'adapt'.
	"""
	fields = seq.split('=', 1)
	if len(fields) > 1:
		name, seq = fields
		name = name.strip()
	else:
		name = None
	seq = seq.strip()
	return name, seq


def gather_adapters(back, anywhere, front):
	"""
	Yield (name, seq, where) tuples from which Adapter instances can be built.
	This generator deals with the notation for anchored 5' adapters and also
	understands the ``file:`` syntax for reading adapters from an external FASTA
	file.
	"""
	for adapter_list, where in ((back, BACK), (anywhere, ANYWHERE), (front, FRONT)):
		for seq in adapter_list:
			if seq.startswith('file:'):
				# read adapter sequences from a file
				path = seq[5:]
				with seqio.FastaReader(path) as fasta:
					for record in fasta:
						name = record.name.split(' ', 1)[0]
						yield (name, record.sequence, where)
			else:
				name, seq = parse_adapter_name(seq)
				w = where
				if w == FRONT and seq.startswith('^'):
					seq = seq[1:]
					w = PREFIX
				elif w == BACK and seq.endswith('$'):
					seq = seq[:-1]
					w = SUFFIX
				yield (name, seq, w)


def trimmed_and_untrimmed_files(
		default_output,
		output_path,
		untrimmed_path,
		discard_trimmed,
		discard_untrimmed
		):
	"""
	Figure out (from command-line parameters) where trimmed and untrimmed reads
	should be written.

	Return a pair (trimmed, untrimmed). The values are either open file-like
	objects or None, in which case no output should be produced. The objects may
	be identical (for example: (sys.stdout, sys.stdout)).

	The parameters are sorted: Later parameters have higher precedence.

	default_output -- If nothing else is specified below, this file-like object
		is returned for both trimmed and untrimmed output.
	output_path -- Path to output file for both trimmed and untrimmed output.
	untrimmed_path -- Path to an output file for untrimmed reads.
	discard_trimmed -- bool, override earlier options.
	discard_untrimmed -- bool, overrides earlier options.
	"""
	
	if discard_trimmed:
		if discard_untrimmed:
			untrimmed = None
		elif untrimmed_path is not None:
			untrimmed = xopen(untrimmed_path, 'w')
		elif output_path is not None:
			untrimmed = xopen(output_path, 'w')
		else:
			untrimmed = default_output
		return (None, untrimmed)
	if discard_untrimmed:
		trimmed = default_output
		if output_path is not None:
			trimmed = xopen(output_path, 'w')
		return (trimmed, None)
	
	trimmed = default_output
	untrimmed = default_output
	if output_path is not None:
		trimmed = untrimmed = xopen(output_path, 'w')
	if untrimmed_path is not None:
		untrimmed = xopen(untrimmed_path, 'w')
	
	return (trimmed, untrimmed)


def get_option_parser():
	parser = CutadaptOptionParser(usage=__doc__, version=__version__)

	parser.add_option("-f", "--format", default=None,
		help="Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. "
			"Ignored when reading csfasta/qual files (default: auto-detect "
			"from file name extension).")

	group = OptionGroup(parser, "Options that influence how the adapters are found",
		description="Each of the following three parameters (-a, -b, -g) can be used "
			"multiple times and in any combination to search for an entire set of "
			"adapters of possibly different types. Only the best matching "
			"adapter will be trimmed from each read (but see the --times option). "
			"Instead of giving an adapter directly, you can also write "
			"file:FILE and the adapter sequences will be read from the given "
			"FILE (which must be in FASTA format).")
	group.add_option("-a", "--adapter", action="append", metavar="ADAPTER", dest="adapters", default=[],
		help="Sequence of an adapter that was ligated to the 3' end. The "
			"adapter itself and anything that follows is trimmed. If the "
			"adapter sequence ends with the '$' character, the adapter is "
			"anchored to the end of the read. It is only found if it is a "
			"suffix of the read.")
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
	group.add_option("--no-indels", action='store_false', dest='indels', default=True,
		help="Do not allow indels in the alignments, that is, allow only "
			"mismatches. This option is currently only supported for anchored "
			"5' adapters ('-g ^ADAPTER') (default: both mismatches and indels "
			"are allowed)")
	group.add_option("-n", "--times", type=int, metavar="COUNT", default=1,
		help="Try to remove adapters at most COUNT times. Useful when an "
			"adapter gets appended multiple times (default: %default).")
	group.add_option("-O", "--overlap", type=int, metavar="LENGTH", default=3,
		help="Minimum overlap length. If the overlap between the read and the "
			"adapter is shorter than LENGTH, the read is not modified. "
			"This reduces the no. of bases trimmed purely due to short random "
			"adapter matches (default: %default).")
	group.add_option("--match-read-wildcards", action="store_true", default=False,
		help="Allow 'N's in the read as matches to the adapter (default: %default).")
	group.add_option("-N", "--no-match-adapter-wildcards", action="store_false",
		default=True, dest='match_adapter_wildcards',
		help=SUPPRESS_HELP)
	parser.add_option_group(group)

	group = OptionGroup(parser, "Options for filtering of processed reads")
	group.add_option("--discard-trimmed", "--discard", action='store_true', default=False,
		help="Discard reads that contain the adapter instead of trimming them. "
			"Also use -O in order to avoid throwing away too many randomly "
			"matching reads!")
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
	group.add_option("--no-trim", dest='trim', action='store_false', default=True,
		help="Match and redirect reads to output/untrimmed-output as usual, "
			"but don't remove the adapters. (Default: trim the adapters)")
	group.add_option("--mask-adapter", dest='mask_adapter', action='store_true', default=False,
		help="Mask with 'N' adapter bases instead of trim (default: False)")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Options that influence what gets output to where")
	group.add_option("--quiet", default=False, action='store_true',
		help="Do not print a report at the end.")
	group.add_option("-o", "--output", default=None, metavar="FILE",
		help="Write the modified sequences to this file instead of standard "
			"output and send the summary report to standard output. "
			"The format is FASTQ if qualities are available, FASTA "
			"otherwise. (default: standard output)")
	group.add_option("-p", "--paired-output", default=None, metavar="FILE",
		help="Write reads from the paired-end input to FILE.")
	group.add_option("--info-file", metavar="FILE",
		help="Write information about each read and its adapter matches into FILE. "
			"See the documentation for the file format.")
	group.add_option("-r", "--rest-file", default=None, metavar="FILE",
		help="When the adapter matches in the middle of a read, write the "
			"rest (after the adapter) into FILE.")
	group.add_option("--wildcard-file", default=None, metavar="FILE",
		help="When the adapter has wildcard bases ('N's), write adapter bases matching wildcard "
			 "positions to FILE. When there are indels in the alignment, this may occasionally "
			 "not be quite accurate.")
	group.add_option("--too-short-output", default=None, metavar="FILE",
		help="Write reads that are too short (according to length specified by -m) to FILE. (default: discard reads)")
	group.add_option("--too-long-output", default=None, metavar="FILE",
		help="Write reads that are too long (according to length specified by -M) to FILE. (default: discard reads)")
	group.add_option("--untrimmed-output", default=None, metavar="FILE",
		help="Write reads that do not contain the adapter to FILE, instead "
			"of writing them to the regular output file. (Default: output "
			"to same file as trimmed reads.)")
	group.add_option("--untrimmed-paired-output", default=None, metavar="FILE",
		help="Write the second read in a pair to this FILE when no adapter "
			"was found in the first read. Use this option together with "
			"--untrimmed-output when trimming paired-end reads. (Default: output "
			"to same file as trimmed reads.)")
	group.add_option
	parser.add_option_group(group)

	group = OptionGroup(parser, "Additional modifications to the reads")
	group.add_option("-u", "--cut", type=int, default=0, metavar="LENGTH",
		help="Remove bases from the beginning or end of each read. "
			"If LENGTH is positive, the bases are removed from the beginning of each read. "
			"If LENGTH is negative, the bases are removed from the end of each read.")
	group.add_option("-q", "--quality-cutoff", type=int, default=0, metavar="CUTOFF",
		help="Trim low-quality bases from 3' ends of reads before adapter "
			"removal. The algorithm is the same as the one used by BWA "
			"(Subtract CUTOFF from all qualities; "
			"compute partial sums from all indices to the end of the "
			"sequence; cut sequence at the index at which the sum "
			"is minimal) (default: %default)")
	group.add_option("--quality-base", type=int, default=33,
		help="Assume that quality values are encoded as ascii(quality + QUALITY_BASE). The default (33) is usually correct, "
			 "except for reads produced by some versions of the Illumina pipeline, where this should be set to 64. (Default: %default)")
	group.add_option("-x", "--prefix", default='',
		help="Add this prefix to read names")
	group.add_option("-y", "--suffix", default='',
		help="Add this suffix to read names")
	group.add_option("--strip-suffix", action='append', default=[],
		help="Remove this suffix from read names if present. Can be given multiple times.")
	group.add_option("-c", "--colorspace", action='store_true', default=False,
		help="Colorspace mode: Also trim the color that is adjacent to the found adapter.")
	group.add_option("-d", "--double-encode", action='store_true', default=False,
		help="When in colorspace, double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).")
	group.add_option("-t", "--trim-primer", action='store_true', default=False,
		help="When in colorspace, trim primer base and the first color "
			"(which is the transition to the first nucleotide)")
	group.add_option("--strip-f3", action='store_true', default=False,
		help="For colorspace: Strip the _F3 suffix of read names")
	group.add_option("--maq", "--bwa", action='store_true', default=False,
		help="MAQ- and BWA-compatible colorspace output. This enables -c, -d, -t, --strip-f3 and -y '/1'.")
	group.add_option("--length-tag", default=None, metavar="TAG",
		help="Search for TAG followed by a decimal number in the name of the read "
			"(description/comment field of the FASTA or FASTQ file). Replace the "
			"decimal number with the correct length of the trimmed read. "
			"For example, use --length-tag 'length=' to correct fields "
			"like 'length=123'.")
	group.add_option("--no-zero-cap", dest='zero_cap', action='store_false',
		help="Do not change negative quality values to zero. Colorspace "
			"quality values of -1 would appear as spaces in the output FASTQ "
			"file. Since many tools have problems with that, negative qualities "
			"are converted to zero when trimming colorspace data. Use this "
			"option to keep negative qualities.")
	group.add_option("--zero-cap", "-z", action='store_true',
		help="Change negative quality values to zero. This is enabled "
		"by default when -c/--colorspace is also enabled. Use the above option "
		"to disable it.")
	parser.set_defaults(zero_cap=None)
	parser.add_option_group(group)

	return parser


def main(cmdlineargs=None, default_outfile=sys.stdout):
	"""
	Main function that evaluates command-line parameters and iterates
	over all reads.

	default_outfile is the default output file to which trimmed reads
	are sent. It can be overriden by using the '-o' parameter.
	"""
	parser = get_option_parser()
	if cmdlineargs is None:
		cmdlineargs = sys.argv[1:]
	options, args = parser.parse_args(args=cmdlineargs)

	if len(args) == 0:
		parser.error("At least one parameter needed: name of a FASTA or FASTQ file.")
	elif len(args) > 2:
		parser.error("Too many parameters.")

	input_filename = args[0]

	# If a second file name was given, then we either have single-end reads
	# provided as a pair of .fasta/.qual files or we have paired-end reads.
	quality_filename = None
	input_paired_filename = None
	if len(args) == 2:
		if args[1].endswith('.qual'):
			quality_filename = args[1]
		else:
			input_paired_filename = args[1]
			if not options.paired_output:
				parser.error('You must use --paired-output when trimming paired-end reads.')

	if len(args) == 1 and options.paired_output:
		parser.error("You specified a --paired-output file, but gave only one input file.")
	if options.paired_output and bool(options.untrimmed_output) != bool(options.untrimmed_paired_output):
		parser.error("When trimming paired-end reads, you must use either none "
			"or both of the --untrimmed-output/--untrimmed-paired-output options.")
	if options.untrimmed_paired_output and not options.paired_output:
		parser.error("Option --untrimmed-paired-output can only be used when "
			"trimming paired-end reads (with option --paired-output).")
	if input_filename.endswith('.qual') and quality_filename.endswith('fasta'):
		parser.error("FASTA and QUAL file given, but the FASTA file must be first.")

	if options.format is not None and options.format.lower() not in ['fasta', 'fastq', 'sra-fastq']:
		parser.error("The input file format must be either 'fasta', 'fastq' or 'sra-fastq' (not '{0}').".format(options.format))

	# TODO should this really be an error?
	if options.format is not None and quality_filename is not None:
		parser.error("If a pair of .fasta and .qual files is given, the -f/--format parameter cannot be used.")

	trimmed_outfile, untrimmed_outfile = trimmed_and_untrimmed_files(
		default_outfile,
		options.output,
		options.untrimmed_output,
		options.discard_trimmed,
		options.discard_untrimmed)

	trimmed_paired_outfile, untrimmed_paired_outfile = trimmed_and_untrimmed_files(
		None,  # applies when not trimming paired-end data
		options.paired_output,
		options.untrimmed_paired_output,
		options.discard_trimmed,
		options.discard_untrimmed)

	too_short_outfile = None  # too short reads go here
	if options.too_short_output is not None:
		too_short_outfile = xopen(options.too_short_output, 'w')
	too_long_outfile = None  # too long reads go here
	if options.too_long_output is not None:
		too_long_outfile = xopen(options.too_long_output, 'w')

	read_writer = ProcessedReadWriter(
		options.minimum_length,
		options.maximum_length,
		too_short_outfile,
		too_long_outfile,
		trimmed_outfile, trimmed_paired_outfile,
		untrimmed_outfile, untrimmed_paired_outfile
	)

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
	adapters = []
	ADAPTER_CLASS = ColorspaceAdapter if options.colorspace else Adapter

	try:
		for name, seq, where in gather_adapters(options.adapters, options.anywhere, options.front):
			if not seq:
				parser.error("The adapter sequence is empty")
			if not options.indels and where not in (PREFIX, SUFFIX):
				parser.error("Not allowing indels is currently supported only for anchored 5' or 3' adapters.")
			adapter = ADAPTER_CLASS(seq, where, options.error_rate,
				options.overlap, options.match_read_wildcards,
				options.match_adapter_wildcards, name=name, indels=options.indels)
			adapters.append(adapter)
	except IOError as e:
		if e.errno == errno.ENOENT:
			print("Error:", e, file=sys.stderr)
			sys.exit(1)
		raise

	start_time = time.clock()

	if input_paired_filename:
		reader = seqio.PairedSequenceReader(input_filename, input_paired_filename,
			colorspace=options.colorspace, fileformat=options.format)
	else:
		reader = read_sequences(input_filename, quality_filename,
			colorspace=options.colorspace, fileformat=options.format)

	# build up list of modifiers
	modifiers = []
	if not adapters and options.quality_cutoff == 0 and options.cut == 0 and \
			options.minimum_length == 0 and options.maximum_length == sys.maxint:
		parser.error("You need to provide at least one adapter sequence.")
	if options.cut:
		modifiers.append(UnconditionalCutter(options.cut))
	if options.quality_cutoff > 0:
		modifiers.append(QualityTrimmer(options.quality_cutoff, options.quality_base))
	if adapters:
		adapter_cutter = RepeatedAdapterCutter(adapters, options.times,
				options.wildcard_file, options.info_file, options.trim,
				rest_writer, options.mask_adapter)
		modifiers.append(adapter_cutter)
	else:
		adapter_cutter = None
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
	if options.zero_cap and reader.delivers_qualities:
		modifiers.append(ZeroCapper(quality_base=options.quality_base))
	if options.trim_primer:
		modifiers.append(PrimerTrimmer)

	try:
		if input_paired_filename:
			stats = process_paired_reads(reader, modifiers, read_writer)
		else:
			stats = process_single_reads(reader, modifiers, read_writer)
	except KeyboardInterrupt as e:
		print("Interrupted", file=sys.stderr)
		sys.exit(1)
	except IOError as e:
		if e.errno == errno.EPIPE:
			sys.exit(1)
		raise
	except seqio.FormatError as e:
		print("Error:", e, file=sys.stderr)
		sys.exit(1)

	# close open files
	for f in [trimmed_outfile, untrimmed_outfile, trimmed_paired_outfile,
			untrimmed_paired_outfile, options.rest_file, options.wildcard_file, options.info_file,
			too_short_outfile, too_long_outfile, options.info_file]:
		if f is not None and f is not sys.stdin and f is not sys.stdout:
			f.close()

	if not options.quiet:
		# send statistics to stderr if result was sent to stdout
		stat_file = sys.stderr if options.output is None else None
		print_statistics(adapters, time.clock() - start_time, stats,
			options.trim, adapter_cutter.reads_matched if adapter_cutter else 0,
			options.error_rate, read_writer.too_short, read_writer.too_long, cmdlineargs, file=stat_file)


if __name__ == '__main__':
	main()
