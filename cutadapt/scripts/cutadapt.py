#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kate: word-wrap off;
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
import gzip
import time
import errno
from os.path import dirname, join, isfile, realpath
if sys.version_info[0] < 3:
	from string import maketrans
else:
	maketrans = bytes.maketrans
	xrange = range
from optparse import OptionParser, OptionGroup
from contextlib import closing
from collections import defaultdict, namedtuple

from .. import align, seqio, __version__
from ..colorspace import encode as colorspace_encode
from ..xopen import xopen
from ..qualtrim import quality_trim_index

# Constants for the find_best_alignment function.
# The function is called with SEQ1 as the adapter, SEQ2 as the read.
BACK = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ1
FRONT = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.START_WITHIN_SEQ1
PREFIX = align.STOP_WITHIN_SEQ2
ANYWHERE = align.SEMIGLOBAL


class HelpfulOptionParser(OptionParser):
	"""An OptionParser that prints full help on errors."""
	def error(self, msg):
		self.print_help(sys.stderr)
		self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


def print_histogram(d, adapter_length, n, error_rate):
	"""
	Print a histogram. Also, print the no. of reads expected to be
	trimmed by chance (assuming a uniform distribution of nucleotides in the reads).
	d -- a dictionary mapping lengths of trimmed sequences to their respective frequency
	adapter_length -- adapter length
	n -- total no. of reads.
	"""
	print("length", "count", "expected", sep="\t")
	h = []
	for length in sorted(d):
		# when length surpasses adapter_length, the
		# probability does not increase anymore
		estimated = n * 0.25 ** min(length, adapter_length)
		h.append( (length, d[length], estimated) )

	for length, count, estimate in h:
		print(length, count, "{0:.1F}".format(estimate), sep="\t")
	print()


class Statistics(object):
	"""Store statistics about reads and adapters"""

	def __init__(self, adapters):
		self.reads_changed = 0
		self.too_short = 0
		self.too_long = 0
		self.n = 0
		self.total_bp = 0
		self._start_time = time.clock()
		self.time = None
		self.adapters = adapters

	def stop_clock(self):
		"""Stop the timer that was automatically started when the class was instantiated."""
		self.time = time.clock() - self._start_time

	def print_statistics(self, error_rate, file=None):
		"""Print summary to file"""
		if self.time is None:
			self.stop_clock()
		old_stdout = sys.stdout
		if file is not None:
			sys.stdout = file
		print("cutadapt version", __version__)
		print("Command line parameters:", " ".join(sys.argv[1:]))
		print("Maximum error rate: %.2f%%" % (error_rate * 100.))
		print("   Processed reads:", self.n)

		trimmed_bp = 0
		for adapter in self.adapters:
			for d in (adapter.lengths_front, adapter.lengths_back):
				trimmed_bp += sum( seqlen*count for (seqlen, count) in d.items() )

		if self.n > 0:
			print("     Trimmed reads:", self.reads_changed, "(%5.1f%%)" % (100. * self.reads_changed / self.n))
			print("   Total basepairs: {0:12} ({1:.1F} Mbp)".format(self.total_bp, self.total_bp/1E6))
			s = " ({0:.2%} of total)".format(float(trimmed_bp)/self.total_bp) if self.total_bp > 0 else ''
			print(" Trimmed basepairs: {0:12} ({1:.1F} Mbp){2}".format(trimmed_bp, trimmed_bp/1E6, s))
			print("   Too short reads:", self.too_short, "(%5.1f%% of processed reads)" % (100. * self.too_short / self.n))
			print("    Too long reads:", self.too_long, "(%5.1f%% of processed reads)" % (100. * self.too_long / self.n))
		print("        Total time: %9.2f s" % self.time)
		if self.n > 0:
			print("     Time per read: %9.2f ms" % (1000. * self.time / self.n))
		print()
		for index, adapter in enumerate(self.adapters):
			total_front = sum(adapter.lengths_front.values())
			total_back = sum(adapter.lengths_back.values())
			total = total_front + total_back
			where = adapter.where
			assert where == ANYWHERE or (where == BACK and total_front == 0) or (where in (FRONT, PREFIX) and total_back == 0)

			print("=" * 3, "Adapter", index+1, "=" * 3)
			print()
			print("Adapter '%s'," % adapter.sequence, "length %d," % len(adapter.sequence), "was trimmed", total, "times.")
			if where == ANYWHERE:
				print(total_front, "times, it overlapped the 5' end of a read")
				print(total_back, "times, it overlapped the 3' end or was within the read")
				print()
				print("Lengths of removed sequences (5')")
				print_histogram(adapter.lengths_front, len(adapter), self.n, error_rate)
				print()
				print("Lengths of removed sequences (3' or within)")
				print_histogram(adapter.lengths_back, len(adapter), self.n, error_rate)
			elif where in (FRONT, PREFIX):
				print()
				print("Lengths of removed sequences")
				print_histogram(adapter.lengths_front, len(adapter), self.n, error_rate)
			else:
				assert where == BACK
				print()
				print("Lengths of removed sequences")
				print_histogram(adapter.lengths_back, len(adapter), self.n, error_rate)

		if self.n == 0:
			print("No reads were read! Either your input file is empty or you used the wrong -f/--format parameter.")
		sys.stdout = old_stdout


AdapterMatch = namedtuple('AdapterMatch', ['astart', 'astop', 'rstart', 'rstop', 'matches', 'errors', 'adapter', 'read'])


class Adapter(object):
	"""
	An adapter knows how to match itself to a read.
	In particular, it knows where it should be within the read and how to interpret
	wildcard characters.

	where --  One of the BACK, FRONT, PREFIX or ANYWHERE constants.
		If the adapter is located in the middle of the read,
		the constant influences which part of the read gets removed.

	sequence -- The adapter sequence as string. Will be converted to uppercase.

	max_error_rate -- Maximum allowed error rate. The error rate is
		the number of errors in the alignment divided by the length
		of the part of the alignment that matches the adapter.

	minimum_overlap -- Minimum length of the part of the alignment
		that matches the adapter.

	match_read_wildcards -- Whether wildcards ('N' characters) in the read
		are matches (at zero cost).

	match_adapter_wildcards -- Whether wildcards in the adapter are allowed
		to match any character in the read (at zero cost).
	"""
	def __init__(self, sequence, where, max_error_rate, min_overlap,
			match_read_wildcards, match_adapter_wildcards,
			rest_file):
		self.sequence = sequence.upper()
		self.where = where
		self.max_error_rate = max_error_rate
		self.min_overlap = min_overlap
		self.wildcard_flags = 0
		if match_read_wildcards:
			self.wildcard_flags |= align.ALLOW_WILDCARD_SEQ2
		if match_adapter_wildcards and 'N' in self.sequence:
			self.wildcard_flags |= align.ALLOW_WILDCARD_SEQ1
		removers = {
			FRONT: self.remove_front,
			PREFIX: self.remove_front,
			BACK: self.remove_back,
			ANYWHERE: self.remove_anywhere
		}
		self.remove = removers[where]
		self.rest_file = rest_file
		# statistics about length of removed sequences
		self.lengths_front = defaultdict(int)
		self.lengths_back = defaultdict(int)

	def __repr__(self):
		return '<Adapter(sequence="{0}", where={1})>'.format(self.sequence, self.where)

	def match(self, read):
		"""
		Try to match this adapter to the given read and return an AdapterMatch instance.

		Return None if the minimum overlap length is not met or the error rate is too high.
		"""
		# try to find an exact match first
		# TODO do not do this when wildcards are allowed!
		read_seq = read.sequence.upper()
		if self.where == PREFIX:
			pos = 0 if read_seq.startswith(self.sequence) else -1
		else:
			pos = read_seq.find(self.sequence)
		if pos >= 0:
			match = AdapterMatch(0, len(self.sequence), pos, pos + len(self.sequence), len(self.sequence), 0, self, read)
		else:
			alignment = align.globalalign_locate(self.sequence, read_seq,
				self.max_error_rate, self.where, self.wildcard_flags)
			match = AdapterMatch(*(alignment + (self, read)))
		length = match.astop - match.astart
		# TODO globalalign_locate should be modified to allow the following
		# assertion.
		# assert length == 0 or match.errors / length <= self.max_error_rate
		if length < self.min_overlap or match.errors / length > self.max_error_rate:
			return None
		return match

	def remove_anywhere(self, read, match):
		if match.astart == 0 and match.rstart > 0:
			return self.remove_back(read, match)
		else:
			return self.remove_front(read, match)

	def remove_front(self, read, match):
		self.lengths_front[match.rstop] += 1
		self._write_rest(read.sequence[:match.rstart], read)
		return read[match.rstop:]

	def remove_back(self, read, match):
		# The adapter is at the end of the read or within the read
		self._write_rest(read.sequence[match.rstop:], read)
		self.lengths_back[len(read) - match.rstart] += 1
		return read[:match.rstart]

	def _write_rest(self, rest, read):
		if len(rest) > 0 and self.rest_file:
			# The adapter is within the read
			print(rest, read.name, file=self.rest_file)

	def __len__(self):
		return len(self.sequence)


class ColorspaceAdapter(Adapter):
	def __init__(self, *args):
		super(ColorspaceAdapter, self).__init__(*args)
		if set(self.sequence) <= set('ACGT'):
			# adapter was given in basespace
			#self.nucleotide_sequence = self.sequence
			self.sequence = colorspace_encode(self.sequence)[1:]

	def remove_front(self, read, match):
		read = Adapter.remove_front(self, read, match)
		read = read[1:]
		# trim one more color
		#rstop = min(rstop + 1, len(read))
		return read

	def remove_back(self, read, match):
		read = Adapter.remove_back(self, read, match)
		read = read[:-1]
		# TODO avoid the copy (previously, this was just an index operation)
		# trim one more color if long enough
		#rstart = max(0, match.rstart - 1)
		return read

	def __repr__(self):
		return '<ColorspaceAdapter(sequence="{0}", where={1})>'.format(self.sequence, self.where)



def matched_wildcards(match, wildcard_char='N'):
	"""
	TODO doc TODO could be an AdapterMatch method
	"""
	wildcards = [ match.read.sequence[match.rstart + i] for i in range(match.astop - match.astart)
		if match.adapter.sequence[match.astart + i] == wildcard_char ]
	return ''.join(wildcards)


def write_read(read, outfile, twoheaders=False):
	"""
	Write read in either FASTA or FASTQ format
	(depending on whether qualities is None or not) to outfile

	If twoheaders is True and the output is FASTQ, then the sequence name
	(description) is also written after the "+" character in the third line.
	"""
	if read.qualities is None:
		# FASTA
		print('>%s\n%s' % (read.name, read.sequence), file=outfile)
	else:
		# FASTQ
		tmp = read.name if twoheaders else ''
		print('@%s\n%s\n+%s\n%s' % (read.name, read.sequence, tmp, read.qualities), file=outfile)


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
		# read from .(CS)FASTA/.QUAL
		return seqio.FastaQualReader(seqfilename, qualityfilename, colorspace)
	else:
		# read from FASTA or FASTQ
		return seqio.SequenceReader(seqfilename, colorspace, fileformat)


class ReadFilter(object):
	"""Filter reads according to length and according to whether any adapter matches."""

	def __init__(self, minimum_length, maximum_length, too_short_outfile, discard_trimmed, statistics):
		self.minimum_length = minimum_length
		self.maximum_length = maximum_length
		self.too_short_outfile = too_short_outfile
		self.statistics = statistics
		self.discard_trimmed = discard_trimmed

	def keep(self, read, trimmed):
		"""
		Return whether to keep the given read.
		"""
		if self.discard_trimmed and trimmed:
			return False
		if len(read.sequence) < self.minimum_length:
			self.statistics.too_short += 1
			if self.too_short_outfile is not None:
				write_read(read, self.too_short_outfile)
			return False
		if len(read.sequence) > self.maximum_length:
			self.statistics.too_long += 1
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
	Remove any suffix from read names.
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


class AdapterCutter(object):
	"""Cut adapters from reads."""

	def __init__(self, adapters, times, rest_file, wildcard_file, info_file):
		"""Initialize this cutter.
		adapters -- list of Adapter objects
		"""
		self.adapters = adapters
		self.times = times
		self.rest_file = rest_file
		self.info_file = info_file
		self.wildcard_file = wildcard_file
		self.stats = Statistics(self.adapters)

	def _best_match(self, read):
		"""
		Find the best matching adapter.

		read -- The read to which each adapter will be aligned

		Return an AdapterMach object or None if there are no matches.
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

	def _write_info(self, read, match):
		"""write one line to the info file"""
		if not self.info_file:
			return
		seq = read.sequence
		if match is None:
			print(read.name, -1, seq, sep='\t', file=self.info_file)
		else:
			print(read.name, match.errors, seq[0:match.rstart], seq[match.rstart:match.rstop], seq[match.rstop:], sep='\t', file=self.info_file)

	def cut(self, read):
		"""
		Cut adapters from a single read. The read will be converted to uppercase
		before it is compared to the adapter sequences.

		read -- a seqio.Sequence object

		Return a tuple (read, trimmed) with the modified read.
		trimmed is True when any adapter was found and trimmed.
		"""
		self.stats.n += 1
		self.stats.total_bp += len(read.sequence)

		if __debug__:
			old_length = len(read.sequence)

		# try at most self.times times to remove an adapter
		any_adapter_matches = False
		for t in xrange(self.times):
			match = self._best_match(read)
			self._write_info(read, match)
			if match is None:
				# nothing found
				break

			if __debug__:
				length = match.astop - match.astart
				assert length > 0
				assert match.errors / length <= match.adapter.max_error_rate
				assert length - match.errors > 0

			any_adapter_matches = True
			if self.wildcard_file:
				print(matched_wildcards(match), read.name, file=self.wildcard_file)

			read = match.adapter.remove(read, match)

		# if an adapter was found, then the read should now be shorter
		assert (not any_adapter_matches) or (len(read.sequence) < old_length)
		if any_adapter_matches: # trimmed: # TODO move to filter class
			self.stats.reads_changed += 1

		return (read, any_adapter_matches)


#class QualityTrimmer:
	#"""Pipeline step that trims qualities"""
	#def __init__(self, cutoff):
		#self.cutoff = cutoff

	#def process(self, read):
		#index = quality_trim_index(qualities, options.quality_cutoff)
		##total_quality_trimmed += len(qualities) - index
		#qualities = qualities[:index]
		#seq = seq[:index]



def main(cmdlineargs=None):
	"""Main function that evaluates command-line parameters and contains the main loop over all reads."""
	parser = HelpfulOptionParser(usage=__doc__, version=__version__)

	parser.add_option("-f", "--format", default=None,
		help="Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. "
		"Ignored when reading csfasta/qual files (default: auto-detect from file name extension).")

	group = OptionGroup(parser, "Options that influence how the adapters are found",
		description="Each of the following three parameters (-a, -b, -g) can be used " +\
			"multiple times and in any combination to search for an entire set of " + \
			"adapters of possibly different types. All of the "+\
			"given adapters will be searched for in each read, but only the best "+\
			"matching one will be trimmed (but see the --times option).")
	group.add_option("-a", "--adapter", action="append", metavar="ADAPTER", dest="adapters", default=[],
		help="Sequence of an adapter that was ligated to the 3' end. The adapter itself and anything that follows is trimmed.")
	group.add_option("-b", "--anywhere", action="append", metavar="ADAPTER", default=[],
		help="Sequence of an adapter that was ligated to the 5' or 3' end. If the adapter is found within the read or overlapping the 3' end of the read, the behavior is the same as for the -a option. If the adapter overlaps the 5' end (beginning of the read), the initial portion of the read matching the adapter is trimmed, but anything that follows is kept.")
	group.add_option("-g", "--front", action="append", metavar="ADAPTER", default=[],
		help="Sequence of an adapter that was ligated to the 5' end. If the " + \
		"adapter sequence starts with the character '^', the adapter is " + \
		"'anchored'. An anchored adapter must appear in its entirety at the " + \
		"5' end of the read (it is a prefix of the read). A non-anchored adapter may " + \
		"appear partially at the 5' end, or it may occur within the read. If it is " + \
		"found within a read, the sequence preceding the adapter is also trimmed. " + \
		"In all cases the adapter itself is trimmed.")
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
		     "positions to FILE.  Use - for standard output.")
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
			"For example, use --length-tag 'length=' to search for fields "
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
	trimmed_outfile = sys.stdout # reads with adapters go here
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
	if options.colorspace and options.front and not options.trim_primer:
		parser.error("Currently, when you want to trim a 5' adapter in colorspace, you must also specify the --trim-primer option")
	if options.anywhere and options.colorspace:
		parser.error("Using --anywhere with color space reads is currently not supported  (if you think this may be useful, contact the author).")
	if not (0 <= options.error_rate <= 1.):
		parser.error("The maximum error rate must be between 0 and 1.")
	if options.overlap < 1:
		parser.error("The overlap must be at least 1.")

	if options.rest_file is not None:
		options.rest_file = xopen(options.rest_file, 'w')
	if options.info_file is not None:
		options.info_file = xopen(options.info_file, 'w')
	if options.wildcard_file is not None:
		options.wildcard_file = xopen(options.wildcard_file, 'w')

	adapters = []

	ADAPTER_CLASS = ColorspaceAdapter if options.colorspace else Adapter
	def append_adapters(adapter_list, where):
		for seq in adapter_list:
			seq = seq.strip()
			w = where
			if w == FRONT and seq.startswith('^'):
				seq = seq[1:]
				w = PREFIX
			adapter = ADAPTER_CLASS(seq, w, options.error_rate,
				options.overlap, options.match_read_wildcards,
				options.match_adapter_wildcards,
				options.rest_file)
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

	#total_bases = 0
	#total_quality_trimmed = 0

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

	cutter = AdapterCutter(adapters, options.times, options.rest_file,
				options.wildcard_file, options.info_file)
	readfilter = ReadFilter(options.minimum_length, options.maximum_length,
		too_short_outfile, options.discard_trimmed, cutter.stats) # TODO stats?
	try:
		twoheaders = None
		reader = read_sequences(input_filename, quality_filename, colorspace=options.colorspace, fileformat=options.format)
		for read in reader:
			# In colorspace, the first character is the last nucleotide of the primer base
			# and the second character encodes the transition from the primer base to the
			# first real base of the read.
			if options.trim_primer:
				read.sequence = read.sequence[2:]
				if read.qualities is not None: # TODO
					read.qualities = read.qualities[1:]
				initial = ''
			elif options.colorspace:
				initial = read.sequence[0]
				read.sequence = read.sequence[1:]
			else:
				initial = ''

			#total_bases += len(qualities)
			if options.quality_cutoff > 0:
				index = quality_trim_index(read.qualities, options.quality_cutoff, options.quality_base)
				read = read[:index]

			read, trimmed = cutter.cut(read)
			for modifier in modifiers:
				read = modifier.apply(read)
			if twoheaders is None:
				try:
					twoheaders = reader.twoheaders
				except AttributeError:
					twoheaders = False
			if not readfilter.keep(read, trimmed):
				continue
			read.sequence = initial + read.sequence
			try:
				write_read(read, trimmed_outfile if trimmed else untrimmed_outfile, twoheaders)
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
	cutter.stats.print_statistics(options.error_rate, file=stat_file)

	return 0


if __name__ == '__main__':
	if len(sys.argv) > 1 and sys.argv[1] == '--profile':
		del sys.argv[1]
		import cProfile as profile
		profile.run('main()', 'cutadapt.prof')
	else:
		sys.exit(main())
