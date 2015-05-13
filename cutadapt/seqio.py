# coding: utf-8
from __future__ import print_function, division, absolute_import
import sys
from os.path import splitext
from cutadapt.xopen import xopen
from cutadapt.compat import zip, basestring

__author__ = "Marcel Martin"


class FormatError(Exception):
	"""
	Raised when an input file (FASTA or FASTQ) is malformatted.
	"""


def _shorten(s, n=100):
	"""Shorten string s to at most n characters, appending "..." if necessary."""
	if s is None:
		return None
	if len(s) > n:
		s = s[:n-3] + '...'
	return s


class Sequence(object):
	"""qualities is a string and it contains the qualities encoded as ascii(qual+33)."""

	def __init__(self, name, sequence, qualities=None, twoheaders=False, match=None):
		"""Set qualities to None if there are no quality values"""
		self.name = name
		self.sequence = sequence
		self.qualities = qualities
		self.twoheaders = twoheaders
		self.match = match
		if qualities is not None:
			if len(qualities) != len(sequence):
				rname = _shorten(name)
				raise FormatError("In read named {0!r}: Length of quality sequence ({1}) and length of read ({2}) do not match".format(
					rname, len(qualities), len(sequence)))

	def __getitem__(self, key):
		"""slicing"""
		return self.__class__(
			self.name,
			self.sequence[key],
			self.qualities[key] if self.qualities is not None else None,
		    self.twoheaders,
		    self.match)

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = ', qualities={0!r}'.format(_shorten(self.qualities))
		return '<Sequence(name={0!r}, sequence={1!r}{2})>'.format(_shorten(self.name), _shorten(self.sequence), qstr)

	def __len__(self):
		return len(self.sequence)

	def __eq__(self, other):
		return self.name == other.name and \
			self.sequence == other.sequence and \
			self.qualities == other.qualities

	def __ne__(self, other):
		return not self.__eq__(other)

	def write(self, outfile):
		if self.qualities is not None:
			s = '@' + self.name + '\n' + self.sequence + '\n+'
			if self.twoheaders:
				s += self.name
			s += '\n' + self.qualities + '\n'
		else:
			s = '>' + self.name + '\n' + self.sequence + '\n'
		outfile.write(s)


try:
	from ._seqio import Sequence
except ImportError:
	pass


class ColorspaceSequence(Sequence):
	def __init__(self, name, sequence, qualities, primer=None, twoheaders=False, match=None):
		# In colorspace, the first character is the last nucleotide of the primer base
		# and the second character encodes the transition from the primer base to the
		# first real base of the read.
		if primer is None:
			self.primer = sequence[0:1]
			sequence = sequence[1:]
		else:
			self.primer = primer
		if qualities is not None and len(sequence) != len(qualities):
			rname = _shorten(name)
			raise FormatError("In read named {0!r}: length of colorspace quality "
				"sequence ({1}) and length of read ({2}) do not match (primer "
				"is: {3!r})".format(rname, len(qualities), len(sequence), self.primer))
		super(ColorspaceSequence, self).__init__(name, sequence, qualities, twoheaders, match)
		if not self.primer in ('A', 'C', 'G', 'T'):
			raise FormatError("Primer base is {0!r} in read {1!r}, but it "
				"should be one of A, C, G, T.".format(
					self.primer, _shorten(name)))

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = ', qualities={0!r}'.format(_shorten(self.qualities))
		return '<ColorspaceSequence(name={0!r}, primer={1!r}, sequence={2!r}{3})>'.format(_shorten(self.name), self.primer, _shorten(self.sequence), qstr)

	def __getitem__(self, key):
		return self.__class__(
			self.name,
			self.sequence[key],
			self.qualities[key] if self.qualities is not None else None,
			self.primer,
			self.twoheaders,
			self.match)

	def write(self, outfile):
		if self.qualities is not None:
			s = '@' + self.name + '\n' + self.primer + self.sequence + '\n+'
			if self.twoheaders:
				s += self.name
			s += '\n' + self.qualities + '\n'
		else:
			s = '>' + self.name + '\n' + self.primer + self.sequence + '\n'
		outfile.write(s)


def sra_colorspace_sequence(name, sequence, qualities, twoheaders):
	"""Factory for an SRA colorspace sequence (which has one quality value too many)"""
	return ColorspaceSequence(name, sequence, qualities[1:], twoheaders=twoheaders)


class FileWithPrependedLine(object):
	"""
	A file-like object that allows to "prepend" a single
	line to an already opened file. That is, further
	reads on the file will return the provided line and
	only then the actual content. This is needed to solve
	the problem of autodetecting input from a stream:
	As soon as the first line has been read, we know
	the file type, but also that line is "gone" and
	unavailable for further processing.
	"""
	def __init__(self, file, line):
		"""
		file is an already opened file-like object.
		line is a single string (newline will be appended if not included)
		"""
		if not line.endswith('\n'):
			line += '\n'
		self.first_line = line
		self.file = file

	def __iter__(self):
		yield self.first_line
		for line in self.file:
			yield line

	def close(self):
		self.file.close()


class FastaReader(object):
	"""
	Reader for FASTA files.
	"""
	def __init__(self, file, keep_linebreaks=False, sequence_class=Sequence):
		"""
		file is a filename or a file-like object.
		If file is a filename, then it is passed to xopen().

		keep_linebreaks -- whether to keep newline characters in the sequence
		"""
		if isinstance(file, basestring):
			file = xopen(file)
			self._file_passed = False
		else:
			self._file_passed = True
		self.fp = file
		self.sequence_class = sequence_class
		self.delivers_qualities = False
		self._delimiter = '\n' if keep_linebreaks else ''

	def __iter__(self):
		"""
		Read next entry from the file (single entry at a time).
		"""
		name = None
		seq = []
		for i, line in enumerate(self.fp):
			# strip() also removes DOS line breaks
			line = line.strip()
			if not line:
				continue
			if line and line[0] == '>':
				if name is not None:
					yield self.sequence_class(name, self._delimiter.join(seq), None)
				name = line[1:]
				seq = []
			elif line and line[0] == '#':
				continue
			elif name is not None:
				seq.append(line)
			else:
				raise FormatError("At line {0}: Expected '>' at beginning of "
					"FASTA record, but got {1!r}.".format(i+1, _shorten(line)))

		if name is not None:
			yield self.sequence_class(name, self._delimiter.join(seq), None)

	def close(self):
		if not self._file_passed and self.fp is not None:
			self.fp.close()
			self.fp = None

	def __enter__(self):
		if self.fp is None:
			raise ValueError("I/O operation on closed FastaReader")
		return self

	def __exit__(self, *args):
		self.close()


class ColorspaceFastaReader(FastaReader):
	def __init__(self, file, keep_linebreaks=False):
		super(ColorspaceFastaReader, self).__init__(file, keep_linebreaks, sequence_class=ColorspaceSequence)


class FastqReader(object):
	"""
	Reader for FASTQ files. Does not support multi-line FASTQ files.
	"""
	def __init__(self, file, sequence_class=Sequence):
		"""
		file is a filename or a file-like object.
		If file is a filename, then .gz files are supported.

		The sequence_class should be a class such as Sequence or
		ColorspaceSequence.
		"""
		if isinstance(file, basestring):
			file = xopen(file)
			self._file_passed = False
		else:
			self._file_passed = True
		self.fp = file
		self.sequence_class = sequence_class
		self.delivers_qualities = True

	def __iter__(self):
		"""
		Return tuples: (name, sequence, qualities).
		qualities is a string and it contains the unmodified, encoded qualities.
		"""
		for i, line in enumerate(self.fp):
			if i % 4 == 0:
				if not line.startswith('@'):
					raise FormatError("At line {0}: Expected a line starting with '@'".format(i+1))
				name = line.strip()[1:]
			elif i % 4 == 1:
				sequence = line.strip()
			elif i % 4 == 2:
				line = line.strip()
				if not line.startswith('+'):
					raise FormatError("At line {0}: Expected a line starting with '+'".format(i+1))
				if len(line) > 1:
					twoheaders = True
					if not line[1:] == name:
						raise FormatError(
							"At line {0}: Sequence descriptions in the FASTQ file do not match "
							"({1!r} != {2!r}).\n"
							"The second sequence description must be either empty "
							"or equal to the first description.".format(
								i+1, name, line.rstrip()[1:]))
				else:
					twoheaders = False
			elif i % 4 == 3:
				qualities = line.rstrip('\n\r')
				yield self.sequence_class(name, sequence, qualities, twoheaders=twoheaders)

	def close(self):
		if not self._file_passed and self.fp is not None:
			self.fp.close()
			self.fp = None

	def __enter__(self):
		if self.fp is None:
			raise ValueError("I/O operation on closed FastqReader")
		return self

	def __exit__(self, *args):
		self.close()


try:
	from ._seqio import FastqReader, FormatError
except ImportError:
	pass


class ColorspaceFastqReader(FastqReader):
	def __init__(self, file):
		super(ColorspaceFastqReader, self).__init__(file, sequence_class=ColorspaceSequence)


class SRAColorspaceFastqReader(FastqReader):
	def __init__(self, file):
		super(SRAColorspaceFastqReader, self).__init__(file, sequence_class=sra_colorspace_sequence)


class FastaQualReader(object):
	"""
	Reader for reads that are stored in .(CS)FASTA and .QUAL files.
	"""
	def __init__(self, fastafile, qualfile, sequence_class=Sequence):
		"""
		fastafile and qualfile are filenames or file-like objects.
		If a filename is used, then .gz files are recognized.

		The objects returned when iteritng over this file are instances of the
		given sequence_class.
		"""
		self.fastareader = FastaReader(fastafile)
		self.qualreader = FastaReader(qualfile, keep_linebreaks=True)
		self.sequence_class = sequence_class
		self.delivers_qualities = True

	def __iter__(self):
		"""
		Yield Sequence objects.
		"""
		# conversion dictionary: maps strings to the appropriate ASCII-encoded character
		conv = dict()
		for i in range(-5, 256 - 33):
			conv[str(i)] = chr(i + 33)
		for fastaread, qualread in zip(self.fastareader, self.qualreader):
			if fastaread.name != qualread.name:
				raise FormatError("The read names in the FASTA and QUAL file do not match ({0!r} != {1!r})".format(fastaread.name, qualread.name))
			try:
				qualities = ''.join([conv[value] for value in qualread.sequence.split()])
			except KeyError as e:
				raise FormatError("Within read named {0!r}: Found invalid quality value {1}".format(fastaread.name, e))
			assert fastaread.name == qualread.name
			yield self.sequence_class(fastaread.name, fastaread.sequence, qualities)

	def close(self):
		self.fastareader.close()
		self.qualreader.close()

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()


class ColorspaceFastaQualReader(FastaQualReader):
	def __init__(self, fastafile, qualfile):
		super(ColorspaceFastaQualReader, self).__init__(fastafile, qualfile, sequence_class=ColorspaceSequence)


class PairedSequenceReader(object):
	"""
	Wrap two SequenceReader instances, making sure that reads are
	properly paired.
	"""
	def __init__(self, file1, file2, colorspace=False, fileformat=None):
		self.reader1 = open(file1, colorspace=colorspace, fileformat=fileformat)
		self.reader2 = open(file2, colorspace=colorspace, fileformat=fileformat)
		self.delivers_qualities = self.reader1.delivers_qualities

	def __iter__(self):
		# Avoid usage of zip() below since it will consume one item too many.
		it1, it2 = iter(self.reader1), iter(self.reader2)
		while True:
			try:
				r1 = next(it1)
			except StopIteration:
				# End of file 1. Make sure that file 2 is also at end.
				try:
					next(it2)
					raise FormatError("Reads are improperly paired. There are more reads in file 2 than in file 1.")
				except StopIteration:
					pass
				break
			try:
				r2 = next(it2)
			except StopIteration:
				raise FormatError("Reads are improperly paired. There are more reads in file 1 than in file 2.")

			name1 = r1.name.split(None, 1)[0]
			name2 = r2.name.split(None, 1)[0]
			if name1[-2:-1] == '/':
				name1 = name1[:-2]
			if name2[-2:-1] == '/':
				name2 = name2[:-2]
			if name1 != name2:
				raise FormatError("Reads are improperly paired. Read name '{0}' in file 1 not equal to '{1}' in file 2.".format(name1, name2))
			yield (r1, r2)


class UnknownFileType(Exception):
	"""
	Raised when open could not autodetect the file type.
	"""


def open(file1, file2=None, qualfile=None, colorspace=False, fileformat=None):
	"""
	Open sequence file in FASTA or FASTQ format. Parameters file1, file2 and
	qualfile will be passed to xopen and can therefore be paths to regular or
	compressed files or file-like objects. If only file1 is provided, a
	FastaReader or FastqReader (for single-end reads) is returned. If file2
	is also provided, a PairedSequenceReader is returned. If qualfile is
	given, a FastaQualReader from file1 and qualfile is returned. One of file2
	and qualfile must always be None (no paired-end data is supported when
	reading qualfiles).

	If the colorspace parameter is set to True, the returned readers are
	ColorspaceFastaReader, ColorspaceFastqReader or ColorspaceFastaQualReader
	instead.

	If possible, file format is autodetected by inspecting the file name:
	.fasta/.fa, .fastq/.fq and some other extensions are allowed. If the
	file name is not available (when reading from standard input), the file is
	read and the file type determined from the content. The autodetection can
	be skipped by setting fileformat to one of 'fasta', 'fastq', 'sra-fastq'.
	Colorspace is not auto-detected and must always be requested explicitly.
	"""
	if file2 is not None and qualfile is not None:
		raise ValueError("Setting both file2 and qualfile is not supported")
	if file2 is not None:
		return PairedSequenceReader(file1, file2, colorspace, fileformat)

	if qualfile is not None:
		if colorspace:
			# read from .(CS)FASTA/.QUAL
			return ColorspaceFastaQualReader(file1, qualfile)
		else:
			return FastaQualReader(file1, qualfile)
	# read from FASTA or FASTQ
	fastq_reader = ColorspaceFastqReader if colorspace else FastqReader
	fasta_reader = ColorspaceFastaReader if colorspace else FastaReader

	if fileformat is not None:
		fileformat = fileformat.lower()
		if fileformat == 'fasta':
			return fasta_reader(file1)
		elif fileformat == 'fastq':
			return fastq_reader(file1)
		elif fileformat == 'sra-fastq' and colorspace:
			return SRAColorspaceFastqReader(file1)
		else:
			raise UnknownFileType("File format {0} is unknown (expected "
				"'sra-fastq' (only for colorspace), 'fasta' or 'fastq').".format(fileformat))

	name = None
	if file1 == "-":
		file1 = sys.stdin
	elif isinstance(file1, basestring):
		name = file1
	elif hasattr(file1, "name"):
		# Assume that 'file1' is an open file1
		name = file1.name

	if name is not None:
		if name.endswith('.gz'):
			name = name[:-3]
		elif name.endswith('.xz'):
			name = name[:-3]
		elif name.endswith('.bz2'):
			name = name[:-4]
		name, ext = splitext(name)
		ext = ext.lower()
		if ext in ['.fasta', '.fa', '.fna', '.csfasta', '.csfa']:
			return fasta_reader(file1)
		elif ext in ['.fastq', '.fq'] or (ext == '.txt' and name.endswith('_sequence')):
			return fastq_reader(file1)
		else:
			raise UnknownFileType("Could not determine whether this is FASTA "
				"or FASTQ: file name extension {0} not recognized".format(ext))

	# No name available.
	# autodetect type by reading from the file
	for line in file1:
		if line.startswith('#'):
			# Skip comment lines (needed for csfasta)
			continue
		if line.startswith('>'):
			return fasta_reader(FileWithPrependedLine(file1, line))
		if line.startswith('@'):
			return fastq_reader(FileWithPrependedLine(file1, line))
	raise UnknownFileType("File is neither FASTQ nor FASTA.")
