# coding: utf-8
"""
Sequence I/O classes: Reading and writing of FASTA and FASTQ files.

TODO

- Sequence.name should be Sequence.description or so (reserve .name for the part
  before the first space)
"""
from __future__ import print_function, division, absolute_import

import sys
from os.path import splitext
from xopen import xopen

from .compat import zip, basestring

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

	def __init__(self, name, sequence, qualities=None, second_header=False, match=None):
		"""Set qualities to None if there are no quality values"""
		self.name = name
		self.sequence = sequence
		self.qualities = qualities
		self.second_header = second_header
		self.match = match
		self.original_length = len(sequence)
		if qualities is not None:
			if len(qualities) != len(sequence):
				rname = _shorten(name)
				raise FormatError("In read named {0!r}: Length of quality sequence ({1}) and "
					"length of read ({2}) do not match".format(rname, len(qualities), len(sequence)))
	
	def __getitem__(self, key):
		"""slicing"""
		return self.__class__(
			self.name,
			self.sequence[key],
			self.qualities[key] if self.qualities is not None else None,
			self.second_header,
			self.match)

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = ', qualities={0!r}'.format(_shorten(self.qualities))
		return '<Sequence(name={0!r}, sequence={1!r}{2})>'.format(
			_shorten(self.name), _shorten(self.sequence), qstr)

	def __len__(self):
		return len(self.sequence)

	def __eq__(self, other):
		return self.name == other.name and \
			self.sequence == other.sequence and \
			self.qualities == other.qualities

	def __ne__(self, other):
		return not self.__eq__(other)


class SequenceReader(object):
	"""Read possibly compressed files containing sequences"""
	_close_on_exit = False

	def __init__(self, file):
		"""
		file is a path or a file-like object. In both cases, the file may
		be compressed (.gz, .bz2, .xz).
		"""
		if isinstance(file, basestring):
			file = xopen(file)
			self._close_on_exit = True
		self._file = file

	def close(self):
		if self._close_on_exit and self._file is not None:
			self._file.close()
			self._file = None

	def __enter__(self):
		if self._file is None:
			raise ValueError("I/O operation on closed SequenceReader")
		return self

	def __exit__(self, *args):
		self.close()


try:
	from ._seqio import Sequence
except ImportError:
	pass


class ColorspaceSequence(Sequence):
	def __init__(self, name, sequence, qualities, primer=None, second_header=False, match=None):
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
		super(ColorspaceSequence, self).__init__(name, sequence, qualities, second_header, match)
		if not self.primer in ('A', 'C', 'G', 'T'):
			raise FormatError("Primer base is {0!r} in read {1!r}, but it "
				"should be one of A, C, G, T.".format(
					self.primer, _shorten(name)))

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = ', qualities={0!r}'.format(_shorten(self.qualities))
		return '<ColorspaceSequence(name={0!r}, primer={1!r}, sequence={2!r}{3})>'.format(
			_shorten(self.name), self.primer, _shorten(self.sequence), qstr)

	def __getitem__(self, key):
		return self.__class__(
			self.name,
			self.sequence[key],
			self.qualities[key] if self.qualities is not None else None,
			self.primer,
			self.second_header,
			self.match)

	def __reduce__(self):
		return (ColorspaceSequence, (self.name, self.sequence, self.qualities, self.primer,
			self.second_header, self.match))


def sra_colorspace_sequence(name, sequence, qualities, second_header):
	"""Factory for an SRA colorspace sequence (which has one quality value too many)"""
	return ColorspaceSequence(name, sequence, qualities[1:], second_header=second_header)


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
		self._file = file

	def __iter__(self):
		yield self.first_line
		for line in self._file:
			yield line

	def close(self):
		self._file.close()


class FastaReader(SequenceReader):
	"""
	Reader for FASTA files.
	"""
	def __init__(self, file, keep_linebreaks=False, sequence_class=Sequence):
		"""
		file is a path or a file-like object. In both cases, the file may
		be compressed (.gz, .bz2, .xz).

		keep_linebreaks -- whether to keep newline characters in the sequence
		"""
		super(FastaReader, self).__init__(file)
		self.sequence_class = sequence_class
		self.delivers_qualities = False
		self._delimiter = '\n' if keep_linebreaks else ''

	def __iter__(self):
		"""
		Read next entry from the file (single entry at a time).
		"""
		name = None
		seq = []
		for i, line in enumerate(self._file):
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


class ColorspaceFastaReader(FastaReader):
	def __init__(self, file, keep_linebreaks=False):
		super(ColorspaceFastaReader, self).__init__(file, keep_linebreaks, sequence_class=ColorspaceSequence)


class FastqReader(SequenceReader):
	"""
	Reader for FASTQ files. Does not support multi-line FASTQ files.
	"""
	def __init__(self, file, sequence_class=Sequence):  # TODO could be a class attribute
		"""
		file is a path or a file-like object. compressed files are supported.

		The sequence_class should be a class such as Sequence or
		ColorspaceSequence.
		"""
		super(FastqReader, self).__init__(file)
		self.sequence_class = sequence_class
		self.delivers_qualities = True

	def __iter__(self):
		"""
		Return tuples: (name, sequence, qualities).
		qualities is a string and it contains the unmodified, encoded qualities.
		"""
		i = 3
		for i, line in enumerate(self._file):
			if i % 4 == 0:
				if not line.startswith('@'):
					raise FormatError("Line {0} in FASTQ file is expected to start with '@', "
						"but found {1!r}".format(i+1, line[:10]))
				name = line.strip()[1:]
			elif i % 4 == 1:
				sequence = line.strip()
			elif i % 4 == 2:
				line = line.strip()
				if not line.startswith('+'):
					raise FormatError("Line {0} in FASTQ file is expected to start with '+', "
						"but found {1!r}".format(i+1, line[:10]))
				if len(line) > 1:
					if line[1:] != name:
						raise FormatError(
							"At line {0}: Sequence descriptions in the FASTQ file do not match "
							"({1!r} != {2!r}).\n"
							"The second sequence description must be either empty "
							"or equal to the first description.".format(
								i+1, name, line[1:].rstrip()))
					second_header = True
				else:
					second_header = False
			elif i % 4 == 3:
				qualities = line.rstrip('\n\r')
				yield self.sequence_class(name, sequence, qualities, second_header=second_header)
		if i % 4 != 3:
			raise FormatError("FASTQ file ended prematurely")


try:
	from ._seqio import FastqReader
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
	delivers_qualities = True

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
				raise FormatError("The read names in the FASTA and QUAL file "
					"do not match ({0!r} != {1!r})".format(fastaread.name, qualread.name))
			try:
				qualities = ''.join([conv[value] for value in qualread.sequence.split()])
			except KeyError as e:
				raise FormatError("Within read named {0!r}: Found invalid quality "
					"value {1}".format(fastaread.name, e))
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


def sequence_names_match(r1, r2):
	"""
	Check whether the sequences r1 and r2 have identical names, ignoring a
	suffix of '1' or '2'. Some old paired-end reads have names that end in '/1'
	and '/2'. Also, the fastq-dump tool (used for converting SRA files to FASTQ)
	appends a .1 and .2 to paired-end reads if option -I is used.
	"""
	name1 = r1.name.split(None, 1)[0]
	name2 = r2.name.split(None, 1)[0]
	if name1[-1:] in '12' and name2[-1:] in '12':
		name1 = name1[:-1]
		name2 = name2[:-1]
	return name1 == name2


class PairedSequenceReader(object):
	"""
	Read paired-end reads from two files.

	Wraps two SequenceReader instances, making sure that reads are properly
	paired.
	"""
	def __init__(self, file1, file2, colorspace=False, fileformat=None):
		self.reader1 = open(file1, colorspace=colorspace, fileformat=fileformat)
		self.reader2 = open(file2, colorspace=colorspace, fileformat=fileformat)
		self.delivers_qualities = self.reader1.delivers_qualities

	def __iter__(self):
		"""
		Iterate over the paired reads. Each item is a pair of Sequence objects.
		"""
		# Avoid usage of zip() below since it will consume one item too many.
		it1, it2 = iter(self.reader1), iter(self.reader2)
		while True:
			try:
				r1 = next(it1)
			except StopIteration:
				# End of file 1. Make sure that file 2 is also at end.
				try:
					next(it2)
					raise FormatError("Reads are improperly paired. There are more reads in "
						"file 2 than in file 1.")
				except StopIteration:
					pass
				break
			try:
				r2 = next(it2)
			except StopIteration:
				raise FormatError("Reads are improperly paired. There are more reads in "
					"file 1 than in file 2.")
			if not sequence_names_match(r1, r2):
				raise FormatError("Reads are improperly paired. Read name '{0}' "
					"in file 1 does not match '{1}' in file 2.".format(r1.name, r2.name))
			yield (r1, r2)

	def close(self):
		self.reader1.close()
		self.reader2.close()

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()


class InterleavedSequenceReader(object):
	"""
	Read paired-end reads from an interleaved FASTQ file.
	"""
	def __init__(self, file, colorspace=False, fileformat=None):
		self.reader = open(file, colorspace=colorspace, fileformat=fileformat)
		self.delivers_qualities = self.reader.delivers_qualities

	def __iter__(self):
		# Avoid usage of zip() below since it will consume one item too many.
		it = iter(self.reader)
		for r1 in it:
			try:
				r2 = next(it)
			except StopIteration:
				raise FormatError("Interleaved input file incomplete: Last record has no partner.")
			if not sequence_names_match(r1, r2):
				raise FormatError("Reads are improperly paired. Name {0!r} "
					"(first) does not match {1!r} (second).".format(r1.name, r2.name))
			yield (r1, r2)

	def close(self):
		self.reader.close()

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()


class FileWriter(object):
	def __init__(self, file):
		if isinstance(file, str):
			self._file = xopen(file, 'w')
			self._close_on_exit = True
		else:
			self._file = file
			self._close_on_exit = False
	
	def close(self):
		if self._close_on_exit:
			self._file.close()
	
	def __enter__(self):
		if self._file.closed:
			raise ValueError("I/O operation on closed file")
		return self

	def __exit__(self, *args):
		self.close()


class SingleRecordWriter(object):
	"""Public interface to single-record files"""
	def write(self, record):
		raise NotImplementedError()


class FastaWriter(FileWriter, SingleRecordWriter):
	"""
	Write FASTA-formatted sequences to a file.
	"""

	def __init__(self, file, line_length=None):
		"""
		If line_length is not None, the lines will
		be wrapped after line_length characters.
		"""
		FileWriter.__init__(self, file)
		self.line_length = line_length if line_length != 0 else None
	
	def write(self, name_or_seq, sequence=None):
		"""Write an entry to the the FASTA file.

		If only one parameter (name_or_seq) is given, it must have
		attributes .name and .sequence, which are then used.
		Otherwise, the first parameter must be the name and the second
		the sequence.

		The effect is that you can write this:
		writer.write("name", "ACCAT")
		or
		writer.write(Sequence("name", "ACCAT"))
		"""
		if sequence is None:
			name = name_or_seq.name
			sequence = name_or_seq.sequence
		else:
			name = name_or_seq
		
		if self.line_length is not None:
			print('>{0}'.format(name), file=self._file)
			for i in range(0, len(sequence), self.line_length):
				print(sequence[i:i+self.line_length], file=self._file)
			if len(sequence) == 0:
				print(file=self._file)
		else:
			print('>{0}'.format(name), sequence, file=self._file, sep='\n')


class ColorspaceFastaWriter(FastaWriter):
	def write(self, record):
		name = record.name
		sequence = record.primer + record.sequence
		super(ColorspaceFastaWriter, self).write(name, sequence)


class FastqWriter(FileWriter, SingleRecordWriter):
	"""
	Write sequences with qualities in FASTQ format.

	FASTQ files are formatted like this:
	@read name
	SEQUENCE
	+
	QUALITIS
	"""
	def write(self, record):
		"""
		Write a Sequence record to the the FASTQ file.

		The record must have attributes .name, .sequence and .qualities.
		"""
		name2 = record.name if record.second_header else ''
		s = ('@' + record.name + '\n' + record.sequence + '\n+' +
				name2 + '\n' + record.qualities + '\n')
		self._file.write(s)

	def writeseq(self, name, sequence, qualities):
		print("@{0:s}\n{1:s}\n+\n{2:s}".format(
			name, sequence, qualities), file=self._file)


class ColorspaceFastqWriter(FastqWriter):
	def write(self, record):
		name = record.name
		sequence = record.primer + record.sequence
		qualities = record.qualities
		super(ColorspaceFastqWriter, self).writeseq(name, sequence, qualities)


class PairRecordWriter(object):
	"""Public interface to paired-record files"""
	def write(self, read1, read2):
		raise NotImplementedError()

	def close(self):
		raise NotImplementedError()
	
	def __enter__(self):
		# TODO do not allow this twice
		return self

	def __exit__(self, *args):
		self.close()


class PairedSequenceWriter(PairRecordWriter):
	def __init__(self, file1, file2, colorspace=False, fileformat='fastq', qualities=None):
		self._writer1 = open(file1, colorspace=colorspace, fileformat=fileformat, mode='w',
			qualities=qualities)
		self._writer2 = open(file2, colorspace=colorspace, fileformat=fileformat, mode='w',
			qualities=qualities)

	def write(self, read1, read2):
		self._writer1.write(read1)
		self._writer2.write(read2)

	def close(self):
		self._writer1.close()
		self._writer2.close()


class InterleavedSequenceWriter(PairRecordWriter):
	"""
	Write paired-end reads to an interleaved FASTA or FASTQ file
	"""
	def __init__(self, file, colorspace=False, fileformat='fastq', qualities=None):
		self._writer = open(
			file, colorspace=colorspace, fileformat=fileformat, mode='w', qualities=qualities)

	def write(self, read1, read2):
		self._writer.write(read1)
		self._writer.write(read2)

	def close(self):
		self._writer.close()


class UnknownFileType(Exception):
	"""
	Raised when open could not autodetect the file type.
	"""


def open(file1, file2=None, qualfile=None, colorspace=False, fileformat=None,
	interleaved=False, mode='r', qualities=None):
	"""
	Open sequence files in FASTA or FASTQ format for reading or writing. This is
	a factory that returns an instance of one of the ...Reader or ...Writer
	classes also defined in this module.

	file1, file2, qualfile -- Paths to regular or compressed files or file-like
		objects. Use file1 if data is single-end. If also file2 is provided,
		sequences are paired. If qualfile is given, then file1 must be a FASTA
		file and sequences are single-end. One of file2 and qualfile must always
		be None (no paired-end data is supported when reading qualfiles).

	mode -- Either 'r' for reading or 'w' for writing.

	interleaved -- If True, then file1 contains interleaved paired-end data.
		file2 and qualfile must be None in this case.

	colorspace -- If True, instances of the Colorspace... classes
		are returned.

	fileformat -- If set to None, file format is autodetected from the file name
		extension. Set to 'fasta', 'fastq', or 'sra-fastq' to not auto-detect.
		Colorspace is not auto-detected and must always be requested explicitly.

	qualities -- When mode is 'w' and fileformat is None, this can be set to
		True or False to specify whether the written sequences will have quality
		values. This is is used in two ways:
		* If the output format cannot be determined (unrecognized extension
		  etc), no exception is raised, but fasta or fastq format is chosen
		  appropriately.
		* When False (no qualities available), an exception is raised when the
		  auto-detected output format is FASTQ.
	"""
	if mode not in ('r', 'w'):
		raise ValueError("Mode must be 'r' or 'w'")
	if interleaved and (file2 is not None or qualfile is not None):
		raise ValueError("When interleaved is set, file2 and qualfile must be None")
	if file2 is not None and qualfile is not None:
		raise ValueError("Setting both file2 and qualfile is not supported")
	if file2 is not None:
		if mode == 'r':
			return PairedSequenceReader(file1, file2, colorspace, fileformat)
		else:
			return PairedSequenceWriter(file1, file2, colorspace, fileformat, qualities)

	if interleaved:
		if mode == 'r':
			return InterleavedSequenceReader(file1, colorspace, fileformat)
		else:
			return InterleavedSequenceWriter(file1, colorspace, fileformat, qualities)

	if qualfile is not None:
		if mode == 'w':
			raise NotImplementedError('Writing to csfasta/qual not supported')
		if colorspace:
			# read from .(CS)FASTA/.QUAL
			return ColorspaceFastaQualReader(file1, qualfile)
		else:
			return FastaQualReader(file1, qualfile)

	# All the multi-file things have been dealt with, delegate rest to the
	# single-file function.
	return _seqopen1(file1, colorspace=colorspace, fileformat=fileformat,
		mode=mode, qualities=qualities)


def _seqopen1(file, colorspace=False, fileformat=None, mode='r', qualities=None):
	"""
	Open a single sequence file. See description above.
	"""
	if mode == 'r':
		fastq_handler = ColorspaceFastqReader if colorspace else FastqReader
		fasta_handler = ColorspaceFastaReader if colorspace else FastaReader
	elif mode == 'w':
		fastq_handler = ColorspaceFastqWriter if colorspace else FastqWriter
		fasta_handler = ColorspaceFastaWriter if colorspace else FastaWriter
	else:
		raise ValueError("Mode must be 'r' or 'w'")

	if fileformat:  # Explict file format given
		fileformat = fileformat.lower()
		if fileformat == 'fasta':
			return fasta_handler(file)
		elif fileformat == 'fastq':
			return fastq_handler(file)
		elif fileformat == 'sra-fastq' and colorspace:
			if mode == 'w':
				raise NotImplementedError('Writing to sra-fastq not supported')
			return SRAColorspaceFastqReader(file)
		else:
			raise UnknownFileType("File format {0!r} is unknown (expected "
				"'sra-fastq' (only for colorspace), 'fasta' or 'fastq').".format(fileformat))

	# Detect file format
	name = None
	if file == "-":
		file = sys.stdin if mode == 'r' else sys.stdout
	elif isinstance(file, basestring):
		name = file
	elif hasattr(file, "name"):  # seems to be an open file-like object
		name = file.name

	if name:
		for ext in ('.gz', '.xz', '.bz2'):
			if name.endswith(ext):
				name = name[:-len(ext)]
				break
		name, ext = splitext(name)
		ext = ext.lower()
		if ext in ['.fasta', '.fa', '.fna', '.csfasta', '.csfa']:
			format = 'fasta'
		elif ext in ['.fastq', '.fq'] or (ext == '.txt' and name.endswith('_sequence')):
			format = 'fastq'
		elif mode == 'w' and qualities is True:
			# Format not recognized, but know we want to write reads with qualities
			format = 'fastq'
		elif mode == 'w' and qualities is False:
			# Same, but we know that we want to write reads without qualities
			format = 'fasta'
		else:
			raise UnknownFileType("Could not determine whether file {0!r} is FASTA "
				"or FASTQ: file name extension {1!r} not recognized".format(file, ext))
		if format == 'fastq' and qualities is False:
			raise ValueError("Output format cannot be FASTQ since no quality "
				"values are available.")
		if format == 'fastq':
			return fastq_handler(file)
		else:
			return fasta_handler(file)

	if mode == 'w':
		if qualities is True:
			return fastq_handler(file)
		elif qualities is False:
			return fasta_handler(file)
		raise UnknownFileType('Cannot determine whether to write in FASTA or '
			'FASTQ format')
	# No name available. Try to autodetect type by reading from the file.
	for line in file:
		if line.startswith('#'):
			# Skip comment lines (needed for csfasta)
			continue
		if line.startswith('>'):
			return fasta_handler(FileWithPrependedLine(file, line))
		if line.startswith('@'):
			return fastq_handler(FileWithPrependedLine(file, line))
	raise UnknownFileType("File is neither FASTQ nor FASTA.")
