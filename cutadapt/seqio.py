# coding: utf-8
__author__ = "Marcel Martin"

import sys
if sys.version_info[0] < 3:
	from itertools import izip as zip
else:
	basestring = str
from .xopen import xopen
from os.path import splitext
import sys


def _shorten(s, n=20):
	"""Shorten string s to at most n characters, appending "..." if necessary."""
	if s is None:
		return None
	if len(s) > n:
		s = s[:n-3] + '...'
	return s


class Sequence:
	"""qualities is a string and it contains the qualities encoded as ascii(qual+33)."""

	def __init__(self, name, sequence, qualities):
		"""Set qualities to None if there are no quality values"""
		self.name = name
		self.sequence = sequence
		self.qualities = qualities

	def __getitem__(self, key):
		"""slicing"""
		return self.__class__(self.name, self.sequence[key], self.qualities[key] if self.qualities is not None else None)

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = '\', qualities="{0}"'.format(_shorten(self.qualities))
		return 'Sequence(name="{0}", sequence="{1}"{2})'.format(_shorten(self.name), _shorten(self.sequence), qstr)

	def __len__(self):
		return len(self.sequence)

	def __eq__(self, other):
		return self.name == other.name and \
			self.sequence == other.sequence and \
			self.qualities == other.qualities

	def __ne__(self, other):
		return not self.__eq__(other)


class FormatError(Exception):
	"""
	Raised when an input file (FASTA or FASTQ) is malformatted.
	"""
	pass


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


class UnknownFileType(Exception):
	"""
	Raised when SequenceReader could not autodetect the file type.
	"""
	pass


def SequenceReader(file, colorspace=False, fileformat=None):
	"""
	Reader for FASTA and FASTQ files that autodetects the file format.
	Returns either an instance of FastaReader or of FastqReader,
	depending on file type.

	The autodetection can be skipped by setting fileformat to the string
	'fasta' or 'fastq'

	file is a filename or a file-like object.
	If file is a filename, then .gz files are supported.
	If the file name is available, the file type is detected
	by looking at the file name.
	If the file name is not available (for example, reading
	from standard input), then the file is read and the file
	type determined from the content.
	"""
	if fileformat is not None:
		fileformat = fileformat.lower()
		if fileformat == 'fasta':
			return FastaReader(file)
		elif fileformat == 'fastq':
			return FastqReader(file)
		else:
			raise UnknownFileType("File format {0} is unknown (expected 'fasta' or 'fastq').".format(fileformat))

	name = None
	if file == "-":
		file = sys.stdin
	elif isinstance(file, basestring):
		name = file
	elif hasattr(file, "name"):
		name = file.name
	if name is not None:
		if name.endswith('.gz'):
			name = name[:-3]
		name, ext = splitext(name)
		ext = ext.lower()
		if ext in ['.fasta', '.fa', '.fna', '.csfasta', '.csfa']:
			return FastaReader(file)
		elif ext in ['.fastq', '.fq'] or (ext == '.txt' and name.endswith('_sequence')):
			return FastqReader(file, colorspace)
		else:
			raise UnknownFileType("Could not determine whether this is FASTA or FASTQ: file name extension {0} not recognized".format(ext))

	# No name available.
	# Assume that 'file' is an open file
	# and autodetect its type by reading from it.
	for line in file:
		if line.startswith('#'):
			# Skip comment lines (needed for csfasta)
			continue
		if line.startswith('>'):
			return FastaReader(FileWithPrependedLine(file, line))
		if line.startswith('@'):
			return FastqReader(FileWithPrependedLine(file, line), colorspace)
	raise UnknownFileType("File is neither FASTQ nor FASTA.")



class FastaReader(object):
	"""
	Reader for FASTA files.
	"""
	def __init__(self, file, wholefile=False, keep_linebreaks=False):
		"""
		file is a filename or a file-like object.
		If file is a filename, then .gz files are supported.
		If wholefile is True, then it is ok to read the entire file
		into memory. This is faster when there are many newlines in
		the file, but may obviously need a lot of memory.
		keep_linebreaks -- whether to keep the newline characters in the sequence
		"""
		if isinstance(file, basestring):
			file = xopen(file, "r")
		self.fp = file
		self.wholefile = wholefile
		self.keep_linebreaks = keep_linebreaks
		assert not (wholefile and keep_linebreaks), "not supported"

	def __iter__(self):
		"""
		Return instances of the Sequence class.
		The qualities attribute is always None.
		"""
		return self._wholefile_iter() if self.wholefile else self._streaming_iter()

	def _streaming_iter(self):
		"""
		Read next entry from the file (single entry at a time).

		# TODO this can be quadratic since += is used for the sequence string
		"""
		name = None
		seq = ""
		appendchar = '\n' if self.keep_linebreaks else ''
		for line in self.fp:
			# strip() should also take care of DOS line breaks
			line = line.strip()
			if line and line[0] == ">":
				if name is not None:
					assert self.keep_linebreaks or seq.find('\n') == -1
					yield Sequence(name, seq, None)
				name = line[1:]
				seq = ""
			else:
				seq += line + appendchar
		if name is not None:
			assert self.keep_linebreaks or seq.find('\n') == -1
			yield Sequence(name, seq, None)

	def _wholefile_iter(self):
		"""
		This reads in the entire file at once, but is faster than the above code when there are lots of newlines.
		The idea comes from the TAMO package (http://fraenkel.mit.edu/TAMO/), module TAMO.seq.Fasta (author is
		David Benjamin Gordon).
		"""
		wholefile = self.fp.read()
		assert '\r' not in wholefile, "Sorry, currently don't know how to deal with files that contain \\r linebreaks"
		assert len(wholefile) == 0 or wholefile[0] == '>', "FASTA file must start with '>'"
		parts = wholefile.split('\n>')
		# first part has '>' in front
		parts[0] = parts[0][1:]
		for part in parts:
			lines = part.split('\n', 1)
			yield Sequence(name=lines[0], sequence=lines[1].replace('\n', ''), qualities=None)

	def __enter__(self):
		if self.fp is None:
			raise ValueError("I/O operation on closed FastaReader")
		return self

	def __exit__(self, *args):
		self.fp.close()


class FastqReader(object):
	"""
	Reader for FASTQ files. Does not support multi-line FASTQ files.
	"""
	def __init__(self, file, colorspace=False):
		"""
		file is a filename or a file-like object.
		If file is a filename, then .gz files are supported.

		colorspace -- Usually (when this is False), there must be n characters in the sequence and
		n quality values. When this is True, there must be n+1 characters in the sequence and n quality values.
		"""
		if isinstance(file, basestring):
			file = xopen(file, "r")
		self.fp = file
		self.colorspace = colorspace
		self.twoheaders = False

	def __iter__(self):
		"""
		Return tuples: (name, sequence, qualities).
		qualities is a string and it contains the unmodified, encoded qualities.
		"""
		lengthdiff = 1 if self.colorspace else 0
		for i, line in enumerate(self.fp):
			if i % 4 == 0:
				if not line.startswith('@'):
					raise FormatError("at line {0}, expected a line starting with '+'".format(i+1))
				name = line.strip()[1:]
			elif i % 4 == 1:
				sequence = line.strip()
			elif i % 4 == 2:
				line = line.strip()
				if not line.startswith('+'):
					raise FormatError("at line {0}, expected a line starting with '+'".format(i+1))
				if len(line) > 1:
					self.twoheaders = True
					if not line[1:] == name:
						raise FormatError(
							"At line {0}: Two sequence descriptions are given in "
							"the FASTQ file, but they don't match "
							"('{1}' != '{2}')".format(i+1, name, line.rstrip()[1:]))
			elif i % 4 == 3:
				qualities = line.rstrip("\n\r")
				if len(qualities) + lengthdiff != len(sequence):
					raise ValueError("Length of quality sequence and length of read do not match (%d+%d!=%d)" % (len(qualities), lengthdiff, len(sequence)))
				yield Sequence(name, sequence, qualities)

	def __enter__(self):
		if self.fp is None:
			raise ValueError("I/O operation on closed FastqReader")
		return self

	def __exit__(self, *args):
		self.fp.close()


def _quality_to_ascii(qualities, base=33):
	"""
	Convert a list containing qualities given as integer to a string of
	ASCII-encoded qualities.

	base -- ASCII code of quality zero (sensible values are 33 and 64).

	>>> _quality_to_ascii([17, 4, 29, 18])
	'2%>3'
	"""
	qualities = ''.join(chr(q+base) for q in qualities)
	return qualities


class FastaQualReader(object):
	"""
	Reader for reads that are stored in .(CS)FASTA and .QUAL files.
	"""
	def __init__(self, fastafile, qualfile, colorspace=False):
		"""
		fastafile and qualfile are filenames file-like objects.
		If file is a filename, then .gz files are supported.

		colorspace -- Usually (when this is False), there must be n characters in the sequence and
		n quality values. When this is True, there must be n+1 characters in the sequence and n quality values.
		"""
		self.fastareader = FastaReader(fastafile)
		self.qualreader = FastaReader(qualfile, keep_linebreaks=True)
		self.colorspace = colorspace

	def __iter__(self):
		"""
		Yield Sequence objects.
		"""
		lengthdiff = 1 if self.colorspace else 0
		for fastaread, qualread in zip(self.fastareader, self.qualreader):
			qualities = _quality_to_ascii(map(int, qualread.sequence.split()))
			assert fastaread.name == qualread.name
			if len(qualities) + lengthdiff != len(fastaread.sequence):
				raise ValueError("Length of quality sequence and length of read do not match (%d+%d!=%d)" % (
					len(qualities), lengthdiff, len(fastaread.sequence)))
			yield Sequence(fastaread.name, fastaread.sequence, qualities)

	def __enter__(self):
		if self.fastafile is None:
			raise ValueError("I/O operation on closed FastaQualReader")
		return self

	def __exit__(self, *args):
		self.fastareader.close()
		self.qualreader.close()

