# kate: syntax Python;
# cython: profile=False, emit_code_comments=False
from __future__ import print_function, division, absolute_import
from xopen import xopen
from .seqio import _shorten, FormatError, SequenceReader

cimport cython

# It would be nice to be able to have the first parameter be a
# unsigned char[:] (memory view), but this fails with a BufferError
# when a bytes object is passed in.
# See <https://stackoverflow.com/questions/28203670/>

ctypedef fused bytes_or_bytearray:
	bytes
	bytearray


def two_fastq_heads(bytes_or_bytearray buf1, bytes_or_bytearray buf2, Py_ssize_t end1, Py_ssize_t end2):
	"""
	Skip forward in the two buffers by multiples of four lines.

	Return a tuple (length1, length2) such that buf1[:length1] and
	buf2[:length2] contain the same number of lines (where the
	line number is divisible by four).
	"""
	cdef:
		Py_ssize_t pos1 = 0, pos2 = 0
		Py_ssize_t linebreaks = 0
		unsigned char* data1 = buf1
		unsigned char* data2 = buf2
		Py_ssize_t record_start1 = 0
		Py_ssize_t record_start2 = 0

	while True:
		while pos1 < end1 and data1[pos1] != '\n':
			pos1 += 1
		if pos1 == end1:
			break
		pos1 += 1
		while pos2 < end2 and data2[pos2] != '\n':
			pos2 += 1
		if pos2 == end2:
			break
		pos2 += 1
		linebreaks += 1
		if linebreaks == 4:
			linebreaks = 0
			record_start1 = pos1
			record_start2 = pos2

	# Hit the end of the data block
	return record_start1, record_start2


cdef class Sequence(object):
	"""
	A record in a FASTQ file. Also used for FASTA (then the qualities attribute
	is None). qualities is a string and it contains the qualities encoded as
	ascii(qual+33).
	"""
	cdef:
		public str name
		public str sequence
		public str qualities
		public bint second_header

	def __init__(self, str name, str sequence, str qualities=None, bint second_header=False):
		"""Set qualities to None if there are no quality values"""
		self.name = name
		self.sequence = sequence
		self.qualities = qualities
		self.second_header = second_header
		if qualities is not None and len(qualities) != len(sequence):
			rname = _shorten(name)
			raise FormatError("In read named {0!r}: length of quality sequence ({1}) and length "
				"of read ({2}) do not match".format(
					rname, len(qualities), len(sequence)))
	
	def __getitem__(self, key):
		"""slicing"""
		return self.__class__(
			self.name,
			self.sequence[key],
			self.qualities[key] if self.qualities is not None else None,
			self.second_header)

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = ', qualities={0!r}'.format(_shorten(self.qualities))
		return '<Sequence(name={0!r}, sequence={1!r}{2})>'.format(_shorten(self.name), _shorten(self.sequence), qstr)

	def __len__(self):
		return len(self.sequence)

	def __richcmp__(self, other, int op):
		if 2 <= op <= 3:
			eq = self.name == other.name and \
				self.sequence == other.sequence and \
				self.qualities == other.qualities
			if op == 2:
				return eq
			else:
				return not eq
		else:
			raise NotImplementedError()

	def __reduce__(self):
		return (Sequence, (self.name, self.sequence, self.qualities, self.second_header))


class FastqReader(SequenceReader):
	"""
	Reader for FASTQ files. Does not support multi-line FASTQ files.
	"""
	def __init__(self, file, sequence_class=Sequence):
		"""
		file is a filename or a file-like object.
		If file is a filename, then .gz files are supported.
		"""
		super().__init__(file)
		self.sequence_class = sequence_class
		self.delivers_qualities = True

	def __iter__(self):
		"""
		Yield Sequence objects
		"""
		cdef int i = 0
		cdef int strip
		cdef str line, name, qualities, sequence, name2
		sequence_class = self.sequence_class

		it = iter(self._file)
		line = next(it)
		if not (line and line[0] == '@'):
			raise FormatError("Line {0} in FASTQ file is expected to start with '@', but found {1!r}".format(i+1, line[:10]))
		strip = -2 if line.endswith('\r\n') else -1
		name = line[1:strip]

		i = 1
		for line in it:
			if i == 0:
				if not (line and line[0] == '@'):
					raise FormatError("Line {0} in FASTQ file is expected to start with '@', but found {1!r}".format(i+1, line[:10]))
				name = line[1:strip]
			elif i == 1:
				sequence = line[:strip]
			elif i == 2:
				if line == '+\n':  # check most common case first
					name2 = ''
				else:
					line = line[:strip]
					if not (line and line[0] == '+'):
						raise FormatError("Line {0} in FASTQ file is expected to start with '+', but found {1!r}".format(i+1, line[:10]))
					if len(line) > 1:
						if not line[1:] == name:
							raise FormatError(
								"At line {0}: Sequence descriptions in the FASTQ file don't match "
								"({1!r} != {2!r}).\n"
								"The second sequence description must be either empty "
								"or equal to the first description.".format(i+1,
									name, line[1:]))
						second_header = True
					else:
						second_header = False
			elif i == 3:
				if len(line) == len(sequence) - strip:
					qualities = line[:strip]
				else:
					qualities = line.rstrip('\r\n')
				yield sequence_class(name, sequence, qualities, second_header=second_header)
			i = (i + 1) % 4
		if i != 0:
			raise FormatError("FASTQ file ended prematurely")
