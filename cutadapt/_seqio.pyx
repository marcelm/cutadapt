# kate: syntax Python;
# cython: profile=False
from __future__ import print_function, division, absolute_import
from .xopen import xopen

# TODO
# the following function and class cannot be imported from seqio.py
# since we would get circular imports

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


cdef class Sequence(object):
	"""
	A record in a FASTQ file. Also used for FASTA (then the qualities attribute
	is None). qualities is a string and it contains the qualities encoded as
	ascii(qual+33).

	If an adapter has been matched to the sequence, the 'match' attribute is
	set to the corresponding AdapterMatch instance.
	"""
	cdef:
		public str name
		public str sequence
		public str qualities
		public object match
		public bint twoheaders

	def __init__(self, str name, str sequence, str qualities=None,
			  bint twoheaders=False, match=None):
		"""Set qualities to None if there are no quality values"""
		self.name = name
		self.sequence = sequence
		self.qualities = qualities
		self.twoheaders = twoheaders
		self.match = match
		if qualities is not None and len(qualities) != len(sequence):
			rname = _shorten(name)
			raise FormatError("In read named {0!r}: length of quality sequence ({1}) and length of read ({2}) do not match".format(
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
		return (Sequence, (self.name, self.sequence, self.qualities, self.twoheaders))

	def write(self, outfile):
		if self.qualities is not None:
			s = '@' + self.name + '\n' + self.sequence + '\n+'
			if self.twoheaders:
				s += self.name
			s += '\n' + self.qualities + '\n'
		else:
			s = '>' + self.name + '\n' + self.sequence + '\n'
		outfile.write(s)


class FastqReader(object):
	"""
	Reader for FASTQ files. Does not support multi-line FASTQ files.
	"""
	def __init__(self, file, sequence_class=Sequence):
		"""
		file is a filename or a file-like object.
		If file is a filename, then .gz files are supported.

		colorspace -- Usually (when this is False), there must be n characters in the sequence and
		n quality values. When this is True, there must be n+1 characters in the sequence and n quality values.
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
		cdef int i = 0
		cdef int strip
		cdef str line, name, qualities, sequence
		cdef bint twoheaders
		sequence_class = self.sequence_class

		it = iter(self.fp)
		line = next(it)
		if not (line and line[0] == '@'):
			raise FormatError("at line {0}, expected a line starting with '@'".format(i+1))
		strip = -2 if line.endswith('\r\n') else -1
		name = line[1:strip]

		i = 1
		for line in it:
			if i == 0:
				if not (line and line[0] == '@'):
					raise FormatError("at line {0}, expected a line starting with '@'".format(i+1))
				name = line[1:strip]
			elif i == 1:
				sequence = line[:strip]
			elif i == 2:
				if line == '+\n':  # check most common case first
					twoheaders = False
				else:
					line = line[:strip]
					if not (line and line[0] == '+'):
						raise FormatError("at line {0}, expected a line starting with '+'".format(i+1))
					if len(line) > 1:
						twoheaders = True
						if not line[1:] == name:
							raise FormatError(
								"At line {0}: Sequence descriptions in the FASTQ file don't match "
								"({1!r} != {2!r}).\n"
								"The second sequence description must be either empty "
								"or equal to the first description.".format(i+1,
									name, line[1:]))
					else:
						twoheaders = False
			elif i == 3:
				if len(line) == len(sequence) - strip:
					qualities = line[:strip]
				else:
					qualities = line.rstrip('\r\n')
				yield sequence_class(name, sequence, qualities, twoheaders=twoheaders)
			i = (i + 1) % 4
		if i != 0:
			raise FormatError("FASTQ file ended prematurely")

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
