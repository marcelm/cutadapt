# kate: syntax Python;
# cython: profile=False
from __future__ import print_function, division, absolute_import
from .xopen import xopen

from cpython.version cimport PY_MAJOR_VERSION
PY3 = PY_MAJOR_VERSION >= 3

cdef bytes str_to_bytes(object s):
	if PY3:
		return s.encode('ascii')
	else:
		return s

cdef str bytes_to_str(object s):
	if PY3:
		return s.decode('ascii')
	else:
		return s

## TODO the following method and class are copied here from seqio.py to avoid circular imports for now

def _shorten(s, n=20):
	"""Shorten string s to at most n characters, appending "..." if necessary."""
	if s is None:
		return None
	if len(s) > n:
		s = s[:n-3] + '...'
	return s


class FormatError(Exception):
	"""
	Raised when an input file (FASTA or FASTQ) is malformatted.
	"""
	pass


cdef class Sequence(object):
	"""qualities is a string and it contains the qualities encoded as ascii(qual+33)."""
	cdef:
		public str name
		public bytes sequence
		public bytes qualities

	def __init__(self, str name, bytes sequence, bytes qualities=None):
		"""Set qualities to None if there are no quality values"""
		self.name = name
		self.sequence = sequence
		self.qualities = qualities
		if qualities is not None:
			if len(qualities) != len(sequence):
				rname = _shorten(name)
				raise ValueError("In read named {0!r}: length of quality sequence and length of read do not match ({1}!={2})".format(
					rname, len(qualities), len(sequence)))

	def __getitem__(self, key):
		"""slicing"""
		return self.__class__(self.name, self.sequence[key], self.qualities[key] if self.qualities is not None else None)

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = '\', qualities={0!r}'.format(_shorten(self.qualities))
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

	def write(self, outfile, bint twoheaders=False):
		cdef bytes name
		name = str_to_bytes(self.name)
		if self.qualities is not None:
			s = b'@' + name + b'\n' + self.sequence + b'\n+'
			if twoheaders:
				s += name
			s += b'\n' + self.qualities + b'\n'
		else:
			s = b'>' + name + b'\n' + self.sequence + b'\n'
		if PY3:
			outfile.buffer.write(s)
		else:
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
			file = xopen(file, 'rb')
		self.fp = file
		self.twoheaders = False
		self.sequence_class = sequence_class

	def __iter__(self):
		"""
		Return tuples: (name, sequence, qualities).
		qualities is a string and it contains the unmodified, encoded qualities.
		"""
		cdef int i = 0
		cdef bytes line, name, qualities, sequence

		for line in self.fp:
			if i % 4 == 0:
				if not line.startswith(b'@'):
					raise FormatError("at line {0}, expected a line starting with '+'".format(i+1))
				name = line.rstrip(b'\r\n')[1:]
			elif i % 4 == 1:
				sequence = line.rstrip(b'\r\n')
			elif i % 4 == 2:
				line = line.rstrip(b'\r\n')
				if not line.startswith(b'+'):
					raise FormatError("at line {0}, expected a line starting with '+'".format(i+1))
				if len(line) > 1:
					self.twoheaders = True
					if not line[1:] == name:
						raise FormatError(
							"At line {0}: Sequence descriptions in the FASTQ file don't match "
							"({1!r} != {2!r}).\n"
							"The second sequence description must be either empty "
							"or equal to the first description.".format(i+1,
								bytes_to_str(name), bytes_to_str(line.rstrip()[1:])))
			elif i % 4 == 3:
				qualities = line.rstrip(b'\r\n')
				yield self.sequence_class(bytes_to_str(name), sequence, qualities)
			i += 1

	def __enter__(self):
		if self.fp is None:
			raise ValueError("I/O operation on closed FastqReader")
		return self

	def __exit__(self, *args):
		self.fp.close()
