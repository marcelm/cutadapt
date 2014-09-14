"""
Open compressed files transparently.
"""
from __future__ import print_function, division, absolute_import

import gzip
import sys
import io
from subprocess import Popen, PIPE
from cutadapt.compat import PY3, basestring

try:
	import bz2
except ImportError:
	bz2 = None

__author__ = 'Marcel Martin'


if sys.version_info < (2, 7):
	buffered_reader = lambda x: x
	buffered_writer = lambda x: x
else:
	buffered_reader = io.BufferedReader
	buffered_writer = io.BufferedWriter


class GzipWriter:
	def __init__(self, path):
		self.outfile = open(path, 'w')
		try:
			self.process = Popen(['gzip'], stdin=PIPE, stdout=self.outfile)
		except IOError as e:
			self.outfile.close()
			raise

	def write(self, arg):
		self.process.stdin.write(arg)

	def close(self):
		self.process.stdin.close()
		c = self.process.wait()
		if c != 0:
			raise IOError("Output gzip process terminated with exit code {}".format(c))


def xopen(filename, mode='r'):
	"""
	Replacement for the "open" function that can also open files that have
	been compressed with gzip or bzip2. If the filename is '-', standard
	output (mode 'w') or input (mode 'r') is returned. If the filename ends
	with .gz, the file is opened with a pipe to the gzip program. If that
	does not work, then gzip.open() is used (the gzip module is slower than
	the pipe to the gzip program). If the filename ends with .bz2, it's
	opened as a bz2.BZ2File. Otherwise, the regular open() is used.
	"""
	assert isinstance(filename, basestring)
	if filename == '-':
		return sys.stdin if 'r' in mode else sys.stdout
	if filename.endswith('.bz2'):
		if bz2 is None:
			raise ImportError("Cannot open bz2 files: The bz2 module is not available")
		if PY3:
			return io.TextIOWrapper(bz2.BZ2File(filename, mode))
		else:
			return bz2.BZ2File(filename, mode)

	elif filename.endswith('.gz'):
		if PY3:
			return io.TextIOWrapper(gzip.open(filename, mode))
		else:
			if 'r' in mode:
				try:
					return Popen(['gzip', '-cd', filename], stdout=PIPE).stdout
				except IOError:
					# gzip not installed
					return buffered_reader(gzip.open(filename, mode))
			else:
				try:
					return GzipWriter(filename)
				except IOError:
					return buffered_writer(gzip.open(filename, mode))
	else:
		return open(filename, mode)
