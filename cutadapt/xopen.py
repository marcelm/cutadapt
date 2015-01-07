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
		retcode = self.process.wait()
		if retcode != 0:
			raise IOError("Output gzip process terminated with exit code {0}".format(retcode))


class GzipReader:
	def __init__(self, path):
		self.process = Popen(['gzip', '-cd', path], stdout=PIPE)

	def close(self):
		retcode = self.process.poll()
		if retcode is None:
			# still running
			self.process.terminate()
		self._raise_if_error()

	def __iter__(self):
		for line in self.process.stdout:
			yield line
		self._raise_if_error()

	def _raise_if_error(self):
		"""
		Raise EOFError if process if process is not running anymore and the
		exit code is nonzero.
		"""
		retcode = self.process.poll()
		if retcode is not None and retcode != 0:
			raise EOFError("gzip process returned non-zero exit code {0}. Is the input file truncated or corrupt?".format(retcode))

	def read(self, *args):
		data = self.process.stdout.read(*args)
		if len(args) == 0:
			# wait for process to terminate until we check the exit code
			self.process.wait()
		self._raise_if_error()


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
					return GzipReader(filename)
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
