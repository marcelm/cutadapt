"""
Open compressed files transparently.
"""
from __future__ import print_function, division, absolute_import
__author__ = 'Marcel Martin'

import gzip
import sys
import io
import os
from subprocess import Popen, PIPE
from .compat import PY3, basestring

try:
	import bz2
except ImportError:
	bz2 = None

try:
	import lzma
except ImportError:
	lzma = None

if sys.version_info < (2, 7):
	buffered_reader = lambda x: x
	buffered_writer = lambda x: x
else:
	buffered_reader = io.BufferedReader
	buffered_writer = io.BufferedWriter


class GzipWriter:
	def __init__(self, path, mode='w'):
		self.outfile = open(path, mode)
		self.devnull = open(os.devnull, 'w')
		try:
			# Setting close_fds to True is necessary due to
			# http://bugs.python.org/issue12786
			self.process = Popen(['gzip'], stdin=PIPE, stdout=self.outfile,
				stderr=self.devnull, close_fds=True)
		except IOError as e:
			self.outfile.close()
			self.devnull.close()
			raise

	def write(self, arg):
		self.process.stdin.write(arg)

	def close(self):
		self.process.stdin.close()
		retcode = self.process.wait()
		self.outfile.close()
		self.devnull.close()
		if retcode != 0:
			raise IOError("Output gzip process terminated with exit code {0}".format(retcode))

	def __enter__(self):
		return self

	def __exit__(self, *exc_info):
		self.close()


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
		self.process.wait()
		self._raise_if_error()

	def _raise_if_error(self):
		"""
		Raise EOFError if process is not running anymore and the
		exit code is nonzero.
		"""
		retcode = self.process.poll()
		if retcode is not None and retcode != 0:
			raise EOFError("gzip process returned non-zero exit code {0}. Is the "
				"input file truncated or corrupt?".format(retcode))

	def read(self, *args):
		data = self.process.stdout.read(*args)
		if len(args) == 0 or args[0] <= 0:
			# wait for process to terminate until we check the exit code
			self.process.wait()
		self._raise_if_error()

	def __enter__(self):
		return self

	def __exit__(self, *exc_info):
		self.close()


def xopen(filename, mode='r'):
	"""
	Replacement for the "open" function that can also open files that have
	been compressed with gzip, bzip2 or xz. If the filename is '-', standard
	output (mode 'w') or input (mode 'r') is returned. If the filename ends
	with .gz, the file is opened with a pipe to the gzip program. If that
	does not work, then gzip.open() is used (the gzip module is slower than
	the pipe to the gzip program). If the filename ends with .bz2, it's
	opened as a bz2.BZ2File. Otherwise, the regular open() is used.

	mode can be: 'rt', 'rb', 'a', 'wt', or 'wb'
	Instead of 'rt' and 'wt', 'r' and 'w' can be used as abbreviations.

	In Python 2, the 't' and 'b' characters are ignored.

	Append mode ('a') is unavailable with BZ2 compression and will raise an error.
	"""
	if mode == 'r':
		mode = 'rt'
	elif mode == 'w':
		mode = 'wt'
	if mode not in ('rt', 'rb', 'wt', 'wb', 'a'):
		raise ValueError("mode '{0}' not supported".format(mode))
	if not PY3:
		mode = mode[0]
	if not isinstance(filename, basestring):
		raise ValueError("the filename must be a string")

	# standard input and standard output handling
	if filename == '-':
		if not PY3:
			return sys.stdin if 'r' in mode else sys.stdout
		return dict(
			rt=sys.stdin,
			wt=sys.stdout,
			rb=sys.stdin.buffer,
			wb=sys.stdout.buffer)[mode]

	if filename.endswith('.bz2'):
		if bz2 is None:
			raise ImportError("Cannot open bz2 files: The bz2 module is not available")
		if PY3:
			if 't' in mode:
				return io.TextIOWrapper(bz2.BZ2File(filename, mode[0]))
			else:
				return bz2.BZ2File(filename, mode)
		else:
			return bz2.BZ2File(filename, mode)
	elif filename.endswith('.xz'):
		if lzma is None:
			raise ImportError("Cannot open xz files: The lzma module is not available "
				"(use Python 3.3 or newer)")
		return lzma.open(filename, mode)
	elif filename.endswith('.gz'):
		if PY3:
			if 't' in mode:
				# gzip.open in Python 3.2 does not support modes 'rt' and 'wt''
				return io.TextIOWrapper(gzip.open(filename, mode[0]))
			else:
				if 'r' in mode:
					return io.BufferedReader(gzip.open(filename, mode))
				else:
					return io.BufferedWriter(gzip.open(filename, mode))
		else:
			# rb/rt are equivalent in Py2
			if 'r' in mode:
				try:
					return GzipReader(filename)
				except IOError:
					# gzip not installed
					return buffered_reader(gzip.open(filename, mode))
			else:
				try:
					return GzipWriter(filename, mode)
				except IOError:
					return buffered_writer(gzip.open(filename, mode))
	else:
		return open(filename, mode)
