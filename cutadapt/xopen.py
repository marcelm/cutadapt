"""
Open compressed files transparently.
"""
from __future__ import print_function, division, absolute_import

import gzip
import sys
import io
from subprocess import Popen, PIPE

try:
	import bz2
except ImportError:
	bz2 = None

__author__ = 'Marcel Martin'

if sys.version_info[0] >= 3:
	basestring = str
	from codecs import getreader, getwriter


if sys.version_info < (2, 7):
	buffered_reader = lambda x: x
	buffered_writer = lambda x: x
else:
	buffered_reader = io.BufferedReader
	buffered_writer = io.BufferedWriter


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
		if sys.version_info[0] < 3:
			return bz2.BZ2File(filename, mode)
		else:
			if 'r' in mode:
				return getreader('ascii')(bz2.BZ2File(filename, mode))
			else:
				return getwriter('ascii')(bz2.BZ2File(filename, mode))
	elif filename.endswith('.gz'):
		if sys.version_info[0] < 3:
			if 'r' in mode:
				try:
					return Popen(['gzip', '-cd', filename], stdout=PIPE).stdout
				except IOError:
					# gzip not installed
					return buffered_reader(gzip.open(filename, mode))
			else:
				try:
					f = open(filename, 'w')
					# TODO
					# f is not closed when the returned
					# file-like object is closed.
					return Popen(['gzip'], stdin=PIPE, stdout=f).stdin
				except IOError:
					return buffered_writer(gzip.open(filename, mode))
		else:
			if 'r' in mode:
				return getreader('ascii')(gzip.open(filename, mode))
			else:
				return getwriter('ascii')(gzip.open(filename, mode))
	else:
		return open(filename, mode)
