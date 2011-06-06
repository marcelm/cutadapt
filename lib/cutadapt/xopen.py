"""
Open compressed files transparently.
"""
import gzip

__author__ = 'Marcel Martin'

import sys
if sys.version_info[0] >= 3:
	basestring = str
	from codecs import getreader, getwriter


def xopen(filename, mode='r'):
	"""
	Replacement for the "open" function that can also open
	files that have been compressed with gzip. If the filename ends with .gz,
	the file is opened with gzip.open(). If it doesn't, the regular open()
	is used. If the filename is '-', standard output (mode 'w') or input
	(mode 'r') is returned.
	"""
	assert isinstance(filename, basestring)
	if filename == '-':
		return sys.stdin if 'r' in mode else sys.stdout
	if filename.endswith('.gz'):
		if sys.version_info[0] < 3:
			return gzip.open(filename, mode)
		else:
			if 'r' in mode:
				return getreader('ascii')(gzip.open(filename, mode))
			else:
				return getwriter('ascii')(gzip.open(filename, mode))
	else:
		return open(filename, mode)
