"""
Open compressed files transparently.
"""
import gzip

__author__ = 'Marcel Martin'


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
		return sys.stdout if 'r' in mode else sys.stdin
	if filename.endswith('.gz'):
		return gzip.open(filename, mode)
	else:
		return open(filename, mode)
