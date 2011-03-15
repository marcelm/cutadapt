"""
Open compressed files transparently.
"""
import gzip
from contextlib import closing

__author__ = 'Marcel Martin'


def xopen(filename, mode='r', is_closing=True):
	"""
	Replacement for the "open" function that can also open
	files that have been compressed with gzip. If the filename ends with .gz,
	the file is opened with gzip.open(). If it doesn't, the regular open()
	is used. If the filename is '-', standard output (mode 'w') or input
	(mode 'r') is returned. TODO disabled for now
	closing -- whether to wrap the returned GzipFile with contextlib.closing
		to make it usable in 'with' statements. Don't use when filename may be '-'.
		(TODO look for a nicer solution)
	"""
	# disabled for now: does not work
	#if filename == '-':
		#return sys.stdout if mode == 'r' else sys.stdin
	if filename.endswith('.gz'):
		if is_closing:
			return closing(gzip.open(filename, mode))
		else:
			return gzip.open(filename, mode)
	else:
		return open(filename, mode)
