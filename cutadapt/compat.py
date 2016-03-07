# coding: utf-8
"""
Minimal Py2/Py3 compatibility library.
"""
from __future__ import print_function, division, absolute_import
import sys
PY3 = sys.version > '3'


if PY3:
	maketrans = str.maketrans
	basestring = str
	zip = zip
	next = next

	def bytes_to_str(s):
		return s.decode('ascii')

	def str_to_bytes(s):
		return s.encode('ascii')

	def force_str(s):
		if isinstance(s, bytes):
			return s.decode('ascii')
		else:
			return s
	from io import StringIO

else:
	def bytes_to_str(s):
		return s

	def str_to_bytes(s):
		return s

	def force_str(s):
		return s

	def next(it):
		return it.next()

	from string import maketrans
	basestring = basestring
	from itertools import izip as zip
	from StringIO import StringIO
