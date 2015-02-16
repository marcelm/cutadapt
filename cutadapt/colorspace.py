# coding: utf-8
"""
Colorspace conversion routines.

Inspired by agapython/util/Dibase.py from Corona lite,
but reimplemented to avoid licensing issues.

Encoding Table

  A C G T
A 0 1 2 3
C 1 0 3 2
G 2 3 0 1
T 3 2 1 0
"""
from __future__ import print_function, division, absolute_import

__author__ = 'Marcel Martin'


def _initialize_dicts():
	"""
	Create the colorspace encoding and decoding dictionaries.
	"""
	enc = {}
	for i, c1 in enumerate("ACGT"):
		enc['N' + c1] = '4'
		enc[c1 + 'N'] = '4'
		enc['.' + c1] = '4'
		enc[c1 + '.'] = '4'
		for j, c2 in enumerate("ACGT"):
			# XOR of nucleotides gives color
			enc[c1 + c2] = chr(ord('0') + (i ^ j))
	enc.update({ 'NN': '4', 'N.': '4', '.N': '4', '..': '4'})

	dec = {}
	for i, c1 in enumerate("ACGT"):
		dec['.' + str(i)] = 'N'
		dec['N' + str(i)] = 'N'
		dec[c1 + '4'] = 'N'
		dec[c1 + '.'] = 'N'
		for j, c2 in enumerate("ACGT"):
			# XOR of nucleotides gives color
			dec[c1 + chr(ord('0') + (i ^ j))] = c2
	dec['N4'] = 'N'

	return (enc, dec)


def encode(s):
	"""
	Given a sequence of nucleotides, convert them to
	colorspace. Only uppercase characters are allowed.
	>>> encode("ACGGTC")
	"A13012"
	"""
	if not s:
		return s
	r = s[0:1]
	for i in range(len(s) - 1):
		r += ENCODE[s[i:i+2]]
	return r


def decode(s):
	"""
	Decode a sequence of colors to nucleotide space.
	The first character in s must be a nucleotide.
	Only uppercase characters are allowed.
	>>> decode("A13012")
	"ACGGTC"
	"""
	if len(s) < 2:
		return s
	x = s[0]
	result = x
	for c in s[1:]:
		x = DECODE[x + c]
		result += x
	return result


(ENCODE, DECODE) = _initialize_dicts()
