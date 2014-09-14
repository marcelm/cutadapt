"""
Quality trimming.
"""
from __future__ import print_function, division, absolute_import

import sys

if sys.version > '3':
	xrange = range


def quality_trim_index(qualities, cutoff, base=33):
	"""
	Find the position at which to trim a low-quality end from a nucleotide sequence.

	Qualities are assumed to be ASCII-encoded as chr(qual + base).

	The algorithm is the same as the one used by BWA within the function
	'bwa_trim_read':
	- Subtract the cutoff value from all qualities.
	- Compute partial sums from all indices to the end of the sequence.
	- Trim sequence at the index at which the sum is minimal.
	"""
	s = 0
	max_qual = 0
	max_i = len(qualities)
	for i in reversed(xrange(max_i)):
		q = ord(qualities[i]) - base
		s += cutoff - q
		if s < 0:
			break
		if s > max_qual:
			max_qual = s
			max_i = i
	return max_i

try:
	from cutadapt._qualtrim import quality_trim_index
except:
	pass
