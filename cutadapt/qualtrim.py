# coding: utf-8
"""
Quality trimming.
"""
from __future__ import print_function, division, absolute_import

import sys

if sys.version > '3':
	xrange = range


def nextseq_trim_index(sequence, cutoff, base=33):
	"""
	Variant of the above quality trimming routine that works on NextSeq data.
	With Illumina NextSeq, bases are encoded with two colors. 'No color' (a
	dark cycle) usually means that a 'G' was sequenced, but that also occurs
	when sequencing falls off the end of the fragment. The read then contains
	a run of high-quality G bases in the end.

	This routine works as the one above, but counts qualities belonging to 'G'
	bases as being equal to cutoff - 1.
	"""
	bases = sequence.sequence
	qualities = sequence.qualities
	s = 0
	max_qual = 0
	max_i = len(qualities)
	for i in reversed(xrange(max_i)):
		q = ord(qualities[i]) - base
		if bases[i] == 'G':
			q = cutoff - 1
		s += cutoff - q
		if s < 0:
			break
		if s > max_qual:
			max_qual = s
			max_i = i
	return max_i


from cutadapt._qualtrim import quality_trim_index, nextseq_trim_index
