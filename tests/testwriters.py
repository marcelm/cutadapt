# coding: utf-8
"""
Tests write output (should it return True or False or write)
"""
from __future__ import print_function, division, absolute_import

from cutadapt.writers import NContentTrimmer
from cutadapt.seqio import Sequence

def test_ncontenttrimmer():
	# third parameter is True if read should be discarded
	params = [
		('AAA', 0, False),
		('AAA', 1, False),
		('AAACCTTGGN', 1, False),
		('AAACNNNCTTGGN', 0.5, False),
		('NNNNNN', 1, True),
		('ANAAAA', 1/6, False),
		('ANAAAA', 0, True)
	]
	for seq, count, expected in params:
		writer = NContentTrimmer(count=count)
		_seq = Sequence('read1', seq, qualities='#'*len(seq))
		assert writer(_seq) == expected

