# coding: utf-8
"""
Tests write output (should it return True or False or write)
"""
from __future__ import print_function, division, absolute_import

from cutadapt.writers import NContentTrimmer, DISCARD, KEEP
from cutadapt.seqio import Sequence

def test_ncontenttrimmer():
	# third parameter is True if read should be discarded
	params = [
		('AAA', 0, KEEP),
		('AAA', 1, KEEP),
		('AAACCTTGGN', 1, KEEP),
		('AAACNNNCTTGGN', 0.5, KEEP),
		('NNNNNN', 1, DISCARD),
		('ANAAAA', 1/6, KEEP),
		('ANAAAA', 0, DISCARD)
	]
	for seq, count, expected in params:
		writer = NContentTrimmer(count=count)
		_seq = Sequence('read1', seq, qualities='#'*len(seq))
		assert writer(_seq) == expected


def test_ncontenttrimmer_paired():
	params = [
		('AAA', 'AAA', 0, KEEP),
		('AAAN', 'AAA', 0, DISCARD),
		('AAA', 'AANA', 0, DISCARD),
		('ANAA', 'AANA', 1, KEEP),
	]
	for seq1, seq2, count, expected in params:
		writer = NContentTrimmer(count=count, check_second=False)
		writer_cs = NContentTrimmer(count=count, check_second=True)
		read1 = Sequence('read1', seq1, qualities='#'*len(seq1))
		read2 = Sequence('read1', seq2, qualities='#'*len(seq2))
		assert writer(read1, read2) == writer(read1)
		# discard entire pair if one of the reads fulfills criteria
		assert writer_cs(read1, read2) == expected
