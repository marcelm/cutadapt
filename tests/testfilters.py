# coding: utf-8
"""
Tests write output (should it return True or False or write)
"""
from __future__ import print_function, division, absolute_import

from cutadapt.filters import NContentFilter, DISCARD, KEEP, LegacyPairedRedirector, PairedRedirector
from cutadapt.seqio import Sequence

def test_ncontentfilter():
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
		filter = NContentFilter(count=count)
		_seq = Sequence('read1', seq, qualities='#'*len(seq))
		assert filter(_seq) == expected


def test_ncontentfilter_paired():
	params = [
		('AAA', 'AAA', 0, KEEP),
		('AAAN', 'AAA', 0, DISCARD),
		('AAA', 'AANA', 0, DISCARD),
		('ANAA', 'AANA', 1, KEEP),
	]
	for seq1, seq2, count, expected in params:
		filter = NContentFilter(count=count)
		filter_legacy = LegacyPairedRedirector(None, filter)
		filter_both = PairedRedirector(None, filter)
		read1 = Sequence('read1', seq1, qualities='#'*len(seq1))
		read2 = Sequence('read1', seq2, qualities='#'*len(seq2))
		assert filter_legacy(read1, read2) == filter(read1)
		# discard entire pair if one of the reads fulfills criteria
		assert filter_both(read1, read2) == expected
