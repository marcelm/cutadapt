"""
Tests write output (should it return True or False or write)
"""
from cutadapt.filters import NContentFilter, DISCARD, KEEP, PairedRedirector
from dnaio import Sequence

from pytest import mark


@mark.parametrize('seq,count,expected', [
	('AAA', 0, KEEP),
	('AAA', 1, KEEP),
	('AAACCTTGGN', 1, KEEP),
	('AAACNNNCTTGGN', 0.5, KEEP),
	('NNNNNN', 1, DISCARD),
	('ANAAAA', 1 / 6, KEEP),
	('ANAAAA', 0, DISCARD),
])
def test_ncontentfilter(seq, count, expected):
	# third parameter is True if read should be discarded
	filter_ = NContentFilter(count=count)
	_seq = Sequence('read1', seq, qualities='#'*len(seq))
	assert filter_(_seq, []) == expected


@mark.parametrize('seq1,seq2,count,expected', [
	('AAA', 'AAA', 0, KEEP),
	('AAAN', 'AAA', 0, DISCARD),
	('AAA', 'AANA', 0, DISCARD),
	('ANAA', 'AANA', 1, KEEP),
])
def test_ncontentfilter_paired(seq1, seq2, count, expected):
	filter_ = NContentFilter(count=count)
	filter_legacy = PairedRedirector(None, filter_, filter_, pair_filter_mode='first')
	filter_any = PairedRedirector(None, filter_, filter_, pair_filter_mode='any')
	read1 = Sequence('read1', seq1, qualities='#'*len(seq1))
	read2 = Sequence('read1', seq2, qualities='#'*len(seq2))
	assert filter_legacy(read1, read2, [], []) == filter_(read1, [])
	# discard entire pair if one of the reads fulfills criteria
	assert filter_any(read1, read2, [], []) == expected
