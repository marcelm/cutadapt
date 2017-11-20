# coding: utf-8
from __future__ import print_function, division, absolute_import
from nose.tools import raises, assert_raises

from cutadapt.seqio import Sequence
from cutadapt.adapters import (Adapter, Match, ColorspaceAdapter, FRONT, BACK,
	parse_braces, LinkedAdapter, AdapterStatistics)


def test_issue_52():
	adapter = Adapter(
		sequence='GAACTCCAGTCACNNNNN',
		where=BACK,
		max_error_rate=0.12,
		min_overlap=5,
		read_wildcards=False,
		adapter_wildcards=True)
	read = Sequence(name="abc", sequence='CCCCAGAACTACAGTCCCGGC')
	am = Match(astart=0, astop=17, rstart=5, rstop=21, matches=15, errors=2,
		remove_before=False, adapter=adapter, read=read)
	assert am.wildcards() == 'GGC'
	"""
	The result above should actually be 'CGGC' since the correct
	alignment is this one:

	adapter         GAACTCCAGTCACNNNNN
	mismatches           X     X
	read       CCCCAGAACTACAGTC-CCGGC

	Since we do not keep the alignment, guessing 'GGC' is the best we
	can currently do.
	"""


def test_issue_80():
	# This issue turned out to not be an actual issue with the alignment
	# algorithm. The following alignment is found because it has more matches
	# than the 'obvious' one:
	#
	# TCGTATGCCGTCTTC
	# =========X==XX=
	# TCGTATGCCCTC--C
	#
	# This is correct, albeit a little surprising, since an alignment without
	# indels would have only two errors.

	adapter = Adapter(
		sequence="TCGTATGCCGTCTTC",
		where=BACK,
		max_error_rate=0.2,
		min_overlap=3,
		read_wildcards=False,
		adapter_wildcards=False)
	read = Sequence(name="seq2", sequence="TCGTATGCCCTCC")
	result = adapter.match_to(read)
	assert result.errors == 3, result
	assert result.astart == 0, result
	assert result.astop == 15, result


def test_str():
	a = Adapter('ACGT', where=BACK, max_error_rate=0.1)
	str(a)
	str(a.match_to(Sequence(name='seq', sequence='TTACGT')))
	ca = ColorspaceAdapter('0123', where=BACK, max_error_rate=0.1)
	str(ca)


@raises(ValueError)
def test_color():
	ColorspaceAdapter('0123', where=FRONT, max_error_rate=0.1)


def test_parse_braces():
	assert parse_braces('') == ''
	assert parse_braces('A') == 'A'
	assert parse_braces('A{0}') == ''
	assert parse_braces('A{1}') == 'A'
	assert parse_braces('A{2}') == 'AA'
	assert parse_braces('A{2}C') == 'AAC'
	assert parse_braces('ACGTN{3}TGACCC') == 'ACGTNNNTGACCC'
	assert parse_braces('ACGTN{10}TGACCC') == 'ACGTNNNNNNNNNNTGACCC'
	assert parse_braces('ACGTN{3}TGA{4}CCC') == 'ACGTNNNTGAAAACCC'
	assert parse_braces('ACGTN{0}TGA{4}CCC') == 'ACGTTGAAAACCC'


def test_parse_braces_fail():
	for expression in ['{', '}', '{}', '{5', '{1}', 'A{-7}', 'A{', 'A{1', 'N{7', 'AN{7', 'A{4{}',
			'A{4}{3}', 'A{b}', 'A{6X}', 'A{X6}']:
		assert_raises(ValueError, lambda: parse_braces(expression))


def test_linked_adapter():
	linked_adapter = LinkedAdapter('AAAA', 'TTTT', min_overlap=4)
	assert linked_adapter.front_adapter.min_overlap == 4
	assert linked_adapter.back_adapter.min_overlap == 4

	sequence = Sequence(name='seq', sequence='AAAACCCCCTTTT')
	trimmed = linked_adapter.match_to(sequence).trimmed()
	assert trimmed.name == 'seq'
	assert trimmed.sequence == 'CCCCC'


def test_info_record():
	adapter = Adapter(
		sequence='GAACTCCAGTCACNNNNN',
		where=BACK,
		max_error_rate=0.12,
		min_overlap=5,
		read_wildcards=False,
		adapter_wildcards=True,
		name="Foo")
	read = Sequence(name="abc", sequence='CCCCAGAACTACAGTCCCGGC')
	am = Match(astart=0, astop=17, rstart=5, rstop=21, matches=15, errors=2, remove_before=False,
		adapter=adapter, read=read)
	assert am.get_info_record() == (
		"abc",
		2,
		5,
		21,
		'CCCCA',
		'GAACTACAGTCCCGGC',
		'',
		'Foo',
		'', 
		'', 
		''
	)


def test_random_match_probabilities():
	a = Adapter('A', where=BACK, max_error_rate=0.1).create_statistics()
	assert a.back.random_match_probabilities(0.5) == [1, 0.25]
	assert a.back.random_match_probabilities(0.2) == [1, 0.4]

	for s in ('ACTG', 'XMWH'):
		a = Adapter(s, where=BACK, max_error_rate=0.1).create_statistics()
		assert a.back.random_match_probabilities(0.5) == [1, 0.25, 0.25**2, 0.25**3, 0.25**4]
		assert a.back.random_match_probabilities(0.2) == [1, 0.4, 0.4*0.1, 0.4*0.1*0.4, 0.4*0.1*0.4*0.1]

	a = Adapter('GTCA', where=FRONT, max_error_rate=0.1).create_statistics()
	assert a.front.random_match_probabilities(0.5) == [1, 0.25, 0.25**2, 0.25**3, 0.25**4]
	assert a.front.random_match_probabilities(0.2) == [1, 0.4, 0.4*0.1, 0.4*0.1*0.4, 0.4*0.1*0.4*0.1]


def test_add_adapter_statistics():
	stats = Adapter('A', name='name', where=BACK, max_error_rate=0.1).create_statistics()
	end_stats = stats.back
	end_stats.adjacent_bases['A'] = 7
	end_stats.adjacent_bases['C'] = 19
	end_stats.adjacent_bases['G'] = 23
	end_stats.adjacent_bases['T'] = 42
	end_stats.adjacent_bases[''] = 45

	end_stats.errors[10][0] = 100
	end_stats.errors[10][1] = 11
	end_stats.errors[10][2] = 3
	end_stats.errors[20][0] = 600
	end_stats.errors[20][1] = 66
	end_stats.errors[20][2] = 6

	stats2 = Adapter('A', name='name', where=BACK, max_error_rate=0.1).create_statistics()
	end_stats2 = stats2.back
	end_stats2.adjacent_bases['A'] = 43
	end_stats2.adjacent_bases['C'] = 31
	end_stats2.adjacent_bases['G'] = 27
	end_stats2.adjacent_bases['T'] = 8
	end_stats2.adjacent_bases[''] = 5
	end_stats2.errors[10][0] = 234
	end_stats2.errors[10][1] = 14
	end_stats2.errors[10][3] = 5
	end_stats2.errors[15][0] = 90
	end_stats2.errors[15][1] = 17
	end_stats2.errors[15][2] = 2

	stats += stats2
	r = stats.back

	assert r.adjacent_bases == {'A': 50, 'C': 50, 'G': 50, 'T': 50, '': 50}
	assert r.errors == {
		10: {0: 334, 1: 25, 2: 3, 3: 5},
		15: {0: 90, 1: 17, 2: 2},
		20: {0: 600, 1: 66, 2: 6},
	}


def test_issue_265():
	"""Crash when accessing the matches property of non-anchored linked adapters"""
	s = Sequence('name', 'AAAATTTT')
	la = LinkedAdapter('GGG', 'TTT', front_anchored=False, back_anchored=False)
	assert la.match_to(s).matches == 3
