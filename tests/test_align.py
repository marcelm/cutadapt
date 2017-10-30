# coding: utf-8
from __future__ import print_function, division, absolute_import

from cutadapt.align import (locate, compare_prefixes, compare_suffixes,
	Aligner)
from cutadapt.adapters import BACK


class TestAligner():
	def test(self):
		reference = 'CTCCAGCTTAGACATATC'
		aligner = Aligner(reference, 0.1, flags=BACK)
		aligner.locate('CC')

	def test_100_percent_error_rate(self):
		reference = 'GCTTAGACATATC'
		aligner = Aligner(reference, 1.0, flags=BACK)
		aligner.locate('CAA')


def test_polya():
	s = 'AAAAAAAAAAAAAAAAA'
	t = 'ACAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
	result = locate(s, t, 0.0, BACK)
	#start_s, stop_s, start_t, stop_t, matches, cost = result
	assert result == (0, len(s), 4, 4 + len(s), len(s), 0)


# Sequences with IUPAC wildcards
# R=A|G, Y=C|T, S=G|C, W=A|T, K=G|T, M=A|C, B=C|G|T, D=A|G|T, H=A|C|T, V=A|C|G,
# N=A|C|G|T, X={}
WILDCARD_SEQUENCES = [
	'CCCATTGATC',  # original sequence without wildcards
	'CCCRTTRATC',  # R=A|G
	'YCCATYGATC',  # Y=C|T
	'CSSATTSATC',  # S=G|C
	'CCCWWWGATC',  # W=A|T
	'CCCATKKATC',  # K=G|T
	'CCMATTGMTC',  # M=A|C
	'BCCATTBABC',  # B=C|G|T
	'BCCATTBABC',  # B
	'CCCDTTDADC',  # D=A|G|T
	'CHCATHGATC',  # H=A|C|T
	'CVCVTTVATC',  # V=A|C|G
	'CCNATNGATC',  # N=A|C|G|T
	'CCCNTTNATC',  # N
#   'CCCXTTXATC',  # X
]


def test_compare_prefixes():
	assert compare_prefixes('AAXAA', 'AAAAATTTTTTTTT') == (0, 5, 0, 5, 4, 1)
	assert compare_prefixes('AANAA', 'AACAATTTTTTTTT', wildcard_ref=True) == (0, 5, 0, 5, 5, 0)
	assert compare_prefixes('AANAA', 'AACAATTTTTTTTT', wildcard_ref=True) == (0, 5, 0, 5, 5, 0)
	assert compare_prefixes('XAAAAA', 'AAAAATTTTTTTTT') == (0, 6, 0, 6, 4, 2)

	a = WILDCARD_SEQUENCES[0]
	for s in WILDCARD_SEQUENCES:
		r = s + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
		result = compare_prefixes(a, r, wildcard_query=True)
		assert result == (0, 10, 0, 10, 10, 0), result

		result = compare_prefixes(r, a, wildcard_ref=True)
		assert result == (0, 10, 0, 10, 10, 0)

	for s in WILDCARD_SEQUENCES:
		for t in WILDCARD_SEQUENCES:
			r = s + 'GCCAGGG'
			result = compare_prefixes(s, r, )
			assert result == (0, 10, 0, 10, 10, 0)

			result = compare_prefixes(r, s, wildcard_ref=True, wildcard_query=True)
			assert result == (0, 10, 0, 10, 10, 0)

	r = WILDCARD_SEQUENCES[0] + 'GCCAGG'
	for wildc_ref in (False, True):
		for wildc_query in (False, True):
			result = compare_prefixes('CCCXTTXATC', r, wildcard_ref=wildc_ref, wildcard_query=wildc_query)
			assert result == (0, 10, 0, 10, 8, 2)


def test_compare_suffixes():
	assert compare_suffixes('AAXAA', 'TTTTTTTAAAAA') == (0, 5, 7, 12, 4, 1)
	assert compare_suffixes('AANAA', 'TTTTTTTAACAA', wildcard_ref=True) == (0, 5, 7, 12, 5, 0)
	assert compare_suffixes('AANAA', 'TTTTTTTAACAA', wildcard_ref=True) == (0, 5, 7, 12, 5, 0)
	assert compare_suffixes('AAAAAX', 'TTTTTTTAAAAA') == (0, 6, 6, 12, 4, 2)


def test_wildcards_in_adapter():
	r = 'CATCTGTCC' + WILDCARD_SEQUENCES[0] + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
	for a in WILDCARD_SEQUENCES:
		result = locate(a, r, 0.0, BACK, wildcard_ref=True)
		assert result == (0, 10, 9, 19, 10, 0), result

	a = 'CCCXTTXATC'
	result = locate(a, r, 0.0, BACK, wildcard_ref=True)
	assert result is None


def test_wildcards_in_read():
	a = WILDCARD_SEQUENCES[0]
	for s in WILDCARD_SEQUENCES:
		r = 'CATCTGTCC' + s + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
		result = locate(a, r, 0.0, BACK, wildcard_query=True)
		if 'X' in s:
			assert result is None
		else:
			assert result == (0, 10, 9, 19, 10, 0), result


def test_wildcards_in_both():
	for a in WILDCARD_SEQUENCES:
		for s in WILDCARD_SEQUENCES:
			if 'X' in s or 'X' in a:
				continue
			r = 'CATCTGTCC' + s + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
			result = locate(a, r, 0.0, BACK, wildcard_ref=True, wildcard_query=True)
			assert result == (0, 10, 9, 19, 10, 0), result


def test_no_match():
	a = locate('CTGATCTGGCCG', 'AAAAGGG', 0.1, BACK)
	assert a is None, a
