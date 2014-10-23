from __future__ import print_function, division, absolute_import

from cutadapt.align import (locate, compare_prefixes, compare_suffixes,
	ALLOW_WILDCARD_SEQ1, ALLOW_WILDCARD_SEQ2)
from cutadapt.adapters import BACK

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
	assert compare_prefixes('AANAA', 'AACAATTTTTTTTT', ALLOW_WILDCARD_SEQ1) == (0, 5, 0, 5, 5, 0)
	assert compare_prefixes('AANAA', 'AACAATTTTTTTTT', ALLOW_WILDCARD_SEQ1) == (0, 5, 0, 5, 5, 0)
	assert compare_prefixes('XAAAAA', 'AAAAATTTTTTTTT') == (0, 6, 0, 6, 4, 2)

	a = WILDCARD_SEQUENCES[0]
	for s in WILDCARD_SEQUENCES:
		r = s + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
		result = compare_prefixes(a, r, degenerate=ALLOW_WILDCARD_SEQ2)
		assert result == (0, 10, 0, 10, 10, 0), result

		result = compare_prefixes(r, a, degenerate=ALLOW_WILDCARD_SEQ1)
		assert result == (0, 10, 0, 10, 10, 0)

	for s in WILDCARD_SEQUENCES:
		for t in WILDCARD_SEQUENCES:
			r = s + 'GCCAGGG'
			result = compare_prefixes(s, r, degenerate=ALLOW_WILDCARD_SEQ1|ALLOW_WILDCARD_SEQ2)
			assert result == (0, 10, 0, 10, 10, 0)

			result = compare_prefixes(r, s, degenerate=ALLOW_WILDCARD_SEQ1|ALLOW_WILDCARD_SEQ2)
			assert result == (0, 10, 0, 10, 10, 0)

	r = WILDCARD_SEQUENCES[0] + 'GCCAGG'
	for deg in 0, ALLOW_WILDCARD_SEQ1, ALLOW_WILDCARD_SEQ2, ALLOW_WILDCARD_SEQ1|ALLOW_WILDCARD_SEQ2:
		result = compare_prefixes('CCCXTTXATC', r, degenerate=deg)
		assert result == (0, 10, 0, 10, 8, 2)


def test_compare_suffixes():
	assert compare_suffixes('AAXAA', 'TTTTTTTAAAAA') == (0, 5, 7, 12, 4, 1)
	assert compare_suffixes('AANAA', 'TTTTTTTAACAA', ALLOW_WILDCARD_SEQ1) == (0, 5, 7, 12, 5, 0)
	assert compare_suffixes('AANAA', 'TTTTTTTAACAA', ALLOW_WILDCARD_SEQ1) == (0, 5, 7, 12, 5, 0)
	assert compare_suffixes('AAAAAX', 'TTTTTTTAAAAA') == (0, 6, 6, 12, 4, 2)


def test_wildcards_in_adapter():
	r = 'CATCTGTCC' + WILDCARD_SEQUENCES[0] + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
	for a in WILDCARD_SEQUENCES:
		result = locate(a, r, 0.0, BACK, degenerate=ALLOW_WILDCARD_SEQ1)
		assert result == (0, 10, 9, 19, 10, 0), result

	a = 'CCCXTTXATC'
	result = locate(a, r, 0.0, BACK, degenerate=ALLOW_WILDCARD_SEQ1)
	assert result is None


def test_wildcards_in_read():
	a = WILDCARD_SEQUENCES[0]
	for s in WILDCARD_SEQUENCES:
		r = 'CATCTGTCC' + s + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
		result = locate(a, r, 0.0, BACK, degenerate=ALLOW_WILDCARD_SEQ2)
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
			result = locate(a, r, 0.0, BACK, degenerate=ALLOW_WILDCARD_SEQ1|ALLOW_WILDCARD_SEQ2)
			assert result == (0, 10, 9, 19, 10, 0), result


def test_no_match():
	a = locate('CTGATCTGGCCG', 'AAAAGGG', 0.1, BACK)
	assert a is None, a
