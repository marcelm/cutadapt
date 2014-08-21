from __future__ import print_function, division, absolute_import

from cutadapt.align import (locate, compare_prefixes,
	ALLOW_WILDCARD_SEQ1, ALLOW_WILDCARD_SEQ1)
from cutadapt.adapters import BACK

def test_polya():
	s = b'AAAAAAAAAAAAAAAAA'
	t = b'ACAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
	result = locate(s, t, 0.0, BACK)
	#start_s, stop_s, start_t, stop_t, matches, cost = result
	assert result == (0, len(s), 4, 4 + len(s), len(s), 0)


def test_compare_prefixes():
	assert compare_prefixes(b'AAXAA', b'AAAAATTTTTTTTT') == (0, 5, 0, 5, 4, 1)
	assert compare_prefixes(b'AANAA', b'AACAATTTTTTTTT', ALLOW_WILDCARD_SEQ1) == (0, 5, 0, 5, 5, 0)
	assert compare_prefixes(b'AANAA', b'AACAATTTTTTTTT', ALLOW_WILDCARD_SEQ1) == (0, 5, 0, 5, 5, 0)
	assert compare_prefixes(b'XAAAAA', b'AAAAATTTTTTTTT') == (0, 6, 0, 6, 4, 2)
