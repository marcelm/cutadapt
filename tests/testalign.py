from cutadapt.align import globalalign_locate
from cutadapt.adapters import BACK

def test_polya():
	s = 'AAAAAAAAAAAAAAAAA'
	t = 'ACAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
	result = globalalign_locate(s, t, 0.0, BACK)
	#start_s, stop_s, start_t, stop_t, matches, cost = result
	assert result == (0, len(s), 4, 4 + len(s), len(s), 0)
