from dnaio import Sequence
from cutadapt.qualtrim import nextseq_trim_index


def test_nextseq_trim():
	s = Sequence('n', '', '')
	assert nextseq_trim_index(s, cutoff=22) == 0
	s = Sequence('n',
		'TCTCGTATGCCGTCTTATGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGNNNNNNNNNNNGGNGG',
		'AA//EAEE//A6///E//A//EA/EEEEEEAEA//EEEEEEEEEEEEEEE###########EE#EA'
	)
	assert nextseq_trim_index(s, cutoff=22) == 33
