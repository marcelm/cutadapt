from __future__ import print_function, division, absolute_import

from cutadapt.seqio import Sequence
from cutadapt.adapters import Adapter, AdapterMatch, BACK

def test_issue_52():
	adapter = Adapter(
		sequence=b"GAACTCCAGTCACNNNNN",
		where=BACK,
		max_error_rate=0.1,
		min_overlap=5,
		match_read_wildcards=False,
		match_adapter_wildcards=True)
	read = Sequence(name="abc", sequence=b'CCCCAGAACTACAGTCCCGGC')
	am = AdapterMatch(astart=0, astop=17, rstart=5, rstop=21, matches=15, errors=2, front=None, adapter=adapter, read=read)
	assert am.wildcards() == b'GGC'
	"""
	The result above should actually be 'CGGC' since the correct
	alignment is this one:

	adapter         GAACTCCAGTCACNNNNN
	mismatches           X     X
	read       CCCCAGAACTACAGTC-CCGGC

	Since we do not keep the alignment, guessing 'GGC' is the best we
	can currently do.
	"""
