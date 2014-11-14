from __future__ import print_function, division, absolute_import
from nose.tools import raises

from cutadapt.seqio import Sequence
from cutadapt.adapters import Adapter, AdapterMatch, ColorspaceAdapter, FRONT, BACK

def test_issue_52():
	adapter = Adapter(
		sequence='GAACTCCAGTCACNNNNN',
		where=BACK,
		max_error_rate=0.1,
		min_overlap=5,
		read_wildcards=False,
		adapter_wildcards=True)
	read = Sequence(name="abc", sequence='CCCCAGAACTACAGTCCCGGC')
	am = AdapterMatch(astart=0, astop=17, rstart=5, rstop=21, matches=15, errors=2, front=None, adapter=adapter, read=read)
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

