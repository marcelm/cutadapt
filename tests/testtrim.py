# coding: utf-8
from __future__ import print_function, division, absolute_import

from cutadapt.seqio import ColorspaceSequence, Sequence
from cutadapt.adapters import Adapter, ColorspaceAdapter, PREFIX, BACK
from cutadapt.modifiers import AdapterCutter


def test_cs_5p():
	read = ColorspaceSequence("name", "0123", "DEFG", "T")
	adapter = ColorspaceAdapter("CG", PREFIX, 0.1)
	cutter = AdapterCutter([adapter])
	trimmed_read = cutter(read)
	# no assertion here, just make sure the above code runs without
	# an exception


def test_statistics():
	read = Sequence('name', 'AAAACCCCAAAA')
	adapters = [Adapter('CCCC', BACK, 0.1)]
	cutter = AdapterCutter(adapters, times=3)
	trimmed_read = cutter(read)
	# TODO make this a lot simpler
	trimmed_bp = 0
	for adapter in adapters:
		for d in (cutter.adapter_statistics[adapter].lengths_front, cutter.adapter_statistics[adapter].lengths_back):
			trimmed_bp += sum(seqlen * count for (seqlen, count) in d.items())
	assert trimmed_bp <= len(read), trimmed_bp


def test_end_trim_with_mismatch():
    """
    Test the not-so-obvious case where an adapter of length 13 is trimmed from
    the end of a sequence with overlap 9 and there is one deletion.
    In this case the algorithm starts with 10 bases of the adapter to get
    the hit and so the match is considered good. An insertion or substitution
    at the same spot is not a match.
    """
    adapter = Adapter('TCGATCGATCGAT', BACK, 0.1)

    read = Sequence('foo1', 'AAAAAAAAAAATCGTCGATC')
    cutter = AdapterCutter([adapter], times=1)
    trimmed_read = cutter(read)

    assert trimmed_read.sequence == 'AAAAAAAAAAA'
    assert cutter.adapter_statistics[adapter].lengths_back == {9: 1}
    # We see 1 error at length 9 even though the number of allowed mismatches at
    # length 9 is 0.
    assert cutter.adapter_statistics[adapter].errors_back[9][1] == 1

    read = Sequence('foo2', 'AAAAAAAAAAATCGAACGA')
    cutter = AdapterCutter([adapter], times=1)
    trimmed_read = cutter(read)

    assert trimmed_read.sequence == read.sequence
    assert cutter.adapter_statistics[adapter].lengths_back == {}
