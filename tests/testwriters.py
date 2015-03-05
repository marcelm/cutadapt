# coding: utf-8
"""
Tests write output (should it return True or False or write)
"""
from __future__ import print_function, division, absolute_import

from cutadapt.writers import NContentTrimmer
from cutadapt.seqio import Sequence

def test_ncontenttrimmer():
	seqs = ['AAACCTTGGN', 'AAACNNNCTTGGN', 'NNNNNN', 'ANAAAA', 'ANAAAA']
	params = [(1, True), (0.5, True), (1, False), (1/6, True), (0, False)]
	for seq, trim_params in zip(seqs, params):
		count, result = trim_params
		writer = NContentTrimmer(count=count)
		_seq = Sequence('read1', seq, qualities='#'*len(seq))
		assert writer(_seq) == result