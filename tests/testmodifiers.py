# coding: utf-8
from __future__ import print_function, division, absolute_import

from cutadapt.seqio import Sequence
from cutadapt.modifiers import UnconditionalCutter, NEndTrimmer

def test_unconditional_cutter():
	uc = UnconditionalCutter(length=5)
	s = 'abcdefg'
	assert UnconditionalCutter(length=2)(s) == 'cdefg'
	assert UnconditionalCutter(length=-2)(s) == 'abcde'
	assert UnconditionalCutter(length=100)(s) == ''
	assert UnconditionalCutter(length=-100)(s) == ''

def test_nend_trimmer():
	trimmer = NEndTrimmer()
	seqs = ['NNNNAAACCTTGGNNN', 'NNNNAAACNNNCTTGGNNN', 'NNNNNN']
	trims = ['AAACCTTGG', 'AAACNNNCTTGG', '']
	for seq, trimmed in zip(seqs, trims):
		_seq = Sequence('read1', seq, qualities='#'*len(seq))
		_trimmed = Sequence('read1', trimmed, qualities='#'*len(trimmed))
		assert trimmer(_seq) == _trimmed