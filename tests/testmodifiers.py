# coding: utf-8
from __future__ import print_function, division, absolute_import

from cutadapt.seqio import Sequence
from cutadapt.modifiers import UnconditionalCutter

def test_unconditional_cutter():
	uc = UnconditionalCutter(length=5)
	s = 'abcdefg'
	assert UnconditionalCutter(length=2)(s) == 'cdefg'
	assert UnconditionalCutter(length=-2)(s) == 'abcde'
	assert UnconditionalCutter(length=100)(s) == ''
	assert UnconditionalCutter(length=-100)(s) == ''
