# coding: utf-8
from __future__ import print_function, division, absolute_import

from cutadapt.seqio import Sequence
from cutadapt.modifiers import Modifiers, UnconditionalCutter, NEndTrimmer, QualityTrimmer

def test_unconditional_cutter():
	uc = UnconditionalCutter(lengths=[5])
	s = 'abcdefg'
	assert UnconditionalCutter(lengths=[2])(s) == 'cdefg'
	assert UnconditionalCutter(lengths=[-2])(s) == 'abcde'
	assert UnconditionalCutter(lengths=[100])(s) == ''
	assert UnconditionalCutter(lengths=[-100])(s) == ''


def test_nend_trimmer():
	trimmer = NEndTrimmer()
	seqs = ['NNNNAAACCTTGGNNN', 'NNNNAAACNNNCTTGGNNN', 'NNNNNN']
	trims = ['AAACCTTGG', 'AAACNNNCTTGG', '']
	for seq, trimmed in zip(seqs, trims):
		_seq = Sequence('read1', seq, qualities='#'*len(seq))
		_trimmed = Sequence('read1', trimmed, qualities='#'*len(trimmed))
		assert trimmer(_seq) == _trimmed


def test_quality_trimmer():
	read = Sequence('read1', 'ACGTTTACGTA', '##456789###')

	qt = QualityTrimmer(10, 10, 33)
	assert qt(read) == Sequence('read1', 'GTTTAC', '456789')

	qt = QualityTrimmer(0, 10, 33)
	assert qt(read) == Sequence('read1', 'ACGTTTAC', '##456789')

	qt = QualityTrimmer(10, 0, 33)
	assert qt(read) == Sequence('read1', 'GTTTACGTA', '456789###')


def test_Modifiers_single():
	m = Modifiers()
	m.add_modifier(UnconditionalCutter, lengths=[5])
	assert len(m.modifiers1) == 1
	assert list(m.modifiers1.keys())[0] == UnconditionalCutter
	assert list(m.modifiers1.values())[0].__class__ == UnconditionalCutter
	assert len(m.modifiers2) == 0
	
	# test single-end
	read = Sequence('read1', 'ACGTTTACGTA', '##456789###')
	mod_read = m.modify(read)
	assert mod_read.sequence == 'TACGTA'
	
	# test legacy paired-end
	read2 = Sequence('read1', 'ACGTTTACGTA', '##456789###')
	mod_read, mod_read2 = m.modify(read, read2)
	assert mod_read.sequence == 'TACGTA'
	assert mod_read2.sequence == 'ACGTTTACGTA'


def test_Modifiers_paired_both():
	m = Modifiers()
	m.add_modifier(UnconditionalCutter, read=1|2, lengths=[5])
	assert len(m.modifiers1) == 1
	assert len(m.modifiers2) == 1
	assert list(m.modifiers1.keys())[0] == list(m.modifiers2.keys())[0] == UnconditionalCutter
	assert list(m.modifiers1.values())[0].__class__ == list(m.modifiers2.values())[0].__class__  == UnconditionalCutter
	read = Sequence('read1', 'ACGTTTACGTA', '##456789###')
	read2 = Sequence('read1', 'ACGTTTACGTA', '##456789###')
	mod_read, mod_read2 = m.modify(read, read2)
	assert mod_read.sequence == 'TACGTA'
	assert mod_read2.sequence == 'TACGTA'
