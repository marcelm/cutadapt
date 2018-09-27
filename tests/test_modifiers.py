from dnaio import Sequence
from cutadapt.modifiers import (UnconditionalCutter, NEndTrimmer, QualityTrimmer,
	Shortener)


def test_unconditional_cutter():
	uc = UnconditionalCutter(length=5)
	s = 'abcdefg'
	assert UnconditionalCutter(length=2)(s, []) == 'cdefg'
	assert UnconditionalCutter(length=-2)(s, []) == 'abcde'
	assert UnconditionalCutter(length=100)(s, []) == ''
	assert UnconditionalCutter(length=-100)(s, []) == ''


def test_nend_trimmer():
	trimmer = NEndTrimmer()
	seqs = ['NNNNAAACCTTGGNNN', 'NNNNAAACNNNCTTGGNNN', 'NNNNNN']
	trims = ['AAACCTTGG', 'AAACNNNCTTGG', '']
	for seq, trimmed in zip(seqs, trims):
		_seq = Sequence('read1', seq, qualities='#'*len(seq))
		_trimmed = Sequence('read1', trimmed, qualities='#'*len(trimmed))
		assert trimmer(_seq, []) == _trimmed


def test_quality_trimmer():
	read = Sequence('read1', 'ACGTTTACGTA', '##456789###')

	qt = QualityTrimmer(10, 10, 33)
	assert qt(read, []) == Sequence('read1', 'GTTTAC', '456789')

	qt = QualityTrimmer(0, 10, 33)
	assert qt(read, []) == Sequence('read1', 'ACGTTTAC', '##456789')

	qt = QualityTrimmer(10, 0, 33)
	assert qt(read, []) == Sequence('read1', 'GTTTACGTA', '456789###')


def test_shortener():
	read = Sequence('read1', 'ACGTTTACGTA', '##456789###')

	shortener = Shortener(0)
	assert shortener(read, []) == Sequence('read1', '', '')

	shortener = Shortener(1)
	assert shortener(read, []) == Sequence('read1', 'A', '#')

	shortener = Shortener(5)
	assert shortener(read, []) == Sequence('read1', 'ACGTT', '##456')

	shortener = Shortener(100)
	assert shortener(read, []) == read
