import pytest

from dnaio import Sequence
from cutadapt.modifiers import (UnconditionalCutter, NEndTrimmer, QualityTrimmer,
    Shortener, AdapterCutter, PairedAdapterCutter, ModificationInfo)


def test_unconditional_cutter():
    UnconditionalCutter(length=5)
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
        assert trimmer(_seq, ModificationInfo(_seq)) == _trimmed


def test_quality_trimmer():
    read = Sequence('read1', 'ACGTTTACGTA', '##456789###')

    qt = QualityTrimmer(10, 10, 33)
    assert qt(read, ModificationInfo(read)) == Sequence('read1', 'GTTTAC', '456789')

    qt = QualityTrimmer(0, 10, 33)
    assert qt(read, ModificationInfo(read)) == Sequence('read1', 'ACGTTTAC', '##456789')

    qt = QualityTrimmer(10, 0, 33)
    assert qt(read, ModificationInfo(read)) == Sequence('read1', 'GTTTACGTA', '456789###')


def test_shortener():
    read = Sequence('read1', 'ACGTTTACGTA', '##456789###')

    shortener = Shortener(0)
    assert shortener(read, ModificationInfo(read)) == Sequence('read1', '', '')

    shortener = Shortener(1)
    assert shortener(read, ModificationInfo(read)) == Sequence('read1', 'A', '#')

    shortener = Shortener(5)
    assert shortener(read, ModificationInfo(read)) == Sequence('read1', 'ACGTT', '##456')

    shortener = Shortener(100)
    assert shortener(read, ModificationInfo(read)) == read


def test_adapter_cutter():
    from cutadapt.adapters import BackAdapter
    a1 = BackAdapter("GTAGTCCCGC")
    a2 = BackAdapter("GTAGTCCCCC")
    match = AdapterCutter.best_match([a1, a2], Sequence("name", "ATACCCCTGTAGTCCCC"))
    assert match.adapter is a2


@pytest.mark.parametrize("action,expected_trimmed1,expected_trimmed2", [
    (None, "CCCCGGTTAACCCC", "TTTTAACCGGTTTT"),
    ("trim", "CCCC", "TTTT"),
    ("lowercase", "CCCCggttaacccc", "TTTTaaccggtttt"),
    ("mask", "CCCCNNNNNNNNNN", "TTTTNNNNNNNNNN")
])
def test_paired_adapter_cutter_actions(action, expected_trimmed1, expected_trimmed2):
    from cutadapt.adapters import BackAdapter
    a1 = BackAdapter("GGTTAA")
    a2 = BackAdapter("AACCGG")
    s1 = Sequence("name", "CCCCGGTTAACCCC")
    s2 = Sequence("name", "TTTTAACCGGTTTT")
    pac = PairedAdapterCutter([a1], [a2], action=action)
    info1 = ModificationInfo(s1)
    info2 = ModificationInfo(s2)
    trimmed1, trimmed2 = pac(s1, s2, info1, info2)
    assert expected_trimmed1 == trimmed1.sequence
    assert expected_trimmed2 == trimmed2.sequence
