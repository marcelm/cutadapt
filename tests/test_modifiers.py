from typing import List

import pytest

from dnaio import Sequence
from cutadapt.adapters import BackAdapter, PrefixAdapter, IndexedPrefixAdapters, LinkedAdapter, \
    FrontAdapter, Adapter
from cutadapt.modifiers import (UnconditionalCutter, NEndTrimmer, QualityTrimmer,
    Shortener, AdapterCutter, PairedAdapterCutter, ModificationInfo, ZeroCapper)


def test_unconditional_cutter():
    UnconditionalCutter(length=5)
    read = Sequence('r1', 'abcdefg')
    info = ModificationInfo(read)
    assert UnconditionalCutter(length=2)(read, info).sequence == 'cdefg'
    assert UnconditionalCutter(length=-2)(read, info).sequence == 'abcde'
    assert UnconditionalCutter(length=100)(read, info).sequence == ''
    assert UnconditionalCutter(length=-100)(read, info).sequence == ''


def test_zero_capper():
    zc = ZeroCapper()
    read = Sequence("r1", "ACGT", "# !%")
    result = zc(read, ModificationInfo(read))
    assert result.sequence == "ACGT"
    assert result.qualities == "#!!%"


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


def test_adapter_cutter_indexing():
    adapters = [
        PrefixAdapter(sequence, max_errors=1, indels=False)
        for sequence in ["ACGAT", "GGAC", "TTTACTTA", "TAACCGGT", "GTTTACGTA", "CGATA"]
    ]
    ac = AdapterCutter(adapters)
    assert len(ac.adapters) == 1
    assert isinstance(ac.adapters[0], IndexedPrefixAdapters)

    ac = AdapterCutter(adapters, index=False)
    assert len(ac.adapters) == len(adapters)


@pytest.mark.parametrize("action,expected_trimmed1,expected_trimmed2", [
    (None, "CCCCGGTTAACCCC", "TTTTAACCGGTTTT"),
    ("trim", "CCCC", "TTTT"),
    ("lowercase", "CCCCggttaacccc", "TTTTaaccggtttt"),
    ("mask", "CCCCNNNNNNNNNN", "TTTTNNNNNNNNNN"),
    ("retain", "CCCCGGTTAA", "TTTTAACCGG"),
])
def test_paired_adapter_cutter_actions(action, expected_trimmed1, expected_trimmed2):
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


def test_retain_times():
    with pytest.raises(ValueError) as e:
        AdapterCutter([BackAdapter("ACGT")], times=2, action="retain")
    assert "cannot be combined with times" in e.value.args[0]


def test_action_retain():
    back = BackAdapter("AACCGG")
    ac = AdapterCutter([back], action="retain")
    seq = Sequence("r1", "ATTGCCAACCGGTATATAT")
    info = ModificationInfo(seq)
    trimmed = ac(seq, info)
    assert "ATTGCCAACCGG" == trimmed.sequence


@pytest.mark.parametrize("s,expected", [
    ("ATTATTggttaaccAAAAAaaccggTATT", "ggttaaccAAAAAaaccgg"),
    ("AAAAAaaccggTATT", "AAAAAaaccgg"),
    ("ATTATTggttaaccAAAAA", "ggttaaccAAAAA"),
    ("ATTATT", "ATTATT"),
])
def test_linked_action_retain(s, expected):
    front = FrontAdapter("GGTTAACC")
    back = BackAdapter("AACCGG")
    adapters: List[Adapter] = [
        LinkedAdapter(front, back, front_required=False, back_required=False, name="linked")
    ]
    ac = AdapterCutter(adapters, action="retain")
    seq = Sequence("r1", s)
    info = ModificationInfo(seq)
    trimmed = ac(seq, info)
    assert expected == trimmed.sequence
