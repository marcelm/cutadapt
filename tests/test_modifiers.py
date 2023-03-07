from typing import List

import pytest
from dnaio import Sequence
from cutadapt.adapters import (
    BackAdapter,
    PrefixAdapter,
    IndexedPrefixAdapters,
    LinkedAdapter,
    FrontAdapter,
    Adapter,
    RemoveBeforeMatch,
    RemoveAfterMatch,
    LinkedMatch,
)
from cutadapt.modifiers import (
    UnconditionalCutter,
    NEndTrimmer,
    QualityTrimmer,
    Shortener,
    AdapterCutter,
    PairedAdapterCutter,
    ModificationInfo,
    ZeroCapper,
    Renamer,
    ReverseComplementer,
    InvalidTemplate,
    PairedEndRenamer,
)


def test_unconditional_cutter():
    UnconditionalCutter(length=5)
    read = Sequence("r1", "abcdefg")

    info = ModificationInfo(read)
    assert UnconditionalCutter(length=2)(read, info).sequence == "cdefg"
    assert info.cut_prefix == "ab"
    assert info.cut_suffix is None

    info = ModificationInfo(read)
    assert UnconditionalCutter(length=-2)(read, info).sequence == "abcde"
    assert info.cut_suffix == "fg"
    assert info.cut_prefix is None

    assert UnconditionalCutter(length=100)(read, info).sequence == ""
    assert UnconditionalCutter(length=-100)(read, info).sequence == ""


def test_reverse_complementer():
    adapters = [
        PrefixAdapter("TTATTTGTCT"),
        PrefixAdapter("TCCGCACTGG"),
    ]
    adapter_cutter = AdapterCutter(adapters, index=False)
    reverse_complementer = ReverseComplementer(adapter_cutter)

    read = Sequence("r", "ttatttgtctCCAGCTTAGACATATCGCCT")
    info = ModificationInfo(read)
    trimmed = reverse_complementer(read, info)
    assert trimmed.sequence == "CCAGCTTAGACATATCGCCT"
    assert not info.is_rc

    read = Sequence("r", "CAACAGGCCACATTAGACATATCGGATGGTagacaaataa")
    info = ModificationInfo(read)
    trimmed = reverse_complementer(read, info)
    assert trimmed.sequence == "ACCATCCGATATGTCTAATGTGGCCTGTTG"
    assert info.is_rc


def test_zero_capper():
    zc = ZeroCapper()
    read = Sequence("r1", "ACGT", "# !%")
    result = zc(read, ModificationInfo(read))
    assert result.sequence == "ACGT"
    assert result.qualities == "#!!%"


def test_nend_trimmer():
    trimmer = NEndTrimmer()
    seqs = ["NNNNAAACCTTGGNNN", "NNNNAAACNNNCTTGGNNN", "NNNNNN"]
    trims = ["AAACCTTGG", "AAACNNNCTTGG", ""]
    for seq, trimmed in zip(seqs, trims):
        _seq = Sequence("read1", seq, qualities="#" * len(seq))
        _trimmed = Sequence("read1", trimmed, qualities="#" * len(trimmed))
        assert trimmer(_seq, ModificationInfo(_seq)) == _trimmed


def test_quality_trimmer():
    read = Sequence("read1", "ACGTTTACGTA", "##456789###")

    qt = QualityTrimmer(10, 10, 33)
    assert qt(read, ModificationInfo(read)) == Sequence("read1", "GTTTAC", "456789")

    qt = QualityTrimmer(0, 10, 33)
    assert qt(read, ModificationInfo(read)) == Sequence("read1", "ACGTTTAC", "##456789")

    qt = QualityTrimmer(10, 0, 33)
    assert qt(read, ModificationInfo(read)) == Sequence(
        "read1", "GTTTACGTA", "456789###"
    )


def test_shortener():
    read = Sequence("read1", "ACGTTTACGTA", "##456789###")

    shortener = Shortener(0)
    assert shortener(read, ModificationInfo(read)) == Sequence("read1", "", "")

    shortener = Shortener(1)
    assert shortener(read, ModificationInfo(read)) == Sequence("read1", "A", "#")

    shortener = Shortener(5)
    assert shortener(read, ModificationInfo(read)) == Sequence(
        "read1", "ACGTT", "##456"
    )

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


class TestPairedAdapterCutter:
    @pytest.mark.parametrize(
        "action,expected_trimmed1,expected_trimmed2",
        [
            (None, "CCCCGGTTAACCCC", "TTTTAACCGGTTTT"),
            ("trim", "CCCC", "TTTT"),
            ("lowercase", "CCCCggttaacccc", "TTTTaaccggtttt"),
            ("mask", "CCCCNNNNNNNNNN", "TTTTNNNNNNNNNN"),
            ("retain", "CCCCGGTTAA", "TTTTAACCGG"),
        ],
    )
    def test_actions(self, action, expected_trimmed1, expected_trimmed2):
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

    def test_multiple_occurrences(self):
        r1_a1 = BackAdapter("AAAAAA")
        r1_a2 = BackAdapter("CCCC")
        r2_a1 = BackAdapter("GGGG")
        r2_a2 = BackAdapter("TTTT")
        s1 = Sequence("name", "TTAAAAAATTCCCCTT")
        s2 = Sequence("name", "ACACTTTTACAC")
        pac = PairedAdapterCutter([r1_a1, r1_a2], [r2_a1, r2_a2], action="lowercase")
        info1 = ModificationInfo(s1)
        info2 = ModificationInfo(s2)
        trimmed1, trimmed2 = pac(s1, s2, info1, info2)
        assert len(info1.matches) == 1 and info1.matches[0].adapter is r1_a2
        assert len(info2.matches) == 1 and info2.matches[0].adapter is r2_a2
        assert "TTAAAAAATTcccctt" == trimmed1.sequence
        assert "ACACttttacac" == trimmed2.sequence


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


@pytest.mark.parametrize(
    "s,expected",
    [
        ("ATTATTggttaaccAAAAAaaccggTATT", "ggttaaccAAAAAaaccgg"),
        ("AAAAAaaccggTATT", "AAAAAaaccgg"),
        ("ATTATTggttaaccAAAAA", "ggttaaccAAAAA"),
        ("ATTATT", "ATTATT"),
    ],
)
def test_linked_action_retain(s, expected):
    front = FrontAdapter("GGTTAACC")
    back = BackAdapter("AACCGG")
    adapters: List[Adapter] = [
        LinkedAdapter(
            front, back, front_required=False, back_required=False, name="linked"
        )
    ]
    ac = AdapterCutter(adapters, action="retain")
    seq = Sequence("r1", s)
    info = ModificationInfo(seq)
    trimmed = ac(seq, info)
    assert expected == trimmed.sequence


class TestRenamer:
    def test_invalid_template_variable(self):
        with pytest.raises(InvalidTemplate):
            Renamer("{id} {invalid}")

    def test_header_template_variable(self):
        renamer = Renamer("{header} extra")
        read = Sequence("theid thecomment", "ACGT")
        info = ModificationInfo(read)
        assert renamer(read, info).name == "theid thecomment extra"

    def test_id_template_variable(self):
        renamer = Renamer("{id} extra")
        read = Sequence("theid thecomment", "ACGT")
        info = ModificationInfo(read)
        assert renamer(read, info).name == "theid extra"

    def test_tab_escape(self):
        renamer = Renamer(r"{id} extra\tand a tab")
        read = Sequence("theid thecomment", "ACGT")
        info = ModificationInfo(read)
        assert renamer(read, info).name == "theid extra\tand a tab"

    def test_comment_template_variable(self):
        renamer = Renamer("{id}_extra {comment}")
        read = Sequence("theid thecomment", "ACGT")
        info = ModificationInfo(read)
        assert renamer(read, info).name == "theid_extra thecomment"

    def test_comment_template_variable_missing_comment(self):
        renamer = Renamer("{id}_extra {comment}")
        read = Sequence("theid", "ACGT")
        info = ModificationInfo(read)
        assert renamer(read, info).name == "theid_extra "

    def test_cut_prefix_template_variable(self):
        renamer = Renamer("{id}_{cut_prefix} {comment}")
        read = Sequence("theid thecomment", "ACGT")
        info = ModificationInfo(read)
        info.cut_prefix = "TTAAGG"
        assert renamer(read, info).name == "theid_TTAAGG thecomment"

    def test_cut_suffix_template_variable(self):
        renamer = Renamer("{id}_{cut_suffix} {comment}")
        read = Sequence("theid thecomment", "ACGT")
        info = ModificationInfo(read)
        info.cut_suffix = "TTAAGG"
        assert renamer(read, info).name == "theid_TTAAGG thecomment"

    def test_rc_template_variable(self):
        renamer = Renamer("{id} rc={rc} {comment}")
        read = Sequence("theid thecomment", "ACGT")
        info = ModificationInfo(read)
        assert renamer(read, info).name == "theid rc= thecomment"

        read = Sequence("theid thecomment", "ACGT")
        info.is_rc = True
        assert renamer(read, info).name == "theid rc=rc thecomment"

    def test_match_sequence(self):
        sequence = "TTTTCCCCACGTGGGG"
        read = Sequence("theid thecomment", sequence)
        adapter = BackAdapter("AGGT")
        info = ModificationInfo(read)
        info.matches.append(
            RemoveBeforeMatch(
                astart=0,
                astop=4,
                rstart=8,
                rstop=12,
                score=3,
                errors=1,
                adapter=adapter,
                sequence=sequence,
            )
        )
        renamer = Renamer("{header} match={match_sequence}")

        renamer(read, info)

        assert read.name == "theid thecomment match=ACGT"

    def test_match_sequence_linked_match(self):
        sequence = "TATTCCCCACGTGGGG"
        read = Sequence("theid thecomment", sequence)
        adapter1 = PrefixAdapter("TTTT")
        adapter2 = BackAdapter("AGGT")
        linked_adapter = LinkedAdapter(
            adapter1,
            adapter2,
            front_required=True,
            back_required=False,
            name="name",
        )
        info = ModificationInfo(read)
        before_match = RemoveBeforeMatch(
            astart=0,
            astop=4,
            rstart=0,
            rstop=4,
            score=3,
            errors=1,
            adapter=adapter1,
            sequence=sequence,
        )
        after_match = RemoveAfterMatch(
            astart=0,
            astop=4,
            rstart=4,
            rstop=8,
            score=3,
            errors=1,
            adapter=adapter2,
            sequence=sequence[4:],
        )
        info.matches.append(LinkedMatch(before_match, after_match, linked_adapter))
        renamer = Renamer("{header} match={match_sequence}")

        renamer(read, info)

        assert read.name == "theid thecomment match=TATT,ACGT"


class TestPairedEndRenamer:
    def test_invalid_template_variable(self):
        with pytest.raises(InvalidTemplate):
            PairedEndRenamer("{id} {invalid}")

    def test_tab_escape(self):
        renamer = PairedEndRenamer(r"{id} {comment}\tand a tab")
        r1 = Sequence("theid comment1", "ACGT")
        r2 = Sequence("theid comment2", "ACGT")
        info1 = ModificationInfo(r1)
        info2 = ModificationInfo(r2)
        renamed1, renamed2 = renamer(r1, r2, info1, info2)
        assert renamed1.name == "theid comment1\tand a tab"
        assert renamed2.name == "theid comment2\tand a tab"

    def test_ids_not_identical(self):
        renamer = PairedEndRenamer("{id} abc {comment} xyz")
        r1 = Sequence("theid_a cmtx", "ACGT")
        r2 = Sequence("theid_b cmty", "ACGT")
        info1 = ModificationInfo(r1)
        info2 = ModificationInfo(r2)
        with pytest.raises(ValueError) as e:
            renamer(r1, r2, info1, info2)
        assert "not identical" in e.value.args[0]

    def test_comment(self):
        renamer = PairedEndRenamer("{id} abc {comment} xyz")
        r1 = Sequence("theid cmtx", "ACGT")
        r2 = Sequence("theid cmty", "ACGT")
        info1 = ModificationInfo(r1)
        info2 = ModificationInfo(r2)
        renamed1, renamed2 = renamer(r1, r2, info1, info2)
        assert renamed1.name == "theid abc cmtx xyz"
        assert renamed2.name == "theid abc cmty xyz"

    def test_r1_comment(self):
        renamer = PairedEndRenamer("{id} abc {r1.comment} xyz")
        r1 = Sequence("theid cmtx", "ACGT")
        r2 = Sequence("theid cmty", "ACGT")
        info1 = ModificationInfo(r1)
        info2 = ModificationInfo(r2)
        renamed1, renamed2 = renamer(r1, r2, info1, info2)
        assert renamed1.name == "theid abc cmtx xyz"
        assert renamed2.name == "theid abc cmtx xyz"

    def test_r2_comment(self):
        renamer = PairedEndRenamer("{id} abc {r2.comment} xyz")
        r1 = Sequence("theid cmtx", "ACGT")
        r2 = Sequence("theid cmty", "ACGT")
        info1 = ModificationInfo(r1)
        info2 = ModificationInfo(r2)
        renamed1, renamed2 = renamer(r1, r2, info1, info2)
        assert renamed1.name == "theid abc cmty xyz"
        assert renamed2.name == "theid abc cmty xyz"

    def test_read_number(self):
        renamer = PairedEndRenamer("{id} read no. is: {rn}")
        r1 = Sequence("theid cmtx", "ACGT")
        r2 = Sequence("theid cmty", "ACGT")
        info1 = ModificationInfo(r1)
        info2 = ModificationInfo(r2)
        renamed1, renamed2 = renamer(r1, r2, info1, info2)
        assert renamed1.name == "theid read no. is: 1"
        assert renamed2.name == "theid read no. is: 2"

    def test_match_sequence(self):
        r1 = Sequence("theid first", "AACC")
        info1 = ModificationInfo(r1)
        info1.matches.append(
            RemoveBeforeMatch(
                astart=2,
                astop=4,
                rstart=1,
                rstop=3,
                score=1,
                errors=1,
                adapter=FrontAdapter("AT"),
                sequence="AACC",
            )
        )
        r2 = Sequence("theid second", "GGTT")
        info2 = ModificationInfo(r2)
        info2.matches.append(
            RemoveBeforeMatch(
                astart=2,
                astop=4,
                rstart=1,
                rstop=3,
                score=1,
                errors=1,
                adapter=FrontAdapter("GA"),
                sequence="GGTT",
            )
        )
        renamer = PairedEndRenamer("{header} s={match_sequence}")

        renamed1, renamed2 = renamer(r1[:], r2[:], info1, info2)
        assert renamed1.name == "theid first s=AC"
        assert renamed2.name == "theid second s=GT"

        renamer = PairedEndRenamer("{header} s={r1.match_sequence}")
        renamed1, renamed2 = renamer(r1[:], r2[:], info1, info2)
        assert renamed1.name == "theid first s=AC"
        assert renamed2.name == "theid second s=AC"

        renamer = PairedEndRenamer("{header} s={r2.match_sequence}")
        renamed1, renamed2 = renamer(r1[:], r2[:], info1, info2)
        assert renamed1.name == "theid first s=GT"
        assert renamed2.name == "theid second s=GT"
