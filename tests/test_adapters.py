import pytest

from dnaio import Sequence
from cutadapt.adapters import (
    RemoveAfterMatch,
    RemoveBeforeMatch,
    FrontAdapter,
    BackAdapter,
    PrefixAdapter,
    SuffixAdapter,
    LinkedAdapter,
    MultipleAdapters,
    IndexedPrefixAdapters,
    IndexedSuffixAdapters,
)


def test_back_adapter_absolute_number_of_errors():
    adapter = BackAdapter(
        sequence="GATCGGAAGA",
        max_errors=1,
        min_overlap=3,
    )
    assert adapter.max_error_rate == 1 / 10


def test_front_adapter_partial_occurrence_in_back():
    adapter = FrontAdapter("CTGAATT", max_errors=0, min_overlap=4)
    assert adapter.match_to("GGGGGCTGAA") is None


def test_back_adapter_partial_occurrence_in_front():
    adapter = BackAdapter("CTGAATT", max_errors=0, min_overlap=4)
    assert adapter.match_to("AATTGGGGGGG") is None


def test_issue_52():
    adapter = BackAdapter(
        sequence='GAACTCCAGTCACNNNNN',
        max_errors=0.12,
        min_overlap=5,
        read_wildcards=False,
        adapter_wildcards=True)
    sequence = "CCCCAGAACTACAGTCCCGGC"
    am = RemoveAfterMatch(astart=0, astop=17, rstart=5, rstop=21, matches=15, errors=2,
        adapter=adapter, sequence=sequence)
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

    adapter = BackAdapter(
        sequence="TCGTATGCCGTCTTC",
        max_errors=0.2,
        min_overlap=3,
        read_wildcards=False,
        adapter_wildcards=False)
    result = adapter.match_to("TCGTATGCCCTCC")
    assert result.errors == 3, result
    assert result.astart == 0, result
    assert result.astop == 15, result


@pytest.mark.xfail(strict=True)
def test_back_adapter_indel_and_exact_occurrence():
    adapter = BackAdapter(
        sequence="GATCGGAAGA",
        max_errors=0.1,
        min_overlap=3,
    )
    match = adapter.match_to("GATCGTGAAGAGATCGGAAGA")
    # We want the leftmost match of these two possible ones:
    # GATCGTGAAGAGATCGGAAGA
    # GATCG-GAAGA
    #            GATCGGAAGA
    assert match.errors == 0
    assert match.matches == 10
    assert match.astart == 0
    assert match.astop == 10
    assert match.rstart == 0
    assert match.rstop == 10


def test_back_adapter_indel_and_mismatch_occurrence():
    adapter = BackAdapter(
        sequence="GATCGGAAGA",
        max_errors=0.1,
        min_overlap=3,
    )
    match = adapter.match_to("CTGGATCGGAGAGCCGTAGATCGGGAGAGGC")
    # CTGGATCGGA-GAGCCGTAGATCGGGAGAGGC
    #    ||||||| ||      ||||||X|||
    #    GATCGGAAGA      GATCGGAAGA
    assert match.errors == 1
    assert match.matches == 9
    assert match.astart == 0
    assert match.astop == 10
    assert match.rstart == 3
    assert match.rstop == 12


def test_str():
    a = BackAdapter('ACGT', max_errors=0.1)
    str(a)
    str(a.match_to("TTACGT"))


def test_prefix_with_indels_one_mismatch():
    a = PrefixAdapter(
        sequence="GCACATCT",
        max_errors=0.15,
        min_overlap=1,
        read_wildcards=False,
        adapter_wildcards=False,
        indels=True,
    )
    result = a.match_to("GCACATCGGAA")
    assert result.errors == 1
    assert result.matches == 7
    assert result.astart == 0
    assert result.astop == 8
    assert result.rstart == 0
    assert result.rstop == 8


def test_prefix_with_indels_two_mismatches():
    a = PrefixAdapter(
        sequence="GCACATTT",
        max_errors=0.3,
        min_overlap=1,
        read_wildcards=False,
        adapter_wildcards=False,
        indels=True,
    )
    result = a.match_to("GCACATCGGAA")
    assert result.errors == 2
    assert result.matches == 6
    assert result.astart == 0
    assert result.astop == 8
    assert result.rstart == 0
    assert result.rstop == 8


def test_linked_adapter():
    front_adapter = PrefixAdapter('AAAA', min_overlap=4)
    back_adapter = BackAdapter('TTTT', min_overlap=3)

    linked_adapter = LinkedAdapter(
        front_adapter, back_adapter, front_required=True, back_required=False, name='name')
    assert linked_adapter.front_adapter.min_overlap == 4
    assert linked_adapter.back_adapter.min_overlap == 3

    read = Sequence(name='seq', sequence='AAAACCCCCTTTT')
    trimmed = linked_adapter.match_to(read.sequence).trimmed(read)
    assert trimmed.name == 'seq'
    assert trimmed.sequence == 'CCCCC'


def test_info_record():
    adapter = BackAdapter(
        sequence='GAACTCCAGTCACNNNNN',
        max_errors=0.12,
        min_overlap=5,
        read_wildcards=False,
        adapter_wildcards=True,
        name="Foo")
    read = Sequence(name="abc", sequence='CCCCAGAACTACAGTCCCGGC')
    am = RemoveAfterMatch(astart=0, astop=17, rstart=5, rstop=21, matches=15, errors=2,
        adapter=adapter, sequence=read.sequence)
    assert am.get_info_records(read) == [[
        "",
        2,
        5,
        21,
        'CCCCA',
        'GAACTACAGTCCCGGC',
        '',
        'Foo',
        '',
        '',
        '',
    ]]


def test_random_match_probabilities():
    a = BackAdapter('A', max_errors=0.1).create_statistics()
    assert a.back.random_match_probabilities(0.5) == [1, 0.25]
    assert a.back.random_match_probabilities(0.2) == [1, 0.4]

    for s in ('ACTG', 'XMWH'):
        a = BackAdapter(s, max_errors=0.1).create_statistics()
        assert a.back.random_match_probabilities(0.5) == [1, 0.25, 0.25**2, 0.25**3, 0.25**4]
        assert a.back.random_match_probabilities(0.2) == [1, 0.4, 0.4*0.1, 0.4*0.1*0.4, 0.4*0.1*0.4*0.1]

    a = FrontAdapter('GTCA', max_errors=0.1).create_statistics()
    assert a.front.random_match_probabilities(0.5) == [1, 0.25, 0.25**2, 0.25**3, 0.25**4]
    assert a.front.random_match_probabilities(0.2) == [1, 0.4, 0.4*0.1, 0.4*0.1*0.4, 0.4*0.1*0.4*0.1]


def test_add_adapter_statistics():
    stats = BackAdapter('A', name='name', max_errors=0.1).create_statistics()
    end_stats = stats.back
    end_stats.adjacent_bases['A'] = 7
    end_stats.adjacent_bases['C'] = 19
    end_stats.adjacent_bases['G'] = 23
    end_stats.adjacent_bases['T'] = 42
    end_stats.adjacent_bases[''] = 45

    end_stats.errors[10][0] = 100
    end_stats.errors[10][1] = 11
    end_stats.errors[10][2] = 3
    end_stats.errors[20][0] = 600
    end_stats.errors[20][1] = 66
    end_stats.errors[20][2] = 6

    stats2 = BackAdapter('A', name='name', max_errors=0.1).create_statistics()
    end_stats2 = stats2.back
    end_stats2.adjacent_bases['A'] = 43
    end_stats2.adjacent_bases['C'] = 31
    end_stats2.adjacent_bases['G'] = 27
    end_stats2.adjacent_bases['T'] = 8
    end_stats2.adjacent_bases[''] = 5
    end_stats2.errors[10][0] = 234
    end_stats2.errors[10][1] = 14
    end_stats2.errors[10][3] = 5
    end_stats2.errors[15][0] = 90
    end_stats2.errors[15][1] = 17
    end_stats2.errors[15][2] = 2

    stats += stats2
    r = stats.back

    assert r.adjacent_bases == {'A': 50, 'C': 50, 'G': 50, 'T': 50, '': 50}
    assert r.errors == {
        10: {0: 334, 1: 25, 2: 3, 3: 5},
        15: {0: 90, 1: 17, 2: 2},
        20: {0: 600, 1: 66, 2: 6},
    }


def test_linked_matches_property():
    """Accessing matches property of non-anchored linked adapters"""
    # Issue #265
    front_adapter = FrontAdapter("GGG")
    back_adapter = BackAdapter("TTT")
    la = LinkedAdapter(front_adapter, back_adapter, front_required=False, back_required=False, name='name')
    assert la.match_to("AAAATTTT").matches == 3


@pytest.mark.parametrize("adapter_class", [PrefixAdapter, SuffixAdapter])
def test_no_indels_empty_read(adapter_class):
    # Issue #376
    adapter = adapter_class("ACGT", indels=False)
    adapter.match_to("")


def test_prefix_match_with_n_wildcard_in_read():
    adapter = PrefixAdapter("NNNACGT", indels=False)
    match = adapter.match_to("TTTACGTAAAA")
    assert match is not None and (0, 7) == (match.rstart, match.rstop)
    match = adapter.match_to("NTTACGTAAAA")
    assert match is not None and (0, 7) == (match.rstart, match.rstop)


def test_suffix_match_with_n_wildcard_in_read():
    adapter = SuffixAdapter("ACGTNNN", indels=False)
    match = adapter.match_to("TTTTACGTTTT")
    assert match is not None and (4, 11) == (match.rstart, match.rstop)
    match = adapter.match_to("TTTTACGTCNC")
    assert match is not None and (4, 11) == (match.rstart, match.rstop)


def test_multiple_adapters():
    a1 = BackAdapter("GTAGTCCCGC")
    a2 = BackAdapter("GTAGTCCCCC")
    ma = MultipleAdapters([a1, a2])
    match = ma.match_to("ATACCCCTGTAGTCCCC")
    assert match.adapter is a2


def test_indexed_prefix_adapters():
    adapters = [
        PrefixAdapter("GAAC", indels=False),
        PrefixAdapter("TGCT", indels=False),
    ]
    ma = IndexedPrefixAdapters(adapters)
    match = ma.match_to("GAACTT")
    assert match.adapter is adapters[0]
    match = ma.match_to("TGCTAA")
    assert match.adapter is adapters[1]
    assert ma.match_to("GGGGGGG") is None


def test_indexed_prefix_adapters_incorrect_type():
    with pytest.raises(ValueError):
        IndexedPrefixAdapters([
            PrefixAdapter("GAAC", indels=False),
            SuffixAdapter("TGCT", indels=False),
        ])


def test_indexed_very_similar(caplog):
    IndexedPrefixAdapters([
        PrefixAdapter("GAAC", max_errors=1, indels=False),
        PrefixAdapter("GAAG", max_errors=1, indels=False),
    ])
    assert "cannot be assigned uniquely" in caplog.text


def test_indexed_too_high_k():
    with pytest.raises(ValueError) as e:
        IndexedPrefixAdapters([
            PrefixAdapter("ACGTACGT", max_errors=3, indels=False),
            PrefixAdapter("AAGGTTCC", max_errors=2, indels=False),
        ])
    assert "Error rate too high" in e.value.args[0]


def test_indexed_suffix_adapters():
    adapters = [
        SuffixAdapter("GAAC", indels=False),
        SuffixAdapter("TGCT", indels=False),
    ]
    ma = IndexedSuffixAdapters(adapters)
    match = ma.match_to("TTGAAC")
    assert match.adapter is adapters[0]
    match = ma.match_to("AATGCT")
    assert match.adapter is adapters[1]


def test_indexed_suffix_adapters_incorrect_type():
    with pytest.raises(ValueError):
        IndexedSuffixAdapters([
            SuffixAdapter("GAAC", indels=False),
            PrefixAdapter("TGCT", indels=False),
        ])


def test_multi_prefix_adapter_with_indels():
    adapters = [
        PrefixAdapter("GTAC", max_errors=1, indels=True),
        PrefixAdapter("TGCT", max_errors=1, indels=True),
    ]
    ma = IndexedPrefixAdapters(adapters)
    match = ma.match_to("GATACGGG")
    assert match.adapter is adapters[0]
    match = ma.match_to("TAGCTAA")
    assert match.adapter is adapters[1]


def test_indexed_prefix_adapters_with_n_wildcard():
    sequence = "GGTCCAGA"
    ma = IndexedPrefixAdapters([PrefixAdapter(sequence, max_errors=1, indels=False)])
    for i in range(len(sequence)):
        t = sequence[:i] + "N" + sequence[i+1:] + "TGCT"
        result = ma.match_to(t)
        assert isinstance(result, RemoveBeforeMatch)
        assert (result.rstart, result.rstop) == (0, 8)
        assert result.errors == 1
        assert result.matches == 7
