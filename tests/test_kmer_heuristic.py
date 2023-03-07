import pytest

from cutadapt.kmer_heuristic import (
    kmer_possibilities,
    minimize_kmer_search_list,
    create_back_overlap_searchsets,
    create_positions_and_kmers,
)


@pytest.mark.parametrize(
    ["sequence", "chunks", "expected"],
    [
        ("ABC", 3, [{"A", "B", "C"}]),
        ("ABCD", 3, [{"A", "B", "CD"}, {"A", "BC", "D"}, {"AB", "C", "D"}]),
    ],
)
def test_kmer_possibilities(sequence, chunks, expected):
    frozen_expected = set(frozenset(s) for s in expected)
    result = kmer_possibilities(sequence, chunks)
    frozen_result = set(frozenset(s) for s in result)
    assert frozen_expected == frozen_result


@pytest.mark.parametrize(
    ["kmer_search_list", "expected"],
    [
        ([("ABC", -33, None), ("ABC", -19, None)], [("ABC", -33, None)]),
        (
            [("ABC", -33, None), ("ABC", -19, None), ("ABC", 0, None)],
            [("ABC", 0, None)],
        ),
        ([("ABC", 0, 10), ("ABC", 0, 20)], [("ABC", 0, 20)]),
        ([("ABC", 0, 10), ("ABC", 0, 20), ("ABC", 0, None)], [("ABC", 0, None)]),
        ([("ABC", 0, 10), ("ABC", -19, None), ("ABC", 0, None)], [("ABC", 0, None)]),
        ([("ABC", 0, 10), ("ABC", -19, None)], [("ABC", 0, 10), ("ABC", -19, None)]),
    ],
)
def test_minimize_kmer_search_list(kmer_search_list, expected):
    result = minimize_kmer_search_list(kmer_search_list)
    assert set(result) == set(expected)


def test_create_back_overlap_searchsets():
    adapter = "ABCDEFGHIJ0123456789"
    searchsets = create_back_overlap_searchsets(adapter, 3, 0.1)
    assert len(searchsets) == 5
    assert (-3, None, [{"ABC"}]) in searchsets
    assert (-4, None, [{"ABCD"}]) in searchsets
    assert (-9, None, [{"ABCDE"}]) in searchsets
    assert (-19, None, kmer_possibilities(adapter[:10], 2)) in searchsets
    assert (-20, None, kmer_possibilities(adapter, 3)) in searchsets


@pytest.mark.parametrize(
    ["kwargs", "expected"],
    [
        (
            dict(back_adapter=True, front_adapter=False, internal=True, min_overlap=3),
            [
                (-3, None, ["ABC"]),
                (-4, None, ["ABCD"]),
                (-19, None, ["ABCDE", "FGHIJ"]),
                (0, None, ["ABCDEFG", "HIJ0123", "456789"]),
            ],
        ),
        (
            dict(back_adapter=True, front_adapter=False, internal=False, min_overlap=3),
            [
                (-3, None, ["ABC"]),
                (-4, None, ["ABCD"]),
                (-19, None, ["ABCDE", "FGHIJ"]),
                (-20, None, ["ABCDEFG", "HIJ0123", "456789"]),
            ],
        ),
        (
            dict(back_adapter=False, front_adapter=True, internal=False, min_overlap=3),
            [
                (0, 3, ["789"]),
                (0, 4, ["6789"]),
                (0, 19, ["01234", "56789"]),
                (0, 20, ["ABCDEF", "GHIJ012", "3456789"]),
            ],
        ),
        (
            dict(back_adapter=True, front_adapter=False, internal=True, min_overlap=20),
            [
                (0, None, ["ABCDEFG", "HIJ0123", "456789"]),
            ],
        ),
        (
            dict(back_adapter=False, front_adapter=False, internal=True, min_overlap=3),
            [
                (0, None, ["ABCDEFG", "HIJ0123", "456789"]),
            ],
        ),
    ],
)
def test_create_kmers_and_positions(kwargs, expected):
    adapter = "ABCDEFGHIJ0123456789"
    result = create_positions_and_kmers(
        adapter,
        error_rate=0.1,
        **kwargs,
    )
    assert {(start, stop): frozenset(kmers) for start, stop, kmers in result} == {
        (start, stop): frozenset(kmers) for start, stop, kmers in expected
    }
