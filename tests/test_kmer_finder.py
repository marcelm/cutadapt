import pytest

from cutadapt.adapters import KmerFinder


KMER_FINDER_TESTS = [
    # kmer, start, stop, ref_wildcards, query_wildcards, sequence, expected
    ("ACGT", 0, None, False, False, "ACGTACG", True),
    ("ACGT", 0, None, False, False, "ACgtACG", True),
    ("acgt", 0, None, False, False, "ACgtACG", True),
    ("ACGT", 0, None, False, False, "acgtacg", True),
    ("ACGT", 0, None, False, False, "gacgact", False),
    ("ACGT", 0, None, False, True, "ACGNACG", True),
    ("ACGT", 0, None, False, False, "ACGNACG", False),
    ("ACGN", 0, None, True, False, "ACGTACG", True),
    ("ACGN", 0, None, True, False, "ACGxACG", True),
    ("ACKN", 0, None, True, False, "ACGTACG", True),
    ("ACKN", 0, None, True, True, "ACWRACG", True),
    ("ACKN", 0, None, True, True, "ACWxACG", False),
]


@pytest.mark.parametrize(
    [
        "kmer",
        "start",
        "stop",
        "ref_wildcards",
        "query_wildcards",
        "sequence",
        "expected",
    ],
    KMER_FINDER_TESTS,
)
def test_kmer_finder(
    kmer: str,
    start: int,
    stop: int,
    ref_wildcards: bool,
    query_wildcards: bool,
    sequence: str,
    expected: bool,
):
    kmer_finder = KmerFinder([(kmer, start, stop)], ref_wildcards, query_wildcards)
    assert kmer_finder.kmers_present(sequence) is expected
