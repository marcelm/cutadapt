import operator

import pytest

from cutadapt._align import _upper_table, _acgt_table, _iupac_table
from cutadapt.adapters import KmerFinder

UPPER_TABLE: bytes = _upper_table()
ACGT_TABLE: bytes = _acgt_table()
IUPAC_TABLE: bytes = _iupac_table()

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


@pytest.mark.parametrize(
    ["ref_table", "query_table", "comp_op", "ref_wildcards", "query_wildcards"],
    [
        (UPPER_TABLE, UPPER_TABLE, operator.eq, False, False),
        (IUPAC_TABLE, ACGT_TABLE, operator.or_, True, False),
        (ACGT_TABLE, IUPAC_TABLE, operator.or_, False, True),
        (IUPAC_TABLE, IUPAC_TABLE, operator.or_, True, True),
    ],
)
def test_kmer_finder_per_char_matching(
    ref_table, query_table, comp_op, ref_wildcards, query_wildcards
):
    iupac_letters = "ACGTURYSWKMBDHVN"
    iupac_letters = iupac_letters + iupac_letters.lower()

    for char in iupac_letters:
        kmer_finder = KmerFinder(
            [(char, 0, None)],
            ref_wildcards=ref_wildcards,
            query_wildcards=query_wildcards,
        )
        for comp_char in iupac_letters:
            should_match = comp_op(ref_table[ord(char)], query_table[ord(comp_char)])
            if kmer_finder.kmers_present(comp_char) is not should_match:
                raise ValueError(f"{char} should match {comp_char}")
