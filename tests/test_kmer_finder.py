import operator
import string

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
    kmer_finder = KmerFinder([(start, stop, [kmer])], ref_wildcards, query_wildcards)
    assert kmer_finder.kmers_present(sequence) is expected


@pytest.mark.parametrize(
    ["ref_table", "query_table", "comp_op", "ref_wildcards", "query_wildcards"],
    [
        (UPPER_TABLE, UPPER_TABLE, operator.eq, False, False),
        (IUPAC_TABLE, ACGT_TABLE, operator.and_, True, False),
        (ACGT_TABLE, IUPAC_TABLE, operator.and_, False, True),
        (IUPAC_TABLE, IUPAC_TABLE, operator.and_, True, True),
    ],
)
def test_kmer_finder_per_char_matching(
    ref_table, query_table, comp_op, ref_wildcards, query_wildcards
):
    for char in string.ascii_letters:
        kmer_finder = KmerFinder(
            [(0, None, [char])],
            ref_wildcards=ref_wildcards,
            query_wildcards=query_wildcards,
        )
        ref_char = ref_table[ord(char)]
        for comp_char in string.ascii_letters:
            query_char = query_table[ord(comp_char)]
            should_match = bool(comp_op(ref_char, query_char))
            if kmer_finder.kmers_present(comp_char) is not should_match:
                raise ValueError(
                    f"{char} should{' ' if should_match else ' not '}match {comp_char}"
                )


def test_kmer_finder_initialize_bigword():
    with pytest.raises(ValueError) as error:
        KmerFinder([(0, None, ["A" * 64])])
    error.match("A" * 64)
    error.match("64")


def test_kmer_finder_initialize_total_greater_than_max():
    kmer_finder = KmerFinder([(0, None, ["A" * 31, "B" * 31, "C" * 31, "D" * 43])])
    assert kmer_finder.kmers_present("X" * 100 + "A" * 31)
    assert kmer_finder.kmers_present("X" * 100 + "B" * 31)
    assert kmer_finder.kmers_present("X" * 100 + "C" * 31)
    assert kmer_finder.kmers_present("X" * 100 + "D" * 43)
    assert not kmer_finder.kmers_present(string.ascii_letters)


def test_kmer_finder_finds_all():
    kmer_finder = KmerFinder([(0, None, ["teenage", "mutant", "ninja", "turtles"])])
    assert kmer_finder.kmers_present("Smells like teenage spirit")
    assert kmer_finder.kmers_present("Everyone with a SNP is technically a mutant.")
    assert kmer_finder.kmers_present("He made a ninja PR that was merged before review")
    assert kmer_finder.kmers_present(
        "Turtles are treated as outgroup, for 'more advanced' reptiles but "
        "molecular evidence suggests they are more close to the dinosaurs than "
        "previously thought."
    )
    assert not kmer_finder.kmers_present(
        "A turtle may be slow, but it also lives for a long time."
    )


def test_kmer_finder_finds_in_region():
    kmer_finder = KmerFinder([(-20, None, ["peace"])])
    # Finding peace, quotes from Mahatma Gandhi
    assert kmer_finder.kmers_present("Each one has to find his peace from within")
    # Peace not found here because outside of the search range.
    assert not kmer_finder.kmers_present(
        "And peace to be real must be unaffected by outside circumstances."
    )
