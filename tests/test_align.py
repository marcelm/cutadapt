from typing import NamedTuple

import pytest

from cutadapt.align import (
    EndSkip,
    Aligner,
    PrefixComparer,
    SuffixComparer,
    hamming_sphere,
    edit_environment,
    edit_distance,
    naive_edit_environment,
    slow_edit_environment,
)
from cutadapt.adapters import Where

from utils import binomial


class AlignmentResult(NamedTuple):
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    score: int
    errors: int


# convenience function (to avoid having to instantiate an Aligner manually)
def locate(
    reference,
    query,
    max_error_rate,
    flags=EndSkip.SEMIGLOBAL,
    wildcard_ref=False,
    wildcard_query=False,
    min_overlap=1,
):
    aligner = Aligner(
        reference,
        max_error_rate,
        flags,
        wildcard_ref,
        wildcard_query,
        min_overlap=min_overlap,
    )
    return aligner.locate(query)


class TestAligner:
    def test(self):
        reference = "CTCCAGCTTAGACATATC"
        aligner = Aligner(reference, 0.1, flags=Where.BACK.value)
        aligner.locate("CC")

    def test_100_percent_error_rate(self):
        reference = "GCTTAGACATATC"
        aligner = Aligner(reference, 1.0, flags=Where.BACK.value)
        aligner.locate("CAA")

    def test_not_only_n_wildcards(self):
        reference = "NNNNN"
        with pytest.raises(ValueError) as info:
            Aligner(reference, 0.1, wildcard_ref=True)
        assert "only N wildcards" in info.value.args[0]

    def test_find_empty_in_empty(self):
        aligner = Aligner("", 0, flags=0, min_overlap=0)
        result = aligner.locate("")
        assert (0, 0, 0, 0, 0, 0) == result

    def test_indels_penalized(self):
        # Alignment in older versions:
        # CCAGTCCTTTCCTGAGAGT
        # CCAGTCCT---CT
        #
        # Should now be:
        # CCAGTCCTTTCCTGAGAGT
        # CCAGTCCTCT
        aligner = Aligner("CCAGTCCTCT", 0.3, flags=Where.PREFIX)
        result = aligner.locate("CCAGTCCTTTCCTGAGAGT")
        assert (0, 10, 0, 10, 9 - 1, 1) == result
        # refstart, refstop, querystart, querystop, score, errors)

        # Alignment:
        # TCGATGC
        # TCGATC
        aligner = Aligner("TCGATC", 1.5 / 6, flags=Where.PREFIX)
        result = aligner.locate("TCGATGC")
        assert (0, 6, 0, 6, 4, 1) == result

    def test_align_illumina(self):
        aligner = Aligner("GCCGAACTTCTTAGACTGCCTTAAGGACGT", 0.1, flags=Where.BACK)
        result = AlignmentResult(
            *aligner.locate("CAAATCACCAGAAGGCGCCTAACTTCTTAGACTGCC")
        )
        #                 GCCGAACTTCTTAGACTGCCTTAAGGACGT (ref)
        #                 |||X||||||||||||||||
        # CAAATCACCAGAAGGCGCCTAACTTCTTAGACTGCC         (query)
        assert result.ref_start == 0
        assert result.ref_end == 20
        assert result.query_start == 16
        assert result.query_end == 36
        assert result.score == 18
        assert result.errors == 1


def test_poly_t():
    aligner = Aligner("TTTT", 0.25, flags=Where.BACK)
    result = AlignmentResult(*aligner.locate("CCTTTT"))
    assert result.ref_start == 0
    assert result.ref_end == 4
    assert result.query_start == 2
    assert result.query_end == 6
    assert result.score == 4
    assert result.errors == 0


def test_poly_t_partial_match():
    aligner = Aligner("TTTTTT", 0.25, flags=Where.BACK)
    result = AlignmentResult(*aligner.locate("CCTTTT"))
    assert result.ref_start == 0
    assert result.ref_end == 4
    assert result.query_start == 2
    assert result.query_end == 6
    assert result.score == 4
    assert result.errors == 0


def test_poly_t_2():
    aligner = Aligner("TTT", 1 / 3, flags=Where.BACK)
    result = AlignmentResult(*aligner.locate("CCTTTT"))
    assert result.ref_start == 0
    assert result.ref_end == 3
    assert result.query_start == 2
    assert result.query_end == 5


def test_poly_a():
    s = "AAAAAAAAAAAAAAAAA"
    t = "ACAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    result = locate(s, t, 0.0, Where.BACK.value)
    # start_s, stop_s, start_t, stop_t, score, cost = result
    assert result == (0, len(s), 4, 4 + len(s), len(s), 0)


# Sequences with IUPAC wildcards
# R=A|G, Y=C|T, S=G|C, W=A|T, K=G|T, M=A|C, B=C|G|T, D=A|G|T, H=A|C|T, V=A|C|G,
# N=A|C|G|T, X={}
WILDCARD_SEQUENCES = [
    "CCCATTGATC",  # original sequence without wildcards
    "CCCRTTRATC",  # R=A|G
    "YCCATYGATC",  # Y=C|T
    "CSSATTSATC",  # S=G|C
    "CCCWWWGATC",  # W=A|T
    "CCCATKKATC",  # K=G|T
    "CCMATTGMTC",  # M=A|C
    "BCCATTBABC",  # B=C|G|T
    "BCCATTBABC",  # B
    "CCCDTTDADC",  # D=A|G|T
    "CHCATHGATC",  # H=A|C|T
    "CVCVTTVATC",  # V=A|C|G
    "CCNATNGATC",  # N=A|C|G|T
    "CCCNTTNATC",  # N
    # 'CCCXTTXATC',  # X
]


def compare_prefixes(ref, query, wildcard_ref=False, wildcard_query=False):
    aligner = PrefixComparer(
        ref,
        max_error_rate=0.9,
        wildcard_ref=wildcard_ref,
        wildcard_query=wildcard_query,
    )
    return aligner.locate(query)


def compare_suffixes(ref, query, wildcard_ref=False, wildcard_query=False):
    aligner = SuffixComparer(
        ref,
        max_error_rate=0.9,
        wildcard_ref=wildcard_ref,
        wildcard_query=wildcard_query,
    )
    return aligner.locate(query)


def test_compare_prefixes():
    assert compare_prefixes("AAXAA", "AAAAATTTTTTTTT") == (0, 5, 0, 5, 3, 1)
    assert compare_prefixes("AANAA", "AACAATTTTTTTTT", wildcard_ref=True) == (
        0,
        5,
        0,
        5,
        5,
        0,
    )
    assert compare_prefixes("AANAA", "AACAATTTTTTTTT", wildcard_ref=True) == (
        0,
        5,
        0,
        5,
        5,
        0,
    )
    assert compare_prefixes("XAAAAA", "AAAAATTTTTTTTT") == (0, 6, 0, 6, 2, 2)

    a = WILDCARD_SEQUENCES[0]
    for s in WILDCARD_SEQUENCES:
        r = s + "GCCAGGGTTGATTCGGCTGATCTGGCCG"
        result = compare_prefixes(a, r, wildcard_query=True)
        assert result == (0, 10, 0, 10, 10, 0), result

        result = compare_prefixes(r, a, wildcard_ref=True)
        assert result == (0, 10, 0, 10, 10, 0)

    for s in WILDCARD_SEQUENCES:
        for t in WILDCARD_SEQUENCES:  # FIXME what is this t doing?
            r = s + "GCCAGGG"
            result = compare_prefixes(
                s,
                r,
            )
            assert result == (0, 10, 0, 10, 10, 0)

            result = compare_prefixes(r, s, wildcard_ref=True, wildcard_query=True)
            assert result == (0, 10, 0, 10, 10, 0)

    r = WILDCARD_SEQUENCES[0] + "GCCAGG"
    for wildc_ref in (False, True):
        for wildc_query in (False, True):
            result = compare_prefixes(
                "CCCXTTXATC", r, wildcard_ref=wildc_ref, wildcard_query=wildc_query
            )
            assert result == (0, 10, 0, 10, 6, 2)


def test_n_wildcard_in_ref_matches_n_wildcard_in_query_prefix():
    # With allowed wildcards in the ref, an N wildcard in the ref should never count as an error,
    # even if matched against an N wildcard in the query while wildcard_query is False
    # Issue #453
    match = compare_prefixes(
        "NNACGT", "NTACGTAA", wildcard_ref=True, wildcard_query=False
    )
    assert match == (0, 6, 0, 6, 6, 0)

    match = compare_prefixes(
        "NNACGT", "YTACGTAA", wildcard_ref=True, wildcard_query=False
    )
    assert match == (0, 6, 0, 6, 6, 0)


def test_n_wildcard_in_ref_matches_n_wildcard_in_query_back():
    aligner = Aligner(
        "NNACGT", max_error_rate=0, wildcard_ref=True, flags=Where.BACK.value
    )
    match = aligner.locate("AAANTACGTAAA")
    assert match == (0, 6, 3, 9, 6, 0)


def test_compare_suffixes():
    assert compare_suffixes("AAXAA", "TTTTTTTAAAAA") == (0, 5, 7, 12, 3, 1)
    assert compare_suffixes("AANAA", "TTTTTTTAACAA", wildcard_ref=True) == (
        0,
        5,
        7,
        12,
        5,
        0,
    )
    assert compare_suffixes("AANAA", "TTTTTTTAACAA", wildcard_ref=True) == (
        0,
        5,
        7,
        12,
        5,
        0,
    )
    assert compare_suffixes("AAAAAX", "TTTTTTTAAAAA") == (0, 6, 6, 12, 2, 2)


@pytest.mark.parametrize("upper", (True, False))
def test_prefix_comparer(upper):
    # only need to test whether None is returned on too many errors, the rest is tested above
    ref = "axcgt"
    if upper:
        ref = ref.upper()
    comparer = PrefixComparer(ref, max_error_rate=0.4)
    repr(comparer)
    assert comparer.locate("TTG") is None
    assert comparer.locate("AGT") is not None
    assert comparer.locate("agt") is not None
    assert comparer.locate("CGT") is None
    assert comparer.locate("TTG") is None


@pytest.mark.parametrize("upper", (True, False))
def test_suffix_comparer(upper):
    # only need to test whether None is returned on too many errors, the rest is tested above
    ref = "axcgt"
    if upper:
        ref = ref.upper()
    comparer = SuffixComparer(ref, max_error_rate=0.4)
    repr(comparer)
    assert comparer.locate("TTG") is None
    assert comparer.locate("AGT") is not None
    assert comparer.locate("agt") is not None
    assert comparer.locate("CGT") is not None
    assert comparer.locate("TTG") is None


@pytest.mark.parametrize("comparer_class", [PrefixComparer, SuffixComparer])
def test_n_wildcards_not_counted_affix(comparer_class):
    # N bases should not contribute to effective adapter length, so only 1 mismatch is allowed
    ref = "CNNNNNNNNGTT"
    assert len(ref) == 12
    comparer = comparer_class(ref, max_error_rate=0.25, wildcard_ref=True)
    assert comparer.locate("CAAAAAAAAGTT") is not None
    assert comparer.locate("CAAAAAAAAGTA") is not None
    assert comparer.locate("CAAAAAAAAGAA") is None  # two mismatches


def test_n_wildcards_not_counted_aligner_back():
    ref = "AGGNNNNNNNNNNNNNNTTC"
    assert len(ref) == 20
    aligner = Aligner(
        ref,
        max_error_rate=0.1,
        wildcard_ref=True,
        flags=Where.BACK.value,
        min_overlap=3,
    )
    assert aligner.effective_length == 6
    assert aligner.locate("TTC") is None
    # adapter start, adapter stop, read start, read stop
    assert aligner.locate("AGG")[:4] == (0, 3, 0, 3)
    assert aligner.locate("AGGCCCCCCC")[:4] == (0, 10, 0, 10)
    assert aligner.locate("ATGCCCCCCC") is None
    assert aligner.locate("AGGCCCCCCCCCCCCCCATC") is None
    assert aligner.locate("CCC" + ref.replace("N", "G") + "AAA") == (
        0,
        20,
        3,
        23,
        20,
        0,
    )


def test_n_wildcards_not_counted_aligner_front():
    ref = "AGGNNNNNNNNNNNNNNTTC"
    assert len(ref) == 20
    aligner = Aligner(
        ref,
        max_error_rate=0.1,
        wildcard_ref=True,
        flags=Where.FRONT.value,
        min_overlap=3,
    )
    assert aligner.effective_length == 6
    # adapter start, adapter stop, read start, read stop
    assert aligner.locate("TTC")[:4] == (17, 20, 0, 3)
    assert aligner.locate("TGC") is None
    assert aligner.locate("CCCCCCCTTC")[:4] == (10, 20, 0, 10)
    assert aligner.locate("CCCCCCCGTC") is None
    assert aligner.locate("CCC" + ref.replace("N", "G") + "AAA") == (
        0,
        20,
        3,
        23,
        20,
        0,
    )


def test_wildcards_in_adapter():
    r = "CATCTGTCC" + WILDCARD_SEQUENCES[0] + "GCCAGGGTTGATTCGGCTGATCTGGCCG"
    for a in WILDCARD_SEQUENCES:
        result = locate(a, r, 0.0, Where.BACK.value, wildcard_ref=True)
        assert result == (0, 10, 9, 19, 10, 0), result

    a = "CCCXTTXATC"
    result = locate(a, r, 0.0, Where.BACK.value, wildcard_ref=True)
    assert result is None


def test_wildcards_in_read():
    a = WILDCARD_SEQUENCES[0]
    for s in WILDCARD_SEQUENCES + ["CCCXTTXATC"]:
        r = "CATCTGTCC" + s + "GCCAGGGTTGATTCGGCTGATCTGGCCG"
        result = locate(a, r, 0.0, Where.BACK.value, wildcard_query=True)
        if "X" in s:
            assert result is None
        else:
            assert result == (0, 10, 9, 19, 10, 0), result


def test_wildcards_in_both():
    for a in WILDCARD_SEQUENCES:
        for s in WILDCARD_SEQUENCES:
            r = "CATCTGTCC" + s + "GCCAGGGTTGATTCGGCTGATCTGGCCG"
            result = locate(
                a, r, 0.0, Where.BACK.value, wildcard_ref=True, wildcard_query=True
            )
            assert result == (0, 10, 9, 19, 10, 0), result


def test_no_match():
    a = locate("CTGATCTGGCCG", "AAAAGGG", 0.1, Where.BACK.value)
    assert a is None, a


def test_hamming_sphere_explicit():
    assert list(hamming_sphere("", 0)) == [""]
    assert list(hamming_sphere("A", 0)) == ["A"]
    assert list(hamming_sphere("A", 1)) == ["C", "G", "T"]
    assert list(hamming_sphere("GTC", 0)) == ["GTC"]
    assert list(hamming_sphere("GTC", 1)) == [
        "ATC",
        "CTC",
        "TTC",
        "GAC",
        "GCC",
        "GGC",
        "GTA",
        "GTG",
        "GTT",
    ]


def hamming_distance(s, t):
    return sum(1 if c != d else 0 for c, d in zip(s, t))


@pytest.mark.parametrize(
    "sk",
    [
        ("", 0),
        ("A", 0),
        ("AAA", 1),
        ("ACC", 2),
        ("TCATTA", 3),
        ("AAAAAAAA", 1),
        ("A" * 15, 2),
    ],
)
def test_hamming_sphere(sk):
    s, k = sk
    result = list(hamming_sphere(s, k))
    result_set = set(result)
    assert len(result) == len(result_set)
    assert len(result) == 3**k * binomial(len(s), k)
    for t in result:
        assert hamming_distance(s, t) == k


@pytest.mark.parametrize(
    "k,s",
    [
        (0, ""),
        (0, "A"),
        (1, "AAA"),
        (1, "TCATTAGA"),
        (2, "ACC"),
        (2, "A" * 10),
        (3, "TCATTA"),
    ],
)
@pytest.mark.parametrize("environment_func", [edit_environment, slow_edit_environment])
def test_edit_environment(k, s, environment_func):
    result = list(environment_func(s, k))
    strings, distances, matches = zip(*result)
    naive = set(naive_edit_environment(s, k))
    assert len(set(strings)) == len(strings)
    assert set(strings) == naive

    error_rate = k / len(s) if s else 0.0
    aligner = Aligner(s, max_error_rate=error_rate, flags=0, min_overlap=len(s))
    for t, dist, m in result:
        result = aligner.locate(t)
        start1, stop1, start2, stop2, score, errors = result
        assert errors == dist
        assert start1 == 0
        assert stop1 == len(s)
        assert start2 == 0
        assert stop2 == len(t)
        assert edit_distance(s, t) == dist
        if environment_func is edit_environment:
            assert m == score
            assert m <= len(s), (s, t, dist)
            assert m <= len(t), (s, t, dist)
