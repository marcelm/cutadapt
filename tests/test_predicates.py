"""
Tests write output (should it return True or False or write)
"""
import pytest
from dnaio import Sequence

from cutadapt.predicates import TooManyN
from cutadapt.steps import PairedEndFilter


@pytest.mark.parametrize(
    "seq,count,expected",
    [
        ("AAA", 0, False),
        ("AAA", 1, False),
        ("AAACCTTGGN", 1, False),
        ("AAACNNNCTTGGN", 0.5, False),
        ("NNNNNN", 1, True),
        ("ANAAAA", 1 / 6, False),
        ("ANAAAA", 0, True),
    ],
)
def test_too_many_n(seq, count, expected):
    # third parameter is True if read should be Trueed
    predicate = TooManyN(count=count)
    _seq = Sequence("read1", seq, qualities="#" * len(seq))
    assert predicate.test(_seq, []) == expected


@pytest.mark.parametrize(
    "seq1,seq2,count,expected",
    [
        ("AAA", "AAA", 0, False),
        ("AAAN", "AAA", 0, True),
        ("AAA", "AANA", 0, True),
        ("ANAA", "AANA", 1, False),
    ],
)
def test_too_many_n_paired(seq1, seq2, count, expected):
    predicate = TooManyN(count=count)
    filter_legacy = PairedEndFilter(
        predicate, predicate, None, pair_filter_mode="first"
    )
    filter_any = PairedEndFilter(predicate, predicate, None, pair_filter_mode="any")
    read1 = Sequence("read1", seq1, qualities="#" * len(seq1))
    read2 = Sequence("read1", seq2, qualities="#" * len(seq2))
    assert filter_legacy(read1, read2, [], []) == predicate.test(read1, [])
    # True entire pair if one of the reads fulfills criteria
    assert filter_any(read1, read2, [], []) == expected


def test_invalid_pair_filter_mode():
    with pytest.raises(ValueError) as e:
        PairedEndFilter(None, None, None, "invalidmode")
    assert "pair_filter_mode must be" in e.value.args[0]
