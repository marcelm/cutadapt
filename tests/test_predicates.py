"""
Tests write output (should it return True or False or write)
"""
import pytest
from dnaio import SequenceRecord

from cutadapt.predicates import TooManyN, TooHighAverageErrorRate
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
    predicate = TooManyN(count=count)
    _seq = SequenceRecord("read1", seq, qualities="#" * len(seq))
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
    read1 = SequenceRecord("read1", seq1, qualities="#" * len(seq1))
    read2 = SequenceRecord("read1", seq2, qualities="#" * len(seq2))
    assert (filter_legacy(read1, read2, [], []) is None) == predicate.test(read1, [])
    # True entire pair if one of the reads fulfills criteria
    assert (filter_any(read1, read2, [], []) is None) == expected


def test_invalid_pair_filter_mode():
    with pytest.raises(ValueError) as e:
        PairedEndFilter(None, None, None, "invalidmode")
    assert "pair_filter_mode must be" in e.value.args[0]


@pytest.mark.parametrize(
    "quals,rate,expected",
    [
        # 3 * 0.1 is larger than 0.3 due to floating point rounding.
        (chr(43) * 3, 0.1, True),
        (chr(43) * 3 + chr(33), 0.1, True),  # 3 * 0.1 + 1
        (chr(43) * 3 + chr(33), 0.33, False),  # 3 * 0.1 + 1
        (chr(43) * 3 + chr(33), 0.32, True),  # 3 * 0.1 + 1
        (chr(126) * 9 + chr(33), 0.1, True),  # 9 * 10^-9.3 + 1
    ],
)
def test_too_high_average_error_rate(quals, rate, expected):
    predicate = TooHighAverageErrorRate(rate)
    _seq = SequenceRecord("read1", "A" * len(quals), qualities=quals)
    assert predicate.test(_seq, []) == expected
