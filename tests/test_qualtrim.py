import pytest

from dnaio import Sequence
from cutadapt.qualtrim import nextseq_trim_index, expected_errors


def test_nextseq_trim():
    s = Sequence("n", "", "")
    assert nextseq_trim_index(s, cutoff=22) == 0
    s = Sequence(
        "n",
        "TCTCGTATGCCGTCTTATGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGNNNNNNNNNNNGGNGG",
        "AA//EAEE//A6///E//A//EA/EEEEEEAEA//EEEEEEEEEEEEEEE###########EE#EA",
    )
    assert nextseq_trim_index(s, cutoff=22) == 33


def test_expected_errors():
    def encode_qualities(quals):
        return "".join(chr(q + 33) for q in quals)

    assert pytest.approx(0.0) == expected_errors("")

    assert pytest.approx(0.1) == expected_errors(encode_qualities([10]))
    assert pytest.approx(0.01) == expected_errors(encode_qualities([20]))
    assert pytest.approx(0.001) == expected_errors(encode_qualities([30]))

    assert pytest.approx(0.2) == expected_errors(encode_qualities([10, 10]))
    assert pytest.approx(0.11) == expected_errors(encode_qualities([10, 20]))
    assert pytest.approx(0.11) == expected_errors(encode_qualities([20, 10]))

    assert pytest.approx(0.3) == expected_errors(encode_qualities([10, 10, 10]))
    assert pytest.approx(0.111) == expected_errors(encode_qualities([10, 20, 30]))
    assert pytest.approx(0.2111) == expected_errors(
        encode_qualities([10, 10, 20, 30, 40])
    )
