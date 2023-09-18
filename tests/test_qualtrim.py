import pytest

from dnaio import SequenceRecord
from cutadapt.qualtrim import nextseq_trim_index, expected_errors, poly_a_trim_index


def test_nextseq_trim():
    s = SequenceRecord("n", "", "")
    assert nextseq_trim_index(s, cutoff=22) == 0
    s = SequenceRecord(
        "n",
        "TCTCGTATGCCGTCTTATGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGNNNNNNNNNNNGGNGG",
        "AA//EAEE//A6///E//A//EA/EEEEEEAEA//EEEEEEEEEEEEEEE###########EE#EA",
    )
    assert nextseq_trim_index(s, cutoff=22) == 33


@pytest.mark.parametrize(
    "sequence,tail",
    [
        ("", ""),
        ("GGGGGGGGAAAGAAGAAGAAGAAGAAGAAG", ""),
        ("TTTAGA", ""),  # shorter than three nucleotides
        ("TTTAGAA", ""),  # shorter than three nucleotides
        ("TTTAG", "AAA"),
        ("TCAAGAAGTCCTTTACCAGCTTTC", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        ("TCAAGAAGTCCTTTACCAGCTTTC", "AAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        ("GCAGATCACCTT", "AAAAAAAAAAAAAAAAAAAAAAAAAAAATAAA"),
        ("GCAGATCACCTT", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAT"),
        ("GCAGATCACCTT", "AAAAAAAAAAAAAAAAAAAAAAAAAAAATCG"),
        ("GCAGATCACCTAT", "AAAACAAAAAAACAAAAAAAACAAAAAA"),
        ("TTTT", "AAATAAAA"),
        ("GGGGGGGGAAAGAAGAAGAAGAAGAAGAAG", "AAA"),
    ],
)
def test_poly_a_trim_index(sequence, tail):
    assert poly_a_trim_index(sequence + tail) == len(sequence)


@pytest.mark.parametrize(
    "head,sequence",
    [
        ("", ""),
        ("", "GGGGGGGGAAAGAAGAAGAAGAAGAAGAAG"),
        ("", "TGTCCC"),
        ("", "TTGTCCC"),
        ("TTT", "GTCCC"),
        (
            "TTTTTTTTTTTTTTTTTTTTT",
            "CAAGAAGTCCCCAGCTTTC",
        ),
        ("TTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "CAAGAAGTCCTTTACCAGCTTTC"),
        ("TTTTTATTTTTTTTTTTTTTTTTTTTTTTTTT", "GCAGATCACCTT"),
        ("ATTTTTTTTTTTTTTTTTTTTTTTTTTTT", "GCAGATCACCTT"),
        ("AGCTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "GCAGATCACCTT"),
        ("TTTTGTTTTTTTGTTTTTTTTGTTTTTT", "GCAGATCACCTAT"),
        ("TTTATTTT", "AAAA"),
        ("TTT", "GGGGGGGGAAAGAAGAAGAAGAAGAAGAAG"),
    ],
)
def test_poly_t_trim_index(head, sequence):
    assert poly_a_trim_index(head + sequence, revcomp=True) == len(head)


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
