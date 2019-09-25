import pytest

from utils import assert_files_equal, FilesDifferent, binomial


def test_files_different():
    with pytest.raises(FilesDifferent):
        assert_files_equal("simple.fasta", "simple.fastq")


def test_binomial():
    assert binomial(0, 0) == 1
    assert binomial(0, 1) == 0
    assert binomial(0, -1) == 0
    assert binomial(1, 0) == 1
    assert binomial(1, 1) == 1
    assert binomial(1, 2) == 0
    assert binomial(10, 5) == 10 * 9 * 8 * 7 * 6 // (2 * 3 * 4 * 5)
