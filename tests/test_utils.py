import dnaio

from cutadapt.utils import (raise_open_files_limit, Progress, DummyProgress, reverse_complement,
    reverse_complemented_sequence)


def test_raise_open_files_limit():
    try:
        raise_open_files_limit(1)
    except ValueError:
        pass


def test_progress():
    p = Progress()
    p.update(100)
    p.update(1000)
    p.stop(100000)


def test_dummy_progress():
    p = DummyProgress()
    p.update(100)
    p.stop(1000)


def test_reverse_complement():
    rc = reverse_complement
    assert rc("") == ""
    assert rc("A") == "T"
    assert rc("C") == "G"
    assert rc("TG") == "CA"
    assert rc("ATG") == "CAT"
    assert rc("N") == "N"
    assert rc("a") == "t"

    assert rc("ACGTUMRWSYKVHDBN") == "NVHDBMRSWYKAACGT"
    assert rc("acgtumrwsykvhdbn") == "nvhdbmrswykaacgt"
    assert rc("ACGTUMRWSYKVHDBNacgtumrwsykvhdbn") == "nvhdbmrswykaacgtNVHDBMRSWYKAACGT"


def test_reverse_complemented_sequence():
    s = dnaio.Sequence("the_name", "ACGTTTGA", "B>%%BB5#")
    assert reverse_complemented_sequence(s) == dnaio.Sequence("the_name", "TCAAACGT", "#5BB%%>B")

    s = dnaio.Sequence("the_name", "ACGTTTGA")
    assert reverse_complemented_sequence(s) == dnaio.Sequence("the_name", "TCAAACGT")
