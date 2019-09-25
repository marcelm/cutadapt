import pytest

from cutadapt.__main__ import main


def test_help():
    with pytest.raises(SystemExit) as e:
        main(["--help"])
    assert e.value.args[0] == 0
