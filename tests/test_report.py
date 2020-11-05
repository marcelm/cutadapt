from cutadapt.report import (
    safe_divide,
    add_if_not_none,
)


def test_safe_divide():
    assert safe_divide(1, 0) == 0
    assert safe_divide(5, 2) == 2.5


def test_add_if_not_none():
    assert add_if_not_none(3, 5) == 8
    assert add_if_not_none(3, None) == 3
    assert add_if_not_none(None, 5) == 5
