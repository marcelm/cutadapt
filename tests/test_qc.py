from cutadapt.qc import equidistant_ranges

import pytest


@pytest.mark.parametrize(
    ["length", "parts", "expected"],
    (
        (10, 3, ((0, 3), (3, 6), (6, 10))),
        (5, 6, ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5))),
        (6, 3, ((0, 2), (2, 4), (4, 6))),
    ),
)
def test_equidistant_ranges(length, parts, expected):
    assert tuple(equidistant_ranges(length, parts)) == expected
