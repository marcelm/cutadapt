from cutadapt.qc import base_weighted_categories, equidistant_ranges

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


@pytest.mark.parametrize(
    ["base_counts", "parts", "expected"],
    (
        ([10, 10, 4, 3, 2, 1], 3, [(0, 1), (1, 2), (2, 6)]),
        ([8, 8, 5, 3, 2, 2], 3, [(0, 2), (2, 3), (3, 6)])
    )
)
def test_base_weighted_categories(base_counts, parts, expected):
    assert list(base_weighted_categories(base_counts, parts)) == expected
