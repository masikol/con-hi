
from typing import Sequence

import pytest

from src.coverage_array import CoverageArray


@pytest.fixture
def some_array() -> Sequence[int]:
    init_list = [1, 2, 3, 4, 5, 6, 7, 8]
    return CoverageArray(init_list)
# end def


class TestCoverageArray:
    # Class for testing class `src.coverage_array.CoverageArray`

    def test_lengths_equality(self, some_array) -> None:
        # Test that lengths of both arrays are equalsome_array

        expected: int = len(some_array)
        assert expected == len(some_array)
    # end def

    def test_getitem(self, some_array) -> None:
        # Test method __getitem__
        expected: int

        expected = None
        assert expected == some_array[-1]
        expected = 1
        assert expected == some_array[0]
        expected = 5
        assert expected == some_array[4]
        expected = 8
        assert expected == some_array[len(some_array)-1]
        expected = None
        assert expected == some_array[len(some_array)]
    # end def

    def test_min(self, some_array) -> None:
        expected = 1
        obtained = some_array._calc_min_coverage()
        assert expected == obtained
    # end def

    def test_avg(self, some_array) -> None:
        expected = 4.5
        obtained = some_array._calc_avg_coverage()
        assert abs(expected - obtained) < 1e-6
    # end def

    def test_median(self, some_array) -> None:
        expected = 4.5
        obtained = some_array._calc_median_coverage()
        assert abs(expected - obtained) < 1e-6
    # end def

    def test_max(self, some_array) -> None:
        expected = 8
        obtained = some_array._calc_max_coverage()
        assert expected == obtained
    # end def
# end class
