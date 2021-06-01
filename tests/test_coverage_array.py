# -*- encoding: utf-8 -*-

from typing import Sequence

import pytest

from src.coverage_array import CoverageArray


@pytest.fixture
def some_int_array() -> Sequence[int]:
    return [1, 2, 3, 4, 5, 6, 7, 8]
# end def some_int_array


class TestCoverageArray:

    def test_lengths_equality(self, some_int_array) -> None:
        # Test that lengths of both arrays are equal
        cov_aray: CoverageArray = CoverageArray(some_int_array)

        expected: int = len(some_int_array)
        assert expected == len(cov_aray)
    # end def test_lengths_equality

    def test_getitem(self, some_int_array) -> None:
        # Test method __getitem__
        cov_aray: CoverageArray = CoverageArray(some_int_array)
        expected: int

        expected = float('inf')
        assert expected == cov_aray[-1]
        expected = 1
        assert expected == cov_aray[0]
        expected = 5
        assert expected == cov_aray[4]
        expected = 8
        assert expected == cov_aray[len(cov_aray)-1]
        expected = float('inf')
        assert expected == cov_aray[len(cov_aray)]
    # end def test_lengths_equality
# end class TestCoverageArray