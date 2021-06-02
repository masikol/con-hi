# -*- encoding: utf-8 -*-
# Version 1.0.a

import pytest

from src.coverage_threshold import CoverageThreshold


class TestCoverageThreshold:
    # Class for testing class `src.coverage_threshold.CoverageThreshold`

    def test_invalid_coverage(self) -> None:
        # Test if proper exception is raised on invalid coverage value passed to constructor
        with pytest.raises(ValueError):
            CoverageThreshold(-9)
        # end with
    # end def test_invalid_coverage

    def test_valid_coverage(self) -> None:
        thr_value: int = 7
        cov_threshold: CoverageThreshold = CoverageThreshold(thr_value)

        assert cov_threshold.get_coverage() == thr_value
        assert cov_threshold.get_label() == f'coverage < {thr_value}'
        assert cov_threshold.test_coverage(thr_value - 1) == True
        assert cov_threshold.test_coverage(thr_value) == False
        assert cov_threshold.test_coverage(thr_value + 1) == False
    # end def test_valid_coverage


    def test_zero_coverage(self) -> None:
        thr_value: int = 0
        cov_threshold: CoverageThreshold = CoverageThreshold(thr_value)

        assert cov_threshold.get_coverage() == thr_value
        assert cov_threshold.get_label() == 'zero coverage'
        assert cov_threshold.test_coverage(thr_value) == True
        assert cov_threshold.test_coverage(thr_value + 1) == False
    # end def test_valid_coverage
# end class TestCoverageThreshold
