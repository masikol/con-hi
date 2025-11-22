
import statistics as sts

import pytest

from src.coverage_thresholds.upper_coverage_threshold import UpperCoverageThreshold, \
                                                             UpperMedianCoverageThreshold

from tests.fixtures import coverage_array_inner
from tests.fixtures import test_bam_fpath, test_coverage_fpath, test_outdir_path
from tests.fixtures import first_test_seq_id


class TestUpperCoverageThreshold:
    # Class for testing class
    #   `src.coverage_thresholds.upper_coverage_threshold.UpperCoverageThreshold`

    def test_invalid_coverage(self) -> None:
        # Test if proper exception is raised on invalid coverage value passed to constructor
        with pytest.raises(ValueError):
            UpperCoverageThreshold(-9)
        # end with
    # end def

    def test_valid_coverage(self) -> None:
        thr_value: int = 7
        cov_threshold: UpperCoverageThreshold = UpperCoverageThreshold(thr_value)

        assert cov_threshold.get_coverage() == thr_value
        assert cov_threshold.get_label() == f'coverage > {thr_value}'
        assert cov_threshold.test_coverage(thr_value - 1) == False
        assert cov_threshold.test_coverage(thr_value) == False
        assert cov_threshold.test_coverage(thr_value + 1) == True
    # end def
# end class


class TestUpperMedianCoverageThreshold:
    # Class for testing class
    #   `src.coverage_thresholds.lower_coverage_threshold.UpperMedianCoverageThreshold`

    times_sign = b'\xc3\x97'.decode('utf-8')

    def test_invalid_coverage(self, coverage_array_inner) -> None:
        # Test if proper exception is raised on invalid coverage value passed to constructor
        with pytest.raises(ValueError):
            UpperMedianCoverageThreshold(-9, coverage_array_inner.median_coverage)
        # end with
    # end def

    def test_valid_coverage(self, coverage_array_inner) -> None:
        thr_coef: float = 1.5
        cov_threshold: UpperMedianCoverageThreshold = UpperMedianCoverageThreshold(
            thr_coef,
            coverage_array_inner.median_coverage
        )

        median_coverage: float = 2.0
        final_thr_value: int = int(median_coverage * thr_coef)

        assert cov_threshold.get_coverage() == final_thr_value
        assert cov_threshold.get_label() == 'coverage > ({}{}median) = {}'.format(
            thr_coef, self.times_sign, final_thr_value
        )
        assert cov_threshold.test_coverage(final_thr_value - 1) == False
        assert cov_threshold.test_coverage(final_thr_value) == False
        assert cov_threshold.test_coverage(final_thr_value + 1) == True
    # end def
# end class
