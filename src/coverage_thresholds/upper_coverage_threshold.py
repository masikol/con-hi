
from src.coverage_thresholds.lower_coverage_threshold import CoverageThreshold


class UpperCoverageThreshold(CoverageThreshold):
    # Class represents coverage thresholds for annotating high-coverage regions.

    def __init__(self, cov_threshold_value: int) -> None:
        if cov_threshold_value < 0:
            raise ValueError(
                'Negative value passed to constructor of class '
                f'`{self.__class__}`: `{cov_threshold_value}`'
            )
        # end if

        self._value: int = cov_threshold_value

        # Zero coverage requires different label
        self._label: str = f'coverage > {self._value}'
    # end def

    def test_coverage(self, cov: int) -> bool:
        # Function checks if current coverage `cov` is below our threshold
        if cov is None:
            return False
        # end if

        return cov > self._value
    # end def

    def __repr__(self):
        return f'<UpperCoverageThreshold: {self._value}, label:`{self._label}`>'
    # end def
# end class


class UpperMedianCoverageThreshold(CoverageThreshold):
    # Class represents coverage thresholds for annotating low-coverage regions.

    def __init__(self, cov_threshold_coef: int, median_coverage: int) -> None:
        if cov_threshold_coef < 0:
            raise ValueError(
                'Negative value passed to constructor of class '
                f'`{self.__class__}`: `{cov_threshold_coef}`'
            )
        # end if

        # Floor product
        self._value: int = int(median_coverage * cov_threshold_coef)

        times_sign: str = b'\xc3\x97'.decode('utf-8')
        self._label: str = 'coverage > ({}{}median) = {}'.format(
            cov_threshold_coef,
            times_sign,
            self._value
        )
    # end def

    def test_coverage(self, cov: int) -> bool:
        # Function checks if current coverage `cov` is above our threshold
        if cov is None:
            return False
        # end if

        return cov > self._value
    # end def

    def __repr__(self):
        return f'<UpperMedianCoverageThreshold: {self._value}, label:`{self._label}`>'
    # end def
# end class
