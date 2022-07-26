
from src.coverage_thresholds.coverage_threshold import CoverageThreshold


class LowerCoverageThreshold(CoverageThreshold):
    # Class represents coverage thresholds for annotating low-coverage regions.

    def __init__(self, cov_threshold_value: int) -> None:
        if cov_threshold_value < 0:
            raise ValueError(
                'Negative value passed to constructor of class '
                f'`{self.__class__}`: `{cov_threshold_value}`'
            )
        # end if

        self._value: int = cov_threshold_value

        # Zero coverage requires different label
        self._label: str
        if self._value == 0:
            self._label = 'zero coverage'
        else:
            self._label = f'coverage < {self._value}'
        # end if
    # end def

    def test_coverage(self, cov: int) -> bool:
        # Function checks if current coverage `cov` is below our threshold

        # Zero coverage requires different behaviour
        if self._value != 0:
            return cov < self._value
        else:
            return cov == 0
        # end if
    # end def

    def __repr__(self):
        return f'<LowerCoverageThreshold: {self._value}, label:`{self._label}`>'
    # end def
# end class
