# Version 2.1.a


class CoverageThreshold:
    # Class represents coverage thresholds for annotating low-coverage regions.

    def __init__(self, cov_threshold_value: int) -> None:

        if cov_threshold_value < 0:
            raise ValueError(f'Negative value passed to constructor of class\
 CoverageThreshold: `{cov_threshold_value}`')
        # end if

        self._value: int = cov_threshold_value

        # Zero coverage requires different label
        self._label: str
        if self._value == 0:
            self._label = 'zero coverage'
        else:
            self._label = f'coverage < {self._value}'
        # end if
    # end def __init__


    def test_coverage(self, cov: int) -> bool:
        # Function checks if current coverage `cov` is below our threshold

        # Zero coverage requires different behaviour
        if self._value != 0:
            return cov < self._value
        else:
            return cov == 0
        # end if
    # end def test_coverage


    def get_label(self) -> str:
        # Getter for label
        return self._label
    # end def get_label


    def get_coverage(self) -> int:
        # Getter for coverage threshold
        return self._value
    # end def get_coverage


    def __repr__(self):
        return f'<CoverageThreshold: {self._value}, label:`{self._label}`>'
    # end def __repr__
# end class CoverageThreshold
