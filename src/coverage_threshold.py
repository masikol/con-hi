# -*- encoding: utf-8 -*-
# Version 1.0.a

class CoverageTheshold:

    def __init__(self, cov_threshold_value: int) -> None:

        if cov_threshold_value < 0:
            raise ValueError(f'Negative value passed to constructor of class\
 CoverageTheshold: `{cov_threshold_value}`')
        # end if

        self._value: int = cov_threshold_value

        self._label: str
        if self._value == 0:
            self._label = 'zero coverage'
        else:
            self._label = f'coverage < {self._value}'
        # end if
    # end def __init__


    def test_coverage(self, cov: int) -> bool:
        if self._value != 0:
            return cov < self._value
        else:
            return cov == 0
        # end if
    # end def test_coverage


    def get_label(self) -> str:
        return self._label
    # end def get_label


    def get_coverage(self) -> int:
        return self._value
    # end def get_coverage


    def __repr__(self):
        return f'<CoverageThreshold: {self._value}, label:`{self._label}`>'
    # end def __repr__
# end class CoverageTheshold
