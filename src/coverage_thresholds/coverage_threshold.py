
class CoverageThreshold:
    # Class represents coverage thresholds for annotating regions of coverage of some type.

    def test_coverage(self, cov: int) -> bool:
        # Function checks if current coverage `cov` is below our threshold
        raise NotImplementedError
    # end def


    def get_label(self) -> str:
        # Getter for label
        return self._label
    # end def


    def get_coverage(self) -> int:
        # Getter for coverage threshold
        return self._value
    # end def
# end class
