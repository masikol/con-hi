
import array
import statistics as sts
from typing import Sequence


class CoverageArray:
    # Class represents array of coverage values.

    def __init__(self, initializer: Sequence[int]) -> None:
        self.coverages = array.array('I', initializer)
        self.min_coverage = self._calc_min_coverage()
        self.max_coverage = self._calc_max_coverage()
        self.avg_coverage = self._calc_avg_coverage()
        self.median_coverage = self._calc_median_coverage()
        self.zero_coverage_bases = self._count(0)
    # end def

    def __getitem__(self, key: int) -> int:
        # Returns coverage at given position `key`

        # We will return infinity if inde is out of bounds
        if key < 0 or key > len(self.coverages) - 1:
            return None
        else:
            return self.coverages[key]
        # end if
    # end def

    def __len__(self) -> int:
        return len(self.coverages)
    # end def

    def _count(self, cov_value: int) -> int:
        return self.coverages.count(cov_value)
    # end def

    def _calc_min_coverage(self):
        return round(min(self.coverages), 2)
    # end def

    def _calc_max_coverage(self):
        return round(max(self.coverages), 2)
    # end def

    def _calc_avg_coverage(self):
        return round(sts.mean(self.coverages), 2)
    # end def

    def _calc_median_coverage(self):
        return round(sts.median(self.coverages), 2)
    # end def
# end class
