# -*- encoding: utf-8 -*-
# Version 1.0.a

import array
import statistics as sts
from typing import Sequence


class CoverageArray:
    # Class represents array of coverage values.

    def __init__(self, initializer: Sequence[int]) -> None:
        self.coverages = array.array('I', initializer)
    # end def __init__

    def __getitem__(self, key: int) -> int:
        # Returns coverage at given position `key`

        # We will return infinity if inde is out of bounds
        if key < 0 or key > len(self.coverages) - 1:
            return float('inf')
        else:
            return self.coverages[key]
        # end if
    # end def __getitem__

    def __len__(self) -> int:
        return len(self.coverages)
    # end def __len__

    def calc_mean_coverage(self):
        return round(sts.mean(self.coverages), 2)
    # end def
# end class CoverageArray
