# -*- encoding: utf-8 -*-
# Version 1.0.a

import array
from typing import Sequence

class CoverageArray:

    def __init__(self, initializer: Sequence[int]) -> None:
        self.coverages = array.array('I', initializer)
        self.infinite_coverage = max(self.coverages) + 1
    # end def __init__

    def __getitem__(self, key: int) -> int:

        if key < 0 or key > len(self.coverages) - 1:
            return self.infinite_coverage
        else:
            return self.coverages[key]
        # end if
    # end def __getitem__

    def __len__(self) -> int:
        return len(self.coverages)
    # end def __len__
# end class CoverageArray
