# -*- encoding: utf-8 -*-
# Version 1.0.a

from typing import Callable, List

from Bio.SeqFeature import SeqFeature, FeatureLocation

from src.coverage_array import CoverageArray


def highlight_coverage_features(
    cov_array: CoverageArray,
    feature_condition: Callable[[int], bool],
    feature_label: str):

    feature_list: List[SeqFeature] = list()

    i: int
    feature_start: int = None
    feature_end: int = None

    for i in range(len(cov_array)):

        if feature_condition(i) and not feature_condition(i-1):
            feature_start = i
        # end if

        if feature_condition(i) and not feature_condition(i+1):

            feature_end = i

            feature_list.append(
                SeqFeature(
                    FeatureLocation(start=feature_start, end=feature_end),
                    type='misc_feature',
                    qualifiers={'label': feature_label}
                )
            )

            feature_start = None
            feature_end = None
        # end if
    # end for

    return feature_list
# end def highlight_coverage_features
