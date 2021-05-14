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
    prev_cov: int
    curr_cov: int
    next_cov: int

    for i in range(len(cov_array)):

        prev_cov = cov_array[i-1]
        curr_cov = cov_array[ i ]
        next_cov = cov_array[i+1]

        if feature_condition(curr_cov) and not feature_condition(prev_cov):
            # Note that the start and end location numbering follow Python’s scheme,
            #   thus a GenBank entry of 123..150 (one based counting)
            #   becomes a location of [122:150] (zero based counting).
            feature_start = i # do not make it 1-based
        # end if

        if feature_condition(curr_cov) and not feature_condition(next_cov):

            feature_end = i + 1 # make it 1-based

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
