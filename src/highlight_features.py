# -*- encoding: utf-8 -*-
# Version 1.0.a

from typing import Callable, List, MutableSequence

from Bio.SeqFeature import SeqFeature, FeatureLocation

from src.coverage_array import CoverageArray
from src.coverage_threshold import CoverageTheshold


def highlight_coverage_features(
        cov_array: CoverageArray,
        coverage_threshold: CoverageTheshold,
        base_feature_note: str
    ) -> MutableSequence[SeqFeature]:

    test_coverage = coverage_threshold.test_coverage
    coverage_label = coverage_threshold.get_label()

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

        if test_coverage(curr_cov) and not test_coverage(prev_cov):
            # Note that the start and end location numbering follow Python’s scheme,
            #   thus a GenBank entry of 123..150 (one based counting)
            #   becomes a location of [122:150] (zero based counting).
            feature_start = i # do not make it 1-based
        # end if

        if test_coverage(curr_cov) and not test_coverage(next_cov):

            feature_end = i + 1 # make it 1-based

            feature_list.append(
                SeqFeature(
                    FeatureLocation(start=feature_start, end=feature_end),
                    type='misc_feature',
                    qualifiers={
                        'label': coverage_label,
                        'note': base_feature_note
                        }
                )
            )

            feature_start = None
            feature_end = None
        # end if
    # end for

    return feature_list
# end def highlight_coverage_features
