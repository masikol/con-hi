# -*- encoding: utf-8 -*-
# Version 1.0.a

from typing import List, MutableSequence

from Bio.SeqFeature import SeqFeature, FeatureLocation

from src.coverage_array import CoverageArray
from src.coverage_threshold import CoverageTheshold


def highlight_coverage_features(
        cov_array: CoverageArray,
        coverage_threshold: CoverageTheshold,
        base_feature_note: str
    ) -> MutableSequence[SeqFeature]:
    # Function annotates low-coverage regions according to `coverage_threshold`
    # :param cov_array: array of coverage values;
    # :param coverage_threshold: coverage thresold to search low-coverage regions;
    # :param base_feature_note: string to `note` qualifier of output GenBank file;

    # Save following objects to local variables
    test_coverage = coverage_threshold.test_coverage
    coverage_label = coverage_threshold.get_label()

    # Init list of features
    feature_list: List[SeqFeature] = list()

    # Variables defined in the loop below
    i: int
    feature_start: int = None
    feature_end: int = None
    prev_cov: int
    curr_cov: int
    next_cov: int

    # Go through the sequence (coverage array) and anotate low-coverage regions
    for i in range(len(cov_array)):

        # Save three coverages to local variables
        prev_cov = cov_array[i-1]
        curr_cov = cov_array[ i ]
        next_cov = cov_array[i+1]

        # If we enter low-coverage region -- set `feature_start`
        if test_coverage(curr_cov) and not test_coverage(prev_cov):
            # Note that the start and end location numbering follow Python’s scheme,
            #   thus a GenBank entry of 123..150 (one based counting)
            #   becomes a location of [122:150] (zero based counting).
            feature_start = i # do not make it 1-based
        # end if

        # If we leave low-coverage region -- set `feature_end`
        if test_coverage(curr_cov) and not test_coverage(next_cov):

            feature_end = i + 1 # make it 1-based

            # Create new feature and append it to the list
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

            # Reset
            feature_start = None
            feature_end = None
        # end if
    # end for

    return feature_list
# end def highlight_coverage_features
