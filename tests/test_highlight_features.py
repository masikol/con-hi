# -*- encoding: utf-8 -*-
# Version 1.0.a

from typing import List

import pytest
from Bio.SeqFeature import SeqFeature

import src.obtain_coverage as oc
import src.highlight_features as hlft
from src.coverage_array import CoverageArray
from src.coverage_threshold import CoverageThreshold

from tests.fixtures import test_outdir_path
from tests.fixtures import test_fasta_fpath, test_bam_fpath, test_coverage_fpath, test_outfpath
from tests.fixtures import first_test_seq_id, second_test_seq_id


@pytest.fixture
def test_nonzero_cov_threshold() -> CoverageThreshold:
    return CoverageThreshold(2)
# end def test_nonzero_cov_threshold

@pytest.fixture
def test_zero_cov_threshold() -> CoverageThreshold:
    return CoverageThreshold(0)
# end def test_zero_cov_threshold

@pytest.fixture
def test_coverage_array_inner(
    test_fasta_fpath,
    test_bam_fpath,
    test_outfpath,
    first_test_seq_id) -> CoverageArray:

    cov_fpath: str = oc.count_cov_for_all_refs(test_bam_fpath, test_outfpath)
    return oc.get_coverage_for_reference(first_test_seq_id, cov_fpath)
# end def test_coverage_array_inner

@pytest.fixture
def test_coverage_array_edge(
    test_fasta_fpath,
    test_bam_fpath,
    test_outfpath,
    second_test_seq_id) -> CoverageArray:

    cov_fpath: str = oc.count_cov_for_all_refs(test_bam_fpath, test_outfpath)
    return oc.get_coverage_for_reference(second_test_seq_id, cov_fpath)
# end def test_coverage_array_edge


class TestHighlightCoverageFeatures:
    # Class for testing fuction `src.highlight_features.highlight_coverage_features`

    def test_highlight_inner_region(self, test_coverage_array_inner, test_nonzero_cov_threshold) -> None:
        # Function tests how `highlight_coverage_features` deals
        #   with single inner low-coverage region

        note: str = 'some note whatever'

        feature_list: List[SeqFeature] = hlft.highlight_coverage_features(
            test_coverage_array_inner,
            test_nonzero_cov_threshold,
            note
        )

        # Test number of features
        expected_feature_num: int = 1
        assert len(feature_list) == expected_feature_num

        ftr: SeqFeature = feature_list[0]

        # Test feature location(s)
        expected_start_pos: int = 23 - 1 # to 0-based (1-based position is 23)
        assert ftr.location.start == expected_start_pos
        expected_end_pos: int = 49
        assert ftr.location.end   == expected_end_pos

        # Test feature type
        expected_feature_type: str = 'misc_feature'
        assert ftr.type == expected_feature_type

        # Test feature label
        expected_label: str = test_nonzero_cov_threshold.get_label()
        assert ftr.qualifiers['label'] == expected_label

        # Test feature note
        expected_note: str = note
        assert ftr.qualifiers['note'] == expected_note
    # end def test_highlight_inner_region


    def test_highlight_inner_zerocov_region(self, test_coverage_array_inner, test_zero_cov_threshold) -> None:
        # Function tests how `highlight_coverage_features` deals
        #   with single inner zero-coverage region

        note: str = 'some note whatever'

        feature_list: List[SeqFeature] = hlft.highlight_coverage_features(
            test_coverage_array_inner,
            test_zero_cov_threshold,
            note
        )

        # Test number of features
        expected_feature_num: int = 1
        assert len(feature_list) == expected_feature_num

        ftr: SeqFeature = feature_list[0]

        # Test feature location(s)
        expected_start_pos: int = 27 - 1 # to 0-based (1-based position is 23)
        assert ftr.location.start == expected_start_pos
        expected_end_pos: int = 40
        assert ftr.location.end   == expected_end_pos

        # Test feature type
        expected_feature_type: str = 'misc_feature'
        assert ftr.type == expected_feature_type

        # Test feature label
        expected_label: str = test_zero_cov_threshold.get_label()
        assert ftr.qualifiers['label'] == expected_label

        # Test feature note
        expected_note: str = note
        assert ftr.qualifiers['note'] == expected_note
    # end def test_highlight_inner_zerocov_region


    def test_highlight_edge_region(self, test_coverage_array_edge, test_nonzero_cov_threshold) -> None:
        # Function tests how `highlight_coverage_features` deals
        #   with single edge low-coverage region

        note: str = 'some note whatever'

        feature_list: List[SeqFeature] = hlft.highlight_coverage_features(
            test_coverage_array_edge,
            test_nonzero_cov_threshold,
            note
        )

        # Test number of features
        expected_feature_num: int = 1
        assert len(feature_list) == expected_feature_num

        ftr: SeqFeature = feature_list[0]

        # Test feature location(s)
        expected_start_pos: int = 35 - 1 # to 0-based (1-based position is 23)
        assert ftr.location.start == expected_start_pos
        expected_end_pos: int = 70
        assert ftr.location.end   == expected_end_pos

        # Test feature type
        expected_feature_type: str = 'misc_feature'
        assert ftr.type == expected_feature_type

        # Test feature label
        expected_label: str = test_nonzero_cov_threshold.get_label()
        assert ftr.qualifiers['label'] == expected_label

        # Test feature note
        expected_note: str = note
        assert ftr.qualifiers['note'] == expected_note
    # end def test_highlight_edge_region
# end class TestHighlightCoverageFeatures
