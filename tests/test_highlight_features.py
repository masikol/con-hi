
from typing import List

import pytest
from Bio.SeqFeature import SeqFeature

import src.obtain_coverage as oc
import src.highlight_features as hlft
from src.coverage_array import CoverageArray
from src.coverage_thresholds.lower_coverage_threshold import LowerCoverageThreshold

from tests.fixtures import test_bam_fpath, test_outfpath, \
                           test_outdir_path, test_coverage_fpath, \
                           first_test_seq_id, second_test_seq_id, \
                           coverage_array_inner,coverage_array_edge, \
                           nonzero_cov_threshold


@pytest.fixture
def zero_cov_threshold() -> LowerCoverageThreshold:
    return LowerCoverageThreshold(0)
# end def


class TestHighlightCoverageFeatures:
    # Class for testing fuction `src.highlight_features.highlight_coverage_features`

    def test_highlight_inner_region(self, coverage_array_inner, nonzero_cov_threshold) -> None:
        # Function tests how `highlight_coverage_features` deals
        #   with single inner low-coverage region

        note: str = 'some note whatever'

        feature_list: List[SeqFeature] = hlft.highlight_coverage_features(
            coverage_array_inner,
            nonzero_cov_threshold,
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
        expected_label: str = nonzero_cov_threshold.get_label()
        assert ftr.qualifiers['label'] == expected_label

        # Test feature note
        expected_note: str = note
        assert ftr.qualifiers['note'] == expected_note
    # end def


    def test_highlight_inner_zerocov_region(self, coverage_array_inner, zero_cov_threshold) -> None:
        # Function tests how `highlight_coverage_features` deals
        #   with single inner zero-coverage region

        note: str = 'some note whatever'

        feature_list: List[SeqFeature] = hlft.highlight_coverage_features(
            coverage_array_inner,
            zero_cov_threshold,
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
        expected_label: str = zero_cov_threshold.get_label()
        assert ftr.qualifiers['label'] == expected_label

        # Test feature note
        expected_note: str = note
        assert ftr.qualifiers['note'] == expected_note
    # end def


    def test_highlight_edge_region(self, coverage_array_edge, nonzero_cov_threshold) -> None:
        # Function tests how `highlight_coverage_features` deals
        #   with single edge low-coverage region

        note: str = 'some note whatever'

        feature_list: List[SeqFeature] = hlft.highlight_coverage_features(
            coverage_array_edge,
            nonzero_cov_threshold,
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
        expected_label: str = nonzero_cov_threshold.get_label()
        assert ftr.qualifiers['label'] == expected_label

        # Test feature note
        expected_note: str = note
        assert ftr.qualifiers['note'] == expected_note
    # end def
# end class
