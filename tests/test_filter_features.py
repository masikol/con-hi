
from typing import Sequence, MutableSequence, NewType

import pytest
from Bio.SeqFeature import SeqFeature, FeatureLocation

from src.filter_features import filter_features

FeatureList = NewType('FeatureList', MutableSequence[SeqFeature])


@pytest.fixture
def empty_feature_list() -> FeatureList:
    return []
# end def

@pytest.fixture
def feature_list_1() -> FeatureList:
    offset: int = 5
    lengths: Sequence = (1, 5, 10)

    return [
        SeqFeature(FeatureLocation(offset, offset + l)) for l in lengths
    ]
# end def


class TestFilterFeatures:
    # Class for testing function `src.filter_features.filter_features`

    def test_empty_feature_list(self, empty_feature_list) -> None:
        # Test how the function hadles empty feature list
        min_feature_len: int = 1
        filtered_features: FeatureList = filter_features(
            empty_feature_list,
            min_feature_len
        )

        assert len(filtered_features) == 0
    # end def

    def test_all_features_fail(self, feature_list_1) -> None:
        # Test how the function hadles empty feature list
        min_feature_len: int = 20
        filtered_features: FeatureList = filter_features(
            feature_list_1,
            min_feature_len
        )

        assert len(filtered_features) == 0
    # end def

    def test_some_features_fail(self, feature_list_1) -> None:
        # Test how the function hadles empty feature list
        min_feature_len: int = 4
        filtered_features: FeatureList = filter_features(
            feature_list_1,
            min_feature_len
        )

        assert len(filtered_features) == 2
    # end def

    def test_all_features_pass(self, feature_list_1) -> None:
        # Test how the function hadles empty feature list
        min_feature_len: int = 1
        filtered_features: FeatureList = filter_features(
            feature_list_1,
            min_feature_len
        )

        assert len(filtered_features) == len(feature_list_1)
    # end def
# end class
