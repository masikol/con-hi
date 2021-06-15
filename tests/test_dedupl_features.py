# -*- encoding: utf-8 -*-
# Version 1.1.a

import copy
from typing import MutableSequence, NewType

import pytest
from Bio.SeqFeature import SeqFeature, FeatureLocation

from src.dedupl_features import dedupl_features


FeatureList = NewType('FeatureList', MutableSequence[SeqFeature])


@pytest.fixture
def empty_feature_list() -> FeatureList:
    return []
# end def empty_feature_list


@pytest.fixture
def extant_feature_list_1() -> FeatureList:
    return [
        SeqFeature(FeatureLocation(5, 10)),
        SeqFeature(FeatureLocation(25, 60))
    ]
# end def extant_feature_list_1


@pytest.fixture
def new_feature_list_1() -> FeatureList:
    return [
        SeqFeature(FeatureLocation(6, 11)),
        SeqFeature(FeatureLocation(38, 90))
    ]
# end def new_feature_list_1


@pytest.fixture
def new_feature_list_2() -> FeatureList:
    return [
        SeqFeature(FeatureLocation(6, 10)),
        SeqFeature(FeatureLocation(38, 90))
    ]
# end def new_feature_list_2


@pytest.fixture
def new_feature_list_3() -> FeatureList:
    return [
        SeqFeature(FeatureLocation(5, 10)),
        SeqFeature(FeatureLocation(38, 90))
    ]
# end def new_feature_list_3


@pytest.fixture
def new_feature_list_4() -> FeatureList:
    return [
        SeqFeature(FeatureLocation(5, 10)),
        SeqFeature(FeatureLocation(25, 60))
    ]
# end def new_feature_list_4


class TestDeduplFeatures:
    # Class for testing function `src.dedupl_features.dedupl_features`

    def test_both_empty(self, empty_feature_list) -> None:
        # Test how the function hadles both empty lists
        new_feature_list: FeatureList = copy.deepcopy(empty_feature_list)
        extant_feature_list: FeatureList = copy.deepcopy(empty_feature_list)

        deduplicated_features: FeatureList = dedupl_features(
            new_feature_list,
            extant_feature_list
        )

        assert len(deduplicated_features) == 0
    # end def test_both_empty

    def test_new_empty(self, empty_feature_list, extant_feature_list_1) -> None:
        # Test how the function hadles case when `new_features` is empty
        new_feature_list: FeatureList = copy.deepcopy(empty_feature_list)
        extant_feature_list: FeatureList = copy.deepcopy(extant_feature_list_1)

        deduplicated_features: FeatureList = dedupl_features(
            new_feature_list,
            extant_feature_list
        )

        assert len(deduplicated_features) == 0
    # end def test_new_empty

    def test_extant_empty(self, new_feature_list_1, empty_feature_list) -> None:
        # Test how the function hadles case when `extant_features` is empty
        new_feature_list: FeatureList = copy.deepcopy(new_feature_list_1)
        extant_feature_list: FeatureList = copy.deepcopy(empty_feature_list)
        int_len_new: int = len(new_feature_list)

        deduplicated_features: FeatureList = dedupl_features(
            new_feature_list,
            extant_feature_list
        )

        assert len(deduplicated_features) == int_len_new
    # end def test_extant_empty

    def test_all_new_unique(self, new_feature_list_1, extant_feature_list_1) -> None:
        # Test how the function hadles case when all `new_features` are unique
        new_feature_list: FeatureList = copy.deepcopy(new_feature_list_1)
        extant_feature_list: FeatureList = copy.deepcopy(extant_feature_list_1)
        int_len_new: int = len(new_feature_list)

        deduplicated_features: FeatureList = dedupl_features(
            new_feature_list,
            extant_feature_list
        )

        assert len(deduplicated_features) == int_len_new
    # end def test_all_new_unique

    def test_one_new_half_unique(self, new_feature_list_2, extant_feature_list_1) -> None:
        # Test how the function hadles case when one of `new_features` has only start
        #   occuring somewhere in `extant_features`
        new_feature_list: FeatureList = copy.deepcopy(new_feature_list_2)
        extant_feature_list: FeatureList = copy.deepcopy(extant_feature_list_1)
        int_len_new: int = len(new_feature_list)

        deduplicated_features: FeatureList = dedupl_features(
            new_feature_list,
            extant_feature_list
        )

        assert len(deduplicated_features) == int_len_new
    # end def test_one_new_half_unique

    def test_one_new_dupl(self, new_feature_list_3, extant_feature_list_1) -> None:
        # Test how the function hadles case when one of `new_features` is duplicated
        new_feature_list: FeatureList = copy.deepcopy(new_feature_list_3)
        extant_feature_list: FeatureList = copy.deepcopy(extant_feature_list_1)
        int_len_new: int = len(new_feature_list)

        deduplicated_features: FeatureList = dedupl_features(
            new_feature_list,
            extant_feature_list
        )

        assert len(deduplicated_features) == int_len_new - 1
    # end def test_one_new_dupl

    def test_both_new_dupl(self, new_feature_list_4, extant_feature_list_1) -> None:
        # Test how the function hadles case when two of `new_features` are duplicated
        new_feature_list: FeatureList = copy.deepcopy(new_feature_list_4)
        extant_feature_list: FeatureList = copy.deepcopy(extant_feature_list_1)
        int_len_new: int = len(new_feature_list)

        deduplicated_features: FeatureList = dedupl_features(
            new_feature_list,
            extant_feature_list
        )

        assert len(deduplicated_features) == int_len_new - 2
    # end def test_both_new_dupl
# end class TestDeduplFeatures
