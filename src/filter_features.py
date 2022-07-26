
from typing import MutableSequence

from Bio.SeqFeature import SeqFeature


def filter_features(features: MutableSequence[SeqFeature],
                    min_feature_len: int):
    len_is_enough = lambda f: len(f) >= min_feature_len
    return list(
        filter(len_is_enough, features)
    )
# end def
