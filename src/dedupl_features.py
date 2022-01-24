# Version 2.1.a

from typing import Sequence, MutableSequence, Tuple

from Bio.SeqFeature import SeqFeature


def dedupl_features(
    new_features: MutableSequence[SeqFeature],
    extant_features: Sequence[SeqFeature]) -> MutableSequence[SeqFeature]:
    # Function absolutely deduplicates feature list.
    # "Absolutely" means that only features starting at the same position
    #   and ending at the same posotion will be removed.
    # Requires that coverage thresholds to be sorted in ascending order.
    # :param new_features: list of recently discovered features;
    # :param extant_features: list of features that existed before `new_features`
    #   was discovered;
    
    # Create the tuple of coordinates in order not to extract these locations each time
    extant_coords: MutableSequence[Tuple[int, int]] = tuple(
        (
            (f.location.start, f.location.end) for f in extant_features
        )
    )

    # Deduplication goes below
    i: int
    start_end: Tuple[int, int]
    for i in range(len(new_features)-1, -1, -1):
        # Tuple of structure (<START_POS>, <END_POS>)
        start_end = (new_features[i].location.start, new_features[i].location.end)
        # If feature at this location exists -- remove it from `extant_features`
        if start_end in extant_coords:
            del new_features[i]
        # end if
    # end for

    return new_features
# end def dedupl_features
