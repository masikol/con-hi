# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
from typing import List

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import src.output as out
import src.highlight_features as hlft
from src.coverage_array import CoverageArray
from src.coverage_threshold import CoverageThreshold

from tests.fixtures import test_fasta_records
from tests.fixtures import nonzero_cov_threshold
from tests.fixtures import test_outfpath, test_outdir_path
from tests.fixtures import first_test_seq_id, second_test_seq_id
from tests.fixtures import coverage_array_inner, coverage_array_edge


@pytest.fixture
def first_annotated_seq_record(test_fasta_records: SeqRecord,
                               coverage_array_inner: CoverageArray,
                               nonzero_cov_threshold: CoverageThreshold) -> str:

    first_record = next(iter(
        filter(lambda r: r.id == first_test_seq_id,
            test_fasta_records)
    ))

    note: str = 'some note whatever'

    coverage_features: List[SeqFeature] = hlft.highlight_coverage_features(
        coverage_array_inner,
        nonzero_cov_threshold,
        note
    )

    first_record.features.extend(coverage_features)

    return first_record
# end def


class TestWriteGenbankOutput:
    # Class for testing function `src.output.write_genbank_output`

    def test_write_genbank(self,
                             first_annotated_seq_record,
                             coverage_array_inner,
                             test_outfpath) -> None:
        # Function tests how `write_genbank_output` writes correct GenBank file

        topology: str = 'linear'
        organism: str = 'Czort lysy'
        outfpath: str = test_outfpath

        record_seq_len: int = len(first_annotated_seq_record)

        out.write_genbank_output(
            first_annotated_seq_record,
            topology,
            organism,
            coverage_array_inner,
            outfpath
        )

        # Check if outfile exists
        assert os.path.exists(outfpath)

        # Parse this file
        parsed_seq_record: SeqRecord = list(SeqIO.parse(outfpath, 'genbank'))[0]

        # Test simple file contents
        assert len(parsed_seq_record) == record_seq_len
        assert first_annotated_seq_record.seq == parsed_seq_record.seq
        assert parsed_seq_record.annotations['molecule_type'] == 'DNA'
        assert parsed_seq_record.annotations['organism'] == organism
        assert parsed_seq_record.annotations['topology'] == topology

        # Test if features are sorted in an ascending order
        for i in range(len(parsed_seq_record.features) - 1):
            curr_pos: int = parsed_seq_record.features[i].location.start
            next_pos: int = parsed_seq_record.features[i+1].location.start
            assert curr_pos <= next_pos
        # end for
    # end def
# end class
