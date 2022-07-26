
import os
from typing import List, Tuple

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import src.output as out
import src.highlight_features as hlft
from src.coverage_array import CoverageArray
from src.coverage_thresholds.lower_coverage_threshold import LowerCoverageThreshold

from tests.fixtures import test_fasta_records, test_bam_fpath
from tests.fixtures import nonzero_cov_threshold, test_coverage_fpath
from tests.fixtures import test_outfpath, test_outdir_path
from tests.fixtures import first_test_seq_id, second_test_seq_id
from tests.fixtures import coverage_array_inner, coverage_array_edge


@pytest.fixture
def annotated_seq_records(test_fasta_records: SeqRecord,
                          coverage_array_inner: CoverageArray,
                          coverage_array_edge: CoverageArray,
                          nonzero_cov_threshold: LowerCoverageThreshold) -> Tuple:

    seq_records = list(test_fasta_records)

    note: str = 'some note whatever'

    coverage_features: List[SeqFeature] = hlft.highlight_coverage_features(
        coverage_array_inner,
        nonzero_cov_threshold,
        note
    )

    seq_records[0].features.extend(coverage_features)

    coverage_features: List[SeqFeature] = hlft.highlight_coverage_features(
        coverage_array_edge,
        nonzero_cov_threshold,
        note
    )

    seq_records[1].features.extend(coverage_features)

    return tuple(seq_records)
# end def


class TestWriteGenbankOutput:
    # Class for testing function `src.output.write_genbank_output`

    def test_write_genbank(self,
                           annotated_seq_records,
                           coverage_array_inner,
                           coverage_array_edge,
                           test_outfpath) -> None:
        # Function tests how `write_genbank_output` writes correct GenBank file
        first_annotated_seq_record, second_annotated_seq_record = \
            annotated_seq_records

        testing_zip = zip(
            (first_annotated_seq_record, second_annotated_seq_record),
            (coverage_array_inner, coverage_array_edge)
        )

        for seq_record, cov_array in testing_zip:
            self.check_genbank_file(seq_record, cov_array, test_outfpath)
        # end def
    # end def

    def check_genbank_file(self, seq_record, cov_array, test_outfpath):
        topology: str = 'linear'
        organism: str = 'Czort lysy'
        outfpath: str = test_outfpath

        record_seq_len: int = len(seq_record)

        out.write_genbank_output(
            seq_record,
            topology,
            organism,
            cov_array,
            outfpath
        )

        # Check if outfile exists
        assert os.path.exists(outfpath)

        # Parse this file
        parsed_seq_record: SeqRecord = list(SeqIO.parse(outfpath, 'genbank'))[0]
        print('\n\n\n{}'.format(seq_record.id))
        print('{}'.format(parsed_seq_record.id))

        # Test simple file contents
        assert len(parsed_seq_record) == record_seq_len
        assert seq_record.seq == parsed_seq_record.seq
        assert parsed_seq_record.annotations['molecule_type'] == 'DNA'
        assert parsed_seq_record.annotations['organism'] == organism
        assert parsed_seq_record.annotations['topology'] == topology

        # Test if features are sorted in an ascending order
        for i in range(len(parsed_seq_record.features) - 1):
            curr_pos: int = parsed_seq_record.features[i].location.start
            next_pos: int = parsed_seq_record.features[i+1].location.start
            assert curr_pos <= next_pos
        # end for

        assert 'structured_comment' in parsed_seq_record.annotations
        assert 'Coverage-Data' in parsed_seq_record.annotations['structured_comment']
        assert 'Minimum Coverage' in parsed_seq_record.annotations['structured_comment']['Coverage-Data']
        assert 'Average Coverage' in parsed_seq_record.annotations['structured_comment']['Coverage-Data']
        assert 'Median Coverage' in parsed_seq_record.annotations['structured_comment']['Coverage-Data']
        assert 'Maximum Coverage' in parsed_seq_record.annotations['structured_comment']['Coverage-Data']

        empty_file(outfpath)
    # end def
# end class


def empty_file(fpath):
    with open(fpath, 'w'):
        pass
    # end with
# end def
