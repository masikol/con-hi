# -*- encoding: utf-8 -*-
# Version 1.0.a

import os

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import src.output as out

from tests.fixtures import test_outfpath, test_outdir_path


@pytest.fixture
def some_seq_record() -> str:
    record: SeqRecord = SeqRecord(Seq('ACTAGCTAGTCTAGTGCAGCGTAGCTGTGTGTA'), id='test_id')
    record.features = [
        SeqFeature(
            FeatureLocation(start=15, end=28),
            type='misc_feature',
            qualifiers={
                'label': 'coverage < 2',
                'note': 'never mind'
            }
        ),
        SeqFeature(
            FeatureLocation(start=7, end=12),
            type='misc_feature',
            qualifiers={
                'label': 'coverage < 2',
                'note': 'never mind'
            }
        )
    ]
    return record
# end def some_seq_record


class TestWriteGenbankOutput:
    # Class for testing function `src.output.write_genbank_output`

    def test_some_seq_record(self, some_seq_record, test_outfpath) -> None:
        # Function tests how `write_genbank_output` writes correct GenBank file

        topology: str = 'linear'
        organism: str = 'Czort lysy'
        outfpath: str = test_outfpath

        record_seq_len: int = len(some_seq_record)

        out.write_genbank_output(some_seq_record, topology, organism, outfpath)

        # Check if outfile exists
        assert os.path.exists(outfpath)

        # Parse this file
        parsed_seq_record: SeqRecord = list(SeqIO.parse(outfpath, 'genbank'))[0]

        # Test simple file contents
        assert len(parsed_seq_record) == record_seq_len
        assert some_seq_record.seq == parsed_seq_record.seq
        assert parsed_seq_record.annotations['molecule_type'] == 'DNA'
        assert parsed_seq_record.annotations['organism'] == organism
        assert parsed_seq_record.annotations['topology'] == topology

        # Test if features are sorted in an ascending order
        for i in range(len(parsed_seq_record.features) - 1):
            curr_pos: int = parsed_seq_record.features[i].location.start
            next_pos: int = parsed_seq_record.features[i+1].location.start
            assert curr_pos <= next_pos
        # end for
    # end def test_some_seq_record
# end class TestWriteGenbankOutput