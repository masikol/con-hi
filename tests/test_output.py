# -*- encoding: utf-8 -*-
# Version 1.0.a

import os

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import src.output as out

from tests.fixtures import test_outdir_path


@pytest.fixture
def test_seq_id() -> str:
    return 'CP000000.1'
# end def test_seq_id

@pytest.fixture
def test_prefix() -> str:
    return 'some_prefix'
# end def test_prefix


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


class TestConfigureOutfilePath:
    # Class for testing function `src.output.configure_outfile_path`

    def test_conf_outfpath_with_prefix(self, test_outdir_path, test_prefix, test_seq_id) -> None:
        # Function for testing how `configure_outfile_path` configures prefix
        #   with specified prefix

        curr_outfpath: str = out.configure_outfile_path(
            test_outdir_path,
            test_prefix,
            test_seq_id
        )

        expeced_fpath: str = os.path.join(test_outdir_path, f'{test_prefix}_{test_seq_id}.gbk')
        assert curr_outfpath == expeced_fpath
    # end def test_conf_outfpath_with_prefix

    def test_conf_outfpath_no_prefix(self, test_outdir_path, test_seq_id) -> None:
        # Function for testing how `configure_outfile_path` configures prefix
        #   without specified prefix

        curr_outfpath: str = out.configure_outfile_path(
            test_outdir_path,
            '',
            test_seq_id
        )

        expeced_fpath: str = os.path.join(test_outdir_path, f'{test_seq_id}.gbk')
        assert curr_outfpath == expeced_fpath
    # end def test_conf_outfpath_no_prefix
# end class TestConfigureOutfilePath


class TestWriteGenbankOutput:
    # Class for testing function `src.output.write_genbank_output`

    def test_some_seq_record(self, some_seq_record, test_outdir_path, test_seq_id) -> None:
        # Function tests how `write_genbank_output` writes correct GenBank file

        topology: str = 'linear'
        organism: str = 'Czort lysy'
        outfpath: str = out.configure_outfile_path(test_outdir_path, '', test_seq_id)

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