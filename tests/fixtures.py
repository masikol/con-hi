# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
from typing import Sequence

import pytest
from Bio.SeqRecord import SeqRecord

import src.parse_fasta_reference as pfr


@pytest.fixture(scope='session')
def first_test_seq_id() -> str:
    return 'test_seq_1'
# end def first_test_seq_id

@pytest.fixture(scope='session')
def first_test_seq_len() -> int:
    return 70
# end def first_test_seq_len


@pytest.fixture(scope='session')
def second_test_seq_id() -> str:
    return 'test_seq_2'
# end def second_test_seq_id

@pytest.fixture(scope='session')
def second_test_seq_len() -> int:
    return 70
# end def second_test_seq_len


@pytest.fixture(scope='session')
def test_fasta_fpath() -> str:
    return os.path.join(os.getcwd(), 'tests', 'data', 'test_reference.fasta')
# end def test_fasta_fpath


@pytest.fixture(scope='session')
def test_bam_fpath() -> str:
    return os.path.join(os.getcwd(), 'tests', 'data', 'test-mapping.sorted.bam')
# end def test_bam_fpath



@pytest.fixture(scope='session')
def test_fasta_records() -> Sequence[SeqRecord]:
    return pfr.parse_fasta_reference(
        os.path.join(os.getcwd(), 'tests', 'data', 'test_reference.fasta')
    )
# end def test_fasta_records


@pytest.fixture(scope='session')
def test_outdir_path(tmpdir_factory) -> str:
    outdir_path: str = tmpdir_factory.mktemp('test_outdir')
    return outdir_path
# end def test_outdir_path

@pytest.fixture(scope='session')
def test_coverage_fpath(test_outdir_path) -> str:
    return os.path.join(test_outdir_path, 'coverages.tsv')
# end def test_coverage_fpath
