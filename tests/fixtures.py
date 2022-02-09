# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
from typing import Sequence

import pytest
from Bio.SeqRecord import SeqRecord

import src.obtain_coverage as oc
import src.parse_fasta_reference as pfr
from src.coverage_array import CoverageArray
from src.coverage_threshold import CoverageThreshold


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
def test_outfpath(test_outdir_path) -> str:
    return os.path.join(test_outdir_path, 'test_annotated_seq.gbk')
# end def test_outfpath


@pytest.fixture(scope='session')
def test_coverage_fpath(test_outdir_path) -> str:
    return os.path.join(test_outdir_path, 'coverages.tsv')
# end def test_coverage_fpath


@pytest.fixture(scope='session')
def nonzero_cov_threshold() -> CoverageThreshold:
    return CoverageThreshold(2)
# end def nonzero_cov_threshold


@pytest.fixture(scope='session')
def coverage_array_inner(test_coverage_fpath: str,
                         first_test_seq_id: str) -> CoverageArray:

    return oc.get_coverage_for_reference(first_test_seq_id, test_coverage_fpath)
# end def coverage_array_inner


@pytest.fixture(scope='session')
def coverage_array_edge(test_coverage_fpath,
                        second_test_seq_id) -> CoverageArray:

    return oc.get_coverage_for_reference(second_test_seq_id, test_coverage_fpath)
# end def coverage_array_edge
