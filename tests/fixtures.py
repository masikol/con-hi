
import os
from typing import Sequence

import pytest
from Bio.SeqRecord import SeqRecord

import src.obtain_coverage as oc
import src.parse_fasta_reference as pfr
from src.coverage_array import CoverageArray
from src.coverage_thresholds.lower_coverage_threshold import LowerCoverageThreshold
from src.coverage_thresholds.upper_coverage_threshold import UpperCoverageThreshold



@pytest.fixture(scope='session')
def first_test_seq_id() -> str:
    return 'test_seq_1'
# end def


@pytest.fixture(scope='session')
def first_test_seq_len() -> int:
    return 70
# end def


@pytest.fixture(scope='session')
def second_test_seq_id() -> str:
    return 'test_seq_2'
# end def


@pytest.fixture(scope='session')
def second_test_seq_len() -> int:
    return 70
# end def


@pytest.fixture(scope='session')
def third_test_seq_id() -> str:
    return 'test_seq_3'
# end def

@pytest.fixture(scope='session')
def third_test_seq_len() -> int:
    return 70
# end def


@pytest.fixture(scope='session')
def test_fasta_fpath() -> str:
    return os.path.join(os.getcwd(), 'tests', 'data', 'test_reference.fasta')
# end def


@pytest.fixture(scope='session')
def test_bam_fpath() -> str:
    return os.path.join(os.getcwd(), 'tests', 'data', 'test-mapping.sorted.bam')
# end def


@pytest.fixture(scope='session')
def test_fasta_records() -> Sequence[SeqRecord]:
    return pfr.parse_fasta_reference(
        os.path.join(os.getcwd(), 'tests', 'data', 'test_reference.fasta')
    )
# end def


@pytest.fixture(scope='session')
def test_outdir_path(tmpdir_factory) -> str:
    outdir_path: str = tmpdir_factory.mktemp('test_outdir')
    return outdir_path
# end def


@pytest.fixture(scope='session')
def test_outfpath(test_outdir_path) -> str:
    return os.path.join(test_outdir_path, 'test_annotated_seq.gbk')
# end def


@pytest.fixture(scope='session')
def test_coverage_fpath(test_outdir_path) -> str:
    return os.path.join(test_outdir_path, 'coverages.tsv')
# end def


@pytest.fixture(scope='session')
def nonzero_cov_threshold() -> LowerCoverageThreshold:
    return LowerCoverageThreshold(2)
# end def


@pytest.fixture(scope='session')
def upper_cov_threshold() -> UpperCoverageThreshold:
    return UpperCoverageThreshold(1, 1)
# end def


@pytest.fixture(scope='session')
def coverage_array_inner(test_bam_fpath: str,
                         test_coverage_fpath: str,
                         first_test_seq_id: str) -> CoverageArray:
    cov_fpath: str = oc.count_coverage(test_bam_fpath, first_test_seq_id, test_coverage_fpath)
    return oc.get_coverage_for_reference(cov_fpath)
# end def


@pytest.fixture(scope='session')
def coverage_array_edge(test_bam_fpath: str,
                        test_coverage_fpath: str,
                        second_test_seq_id: str) -> CoverageArray:
    cov_fpath: str = oc.count_coverage(test_bam_fpath, second_test_seq_id, test_coverage_fpath)
    return oc.get_coverage_for_reference(cov_fpath)
# end def


@pytest.fixture(scope='session')
def coverage_array_inner_upper(test_bam_fpath: str,
                               test_coverage_fpath: str,
                               third_test_seq_id: str) -> CoverageArray:
    cov_fpath: str = oc.count_coverage(test_bam_fpath, third_test_seq_id, test_coverage_fpath)
    return oc.get_coverage_for_reference(cov_fpath)
# end def
