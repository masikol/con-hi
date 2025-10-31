# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
import re
from typing import Sequence

import src.obtain_coverage as oc
from src.coverage_array import CoverageArray

from tests.fixtures import test_outdir_path
from tests.fixtures import test_fasta_fpath, test_bam_fpath, test_coverage_fpath, test_outfpath
from tests.fixtures import first_test_seq_id, first_test_seq_len, \
                           second_test_seq_id, second_test_seq_len, \
                           third_test_seq_id, third_test_seq_len


class TestConfSamtoolsDepthCmd:
    # Class for testing function `src.obtain_coverage._conf_samtools_depth_cmd`

    def test_samtools_depth_cmd(
            self,
            test_coverage_fpath,
            first_test_seq_id) -> None:
        # Function tests that function `_conf_samtools_depth_cmd` returns correct command
        ref_fasta_fpath: str = 'some_file.fasta'
        bam_fpath: str = 'mapping.sorted.bam'

        expected: str = f'samtools depth -r {first_test_seq_id} -a -J -g 256 -o {test_coverage_fpath} {bam_fpath}'
        observed: str = oc._conf_samtools_depth_cmd(
            bam_fpath,
            first_test_seq_id,
            test_coverage_fpath
        )

        assert observed == expected
    # end def
# end class


class TestCountCovForAllRefs:
    # Class for testing function `src.obtain_coverage.count_coverage`

    def test_coverage_counting(
        self,
        test_fasta_fpath,
        test_bam_fpath,
        test_coverage_fpath,
        first_test_seq_id,
        first_test_seq_len,
        second_test_seq_id,
        second_test_seq_len,
        third_test_seq_id,
        third_test_seq_len) -> None:

        fixture_zip = zip(
            (first_test_seq_id,  second_test_seq_id,  third_test_seq_id),
            (first_test_seq_len, second_test_seq_len, third_test_seq_len)
        )
        for seq_id, seq_len in fixture_zip:
            # Function tests that `count_coverage` counts coverage in a plausable way
            cov_fpath: str = oc.count_coverage(
                test_bam_fpath,
                seq_id,
                test_coverage_fpath
            )

            # Check if file exists
            assert os.path.exists(cov_fpath) == True

            # Check size (in lines) of file's contents
            cov_lines: Sequence[str]
            with open(cov_fpath, 'r') as cov_file:
                cov_lines = tuple(map(str.strip, cov_file.readlines()))
            # end with
            assert len(cov_lines) == seq_len

            # Check format of all lines
            cov_line_pattern: str = '.+\t[0-9]+\t[0-9]+'
            def conforms_cov_line_pattern(string: str) -> bool:
                return not re.match(cov_line_pattern, string) is None
            # end def
            assert all(map(conforms_cov_line_pattern, cov_lines))
        # end if
    # end def
# end class


class TestGetCoverageForReference:
    # Class for testing function `src.obtain_coverage.get_coverage_for_reference`

    def test_get_covs_seq_1(
        self,
        first_test_seq_id,
        test_bam_fpath,
        test_coverage_fpath,
        first_test_seq_len) -> None:

        # Funciton tests that `get_coverage_for_reference` retrieves coverage correctly
        #   for the first test sequence

        cov_fpath: str = oc.count_coverage(
            test_bam_fpath,
            first_test_seq_id,
            test_coverage_fpath
        )

        cov_array: CoverageArray = oc.get_coverage_for_reference(
            test_coverage_fpath
        )

        assert len(cov_array) == first_test_seq_len
    # end def

    def test_get_covs_seq_2(
        self,
        second_test_seq_id,
        test_bam_fpath,
        test_coverage_fpath,
        second_test_seq_len) -> None:

        # Funciton tests that `get_coverage_for_reference` retrieves coverage correctly
        #   for the second test sequence

        cov_fpath: str = oc.count_coverage(
            test_bam_fpath,
            second_test_seq_id,
            test_coverage_fpath
        )

        cov_array: CoverageArray = oc.get_coverage_for_reference(
            test_coverage_fpath
        )

        assert len(cov_array) == second_test_seq_len
    # end def
# end class
