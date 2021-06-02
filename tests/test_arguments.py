# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
from typing import List

import pytest

import src.arguments as args
from src.coverage_threshold import CoverageThreshold

from tests.fixtures import test_fasta_fpath, test_bam_fpath, test_coverage_fpath, test_outdir_path


@pytest.fixture
def some_params() -> args.HighlighterParams:
    return args.HighlighterParams(
        target_fasta_fpath=None,
        bam_fpath=None,
        outdir_path=os.path.join(os.getcwd(), 'consensus-highlighter-output'),
        coverage_thresholds=(10, 15, 20),
        suppress_zero_cov_output=False,
        topology='linear',
        organism='.',
        outfile_prefix=''
    )
# end def some_params


@pytest.fixture
def opts_all_is_ok_short_options(test_fasta_fpath, test_bam_fpath, test_outdir_path) -> List[List[str]]:
    return [
        ['-f', test_fasta_fpath],
        ['-b', test_bam_fpath],
        ['-o', test_outdir_path],
        ['-c', '10,20'],
        ['-n', ''],
    ]
# end def opts_all_is_ok_short_options

@pytest.fixture
def opts_all_is_ok_long_options(test_fasta_fpath, test_bam_fpath, test_outdir_path) -> List[List[str]]:
    return [
        ['--target-fasta', test_fasta_fpath],
        ['--bam', test_bam_fpath],
        ['--outdir', test_outdir_path],
        ['--coverage-thresholds', '10,20'],
        ['--no-zero-output', ''],
        ['--circular', ''],
        ['--organism', 'Czort lysy'],
        ['--prefix', 'some_prefix'],
    ]
# end def opts_all_is_ok_long_options


def test_is_fasta(test_fasta_fpath, test_bam_fpath, test_coverage_fpath) -> None:
    # Function for testing function `src.arguments._is_fasta`
    assert args._is_fasta(test_fasta_fpath)
    assert not args._is_fasta(test_bam_fpath)
    assert not args._is_fasta(test_coverage_fpath)
# end def test_is_fasta


def test_is_bam(test_fasta_fpath, test_bam_fpath, test_coverage_fpath) -> None:
    # Function for testing function `src.arguments._is_bam`
    assert args._is_bam(test_bam_fpath)
    assert not args._is_bam(test_fasta_fpath)
    assert not args._is_bam(test_coverage_fpath)
# end def test_is_fasta


def test_coverage_not_parsable() -> None:
    # Function for testing function `src.arguments._coverage_not_parsable`
    assert args._coverage_not_parsable('-1')
    assert args._coverage_not_parsable('0')
    assert not args._coverage_not_parsable('1')
    assert not args._coverage_not_parsable('12')
# end def test_coverage_not_parsable


def test_add_zero_theshold(some_params) -> None:
    # Function for testing function `src.arguments._add_zero_theshold`
    args._add_zero_theshold(some_params)

    expected_zero_coverage: int = 0
    assert some_params.coverage_thresholds[0].get_coverage() == expected_zero_coverage
# end def test_add_zero_theshold


class TestParseOptions:
    # Class for testing function `src.arguments._parse_options`

    def test_opts_all_is_ok_short_options(
        self,
        opts_all_is_ok_short_options) -> None:
        # Function tests how `_parse_options` parses correct short options

        params: args.HighlighterParams = args._parse_options(opts_all_is_ok_short_options)

        assert params.target_fasta_fpath == opts_all_is_ok_short_options[0][1]
        assert params.bam_fpath == opts_all_is_ok_short_options[1][1]
        assert params.outdir_path == opts_all_is_ok_short_options[2][1]
        assert params.suppress_zero_cov_output == True
        assert params.topology == 'linear'
        assert params.organism == '.'
        assert params.outfile_prefix == ''

        expected_threshold_repr: str = ','.join(
            map(
                lambda covthr: str(covthr.get_coverage()),
                params.coverage_thresholds
            )
        )
        assert expected_threshold_repr == opts_all_is_ok_short_options[3][1]
    # end def test_opts_all_is_ok_short_options


    def test_opts_all_is_ok_long_options(
        self,
        opts_all_is_ok_long_options) -> None:
        # Function tests how `_parse_options` parses correct long options

        params: args.HighlighterParams = args._parse_options(opts_all_is_ok_long_options)

        assert params.target_fasta_fpath == opts_all_is_ok_long_options[0][1]
        assert params.bam_fpath == opts_all_is_ok_long_options[1][1]
        assert params.outdir_path == opts_all_is_ok_long_options[2][1]
        assert params.suppress_zero_cov_output == True
        assert params.topology == 'circular'
        assert params.organism == opts_all_is_ok_long_options[6][1]
        assert params.outfile_prefix == opts_all_is_ok_long_options[7][1]

        expected_threshold_repr: str = ','.join(
            map(
                lambda covthr: str(covthr.get_coverage()),
                params.coverage_thresholds
            )
        )
        assert expected_threshold_repr == opts_all_is_ok_long_options[3][1]
    # end def opts_all_is_ok_long_options
# end class TestParseOptions
