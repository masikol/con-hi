
import os
import sys
from typing import List, Sequence

import pytest

import src.arguments as args

from tests.fixtures import test_outdir_path
from tests.fixtures import test_fasta_fpath, test_bam_fpath, test_coverage_fpath, test_outfpath


@pytest.fixture
def some_params() -> args.HighlighterParams:
    return args.HighlighterParams(
        target_fasta_fpath=None,
        bam_fpath=None,
        outfpath=os.path.join(os.getcwd(), 'highlighted_sequence.gbk'),
        lower_coverage_thresholds=(10, 15, 20),
        upper_coverage_coefficients=(1.5, 2.0),
        suppress_zero_cov_output=False,
        min_feature_len=5,
        topology='linear',
        organism='.'
    )
# end def


# === Fixtures for function `src.arguments._parse_options` ===

@pytest.fixture
def opts_all_is_ok_short_options(test_fasta_fpath, test_bam_fpath, test_outfpath) -> List[List[str]]:
    return [
        ['-f', test_fasta_fpath],
        ['-b', test_bam_fpath],
        ['-o', test_outfpath],
        ['-c', '10,20'],
        ['-C', '1.5,2.0'],
        ['-n', ''],
    ]
# end def

@pytest.fixture
def opts_all_is_ok_long_options(test_fasta_fpath, test_bam_fpath, test_outfpath) -> List[List[str]]:
    return [
        ['--target-fasta', test_fasta_fpath],
        ['--bam', test_bam_fpath],
        ['--outfile', test_outfpath],
        ['--lower-coverage-thresholds', '10,20'],
        ['--upper-coverage-coefficients', '1.5,2.0'],
        ['--no-zero-output', ''],
        ['--circular', ''],
        ['--organism', 'Czort lysy']
    ]
# end def

@pytest.fixture
def opts_invalid_fasta(test_fasta_fpath, test_bam_fpath) -> List[List[str]]:
    return [
        ['-f', 'TRASH' + test_fasta_fpath], # <- non-extant file
        ['-b', test_bam_fpath],
        ['-c', '10,20'],
    ]
# end def

@pytest.fixture
def opts_fasta_not_fasta(test_bam_fpath) -> List[List[str]]:
    return [
        ['-f', test_bam_fpath], # <- path to bam, not fasta
        ['-b', test_bam_fpath],
        ['-c', '10,20'],
    ]
# end def

@pytest.fixture
def opts_invalid_bam(test_fasta_fpath, test_bam_fpath) -> List[List[str]]:
    return [
        ['-f', test_fasta_fpath],
        ['-b', 'TRASH' + test_bam_fpath], # <- non-extant file
        ['-c', '10,20'],
    ]
# end def

@pytest.fixture
def opts_bam_not_bam(test_fasta_fpath) -> List[List[str]]:
    return [
        ['-f', test_fasta_fpath], # <- path to fasta, not bam
        ['-b', test_fasta_fpath],
        ['-c', '10,20'],
    ]
# end def

@pytest.fixture
def opts_NAN_threshold(test_fasta_fpath, test_bam_fpath) -> List[List[str]]:
    return [
        ['-f', test_fasta_fpath],
        ['-b', test_bam_fpath],
        ['-c', '10,2X'], # <- non-numeric threshold
    ]
# end def

@pytest.fixture
def opts_zero_threshold(test_fasta_fpath, test_bam_fpath) -> List[List[str]]:
    return [
        ['-f', test_fasta_fpath],
        ['-b', test_bam_fpath],
        ['-c', '0,20'], # <- zero threshold
    ]
# end def

@pytest.fixture
def opts_negative_threshold(test_fasta_fpath, test_bam_fpath) -> List[List[str]]:
    return [
        ['-f', test_fasta_fpath],
        ['-b', test_bam_fpath],
        ['-c', '-1,20'], # <- negative threshold
    ]
# end def

@pytest.fixture
def opts_keep_zero_thresahold(test_fasta_fpath, test_bam_fpath) -> List[List[str]]:
    return [
        ['-f', test_fasta_fpath],
        ['-b', test_bam_fpath],
        ['-c', '10,20'],
    ]
# end def

@pytest.fixture
def opts_suppress_zero_thresahold(test_fasta_fpath, test_bam_fpath) -> List[List[str]]:
    return [
        ['-f', test_fasta_fpath],
        ['-b', test_bam_fpath],
        ['-c', '10,20'],
        ['-n', '']
    ]
# end def


# === Fixtures for function `src.arguments.parse_arguments` ===

@pytest.fixture
def cmd_all_ok_short_options(test_fasta_fpath, test_bam_fpath, test_outfpath) -> Sequence[str]:
    return [
        'consensus-highlighter',
        '-f', test_fasta_fpath,
        '-b', test_bam_fpath,
        '-o', test_outfpath,
        '-c', '10,20',
        '-C', '1.5,2.0',
        '-n',
    ]
# end def

@pytest.fixture
def cmd_all_ok_long_options(test_fasta_fpath, test_bam_fpath, test_outfpath) -> Sequence[str]:
    return [
        'consensus-highlighter',
        '--target-fasta', test_fasta_fpath,
        '--bam', test_bam_fpath,
        '--outfile', test_outfpath,
        '--lower-coverage-thresholds', '10,20',
        '--upper-coverage-coefficients', '1.5,2.0',
        '--no-zero-output',
        '--circular',
        '--organism', 'Czort lysy'
    ]
# end def

@pytest.fixture
def cmd_positional_args(test_fasta_fpath, test_bam_fpath) -> Sequence[str]:
    return [
        'consensus-highlighter',
        '-f', test_fasta_fpath,
        '-b', test_bam_fpath,
        '-c', '10,', '20',
        #          ~~~
        #          imitates space between threshold values
    ]
# end def

@pytest.fixture
def cmd_fasta_missing(test_bam_fpath) -> Sequence[str]:
    return [
        'consensus-highlighter',
        '-b', test_bam_fpath,
    ]
# end def

@pytest.fixture
def cmd_bam_missing(test_fasta_fpath) -> Sequence[str]:
    return [
        'consensus-highlighter',
        '-f', test_fasta_fpath,
    ]
# end def


def test_is_fasta(test_fasta_fpath, test_bam_fpath, test_coverage_fpath) -> None:
    # Function for testing function `src.arguments._is_fasta`
    assert args._is_fasta(test_fasta_fpath)
    assert not args._is_fasta(test_bam_fpath)
    assert not args._is_fasta(test_coverage_fpath)
# end def


def test_is_bam(test_fasta_fpath, test_bam_fpath, test_coverage_fpath) -> None:
    # Function for testing function `src.arguments._is_bam`
    assert args._is_bam(test_bam_fpath)
    assert not args._is_bam(test_fasta_fpath)
    assert not args._is_bam(test_coverage_fpath)
# end def


def test_coverage_not_parsable() -> None:
    # Function for testing function `src.arguments._coverage_not_parsable`
    assert args._coverage_not_parsable('-1')
    assert args._coverage_not_parsable('0')
    assert not args._coverage_not_parsable('1')
    assert not args._coverage_not_parsable('12')
# end def


def test_add_zero_theshold(some_params) -> None:
    # Function for testing function `src.arguments._add_zero_theshold`
    args._add_zero_theshold(some_params)

    expected_zero_coverage: int = 0
    assert some_params.lower_coverage_thresholds[0] == expected_zero_coverage
# end def


class TestParseOptions:
    # Class for testing function `src.arguments._parse_options`

    def test_opts_all_is_ok_short_options(self, opts_all_is_ok_short_options) -> None:
        # Function tests how `_parse_options` parses correct short options

        params: args.HighlighterParams = args._parse_options(opts_all_is_ok_short_options)

        assert params.target_fasta_fpath == opts_all_is_ok_short_options[0][1]
        assert params.bam_fpath == opts_all_is_ok_short_options[1][1]
        assert params.outfpath == opts_all_is_ok_short_options[2][1]
        assert params.suppress_zero_cov_output == True
        assert params.topology == 'linear'
        assert params.organism == '.'

        obtained_threshold_repr: str = ','.join(
            map(str, params.lower_coverage_thresholds)
        )
        expected_threshold_repr: str = opts_all_is_ok_short_options[3][1]
        assert obtained_threshold_repr == expected_threshold_repr

        obtained_coef_repr: str = ','.join(
            map(str, params.upper_coverage_coefficients)
        )
        expected_coef_repr: str = opts_all_is_ok_short_options[4][1]
        assert obtained_coef_repr == expected_coef_repr
    # end def


    def test_opts_all_is_ok_long_options(self, opts_all_is_ok_long_options) -> None:
        # Function tests how `_parse_options` parses correct long options

        params: args.HighlighterParams = args._parse_options(opts_all_is_ok_long_options)

        assert params.target_fasta_fpath == opts_all_is_ok_long_options[0][1]
        assert params.bam_fpath == opts_all_is_ok_long_options[1][1]
        assert params.outfpath == opts_all_is_ok_long_options[2][1]
        assert params.suppress_zero_cov_output == True
        assert params.topology == 'circular'
        assert params.organism == opts_all_is_ok_long_options[7][1]

        obtained_threshold_repr: str = ','.join(
            map(str, params.lower_coverage_thresholds)
        )
        expected_threshold_repr: str = opts_all_is_ok_long_options[3][1]
        assert obtained_threshold_repr == expected_threshold_repr

        obtained_coef_repr: str = ','.join(
            map(str, params.upper_coverage_coefficients)
        )
        expected_coef_repr: str = opts_all_is_ok_long_options[4][1]
        assert obtained_coef_repr == expected_coef_repr
    # end def

    def test_opts_invalid_fasta(self, opts_invalid_fasta) -> None:
        # Function tests how `_parse_options` parses options
        #   where target fasta file does not exist

        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args._parse_options(opts_invalid_fasta)
        # end with
    # end def

    def test_opts_fasta_not_fasta(self, opts_fasta_not_fasta) -> None:
        # Function tests how `_parse_options` parses options
        #   where non-fasta file is specified as target fasta file

        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args._parse_options(opts_fasta_not_fasta)
        # end with
    # end def

    def test_opts_invalid_bam(self, opts_invalid_bam) -> None:
        # Function tests how `_parse_options` parses options
        #   where bam file does not exist

        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args._parse_options(opts_invalid_bam)
        # end with
    # end def

    def test_opts_bam_not_bam(self, opts_bam_not_bam) -> None:
        # Function tests how `_parse_options` parses options
        #   where non-bam file is specified as bam file

        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args._parse_options(opts_bam_not_bam)
        # end with
    # end def

    def test_opts_NAN_threshold(self, opts_NAN_threshold) -> None:
        # Function tests how `_parse_options` parses options
        #   where non-numeric threshold is specified

        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args._parse_options(opts_NAN_threshold)
        # end with
    # end def

    def test_opts_zero_threshold(self, opts_zero_threshold) -> None:
        # Function tests how `_parse_options` parses options
        #   where zero threshold is specified

        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args._parse_options(opts_zero_threshold)
        # end with
    # end def

    def test_opts_negative_threshold(self, opts_negative_threshold) -> None:
        # Function tests how `_parse_options` parses options
        #   where negative threshold is specified

        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args._parse_options(opts_negative_threshold)
        # end with
    # end def

    def test_opts_keep_zero_thresahold(self, opts_keep_zero_thresahold) -> None:
        # Function tests how `_parse_options` adds zero threshold

        params: args.HighlighterParams = args._parse_options(opts_keep_zero_thresahold)

        assert params.suppress_zero_cov_output == False

        zero_coverage: int = 0
        # Zero threshold should be first
        first_coverage_value: int = params.lower_coverage_thresholds[0]

        assert zero_coverage == first_coverage_value
    # end def

    def test_opts_suppress_zero_thresahold(self, opts_suppress_zero_thresahold) -> None:
        # Function tests how `_parse_options` suppresses zero threshold

        params: args.HighlighterParams = args._parse_options(opts_suppress_zero_thresahold)

        assert params.suppress_zero_cov_output == True

        zero_coverage: int = 0
        assert not zero_coverage in params.lower_coverage_thresholds
    # end def
# end class


class TestParseArguments:
    # Class for testing function `src.arguments.parse_arguments`

    def test_cmd_all_ok_short_options(self, cmd_all_ok_short_options) -> None:
        # Function tests how `parse_arguments` parses correct short options

        # Backup argv
        buff_argv: List[str] = sys.argv
        # Set test argv
        sys.argv = cmd_all_ok_short_options

        # Parse arguments
        params: args.HighlighterParams = args.parse_arguments()

        # Restore argv
        argv = buff_argv

        assert params.target_fasta_fpath == cmd_all_ok_short_options[2]
        assert params.bam_fpath == cmd_all_ok_short_options[4]
        assert params.outfpath == cmd_all_ok_short_options[6]
        assert params.suppress_zero_cov_output == True
        assert params.topology == 'linear'
        assert params.organism == '.'

        obtained_threshold_repr: str = ','.join(
            map(str, params.lower_coverage_thresholds)
        )
        expected_threshold_repr: str = cmd_all_ok_short_options[8]
        assert obtained_threshold_repr == expected_threshold_repr

        obtained_coef_repr: str = ','.join(
            map(str, params.upper_coverage_coefficients)
        )
        expected_coef_repr: str = cmd_all_ok_short_options[10]
        assert obtained_coef_repr == expected_coef_repr
    # end def

    def test_cmd_all_ok_long_options(self, cmd_all_ok_long_options) -> None:
        # Function tests how `parse_arguments` parses correct long options

        # Backup argv
        buff_argv: List[str] = sys.argv
        # Set test argv
        sys.argv = cmd_all_ok_long_options

        # Parse arguments
        params: args.HighlighterParams = args.parse_arguments()

        # Restore argv
        argv = buff_argv

        assert params.target_fasta_fpath == cmd_all_ok_long_options[2]
        assert params.bam_fpath == cmd_all_ok_long_options[4]
        assert params.outfpath == cmd_all_ok_long_options[6]
        assert params.suppress_zero_cov_output == True
        assert params.topology == 'circular'
        assert params.organism == cmd_all_ok_long_options[14]

        obtained_threshold_repr: str = ','.join(
            map(str, params.lower_coverage_thresholds)
        )
        expected_threshold_repr: str = cmd_all_ok_long_options[8]
        assert obtained_threshold_repr == expected_threshold_repr

        obtained_coef_repr: str = ','.join(
            map(str, params.upper_coverage_coefficients)
        )
        expected_coef_repr: str = cmd_all_ok_long_options[10]
        assert obtained_coef_repr == expected_coef_repr
    # end def

    def test_cmd_positional_args(self, cmd_positional_args) -> None:
        # Function tests how `parse_arguments` parses command line
        #    with positional arguments

        # Backup argv
        buff_argv: List[str] = sys.argv
        # Set test argv
        sys.argv = cmd_positional_args

        # Parse arguments
        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args.parse_arguments()
        # end with

        # Restore argv
        argv = buff_argv
    # end def

    def test_cmd_fasta_missing(self, cmd_fasta_missing) -> None:
        # Function tests how `parse_arguments` parses command line
        #    with positional arguments

        # Backup argv
        buff_argv: List[str] = sys.argv
        # Set test argv
        sys.argv = cmd_fasta_missing

        # Parse arguments
        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args.parse_arguments()
        # end with

        # Restore argv
        argv = buff_argv
    # end def

    def test_cmd_bam_missing(self, cmd_bam_missing) -> None:
        # Function tests how `parse_arguments` parses command line
        #    with positional arguments

        # Backup argv
        buff_argv: List[str] = sys.argv
        # Set test argv
        sys.argv = cmd_bam_missing

        # Parse arguments
        with pytest.raises(SystemExit):
            params: args.HighlighterParams = args.parse_arguments()
        # end with

        # Restore argv
        argv = buff_argv
    # end def
# end class
