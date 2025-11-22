
import os
import re
import sys
import getopt
import logging
from typing import List, Sequence, Iterable, Set

from src.platform import platf_depend_exit


_ARG_SEP = ','


class HighlighterArgs:
    # Class-container for storing run parameters

    def __init__(
        self,
        target_fasta_fpath: str,
        target_seq_ids: Set[str],
        bam_fpath: str,
        outfpath: str,
        lower_coverage_thresholds: Sequence[int],
        upper_coverage_thresholds: Sequence[int],
        upper_coverage_coefficients: Sequence[float],
        disable_zero_cov_output: bool,
        min_feature_len: int,
        topology: str,
        organism: str,
        keep_tmp_cov_file: bool) -> None:

        self.target_fasta_fpath: str = target_fasta_fpath
        self.target_seq_ids: str = target_seq_ids
        self.bam_fpath: str = bam_fpath
        self.outfpath: str = outfpath

        self.disable_zero_cov_output = disable_zero_cov_output
        self.min_feature_len = min_feature_len

        self.lower_coverage_thresholds = lower_coverage_thresholds
        self.upper_coverage_thresholds = upper_coverage_thresholds
        self.upper_coverage_coefficients = upper_coverage_coefficients

        self.topology = topology
        self.organism: str = organism
        self.keep_tmp_cov_file = keep_tmp_cov_file
    # end def


    def __repr__(self) -> str:
        return f"""HighlighterArgs(
  target_fasta_fpath: `{self.target_fasta_fpath}`
  target_seq_ids: `{self.target_seq_ids}`
  bam_fpath: `{self.bam_fpath}`
  outfpath: `{self.outfpath}`
  lower_coverage_thresholds: {self.lower_coverage_thresholds}
  upper_coverage_thresholds: {self.upper_coverage_thresholds}
  upper_coverage_coefficients: {self.upper_coverage_coefficients}
  disable_zero_cov_output: {self.disable_zero_cov_output}
  topology: {self.topology}
  organism: `{self.organism}`
  keep_tmp_cov_file: {self.keep_tmp_cov_file}
)"""
    # end def

# end class


def parse_arguments() -> HighlighterArgs:
    # Function parses command line arguments and produces
    #   container with program parameters in it.

    # Parse arguments with getopt
    opts: List[List[str]]
    args: List[str]
    try:
        opts, args = getopt.gnu_getopt(
            sys.argv[1:],
            'hvf:r:b:o:c:C:X:l:nk',
            [
                'help',
                'version',
                'target-fasta=',
                'target-seq-ids=',
                'bam=',
                'outfile=',
                'lower-coverage-thresholds=',
                'upper-coverage-thresholds=',
                'upper-coverage-coefficients=',
                'min-feature-len=',
                'no-zero-output',
                'circular',
                'organism=',
                'keep-temp-cov-file',
            ]
        )
    except getopt.GetoptError as err:
        logging.error(str(err))
        platf_depend_exit(2)
    # end try

    # Check positional arguments: their existance is an error signal
    if len(args) != 0:
        logging.error('Error: con-hi.py does not take any positional arguments.')
        logging.error('You passed following positional argument(s):')
        arg: str
        for arg in args:
            logging.error(f'`{arg}`')
        # end for
        platf_depend_exit(2)
    # end if

    # Parse options
    args: HighlighterArgs = _parse_options(opts)

    # Check mandatory options
    if args.target_fasta_fpath is None:
        logging.error('Error: option `-f` (`--target-fasta`) is mandatory.')
        platf_depend_exit(2)
    # end if

    if args.bam_fpath is None:
        logging.error('Error: option `-b` (`--bam`) is mandatory.')
        platf_depend_exit(2)
    # end if

    return args
# end def


def _parse_options(opts: List[List[str]]) -> HighlighterArgs:
    # Function parses program options

    # Initialize run parameters with default values
    args: HighlighterArgs = HighlighterArgs(
        target_fasta_fpath=None,
        target_seq_ids=set(),
        bam_fpath=None,
        outfpath=os.path.join(os.getcwd(), 'highlighted_sequence.gbk'),
        lower_coverage_thresholds=(10,),
        upper_coverage_thresholds=(500,),
        upper_coverage_coefficients=(2.0,),
        disable_zero_cov_output=False,
        min_feature_len=5,
        topology='linear',
        organism='.',
        keep_tmp_cov_file = False
    )

    # Parse command line options
    opt: str
    arg: str
    for opt, arg in opts:

        # Target fasta file
        if opt in ('-f', '--target-fasta'):

            if not os.path.isfile(arg):
                logging.error(f'\aError: file {arg}` does not exist.')
                platf_depend_exit(2)
            # end if

            if not _is_fasta(arg):
                logging.error('\aError: only plain fasta or gzipped fasta are supported.')
                logging.error(f'Erroneous file: `{arg}`.')
                logging.error('Allowed extentions: `.fasta`, `.fa`, `.fasta.gz`, `.fa.gz`.')
                platf_depend_exit(2)
            # end if

            args.target_fasta_fpath = arg

        # Target reference sequence IDs
        elif opt in ('-r', '--target-seq-ids'):
            args.target_seq_ids = set(arg.split(_ARG_SEP))

        # BAM file
        elif opt in ('-b', '--bam'):

            if not os.path.isfile(arg):
                logging.error(f'\aError: file {arg}` does not exist.')
                platf_depend_exit(2)
            # end if

            if not _is_bam(arg):
                logging.error(f'\aError: file `{arg}` does not seem like a BAM file.')
                logging.error('The program only accepts BAM files having `.bam` extentions')
                platf_depend_exit(2)
            # end if

            args.bam_fpath = arg

        # Output directory
        elif opt in ('-o', '--outfile'):
            args.outfpath = os.path.abspath(arg)

        # List of lower coverage thesholds
        elif opt in ('-c', '--lower-coverage-thresholds'):
            if arg == 'off':
                args.lower_coverage_thresholds = list()
                continue
            # end if

            cov_strings: Sequence[str] = arg.split(_ARG_SEP)

            if any(map(_coverage_not_parsable, cov_strings)):
                invalid_strings: Iterable[str] = filter(_coverage_not_parsable, cov_strings)
                logging.error(f'\aError: invalid coverage thresholds in `{arg}`:')
                for s in invalid_strings:
                    logging.error(f'  `{s}`')
                # end for
                logging.error('Coverage thesholds must be positive integer numbers.')
                platf_depend_exit(2)
            # end if

            # Now cast coverages to int and sort them in ascending order
            coverage_thresholds: Sequence[int] = tuple(
                (int(cov_str) for cov_str in sorted(cov_strings, key=int))
            )

            args.lower_coverage_thresholds = coverage_thresholds

        # List of upper coverage thesholds
        elif opt in ('-C', '--upper-coverage-thresholds'):
            if arg == 'off':
                args.upper_coverage_thresholds = list()
                continue
            # end if

            cov_strings: Sequence[str] = arg.split(_ARG_SEP)

            if any(map(_coverage_not_parsable, cov_strings)):
                invalid_strings: Iterable[str] = filter(_coverage_not_parsable, cov_strings)
                logging.error(f'\aError: invalid coverage thresholds in `{arg}`:')
                for s in invalid_strings:
                    logging.error(f'  `{s}`')
                # end for
                logging.error('Coverage thesholds must be positive integer numbers.')
                platf_depend_exit(2)
            # end if

            # Now cast coverages to int and sort them in ascending order
            coverage_thresholds: Sequence[int] = tuple(
                (int(cov_str) for cov_str in sorted(cov_strings, key=int))
            )

            args.upper_coverage_thresholds = coverage_thresholds

        # List of upper coverage coefficients
        elif opt in ('-X', '--upper-coverage-coefficients'):
            if arg == 'off':
                args.upper_coverage_coefficients = list()
                continue
            # end if

            coef_strings: Sequence[str] = arg.split(_ARG_SEP)

            if any(map(_coefficient_not_parsable, coef_strings)):
                invalid_strings: Iterable[str] = filter(_coefficient_not_parsable, coef_strings)
                logging.error(f'\aError: invalid coverage coefficient in `{arg}`:')
                for s in invalid_strings:
                    logging.error(f'`{s}`')
                # end for
                logging.error('Coverage coefficients must be positive numbers.')
                platf_depend_exit(2)
            # end if

            # Now cast coverages to int and sort them in ascending order
            coverage_coefficients: Sequence[float] = tuple(
                (float(coef_str) for coef_str in sorted(coef_strings, key=float))
            )

            args.upper_coverage_coefficients = coverage_coefficients

        # Repress zero output
        elif opt in ('-n', '--no-zero-output'):
            args.disable_zero_cov_output = True

        # Minimum output feature length
        elif opt in ('-l', '--min-feature-len'):
            try:
                min_feature_len = int(arg)
                if min_feature_len < 0:
                    raise ValueError
                # end if
            except ValueError:
                logging.error(f'\aError: invalid minimum feature length: `{arg}`')
                logging.error('It must be non-negative negative number.')
                platf_depend_exit(2)
            # end try
            args.min_feature_len = min_feature_len

        # Molecule topology for annotation
        elif opt == '--circular':
            args.topology = 'circular'

        # Organism name for annotation
        elif opt == '--organism':
            args.organism = arg

        elif opt in ('-k', '--keep-temp-cov-file'):
            args.keep_tmp_cov_file = True
    # end for

    # Add zero coverage threshold, if no disabling is specified
    if not args.disable_zero_cov_output and len(args.lower_coverage_thresholds) != 0:
        _add_zero_theshold(args)
    # end if

    return args
# end def


def _is_fasta(fpath: str) -> bool:
    fasta_pattern: str = r'.+\.f(ast)?a(\.gz)?'
    return not re.match(fasta_pattern, fpath) is None
# end def


def _is_bam(fpath: str) -> bool:
    bam_pattern: str = r'.+\.bam'
    return not re.match(bam_pattern, fpath) is None
# end def


def _add_zero_theshold(args: HighlighterArgs) -> None:
    # Function adds zero threshold to lost of coverage thresholds
    # :param args: program parameters;

    coverage_thresholds_with_zero: Sequence[int] = (0,) + args.lower_coverage_thresholds

    args.lower_coverage_thresholds = coverage_thresholds_with_zero
# end def


def _coverage_not_parsable(string: str) -> bool:
    # Function checks validity of coverage threshold string representation
    # :param string: string to validate;

    cov_is_not_parsable: bool = True

    try:
        cov: int = int(string)
        if cov < 1:
            raise ValueError
        # end if
    except ValueError:
        cov_is_not_parsable = True
    else:
        cov_is_not_parsable = False
    # end try

    return cov_is_not_parsable
# end def


def _coefficient_not_parsable(string: str) -> bool:
    # Function checks validity of coverage coefficient string representation
    # :param string: string to validate;

    coef_is_not_parsable: bool = True

    try:
        ceof: float = float(string)
        if ceof < 1:
            raise ValueError
        # end if
    except ValueError:
        coef_is_not_parsable = True
    else:
        coef_is_not_parsable = False
    # end try

    return coef_is_not_parsable
# end def
