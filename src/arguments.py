
import os
import re
import sys
import getopt
from typing import List, Sequence, Iterable

from src.printing import print_err
from src.platform import platf_depend_exit
from src.coverage_threshold import CoverageThreshold


class HighlighterParams:
    # Class-container for storing run parameters

    def __init__(
        self,
        target_fasta_fpath: str,
        bam_fpath: str,
        outfpath: str,
        coverage_thresholds: Sequence[CoverageThreshold],
        suppress_zero_cov_output: bool,
        min_feature_len: int,
        topology: str,
        organism: str) -> None:

        self.target_fasta_fpath: str = target_fasta_fpath
        self.bam_fpath: str = bam_fpath
        self.outfpath: str = outfpath

        self.suppress_zero_cov_output = suppress_zero_cov_output
        self.min_feature_len = min_feature_len

        self.set_coverage_thresholds(coverage_thresholds)

        self.topology = topology
        self.organism: str = organism
    # end def __init__


    def set_coverage_thresholds(self, coverage_thresholds: Sequence[int]) -> None:
        # Function sets coverage thresholds for the instance.
        # :param coverage_thresholds: collection of coverage thresholds;

        self.coverage_thresholds: Sequence[CoverageThreshold]

        # Set the thresholds
        self.coverage_thresholds = tuple(
            (CoverageThreshold(cov) for cov in coverage_thresholds)
        )
    # end def set_coverage_thresholds


    def __repr__(self) -> str:
        return f"""HighlighterParams(
  target_fasta_fpath: `{self.target_fasta_fpath}`
  bam_fpath: `{self.bam_fpath}`
  outfpath: `{self.outfpath}`
  coverage_thresholds: {self.coverage_thresholds}
  suppress_zero_cov_output: {self.suppress_zero_cov_output}
  topology: {self.topology}
  organism: `{self.organism}`
)"""
    # end def __repr__

# end class HighlighterParams


def parse_arguments() -> HighlighterParams:
    # Function parses command line arguments and produces
    #   container with program parameters in it.

    # Parse arguments with getopt
    opts: List[List[str]]
    args: List[str]
    try:
        opts, args = getopt.gnu_getopt(
            sys.argv[1:],
            'hvf:b:o:c:l:n',
            [
                'help',
                'version',
                'target-fasta=',
                'bam=',
                'outfile=',
                'coverage-thresholds=',
                'min-feature-len=',
                'no-zero-output',
                'circular',
                'organism='
            ]
        )
    except getopt.GetoptError as err:
        print_err(str(err))
        platf_depend_exit(2)
    # end try

    # Check positional arguments: their existance is an error signal
    if len(args) != 0:
        print_err('Error: consensus-highlighter.py does not take any positional arguments.')
        print_err('You passed following positional argument(s):')
        arg: str
        for arg in args:
            print_err(f'  `{arg}`')
        # end for
        platf_depend_exit(2)
    # end if

    # Parse options
    params: HighlighterParams = _parse_options(opts)

    # Check mandatory options
    if params.target_fasta_fpath is None:
        print_err('Error: option `-f` (`--target-fasta`) is mandatory.')
        platf_depend_exit(2)
    # end if

    if params.bam_fpath is None:
        print_err('Error: option `-b` (`--bam`) is mandatory.')
        platf_depend_exit(2)
    # end if

    return params
# end def parse_arguments


def _parse_options(opts: List[List[str]]) -> HighlighterParams:
    # Function parses program options

    # Initialize run parameters with default values
    params: HighlighterParams = HighlighterParams(
        target_fasta_fpath=None,
        bam_fpath=None,
        outfpath=os.path.join(os.getcwd(), 'highlighted_sequence.gbk'),
        coverage_thresholds=(10,),
        suppress_zero_cov_output=False,
        min_feature_len=5,
        topology='linear',
        organism='.'
    )

    # Parse command line options
    opt: str
    arg: str
    for opt, arg in opts:

        # Target fasta file
        if opt in ('-f', '--target-fasta'):

            if not os.path.isfile(arg):
                print_err(f'\aError: file {arg}` does not exist.')
                platf_depend_exit(2)
            # end if

            if not _is_fasta(arg):
                print_err('\aError: only plain fasta or gzipped fasta are supported.')
                print_err(f'Erroneous file: `{arg}`.')
                print_err('Allowed extentions: `.fasta`, `.fa`, `.fasta.gz`, `.fa.gz`.')
                platf_depend_exit(2)
            # end if

            params.target_fasta_fpath = arg

        # BAM file
        elif opt in ('-b', '--bam'):

            if not os.path.isfile(arg):
                print_err(f'\aError: file {arg}` does not exist.')
                platf_depend_exit(2)
            # end if

            if not _is_bam(arg):
                print_err(f'\aError: file `{arg}` does not seem like a BAM file.')
                print_err('The program only accepts BAM files having `.bam` extentions')
                platf_depend_exit(2)
            # end if

            params.bam_fpath = arg

        # Output directory
        elif opt in ('-o', '--outfile'):
            params.outfpath = os.path.abspath(arg)

        # List of coverage thesholds
        elif opt in ('-c', '--coverage-thresholds'):

            cov_strings: Sequence[str] = arg.split(',')

            if any(map(_coverage_not_parsable, cov_strings)):

                invalid_strings: Iterable[str] = filter(_coverage_not_parsable, cov_strings)

                print_err(f'\aError: invalid coverage thresholds in `{arg}`:')
                for s in invalid_strings:
                    print_err(f'  `{s}`')
                # end for
                print_err('Coverage thesholds must be positive integer numbers.')
                platf_depend_exit(2)
            # end if

            # Now cast coverages to int and sort them in ascending order
            coverage_thresholds: Sequence[int] = tuple(
                (int(cov_str) for cov_str in sorted(cov_strings, key=int))
            )

            params.set_coverage_thresholds(coverage_thresholds)

        # Repress zero output
        elif opt in ('-n', '--no-zero-output'):
            params.suppress_zero_cov_output = True

        elif opt in ('-l', '--min-feature-len'):
            try:
                min_feature_len = int(arg)
                if min_feature_len < 0:
                    raise ValueError
                # end if
            except ValueError:
                print_err(f'\aError: invalid minimum feature length: `{arg}`')
                print_err('It must be non-negative negative number.')
                platf_depend_exit(2)
            # end try
            params.min_feature_len = min_feature_len

        # Molecule topology for annotation
        elif opt == '--circular':
            params.topology = 'circular'

        # Organism name for annotation
        elif opt == '--organism':
            params.organism = arg
    # end for

    # Add zero coverage threshold, if no suppression is specified
    if not params.suppress_zero_cov_output:
        _add_zero_theshold(params)
    # end if

    return params
# end def _parse_options


def _is_fasta(fpath: str) -> bool:
    fasta_pattern: str = r'.+\.f(ast)?a(\.gz)?'
    return not re.match(fasta_pattern, fpath) is None
# end def _is_fasta


def _is_bam(fpath: str) -> bool:
    bam_pattern: str = r'.+\.bam'
    return not re.match(bam_pattern, fpath) is None
# end def _is_bam


def _add_zero_theshold(params: HighlighterParams) -> None:
    # Function adds zero threshold to lost of coverage thresholds
    # :param params: program parameters;

    curr_coverage_thresholds: Sequence[int] = tuple(
        map(
            lambda x: x.get_coverage(),
            params.coverage_thresholds
        )
    )

    coverage_thresholds_with_zero: Sequence[int] = (0,) + curr_coverage_thresholds

    params.set_coverage_thresholds(coverage_thresholds_with_zero)
# end def _add_zero_theshold


def _coverage_not_parsable(string: str) -> bool:
    # Function checks validity of coverage threshold strign representation
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
# end def _coverage_not_parsable
