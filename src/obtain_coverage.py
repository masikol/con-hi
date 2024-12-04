
import os
import subprocess as sp
from typing import Tuple, Iterator

from src.printing import print_err
from src.platform import platf_depend_exit
from src.coverage_array import CoverageArray


def count_coverage(bam_fpath: str,
                   target_seq_id: str,
                   coverage_fpath: str) -> str:
    # TODO: doc string
    # Function counts coverages for all sequences from input fasta file.
    # :param bam_fpath: path to bam file;
    # :param coverage_fpath: path to coverage file;
    # Returns path to file that contains coverage:
    #   first column -- sequence id, second column -- position, third column -- coverage.

    # Obtain command to run
    samtools_depth_cmd: str = _conf_samtools_depth_cmd(
        bam_fpath,
        target_seq_id,
        coverage_fpath,
    )

    # Run the command
    pipe: sp.Popen = sp.Popen(samtools_depth_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr: Tuple[bytes, bytes] = pipe.communicate()

    # Print error if it occures
    if pipe.returncode != 0:
        print_err('\nError occured while running `samtools depth`')
        print_err(stdout_stderr[1].decode('utf-8'))
        platf_depend_exit(pipe.returncode)
    # end if

    return coverage_fpath
# end def


def get_coverage_for_reference(coverage_fpath: str) -> CoverageArray:
    # Function retrieves coverage array for a single sequence from input fasta file.
    # :param coverage_fpath: path to coverage file;

    coverages: Iterator[int]

    # Read the coverage file
    with open(coverage_fpath, 'r') as cov_file:
        lines: Iterator[str] = tuple(
            map(
                lambda line: line.strip(),
                cov_file.readlines()
            )
        )
    # end with

    if len(lines) == 0:
        raise MissingCoveragesError()
    # end if

    # Parse coverage
    coverages = map(
        lambda x: int(x.strip().split('\t')[2]),
        lines
    )

    # Convert `coverages` to CoverageArray
    cov_array: CoverageArray = CoverageArray(coverages)

    return cov_array
# end def


def _conf_samtools_depth_cmd(bam_fpath: str,
                             target_seq_id: str,
                             coverage_fpath: str) -> str:
    # TODO: doc string
    # Function configures command to count coverages.
    # :param bam_fpath: path to bam file;
    # :param coverage_fpath: path to coverage file;

    return ' '.join(
        [
            'samtools',
            'depth',
            f'-r {target_seq_id}',
            '-a',
            '-J',
            f'-o {coverage_fpath}',
            bam_fpath,
        ]
    )
# end def


class MissingCoveragesError(Exception):
    pass
# end class
