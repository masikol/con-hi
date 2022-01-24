# Version 2.1.a

import os
import subprocess as sp
from typing import Tuple, Iterator

from src.printing import print_err
from src.platform import platf_depend_exit
from src.coverage_array import CoverageArray


def count_cov_for_all_refs(bam_fpath: str, coverage_fpath: str) -> str:
    # Function counts coverages for all sequences from input fasta file.
    # :param bam_fpath: path to bam file;
    # :param coverage_fpath: path to coverage file;
    # Returns path to file that contains coverage:
    #   first column -- sequence id, second column -- position, third column -- coverage.

    # Obtain command to run
    samtools_depth_cmd: str = _conf_samtools_depth_cmd(bam_fpath, coverage_fpath)

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
# end def count_cov_for_all_refs


def get_coverage_for_reference(sequence_id: str, coverage_fpath: str) -> CoverageArray:
    # Function retrieves coverage array for a single sequence from input fasta file.
    # :param sequence_id: id of sequence to retrieve coverages for;
    # :param coverage_fpath: path to coverage file;

    coverages: Iterator[int]

    # Read the coverage file
    with open(coverage_fpath, 'r') as cov_file:

        lines_of_curr_ref: Iterator[str]

        # Read lines and filter them: we need coverage only for current reference sequence
        lines_of_curr_ref = filter(
            lambda x: x.split('\t')[0] == sequence_id,
            cov_file.readlines()
        )
    # end with

    # Parse coverage
    coverages = map(
        lambda x: int(x.strip().split('\t')[2]),
        lines_of_curr_ref
    )

    # Convert `coverages` to CoverageArray
    cov_array: CoverageArray = CoverageArray(coverages)

    return cov_array
# end def get_coverage_for_reference


def _conf_samtools_depth_cmd(bam_fpath: str, coverage_fpath: str) -> str:
    # Function configures command to count coverages.
    # :param bam_fpath: path to bam file;
    # :param coverage_fpath: path to coverage file;
    return 'samtools depth -aa -J -o {} {}'\
        .format(coverage_fpath, bam_fpath)
# end def _conf_samtools_depth_cmd

