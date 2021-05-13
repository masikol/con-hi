# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
import subprocess as sp
from typing import Tuple, Iterator

from src.platform import platf_depend_exit
from src.coverage_array import CoverageArray

def count_cov_for_all_refs(ref_fasta_fpath: str, bam_fpath: str, outdir: str) -> str:

    coverage_fpath: str = _conf_path_to_depth_file(outdir)
    samtools_depth_cmd: str = _conf_samtools_depth_cmd(ref_fasta_fpath, bam_fpath, coverage_fpath)

    pipe: sp.Popen = sp.Popen(samtools_depth_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr: Tuple[bytes, bytes] = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError occured while running samtools depth')
        print(stdout_stderr[1].decode('utf-8'))
        platf_depend_exit(pipe.returncode)
    # end if

    return coverage_fpath
# end def count_cov_for_all_refs


def get_coverage_for_reference(sequence_id: str, coverage_fpath: str) -> CoverageArray:

    coverages: Iterator[int]

    with open(coverage_fpath, 'r') as cov_file:

        lines_of_curr_ref: Iterator[str]

        # Read lines and filter them: we need coverage only for current reference sequence
        lines_of_curr_ref = filter(
            lambda x: x.split('\t')[0] == sequence_id,
            cov_file.readlines()
        )

        coverages = map(
            lambda x: int(x.strip().split('\t')[2]),
            lines_of_curr_ref
        )
    # end with

    cov_array: CoverageArray = CoverageArray(coverages)

    return cov_array
# end def get_coverage_for_reference


def _conf_samtools_depth_cmd(ref_fasta_fpath: str, bam_fpath: str, coverage_fpath: str) -> str:
    return 'samtools depth -a -J --reference {} -o {} {}'\
        .format(ref_fasta_fpath, coverage_fpath, bam_fpath)
# end def _conf_samtools_depth_cmd


def _conf_path_to_depth_file(outdir: str) -> str:
    return os.path.join(outdir, 'coverages.tsv')
# end def _conf_path_to_depth_file
