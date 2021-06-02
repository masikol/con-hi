# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
import sys
import glob
from typing import Sequence, MutableSequence, List

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

import src.output as out
import src.obtain_coverage as oc
import src.highlight_features as hlft
import src.parse_fasta_reference as pfr
from src.platform import platf_depend_exit
from src.coverage_array import CoverageArray
from src.coverage_threshold import CoverageTheshold
from src.arguments import parse_arguments, HighlighterParams


def main(version: str, last_update_date: str) -> None:

    # Parse arguments
    params: HighlighterParams = parse_arguments()

    # This string will be used for annotation of result GenBank file
    base_feature_note: str = f'generated with consensus-highlighter v{version}'

    # Read fasta records from input file
    print('Importing fasta from `{}`...'.format(params.target_fasta_fpath), end=' ')
    sys.stdout.flush()
    fasta_records: Sequence[SeqRecord] = pfr.parse_fasta_reference(params.target_fasta_fpath)
    print('done')

    # Create ouput directory
    _create_outdir(params.outdir_path)

    # Count coverages with samtools depth
    print('Silently counting coverages with `samtools depth`...', end=' ')
    sys.stdout.flush()
    cov_fpath: str = oc.count_cov_for_all_refs(
        params.target_fasta_fpath,
        params.bam_fpath,
        params.outdir_path
    )
    print('done\n')

    # Proceed with annotation
    rec: SeqRecord
    for rec in fasta_records:

        print(f'Processing sequence `{rec.description}`')

        # Obtain coverages for current sequence
        cov_array: CoverageArray = oc.get_coverage_for_reference(rec.id, cov_fpath)

        cov_threshold: CoverageTheshold
        coverage_features: MutableSequence[SeqFeature]

        # Detect all necessary coverage features
        for cov_threshold in params.coverage_thresholds:

            print(f'Screening the sequence for regions with {cov_threshold.get_label()}...', end=' ')
            sys.stdout.flush()

            # Get coverage features
            coverage_features = hlft.highlight_coverage_features(
                cov_array,
                cov_threshold,
                base_feature_note
            )

            # Append features to list
            rec.features.extend(coverage_features)
            print('done')
        # end for

        # Configure path to output file
        outfpath: str = out.configure_outfile_path(
            params.outdir_path,
            params.outfile_prefix,
            rec.id
        )

        print(f'Writing annotated sequence to `{outfpath}`...', end=' ')
        sys.stdout.flush()

        # Write result GanBank record
        out.write_genbank_output(rec, params.topology, params.organism, outfpath)
        print('done')

        print('=' * 10)
    # end for

    print('Completed!')
# end def main


def _create_outdir(outdpath: str) -> None:
    # Function creates output directory
    # :param outdpath: path to directory to create;

    # Create directory if it does not exist
    if not os.path.isdir(outdpath):
        try:
            os.makedirs(outdpath)
        except OSError as err:
            print(f'Error! Cannot create output directory `{outdpath}`.')
            print(str(err))
            platf_depend_exit(1)
        # end try
    # end if

    dir_content_gbk: List[str] = glob.glob(os.path.join(outdpath, '*.gbk'))

    # Check if there are some .gbk files in the output directory
    #   and warn a user if there are some
    if len(dir_content_gbk) != 0:

        print(f'\nOutput directory `{outdpath}` contain some .gbk files.')
        print('They might be further overwritten by the program.')

        error: bool = True
        while error:

            print('Please, press ENTER to remove all .gbk files in this directory and proceed')
            reply: str = input('  OR enter `a` to accept this risk and just proceed: > ')

            if reply == '':
                # Remove all .gbk files from the directory
                fpath: str
                for fpath in dir_content_gbk:
                    print(f'Removing `{fpath}`')
                    try:
                        os.unlink(fpath)
                    except OSError as err:
                        print(f'Error: cannot remove `{fpath}`: {str(err)}.')
                    # end try
                # end for
                print()
                error = False
            elif reply.upper() == 'A':
                # Just proceed
                print('Proceeding\n')
                error = False
            else:
                print(f'Invalid reply: `{reply}`.\n')
            # end if
        # end while
    # end if
# end def _create_outdir
