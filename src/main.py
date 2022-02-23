# Version 2.2.b

import os
import sys
import glob
from typing import Sequence, MutableSequence, List

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from src.printing import print_err
import src.output as out
import src.obtain_coverage as oc
import src.dedupl_features as ddf
import src.highlight_features as hlft
import src.parse_fasta_reference as pfr
from src.platform import platf_depend_exit
from src.coverage_array import CoverageArray
from src.coverage_threshold import CoverageThreshold
from src.arguments import parse_arguments, HighlighterParams


def main(version: str, last_update_date: str) -> None:

    # Parse arguments
    params: HighlighterParams = parse_arguments()

    # This string will be used for annotation of result GenBank file
    base_feature_note: str = f'generated by consensus-highlighter v{version}'

    # String for storing info about warnings
    with_warnings: str = ''

    # Read fasta records from input file
    print('Importing fasta from `{}`...'.format(params.target_fasta_fpath), end=' ')
    sys.stdout.flush()
    fasta_records: Sequence[SeqRecord] = pfr.parse_fasta_reference(params.target_fasta_fpath)
    print('done')

    # Create ouput directory
    _create_outdir_from_outfile(params.outfpath)
    out.create_or_emply_file(params.outfpath)

    # Obtain path to coverage file
    coverage_fpath: str = out.conf_path_to_depth_file(params.outfpath)

    # Count coverages with samtools depth
    print('Silently counting coverages with `samtools depth`...', end=' ')
    sys.stdout.flush()
    cov_fpath: str = oc.count_cov_for_all_refs(
        params.bam_fpath,
        coverage_fpath
    )
    print('done\n')

    # Proceed with annotation
    rec: SeqRecord
    for rec in fasta_records:

        print(f'Processing sequence `{rec.description}`')

        # Obtain coverages for current sequence
        cov_array: CoverageArray = oc.get_coverage_for_reference(rec.id, cov_fpath)

        # Check length of the coverage array
        if len(cov_array) == 0:
            print_err(f'!  Warning: no coverage information found for sequence `{rec.id}`.')
            print_err(f"""!  Please, make sure that field `RNAME` (3-rd column) in your BAM file contains
!    id of this sequence specified in fasta header (i.e. `{rec.id}`).""")
            print_err('! Omitting this sequence.')
            print_err('=' * 10)
            with_warnings = ' with warnings'
            continue
        # end if

        if len(cov_array) != len(rec.seq):
            print_err(f"""!  Warning: length of sequence `{rec.id}` ({len(rec.seq)} bp)
!    is not equal to number of coverage positions ({len(cov_array)}) reported by `samtools depth`
!    and stored in coverage file `{cov_fpath}`.""")
            print_err('!  Re-creating the bam file might be the solution of this issue.')
            print_err('!  Omitting this sequence.')
            print_err('=' * 10)
            with_warnings = ' with warnings'
            continue
        # end if

        cov_threshold: CoverageThreshold
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

            coverage_features = ddf.dedupl_features(coverage_features, rec.features)

            # Append features to list
            rec.features.extend(coverage_features)
            print('done')
        # end for

        print(f'Writing annotated sequence to `{params.outfpath}`...', end=' ')
        sys.stdout.flush()

        # Write result GanBank record
        out.write_genbank_output(
            rec,
            params.topology,
            params.organism,
            cov_array,
            params.outfpath
        )
        print('done')

        print('=' * 10)
    # end for

    _try_rm_temp_file(coverage_fpath)
    print(f'Completed{with_warnings}!')
# end def main


def _create_outdir_from_outfile(outfpath: str) -> None:
    # Function creates output directory
    # :param outfpath: path to output file;

    outdpath = os.path.dirname(outfpath)

    # Create directory if it does not exist
    if not os.path.isdir(outdpath):
        try:
            os.makedirs(outdpath)
        except OSError as err:
            print_err(f'Error! Cannot create output directory `{outdpath}`.')
            print_err(str(err))
            platf_depend_exit(1)
        # end try
    # end if
# end def _create_outdir_from_outfile


def _try_rm_temp_file(file_path):
    try:
        os.unlink(file_path)
    except OSError as err:
        print_err('Warning: cannot remove temporary file `{}`'.format(file_path))
        print_err(str(err))
    # end try
# end def
