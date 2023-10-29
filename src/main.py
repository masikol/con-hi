
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
from src.filter_features import filter_features
from src.arguments import parse_arguments, HighlighterArgs
from src.coverage_thresholds.lower_coverage_threshold import LowerCoverageThreshold
from src.coverage_thresholds.upper_coverage_threshold import UpperCoverageThreshold


def main(version: str, last_update_date: str) -> None:

    # Parse arguments
    args: HighlighterArgs = parse_arguments()

    # This string will be used for annotation of result GenBank file
    base_feature_note: str = f'generated by consensus-highlighter v{version}'

    # String for storing info about warnings
    with_warnings: str = ''

    # Read fasta records from input file
    print('Importing fasta from `{}`...'.format(args.target_fasta_fpath), end=' ')
    sys.stdout.flush()
    fasta_records: Sequence[SeqRecord] = pfr.parse_fasta_reference(args.target_fasta_fpath)
    print('done')

    # Create ouput directory
    _create_outdir_from_outfile(args.outfpath)
    out.create_or_emply_file(args.outfpath)

    # Obtain path to coverage file
    coverage_fpath: str = out.conf_path_to_depth_file(args.outfpath)

    # Count coverages with samtools depth
    print('Silently counting coverages with `samtools depth`...', end=' ')
    sys.stdout.flush()
    cov_fpath: str = oc.count_cov_for_all_refs(
        args.bam_fpath,
        coverage_fpath
    )
    print('done\n')

    # Proceed with annotation
    rec: SeqRecord
    for rec in fasta_records:

        print(f'Processing sequence `{rec.description}`')

        # Obtain coverages for current sequence
        try:
            cov_array: CoverageArray = oc.get_coverage_for_reference(rec.id, cov_fpath)
        except oc.MissingCoveragesError:
            print_err('Warning: coverage values for this sequence are not found\n')
            continue
        # end try

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
        cov_thresholds: Sequence[cov_threshold: CoverageThreshold]
        coverage_features: MutableSequence[SeqFeature]

        cov_thresholds = _make_coverage_thresholds(args, cov_array)

        # Detect all necessary coverage features
        for cov_threshold in cov_thresholds:
            sys.stdout.write(f'Screening the sequence for regions with {cov_threshold.get_label()}...')
            sys.stdout.flush()

            # Get coverage features
            coverage_features = hlft.highlight_coverage_features(
                cov_array,
                cov_threshold,
                base_feature_note
            )

            coverage_features = filter_features(coverage_features, args.min_feature_len)
            coverage_features = ddf.dedupl_features(coverage_features, rec.features)

            # Append features to list
            rec.features.extend(coverage_features)
            print(' done')
        # end for

        sys.stdout.write(f'Writing annotated sequence to `{args.outfpath}`...')
        sys.stdout.flush()

        # Write result GanBank record
        out.write_genbank_output(
            rec,
            args.topology,
            args.organism,
            cov_array,
            args.outfpath
        )
        print('done')

        print('=' * 10)
    # end for

    if not args.keep_tmp_cov_file:
        _try_rm_temp_file(coverage_fpath)
    # end if
    print(f'Completed{with_warnings}!')
# end def


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
# end def


def _try_rm_temp_file(file_path):
    try:
        os.unlink(file_path)
    except OSError as err:
        print_err('Warning: cannot remove temporary file `{}`'.format(file_path))
        print_err(str(err))
    # end try
# end def


def _make_coverage_thresholds(args, cov_array):
    lower_coverage_thresholds = [
        LowerCoverageThreshold(thr_int)
            for thr_int in args.lower_coverage_thresholds
    ]

    median_coverage = cov_array.median_coverage

    upper_coverage_thresholds = [
        UpperCoverageThreshold(coef_float, median_coverage)
            for coef_float in args.upper_coverage_coefficients
    ]

    return lower_coverage_thresholds + upper_coverage_thresholds
# end def
