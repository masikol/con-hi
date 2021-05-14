# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
import sys

from Bio.SeqRecord import SeqRecord

import src.output as out
import src.obtain_coverage as oc
import src.highlight_features as hlft
import src.parse_fasta_reference as pfr
from src.coverage_array import CoverageArray
from src.coverage_threshold import CoverageTheshold
from src.arguments import parse_arguments, HighlighterParams

def main(version, last_update_date) -> None:

    # ref_fasta_fpath = '/mnt/1.5_drive_0/kromsatel-dev/reference/Wuhan-Hu-1-compele-genome.fasta'
    # bam_fpath = '/mnt/1.5_drive_0/kromsatel-dev/sam/Wuhan-Hu-1_cleaned.sorted.bam'
    # outdir = '/mnt/1.5_drive_0/kromsatel-dev/highlighter-outdir-test'

    params: HighlighterParams = parse_arguments(version, last_update_date)

    base_feature_note: str = f'generated with consensus-highlighter.py v{version}'

    # Read fasta records from input file
    print('Importing fasta from `{}`...'.format(params.target_fasta_fpath), end=' ')
    sys.stdout.flush()
    fasta_records: Sequence[SeqRecord] = pfr.parse_fasta_reference(params.target_fasta_fpath)
    print('done')

    # Count coverages with samtools depth
    print('Counting coverages with `samtools depth`...', end=' ')
    sys.stdout.flush()
    cov_fpath: str = oc.count_cov_for_all_refs(
        params.target_fasta_fpath,
        params.bam_fpath,
        params.outdir_path
    )
    print('done\n')

    rec: SeqRecord
    for rec in fasta_records:

        print(f'Processing sequence `{rec.description}`')

        cov_array: CoverageArray = oc.get_coverage_for_reference(rec.id, cov_fpath)

        cov_threshold: CoverageTheshold
        coverage_features: MutableSequence[SeqFeature]

        for cov_threshold in params.coverage_thresholds:

            print(f'Screening the sequence for segments with {cov_threshold.get_label()}...', end=' ')
            sys.stdout.flush()

            coverage_features = hlft.highlight_coverage_features(
                cov_array,
                cov_threshold,
                base_feature_note
            )

            rec.features.extend(coverage_features)
            print('done')
        # end for

        outfpath: str = out.configure_outfile_path(
            params.outdir_path,
            params.outfile_prefix,
            rec.id
        )

        print(f'Writing annotated sequence to `{outfpath}`...', end=' ')
        sys.stdout.flush()

        out.write_genbank_output(rec, params.topology, params.organism, outfpath)
        print('done')

        print('=' * 10)
    # end for

    print('Completed!')
# end def main
