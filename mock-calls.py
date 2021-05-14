#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os

ref_fasta_fpath = '/mnt/1.5_drive_0/kromsatel-dev/reference/Wuhan-Hu-1-compele-genome.fasta'
bam_fpath = '/mnt/1.5_drive_0/kromsatel-dev/sam/Wuhan-Hu-1_cleaned.sorted.bam'
outdir = '/mnt/1.5_drive_0/kromsatel-dev/highlighter-outdir-test'

import src.obtain_coverage as oc
import src.parse_fasta_reference as pfr
import src.highlight_features as hlft
import src.output as out


fasta_records = pfr.parse_fasta_reference(ref_fasta_fpath)

cov_fpath = oc.count_cov_for_all_refs(ref_fasta_fpath, bam_fpath, outdir)

for rec in fasta_records:

    cov_array = oc.get_coverage_for_reference(rec.id, cov_fpath)
    rec.features = hlft.highlight_coverage_features(cov_array, lambda x: x < 10, 'sabaka')

    for f in rec.features:
        print(f)

    outfpath = os.path.join(outdir, 'output.gbk')
    out.write_genbank_output(rec, outfpath)
    print('Done')
# end for
