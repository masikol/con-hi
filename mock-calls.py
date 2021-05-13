#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

ref_fasta_fpath = '/mnt/1.5_drive_0/kromsatel-dev/reference/Wuhan-Hu-1-compele-genome.fasta'
bam_fpath = '/mnt/1.5_drive_0/kromsatel-dev/sam/Wuhan-Hu-1_cleaned.sorted.bam'
outdir = '/mnt/1.5_drive_0/kromsatel-dev/highlighter-outdir-test'

import src.obtain_coverage as oc
import src.parse_fasta_reference as pfr

# seq_id = 'NC_045512.2'
# cov_fpath = oc.count_cov_for_all_refs(ref_fasta_fpath, bam_fpath, outdir)
# cov_array = oc.get_coverage_for_reference(seq_id, cov_fpath)


fasta_records = pfr.parse_fasta_reference(ref_fasta_fpath)

print(type(fasta_records[0]))

for rec in fasta_records:
    print(rec.seq)

