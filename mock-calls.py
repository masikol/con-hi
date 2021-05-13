#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

ref_fasta_fpath = '/mnt/1.5_drive_0/kromsatel-dev/reference/Wuhan-Hu-1-compele-genome.fasta'
bam_fpath = '/mnt/1.5_drive_0/kromsatel-dev/sam/Wuhan-Hu-1_cleaned.sorted.bam'
outdir = '/mnt/1.5_drive_0/kromsatel-dev/highlighter-outdir-test'

import src.obtain_coverage as oc

seq_id = 'NC_045512.2'

cov_fpath = oc.count_cov_for_all_refs(ref_fasta_fpath, bam_fpath, outdir)

cov_array = oc.get_coverage_for_reference(seq_id, cov_fpath)

print(cov_array[-1])
print(cov_array[29903])
print(len(cov_array))
