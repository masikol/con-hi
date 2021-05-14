# -*- encoding: utf-8 -*-
# Version 1.0.a

from Bio import SeqIO

def write_genbank_output(seq_record, outfpath):

    seq_record.annotations = {'molecule_type': 'DNA'}

    with open(outfpath, 'w') as outfile:
        SeqIO.write(seq_record, outfile, 'genbank')
    # end with
# end def write_genbank_output
