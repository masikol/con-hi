# -*- encoding: utf-8 -*-
# Version 1.0.a

from typing import Sequence

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_fasta_reference(ref_fasta_fpath: str) -> Sequence[SeqRecord]:

    fasta_records: Sequence[SeqRecord]

    with open(ref_fasta_fpath, 'r') as infpath:
        fasta_records = tuple(SeqIO.parse(infpath, 'fasta'))
    # end with

    return fasta_records
# end def parse_fasta_reference
