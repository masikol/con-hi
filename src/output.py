# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def write_genbank_output(seq_record: SeqRecord, topology: str, organism: str, outfpath: str) -> None:
    # Function writes annotated sequence to output GenBank file.
    # :param seq_record: sequence record to output;
    # :param topology: string for `topology` annotation;
    # :param organism: string for `organism` annotation;
    # :param outfpath: path to output file;

    # Set annotations
    seq_record.annotations = {
        'molecule_type': 'DNA',
        'organism': organism,
        'date': _get_date(),
        'topology': topology
    }

    # Sort features by their location if ascending order
    seq_record.features = sorted(
        seq_record.features,
        key = lambda feature: feature.location.start
    )

    # Write output file
    with open(outfpath, 'w') as outfile:
        SeqIO.write(seq_record, outfile, 'genbank')
    # end with
# end def write_genbank_output


def configure_outfile_path(outdir_path: str, prefix: str, seq_id: str) -> str:
    # Function configures path to output GenBank file
    # :param outdir_path: path to output directory;
    # :param prefix: prefix for output files;
    # :param seq_id: id of sequence;

    outfpath: str

    # The difference between these two branches is in the underscrore `_`
    if prefix == '':
        outfpath = os.path.join(outdir_path, f'{seq_id}.gbk')
    else:
        outfpath = os.path.join(outdir_path, f'{prefix}_{seq_id}.gbk')
    # end if

    return outfpath
# end def configure_outfile_path


def _get_date() -> str:
    # Function returns formatted current date.
    now: datetime.datetime = datetime.datetime.now()
    return now.strftime('%d-%b-%Y').upper()
# end def _get_date
