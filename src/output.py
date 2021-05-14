# -*- encoding: utf-8 -*-
# Version 1.0.a

import os
import datetime

from Bio import SeqIO

def write_genbank_output(seq_record, topology, organism, outfpath):

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

    with open(outfpath, 'w') as outfile:
        SeqIO.write(seq_record, outfile, 'genbank')
    # end with
# end def write_genbank_output


def configure_outfile_path(outdir_path: str, prefix: str, seq_id: str) -> str:

    outfpath: str

    if prefix == '':
        outfpath = os.path.join(outdir_path, f'{seq_id}.gbk')
    else:
        outfpath = os.path.join(outdir_path, f'{prefix}_{seq_id}.gbk')
    # end if

    return outfpath
# end def configure_outfile_path


def _get_date() -> str:
    now: datetime.datetime = datetime.datetime.now()
    return now.strftime('%d-%b-%Y').upper()
# end def _get_date
