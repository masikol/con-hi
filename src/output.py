# Version 2.0.a

import os
import sys
import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.platform import platf_depend_exit


def create_or_emply_file(file_path):
    try:
        with open(file_path, 'wt') as _:
            pass
        # end with
    except OSError as err:
        print(f'\nError: cannot create file {file_path}')
        print(str(err))
        platf_depend_exit(1)
    # end try
# end def create_or_emply_file


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
    with open(outfpath, 'a') as outfile:
        SeqIO.write(seq_record, outfile, 'genbank')
    # end with
# end def write_genbank_output


def _get_date() -> str:
    # Function returns formatted current date.
    now: datetime.datetime = datetime.datetime.now()
    return now.strftime('%d-%b-%Y').upper()
# end def _get_date


def conf_path_to_depth_file(outfpath: str) -> str:
    # Function configures path to coverage file.
    # :param outdir: path to output directory;
    return os.path.join(
        os.path.dirname(outfpath),
        'coverages.tsv'
    )
# end def conf_path_to_depth_file
