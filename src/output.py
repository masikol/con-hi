
import os
import sys
import datetime
import warnings

from Bio import SeqIO
from Bio import BiopythonWarning
from Bio.SeqRecord import SeqRecord

from src.printing import print_err
from src.platform import platf_depend_exit
from src.coverage_array import CoverageArray


def create_or_emply_file(file_path):
    try:
        with open(file_path, 'wt') as _:
            pass
        # end with
    except OSError as err:
        print_err(f'\nError: cannot create file {file_path}')
        print_err(str(err))
        platf_depend_exit(1)
    # end try
# end def


def write_genbank_output(seq_record: SeqRecord,
                         topology: str,
                         organism: str,
                         cov_array: CoverageArray,
                         outfpath: str) -> None:
    # Function writes annotated sequence to output GenBank file.
    # :param seq_record: sequence record to output;
    # :param topology: string for `topology` annotation;
    # :param organism: string for `organism` annotation;
    # :param outfpath: path to output file;

    annotations = _create_annotations(
        seq_record,
        organism,
        topology,
        cov_array
    )
    seq_record.annotations = annotations

    # Sort features by their location if ascending order
    seq_record.features = sorted(
        seq_record.features,
        key = lambda feature: feature.location.start
    )

    # As written in Biopython
    #   (file `InsdcIO.py`, class `GenBankWriter`, method `_write_the_first_line`),
    # *quote* Locus name and record length to[o] long to squeeze in.
    # *quote* Per updated GenBank standard (Dec 15, 2018) 229.0
    # *quote* the Locus identifier can be any length, and a space
    # *quote* is added after the identifier to keep the identifier
    # *quote* and length fields separated
    locus_str_len: int = len(seq_record.name)
    seq_len_str_len: int = len(str(len(seq_record)))
    if locus_str_len + 1 + seq_len_str_len > 28:
        print_err('\n! Warning: Locus name and record length to long to squeeze in,\n')
        print_err('   and thus GenBank representation will not be "pretty".\n')
        print_err('The GenBank file is still valid, though.\n')
    # end if

    # Catch and ignore BiopythonWarning, for we have a better one printed already
    warnings.filterwarnings('error')
    try:
        # Write output file
        with open(outfpath, 'a') as outfile:
            SeqIO.write(seq_record, outfile, 'genbank')
        # end with
    except BiopythonWarning:
        pass
    finally:
        warnings.resetwarnings()
    # end try
# end def


def _create_annotations(seq_record: SeqRecord,
                        organism: str,
                        topology: str,
                        cov_array: CoverageArray) -> dict:
    annotations = {
        'molecule_type': 'DNA',
        'organism': organism,
        'date': _get_date(),
        'topology': topology,
        'structured_comment': {
            'Coverage-Data': {
                'Minimum Coverage': cov_array.min_coverage,
                'Average Coverage': cov_array.avg_coverage,
                'Median Coverage' : cov_array.median_coverage,
                'Maximum Coverage': cov_array.max_coverage,
                'Zero-coverage bases': '{} bp'.format(cov_array.zero_coverage_bases),
            }
        }
    }

    return annotations
# end def


def _get_date() -> str:
    # Function returns formatted current date.
    now: datetime.datetime = datetime.datetime.now()
    return now.strftime('%d-%b-%Y').upper()
# end def


def conf_path_to_depth_file(outfpath: str) -> str:
    # Function configures path to coverage file.
    # :param outdir: path to output directory;
    return os.path.join(
        os.path.dirname(outfpath),
        'coverages.tsv'
    )
# end def
