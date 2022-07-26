
import os
import sys
import datetime

from Bio import SeqIO
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
# end def create_or_emply_file


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

    # Write output file
    with open(outfpath, 'a') as outfile:
        SeqIO.write(seq_record, outfile, 'genbank')
    # end with
# end def write_genbank_output


def _create_annotations(seq_record: SeqRecord,
                        organism: str,
                        topology: str,
                        cov_array: CoverageArray) -> dict:

    zero_coverage_bases = cov_array.count(0)

    annotations = {
        'molecule_type': 'DNA',
        'organism': organism,
        'date': _get_date(),
        'topology': topology,
        'structured_comment': {
            'Coverage-Data': {
                'Minimum Coverage': cov_array.calc_min_coverage(),
                'Average Coverage': cov_array.calc_avg_coverage(),
                'Median Coverage' : cov_array.calc_median_coverage(),
                'Maximum Coverage': cov_array.calc_max_coverage(),
                'Zero-coverage bases': '{} bp'.format(zero_coverage_bases),
            }
        }
    }

    return annotations
# end def


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
