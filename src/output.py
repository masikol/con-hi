
import os
import logging
import datetime
import warnings

from Bio import SeqIO
from Bio import BiopythonWarning
from Bio.SeqRecord import SeqRecord

from src.platform import platf_depend_exit
from src.coverage_array import CoverageArray


def create_or_emply_file(file_path):
    try:
        with open(file_path, 'wt') as _:
            pass
        # end with
    except OSError as err:
        logging.error(f'Error: cannot create file {file_path}')
        logging.error(str(err))
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

    _write_gb_warning_safe(seq_record, outfpath, mode='a')
# end def


def _write_gb_warning_safe(seq_record: SeqRecord, outfpath: str, mode: str = 'a'):
    # Generic Python warnings get printed with one extra line, e.g. "warnings.warn("
    # So, we will catch the warning, print its message by ourselves,
    #   and write the record ignoring warnings ignore warnings

    records_written: int = 0

    # Write output file
    # Catch and ignore BiopythonWarning
    warnings.filterwarnings('error')
    try:
        records_written = _write_gb(seq_record, outfpath, mode=mode)
    except BiopythonWarning as w:
        logging.warning(str(w))
        # We will write the record (again) only if it hasn't been already written
        if records_written == 0:
            warnings.filterwarnings('ignore')
            _write_gb(seq_record, outfpath, mode=mode)
        # end if
    finally:
        warnings.resetwarnings()
    # end try
# end def


def _write_gb(seq_record: SeqRecord, outfpath: str, mode: str = 'a') -> int:
    records_written: int = 0
    with open(outfpath, mode) as outfile:
        records_written = SeqIO.write(seq_record, outfile, 'genbank')
    # end with
    return records_written
# end def


def _create_annotations(organism: str,
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


def conf_path_to_depth_file(outfpath: str, target_seq_id: str) -> str:
    # Function configures path to coverage file.
    # :param outfpath: path to output file;
    # :param target_seq_id: target sequence id;
    return os.path.join(
        os.path.dirname(outfpath),
        f'coverage_{target_seq_id}.tsv'
    )
# end def
