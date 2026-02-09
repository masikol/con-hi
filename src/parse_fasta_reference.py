
import re
import gzip
import functools
from typing import Sequence, Callable, TextIO, List, Set

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_fasta_reference(ref_fasta_fpath: str,
                          ref_seq_ids: Set[str] = set()) -> Sequence[SeqRecord]:
    # Function parses fasta file and returns collection of fasta records stored in this file.
    # :param ref_fasta_fpath: path to target fasta file;
    # :param ref_seq_ids: ids of target sequences to parse;

    # Choose open function for the file
    open_func: Callable[[str], TextIO]

    if _is_gzipped(ref_fasta_fpath):
        open_func = functools.partial(gzip.open, mode='rt')
    elif _is_plain_text(ref_fasta_fpath):
        open_func = functools.partial(open, mode='rt')
    else:
        raise _InvalidFileError(f'Format of file `{ref_fasta_fpath}` is not supported.\
 Only plain fasta or gzipped fasta are supported.\
 Allowed extentions: `.fasta(.gz)`, `.fa(.gz)`, `.fna(.gz)`.')
    # end if

    # Parse fasta data
    fasta_records: Sequence[SeqRecord]
    with open_func(ref_fasta_fpath) as infile:
        fasta_records = tuple(
            SeqIO.parse(infile, 'fasta')
        )
    # end with

    # Make all sequences uppercase
    record: SeqRecord
    for record in fasta_records:
        record.seq = record.seq.upper()
    # end for

    # If no ref_seq_ids are provided, use all present in the target fasta
    # Otherwise, subset only sequences having provided target ids
    if len(ref_seq_ids) != 0:
        _validate_user_target_seq_ids(ref_seq_ids, fasta_records)
        fasta_records = tuple(
            filter(
                lambda sr: sr.id in ref_seq_ids,
                fasta_records
            )
        )
    # end if

    # Validate parsed fasta data
    _validate_fasta_records(fasta_records)

    return fasta_records
# end def


def _is_gzipped(fpath: str) -> bool:
    # Function, which returs True if passed file is a file compressed with gzip.
    # :param fpath: path to file ot check;

    if not fpath.endswith('.gz'):
        return False
    # end if

    # Test this file
    try:
        with gzip.open(fpath, 'rt', encoding='utf-8') as gzfile:
            gzfile.readline() # this will fail if file is not valid
        # end with
    except (gzip.BadGzipFile, OSError, UnicodeDecodeError) as err:
        raise _InvalidFileError(f'Error: file `{fpath}` is not a valid gzip file: {err}')
    else:
        return True
    # end try
# end def


def _is_plain_text(fpath: str) -> bool:
    # Function, which returs True if passed file is a plain text file.
    # :param fpath: path to file ot check;
    # :type fpath: str;

    # Test this file
    try:
        with open(fpath, 'r', encoding='utf-8') as file:
            file.readline() # this will fail if file is not valid
        # end with
    except (OSError, UnicodeDecodeError) as err:
        raise _InvalidFileError(f'Error: cannot read file `{fpath}`: {err}.')
    else:
        return True
    # end try
# end def


def _validate_fasta_records(fasta_records: Sequence[SeqRecord]) -> None:
    # Function validates input fasta file.
    # :param fasta_records: list of fasta records to validate;

    # Check if passed list is empty
    if len(fasta_records) == 0:
        raise ValueError("""
  Error: your fasta file contain no sequences.""")
    # end if

    # Check if fasta records contain non-IUPAC characters
    iupac_codes: str = 'AGCTRYSWKMBDHVN'
    invalid_records: List = list()

    for record in fasta_records:
        if not re.search(r'[^{}]'.format(iupac_codes), str(record.seq)) is None:
            invalid_records.append(record.id)
        # end if
    # end for

    if len(invalid_records) != 0:
        raise ValueError(f"""
  Error: some fasta records contain invalid DNA characters.
  Here are are ids of these sequences:
  `{'`, `'.join(invalid_records)}`.
  Invalid characters are not {iupac_codes}.""")
    # end if
# end def

def _validate_user_target_seq_ids(ref_seq_ids: Set[str],
                                  fasta_records: Sequence[SeqRecord]):
    # The function checks if all ref_seq_ids are in fasta_records
    # :param ref_seq_ids: target sequence ids;
    # :param fasta_records: parsed fasta records;
    seq_ids_in_fasta = set(
        map(lambda sr: sr.id, fasta_records)
    )
    missing_target_seq_ids = set(ref_seq_ids) - seq_ids_in_fasta
    if len(missing_target_seq_ids) != 0:
        missing_target_seq_ids = ','.join(missing_target_seq_ids)
        raise ValueError(f"""
  Error: these target seq IDs are not in the target fasta file:
  {missing_target_seq_ids}""")
    # end if
# end def


class _InvalidFileError(OSError):
    # Custom exception for indicating that file compression type is not supported.
    pass
# end class
