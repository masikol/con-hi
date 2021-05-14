# -*- encoding: utf-8 -*-
# Version 1.0.a

import gzip
import functools
from typing import Sequence, Callable, TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_fasta_reference(ref_fasta_fpath: str) -> Sequence[SeqRecord]:

    open_func: Callable[[str], TextIO]

    if _is_plain_text(ref_fasta_fpath):
        open_func = functools.partial(open, mode='rt')
    elif _is_gzipped(ref_fasta_fpath):
        open_func = functools.partial(gzip.open, mode='rt')
    else:
        raise _InvalidFileError(f'Format of file `{ref_fasta_fpath}` is not supported.\
 Only plain fasta or gzipped fasta are supported.\
 Allowed extentions: `.fasta`, `.fa`, `.fasta.gz`, `.fa.gz`.')
    # end if

    fasta_records: Sequence[SeqRecord]

    with open_func(ref_fasta_fpath) as infpath:
        fasta_records = tuple(SeqIO.parse(infpath, 'fasta'))
    # end with

    return fasta_records
# end def parse_fasta_reference


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
# end def _is_gzipped


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
# end def _is_plain_text


class _InvalidFileError(OSError):
    # Custom exception for indicating that file compression type is not supported.
    pass
# end class _InvalidFileError
