
import os
from typing import Sequence, Set

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import src.parse_fasta_reference as pfr

from tests.fixtures import first_test_seq_len, \
                           second_test_seq_len, \
                           third_test_seq_len


# == Fixtures for testing fucntions `src.parse_fasta_reference._is_plain_text`
#    and `src.parse_fasta_reference._is_gzipped`
#    and `src.parse_fasta_reference.parse_fasta_reference` ==

@pytest.fixture
def plain_fasta() -> str:
    return os.path.join(os.getcwd(), 'tests', 'data', 'test_reference.fasta')
# end def

@pytest.fixture
def gzipped_fasta() -> str:
    return os.path.join(os.getcwd(), 'tests', 'data', 'test_reference.fasta.gz')
# end def

@pytest.fixture
def target_seq_ids() -> Set[str]:
    return {
        'test_seq_1',
        'test_seq_3',
    }
# end def


# == Fixtures for testing fucntion `src.parse_fasta_reference._validate_fasta_records` ==

@pytest.fixture
def empty_seq_list() -> Sequence[SeqRecord]:
    return []
# end def

@pytest.fixture
def unambig_dna_seq_records() -> Sequence[SeqRecord]:
    return [
        SeqRecord(Seq('ATTAAAGGTTTATACCTTCCCAGGTA')),
        SeqRecord(Seq('CTTGTAGATCTGTTCTCTAAA')),
        SeqRecord(Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACTTTGGTC')),
    ]
# end def

@pytest.fixture
def ambig_dna_seq_records() -> Sequence[SeqRecord]:
    return [
        SeqRecord(Seq('ATTAAAGRYSWTACCTTCCCAGGTA')),
        #                 ~~~~
        SeqRecord(Seq('CTTGKMBDATCTGTTCTCTAAA')),
        #              ~~~~
        SeqRecord(Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTBDHVNGGAGGATCTGTTTACTTTGGTC')),
        #                                           ~~~~~
    ]
# end def

@pytest.fixture
def invalid_dna_seq_records() -> Sequence[SeqRecord]:
    return [
        SeqRecord(Seq('ATTAAAGRYSWTACC45674CAGGTA')),
        #                         ~~~~~
        SeqRecord(Seq('CTTGKMBDAT546456TCTCTAAA')),
        #                    ~~~~~~
        SeqRecord(Seq('CATTGTTGAGATCACAT888GAGTTBDHVNGXXGGATCTGTTTACTTTGGTC')),
        #                           ~~~           ~~
    ]
# end def


class TestIsPlainText:
    # Class for testing function `src.parse_fasta_reference._is_plain_text`

    def test_true_plain_fasta(self, plain_fasta) -> None:
        # Function testes how `_is_plain_text` recognizes plain files
        assert pfr._is_plain_text(plain_fasta) == True
    # end def

    def test_gzipped_fasta(self, gzipped_fasta) -> None:
        # Function testes how `_is_plain_text` recognizes gzipped files
        with pytest.raises(pfr._InvalidFileError):
            pfr._is_plain_text(gzipped_fasta)
        # end with
    # end def
# end class


class TestIsGzipped:
    # Class for testing function `src.parse_fasta_reference._is_gzipped`

    def test_gzipped_fasta(self, gzipped_fasta) -> None:
        # Function testes how `_is_gzipped` recognizes gzipped files
        assert pfr._is_gzipped(gzipped_fasta) == True
    # end def

    def test_true_plain_fasta(self, plain_fasta) -> None:
        # Function testes how `_is_gzipped` recognizes plain files
        assert pfr._is_gzipped(plain_fasta) == False
    # end def
# end class


class TestValidateFastaRecords:
    # Class for testing function `src.parse_fasta_reference._validate_fasta_records`

    def test_empty_seq_set(self, empty_seq_list) -> None:
        # Function tests how `_validate_fasta_records` handles empty sequence set

        # Should raise a ValueError
        with pytest.raises(ValueError):
            pfr._validate_fasta_records(empty_seq_list)
        # end with
    # end def

    def test_unambig_seq_records(self, unambig_dna_seq_records) -> None:
        # Function tests how `_validate_fasta_records` handles set
        #    of unambiguous DNA records

        # Should not raise any exception
        assert pfr._validate_fasta_records(unambig_dna_seq_records) == None
    # end def

    def test_ambig_seq_records(self, ambig_dna_seq_records) -> None:
        # Function tests how `_validate_fasta_records` handles set
        #    of ambiguous DNA records

        # Should not raise any exception
        assert pfr._validate_fasta_records(ambig_dna_seq_records) == None
    # end def

    def test_invalid_dna_records(self, invalid_dna_seq_records) -> None:
        # Function tests how `_validate_fasta_records` handles set
        #    of invalid DNA records

        # Should raise a ValueError
        with pytest.raises(ValueError):
            pfr._validate_fasta_records(invalid_dna_seq_records)
        # end with
    # end def
# end class


class TestParseFastaReference:
    # Class for testing function `src.parse_fasta_reference.parse_fasta_reference`

    def test_parse_plain_text_fasta(
        self,
        plain_fasta: str,
        first_test_seq_len: int,
        second_test_seq_len: int,
        third_test_seq_len: int) -> None:

        # Function tests how `parse_fasta_reference` parses plain text fasta file
        fasta_records: Sequence[SeqRecord] = pfr.parse_fasta_reference(plain_fasta)

        expected_num_seqs: int = 3
        assert len(fasta_records) == expected_num_seqs

        expected_seq_len: int

        expected_seq_len = first_test_seq_len
        assert len(fasta_records[0]) == expected_seq_len

        expected_seq_len = second_test_seq_len
        assert len(fasta_records[1]) == expected_seq_len

        expected_seq_len = third_test_seq_len
        assert len(fasta_records[2]) == expected_seq_len
    # end def

    def test_parse_gzipped_fasta(
        self,
        gzipped_fasta: str,
        first_test_seq_len: int,
        second_test_seq_len: int,
        third_test_seq_len: int) -> None:

        # Function tests how `parse_fasta_reference` parses gzipped fasta file
        fasta_records: Sequence[SeqRecord] = pfr.parse_fasta_reference(gzipped_fasta)

        expected_num_seqs: int = 3
        assert len(fasta_records) == expected_num_seqs

        expected_seq_len: int

        expected_seq_len = first_test_seq_len
        assert len(fasta_records[0]) == expected_seq_len

        expected_seq_len = second_test_seq_len
        assert len(fasta_records[1]) == expected_seq_len

        expected_seq_len = third_test_seq_len
        assert len(fasta_records[2]) == expected_seq_len
    # end def

    def test_subset_target_seq_ids(self,
                                   plain_fasta: str,
                                   target_seq_ids: Set[str]) -> None:
        # Function tests how `parse_fasta_reference` subsets fasta records
        #   using its `ref_seq_ids` parameter
        fasta_records: Sequence[SeqRecord] = pfr.parse_fasta_reference(
            plain_fasta,
            target_seq_ids
        )

        expected: Set[str] = set(target_seq_ids)
        observed: Set[str] = set(
            map(lambda sr: sr.id, fasta_records)
        )
        assert observed == expected
    # end def
# end class
