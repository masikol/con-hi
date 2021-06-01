# -*- encoding: utf-8 -*-

import os

import pytest

import src.parse_fasta_reference as pfr


@pytest.fixture
def plain_fasta() -> str:
    return os.path.join(os.getcwd(), 'tests', 'data', 'test_reference.fasta')
# end def plain_fasta

@pytest.fixture
def gzipped_fasta() -> str:
    return os.path.join(os.getcwd(), 'tests', 'data', 'test_reference.fasta.gz')
# end def gzipped_fasta


class TestIsPlainText:
    # Class for testing function `src.parse_fasta_reference._is_plain_text`

    def test_true_plain_fasta(self, plain_fasta):
        # Function testes how `_is_plain_text` recognizes plain files
        assert pfr._is_plain_text(plain_fasta) == True
    # end def test_true_plain_fasta

    def test_gzipped_fasta(self, gzipped_fasta):
        # Function testes how `_is_plain_text` recognizes gzipped files
        with pytest.raises(pfr._InvalidFileError):
            pfr._is_plain_text(gzipped_fasta)
        # end with
    # end def test_gzipped_fasta
# end class TestIsPlainText

class TestIsGzipped:
    # Class for testing function `src.parse_fasta_reference._is_gzipped`

    def test_gzipped_fasta(self, gzipped_fasta):
        # Function testes how `_is_gzipped` recognizes gzipped files
        assert pfr._is_gzipped(gzipped_fasta) == True
    # end def test_gzipped_fasta

    def test_true_plain_fasta(self, plain_fasta):
        # Function testes how `_is_gzipped` recognizes plain files
        assert pfr._is_gzipped(plain_fasta) == False
    # end def test_true_plain_fasta
# end class TestIsPlainText
