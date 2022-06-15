from sqlalchemy.sql import text
from sqlalchemy.exc import DataError
from db_tests.engine import BioseqdbTestDatabase
from db_tests.utils import Alphabet
import pytest
import string

def test_basic_types():
    with BioseqdbTestDatabase() as db:
        db(
            f"""
            SELECT * FROM nuclseq_len(nuclseq_in(:seq))
            """,
            [(5,)],
            args=dict(seq="AACTG")
        )
        db(
            f"""
            SELECT * FROM nuclseq_reverse(nuclseq_in(:seq))
            """,
            [("GTCAA",)],
            args=dict(seq="AACTG")
        )
        for letter in string.printable:
            seq = "AACTG"
            db(
                f"""
                SELECT * FROM nuclseq_content(nuclseq_in(:seq), :letter)
                """,
                [(seq.count(letter)/len(seq),)] if letter in Alphabet.NUCL_SEQ.value else DataError,
                args=dict(seq=seq, letter=letter)
            )