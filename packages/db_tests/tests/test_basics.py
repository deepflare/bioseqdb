from sqlalchemy.sql import text
from sqlalchemy.exc import DataError
from db_tests.storage.engine import BioseqdbTestEngine
from db_tests.utils import Alphabet, SequenceGenerator
from db_tests.models import BWAAlignment
from db_tests.storage.records_list import RecordsList
import string
import pytest

class TestBasics:
    """
        Testing sequence basic operations
    """

    def test_bwa_alignments(
        self,
        seq_generator: SequenceGenerator,
        db: BioseqdbTestEngine,
    ):

        assert db(
            f"""
            SELECT * FROM nuclseq_len(nuclseq_in(:seq))
            """,
            args=dict(seq="AACTG"),
        ) == RecordsList([(5,)])

        assert db(
            f"""
            SELECT * FROM nuclseq_reverse(nuclseq_in(:seq))
            """,
            args=dict(seq="AACTG"),
        ) == RecordsList([("GTCAA",)])

        # for letter in string.printable:
        #     seq = "AACTG"
        #     assert db(
        #         f"""
        #         SELECT * FROM nuclseq_content(nuclseq_in(:seq), :letter)
        #         """,
        #         args=dict(seq=seq, letter=letter),
        #     ) == (RecordsList([(seq.count(letter)/len(seq),)]) if letter in Alphabet.NUCL_SEQ.value else DataError)

