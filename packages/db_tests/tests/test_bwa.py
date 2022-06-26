from sqlalchemy.sql import text
from sqlalchemy.exc import DataError
from db_tests.storage.engine import BioseqdbTestEngine
from db_tests.utils import Alphabet, SequenceGenerator
from db_tests.models import BWAAlignment
from db_tests.storage.records_list import RecordsList
import pytest

@pytest.mark.parametrize(
    "type_name,input_alphabet",
    [
        ("nucl_seq", Alphabet.NUCL_SEQ),
    ],
)
class TestBWA:
    """
        Testing database BWA algorithm
    """

    def test_bwa_alignments(
        self,
        seq_generator: SequenceGenerator,
        db: BioseqdbTestEngine,
        type_name: str,
        input_alphabet: Alphabet,
    ):
        input_sequences = seq_generator(
            "Input sequences",
            input_alphabet,
            min_length=250,
        )

        assert db(
            f"""
            CREATE TABLE dataset (
                id SERIAL PRIMARY KEY NOT NULL,
                seq {type_name} NOT NULL
            );
            """,
        ) == None, f"Create new dataset table with sequences of type {type_name}"

        assert db(
            f"INSERT INTO dataset (id, seq) VALUES (:id, :seq);",
            args=[dict(id=index, seq=seq) for index, seq in enumerate(input_sequences)],
        ) == None, f"Insert randomized {type_name} sequences into the dataset table."
        
        assert db(
            "SELECT COUNT(*) FROM dataset;",
        ) == RecordsList([(len(input_sequences),)]), "Count all records in the dataset table"

        first_seq_len = len(input_sequences[0])
        assert db(
            "SELECT * FROM nuclseq_search_bwa(nuclseq_in(:searched_seq), 'SELECT id, seq FROM dataset');",
            args=[dict(searched_seq=input_sequences[0])],
            model=BWAAlignment,
        ) == RecordsList([
            BWAAlignment(
                ref_id=0,
                ref_subseq=input_sequences[0],
                ref_match_start=0,
                ref_match_end=first_seq_len,
                ref_match_len=first_seq_len,
                query_id=None,
                query_subseq=input_sequences[0],
                query_match_start=0,
                query_match_end=first_seq_len,
                query_match_len=first_seq_len,
                is_primary=True,
                is_secondary=False,
                is_reverse=False,
                cigar=f'{first_seq_len}M',
                score=first_seq_len,
            ),
        ]), "Perform BWA alignment sequence-to-table and search first sequence among the table"

        assert db(
            "SELECT * FROM nuclseq_multi_search_bwa('SELECT id, seq FROM dataset', 'SELECT id, seq FROM dataset');",
            model=BWAAlignment,
            ignore_order=True,
        ) == RecordsList([
            BWAAlignment(
                ref_id=0,
                ref_subseq=seq,
                ref_match_start=0,
                ref_match_end=len(seq),
                ref_match_len=len(seq),
                query_id=None,
                query_subseq=seq,
                query_match_start=0,
                query_match_end=len(seq),
                query_match_len=len(seq),
                is_primary=True,
                is_secondary=False,
                is_reverse=False,
                cigar=f'{len(seq)}M',
                score=len(seq),
            ) for seq in input_sequences
        ]), "Perform BWA alignment table-to-table for long sequences"

