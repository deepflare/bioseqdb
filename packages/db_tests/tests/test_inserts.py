from sqlalchemy.sql import text
from sqlalchemy.exc import DataError
from db_tests.engine import BioseqdbTestDatabase
from db_tests.utils import Alphabet
import pytest

@pytest.mark.parametrize("type_name, input_alphabet, is_correct", [
    ("nucl_seq", Alphabet.NUCL_SEQ, True),
    ("nucl_seq", Alphabet.AMB_AA_SEQ, False),
    ("nucl_seq", Alphabet.AA_SEQ, False),
])
def test_basic_types(type_name: str, input_alphabet: Alphabet, is_correct):
    RAND_SEQUENCES = input_alphabet.random_list(
        min_length=0,
        max_length=200,
        min_count=10,
        max_count=20,
    )

    with BioseqdbTestDatabase() as db:
        db(
            f"""
            CREATE TABLE dataset (
                id serial primary key not null,
                seq {type_name} not null
            );
            """,
            None,
        )
        db(
            f"INSERT INTO dataset (id, seq) VALUES (:id, :seq);",
            None if is_correct else DataError,
            args=[dict(id=index, seq=seq) for index, seq in enumerate(RAND_SEQUENCES)],
        )
        if is_correct:
            db(
                "SELECT * FROM dataset;",
                list(enumerate(RAND_SEQUENCES))
            )


    