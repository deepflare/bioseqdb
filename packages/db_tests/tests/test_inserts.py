from sqlalchemy.sql import text
from sqlalchemy.exc import DataError
from db_tests.storage.engine import BioseqdbTestEngine
from db_tests.utils import Alphabet, SequenceGenerator
from db_tests.models import BWAAlignment
from db_tests.storage.records_list import RecordsList
import pytest

@pytest.mark.parametrize("type_name,input_alphabet,is_correct", [
    ("nucl_seq", Alphabet.NUCL_SEQ, True),
    ("nucl_seq", Alphabet.AMB_AA_SEQ, False),
    ("nucl_seq", Alphabet.AA_SEQ, False),
])
class TestInserts:
    """
        Testing database sequence inserts
    """

    def test_inserts(
        self,
        seq_generator: SequenceGenerator,
        db: BioseqdbTestEngine,
        type_name: str,
        input_alphabet: Alphabet,
        is_correct: bool,
    ):
        input_sequences = seq_generator(
            "Input sequences",
            input_alphabet,
            min_length=0,
            max_length=200,
            min_count=10,
            max_count=20,
        )

        assert db(
            f"""
            CREATE TABLE dataset (
                id serial primary key not null,
                seq {type_name} not null
            );
            """,
        ) == None, f"Create new dataset table with sequences of type {type_name}"

        assert db(
            f"INSERT INTO dataset (id, seq) VALUES (:id, :seq);",
            args=[dict(id=index, seq=seq) for index, seq in enumerate(input_sequences)],
        ) == (None if is_correct else DataError), f"Insert sequences into the table of type {type_name} from alphabet {input_alphabet}"

        if is_correct:
            assert db(
                "SELECT * FROM dataset;",
            ) == RecordsList(enumerate(input_sequences)), "Validate the table has all created sequences"

