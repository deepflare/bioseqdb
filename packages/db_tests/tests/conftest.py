from typing import List
from db_tests.storage.engine import BioseqdbTestDatabase
from db_tests.utils import SequenceGenerator
import pytest

@pytest.fixture()
def db():
    with BioseqdbTestDatabase() as db:
        yield db
    # Teardown

@pytest.fixture()
def seq_generator():
    yield SequenceGenerator()