from sqlalchemy import TypeDecorator, String
import sqlalchemy.types as types
from sqlalchemy.util import generic_repr

class NuclSeq(types.String):
    cache_ok = True

    def __init__(self, **kw):
        super().__init__(**kw)

    def get_col_spec(self, **kw):
        return "NUCL_SEQ"

    def bind_processor(self, dialect):
        def process(value):
            return value
        return process

    def result_processor(self, dialect, coltype):
        def process(value):
            return value
        return process

    def __repr__(self):
        s = generic_repr(
            self, to_inspect=[NuclSeq, types.String]
        )
        return f"NUCL({s})"

