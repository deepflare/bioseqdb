from sqlalchemy import create_engine as create_engine_raw
from sqlalchemy.sql import text
import testing.postgresql
from testing.common.database import Database
from sqlalchemy.engine import Engine
from pgsanity.pgsanity import check_string as validate_pg_sql

from sqlalchemy.orm import Session
from sqlalchemy.exc import ResourceClosedError
from typing import Optional, List, Tuple, Any, Union
from sqlalchemy.orm import sessionmaker


AnyRowList = Optional[Tuple[Any, ...]]
SQLQuery = str

class BioseqdbTestEngine():
    engine: Engine
    session_maker: sessionmaker

    def __init__(self, engine: Engine):
        self.engine = engine
        self.session_maker = sessionmaker(bind=engine)

    def assert_sql(
        self,
        sql_content: SQLQuery,
        expected: AnyRowList,
        args = None,
    ) -> AnyRowList:
        #success, msg = validate_pg_sql(sql_content)
        #print(msg)
        #if not success:
        #    raise ValueError(f"Invalid SQL query: {sql_content}")
        if not args:
            args = dict()
        #with self.engine.connect() as con:
        with self.session_maker() as con:
        #with Session(self.engine) as session:
            try:
                statement = text(sql_content)
                #for name in args.keys():
                #    statement = statement.bindparams(bindparam(name, expanding=True))
                print(args)
                r = con.execute(statement, args)
                print(f"EXPECTED: {expected}")
                print(f"EXP?1 => {(expected is None)}")
                print(f"EXP?2 => {(not isinstance(expected, list))}")
                con.commit()
                if (expected is None) or (not isinstance(expected, list)):
                    return None
                actual = list(r)
            except ResourceClosedError:
                actual = None
            except Exception as e:
                if expected and not (isinstance(expected, list) or expected is None):
                    expected_exception_class = expected
                    actual_expection_class = e.__class__
                    print(e, expected_exception_class, actual_expection_class)
                    assert expected_exception_class == actual_expection_class
                else:
                    raise e
                return None
            print(actual)
            assert actual == expected
            return actual

    def __call__(
        self,
        sql_content: SQLQuery,
        expected: AnyRowList,
        args = None,
    ) -> AnyRowList:
         return self.assert_sql(
             sql_content,
             expected,
             args,
         )

class BioseqdbTestDatabase():

    postgresql_instance: Optional[Database]
    engine: Optional[Engine]

    def __init__(
        self,
    ):
        self.postgresql_instance = None
        self.engine = None

    def install_bioseqdb(self, engine: Engine):
        with engine.connect() as con:
            sql_statement = text("CREATE EXTENSION IF NOT EXISTS bioseqdb;")
            con.execute(sql_statement)

    def create_engine(self) -> Engine:
        engine = create_engine_raw(self.postgresql_instance.url())
        self.install_bioseqdb(engine)
        return engine

    def __enter__(self) -> BioseqdbTestEngine:
        if self.engine is not None:
            return self.engine
        if self.postgresql_instance is None:
            self.postgresql_instance = testing.postgresql.Postgresql()
        self.postgresql_instance.__enter__()
        self.engine = self.create_engine()
        return BioseqdbTestEngine(self.engine)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.postgresql_instance.__exit__(self, exc_type, exc_val, exc_tb)
        self.postgresql_instance = None
        self.engine = None