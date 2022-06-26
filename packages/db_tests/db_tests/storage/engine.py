from sqlalchemy import create_engine as create_engine_raw
from sqlalchemy.sql import text
from db_tests.storage.records_list import RecordsList
import testing.postgresql
from testing.common.database import Database
from sqlalchemy.engine import Engine
from pgsanity.pgsanity import check_string as validate_pg_sql

from sqlalchemy.orm import Session
from sqlalchemy.exc import ResourceClosedError
from typing import Callable, Optional, List, Tuple, Any, Union
from sqlalchemy.orm import sessionmaker


AnyRowList = Optional[Tuple[Any, ...]]
SQLQuery = str

class BioseqdbTestEngine():
    engine: Engine
    session_maker: sessionmaker

    def __init__(self, engine: Engine):
        self.engine = engine
        self.session_maker = sessionmaker(bind=engine)

    def query(
        self,
        sql_content: SQLQuery,
        # expected: AnyRowList,
        args = None,
        model = None,
        ignore_order: bool = False,
    ) -> RecordsList:
        #success, msg = validate_pg_sql(sql_content)
        #print(msg)
        #if not success:
        #    raise ValueError(f"Invalid SQL query: {sql_content}")
        if not args:
            args = dict()
        with self.session_maker() as con:
            try:
                statement = text(sql_content)
                r = con.execute(statement, args)
                con.commit()
                actual = list(r)
                if model:
                    col_keys = model.__table__.columns.keys()
                    actual = [model(**dict(zip(col_keys, row))) for row in actual]
                actual = RecordsList(content=actual, ignore_order=ignore_order)
                #for x in actual.content:
                #    print(x)
            except ResourceClosedError:
                actual = None
            except Exception as e:
                return e.__class__
            return actual

    def __call__(
        self,
        sql_content: SQLQuery,
        args = None,
        model = None,
        ignore_order: bool = False,
    ) -> AnyRowList:
         return self.query(
             sql_content,
             args=args,
             model=model,
             ignore_order=ignore_order,
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