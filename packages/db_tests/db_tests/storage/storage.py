from sqlalchemy import Column, Integer, String, Date, create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.ext.mutable import MutableDict
from sqlalchemy import ForeignKey
from db_tests.utils import get_postgresql_env_variables
from db_tests.storage.query import ExtendedQuery
from common_api_utils import LOGS

from typing import Optional, Type

#pg_user, pg_pwd, pg_external_host, pg_db = get_postgresql_env_variables()


class StorageProvider:
    _engine: Optional[Engine]
    _sessionmaker: Optional[sessionmaker]
    _session: Optional[Session]

    def __init__(self) -> None:
        self._engine = None
        self._sessionmaker = None
        self._session = None

    def _get_session(self) -> Session:
        if self._engine is None or self._sessionmaker is None:
            # an Engine, which the Session will use for connection resources
            self._engine = create_engine(
                "sqlite:///foo.db",
                #f"postgresql+psycopg2://{pg_user}:{pg_pwd}@{pg_external_host}/{pg_db}"
            )
            # Create models
            self._sessionmaker = sessionmaker(bind=self._engine, query_cls=ExtendedQuery)
        return self._sessionmaker()

    def _init_db(self, base_model_cls):
        _ = self._get_session()
        base_model_cls.metadata.create_all(self._engine)

    @property
    def db(self):
        if self._session is None:
            self._session = self._get_session()
        return self._session


STORAGE = StorageProvider()