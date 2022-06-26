from sqlalchemy import Column, Integer, String, Date, Enum
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.dialects.sqlite import JSON
from sqlalchemy.ext.mutable import MutableDict
from sqlalchemy import ForeignKey
from db_tests.logs import LOGS


from datetime import datetime, date

from typing import Optional, Type, Any

class OrmFieldsHelper:
    use_nullable: bool

    def __init__(
        self,
        use_nullable: bool = False,
    ):
        self.use_nullable = use_nullable

    @property
    def optional(self):
        return OrmFieldsHelper(use_nullable=True)

    def column(
        self,
        column_type,
        description: str = "",
        regex: Optional[str] = None,
        type: Type[Any] = None,
        gt: Optional[Any] = None,
        lt: Optional[Any] = None,
        **kwargs,
    ) -> Column:
        metadata=dict(
            description=description,
            regex=regex,
            gt=gt,
            lt=lt,
        )
        if type:
            metadata["type"] = type
        col = Column(column_type, **kwargs)
        setattr(col, "_metadata", metadata)
        return col

    def _column(self, column_type, metadata=None, **kwargs) -> Column:
        if self.use_nullable:
            opts = dict(
                nullable=True,
                default=None,
            )
        else:
            opts = dict(
                nullable=False,
            )
        col = Column(column_type, **{**opts, **kwargs})
        setattr(col, "_metadata", metadata)
        return col
    
    def primary_id(
        self,
        description: str,
    ):
        if self.use_nullable:
            raise ValueError(f"Cannot use .optional.primary_id(). Primary ID cannot be optional (nullable)")
        return self._column(
            String(100),
            primary_key=True,
            nullable=False,
            metadata=dict(
                description=description,
            ),
        )

    def str(
        self,
        description: str,
        length: int = 255,
        regex: Optional[str] = None,
        unique: bool = False,
    ):
        return self._column(
            String(length),
            unique=unique,
            metadata=dict(
                regex=regex,
                description=description,
            ),
        )

    def enum(
        self,
        description: str,
        t: Type[Any],
    ):
        return self._column(
            Enum(t),
            metadata=dict(
                description=description,
                validator_post=lambda value, values, config, field: (self.optional and value is None) or isinstance(value, t)
            ),
        )

    def integer(
        self,
        description: str,
        min_value: Optional[int] = None,
        max_value: Optional[int] = None,
    ):
        return self._column(
            Integer,
            metadata=dict(
                description=description,
                gt=min_value,
                lt=max_value,
            ),
        )
    
    def date(
        self,
        description: str,
        format: str,
    ):
        def _validate_date(value):
            if self.optional and value is None:
                return None
            if isinstance(value, date):
                return value
            return datetime.strptime(value, format).date()
        return self._column(
            Date,
            metadata=dict(
                description=description,
                validator=lambda value, values, config, field: _validate_date(value)
            ),
        )
    
    def dict(
        self,
        description: str,
        t: Type[Any],
    ):
        return self._column(
            MutableDict.as_mutable(JSON),
            metadata=dict(
                description=description,
                type=t,
            ),
        )