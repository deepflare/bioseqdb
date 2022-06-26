from db_tests.storage.base_model import OrmBaseModel, BaseModel, FIELDS

import pandas as pd

from typing import List, Dict, Any

from sqlalchemy import Text, BigInteger, Integer, Boolean

from db_tests.storage.types import NuclSeq

# ref_id BIGINT,
#     ref_subseq nucl_seq,
#     ref_match_start INTEGER,
#     ref_match_end INTEGER,
#     ref_match_len INTEGER,
#     query_id BIGINT,
#     query_subseq nucl_seq,
#     query_match_start INTEGER,
#     query_match_end INTEGER,
#     query_match_len INTEGER,
#     is_primary BOOLEAN,
#     is_secondary BOOLEAN,
#     is_reverse BOOLEAN,
#     cigar TEXT,
#     score INTEGER

class BWAAlignment(OrmBaseModel):
    """
        XYZ
    """
    ### Target
    __ignore_attrs_eq__ = ["query_id", "ref_id"]
    ref_id: int = FIELDS.column(
        BigInteger(),
        primary_key=True,
        nullable=False,
        description="Target sequence primary key",
    )
    ref_subseq: NuclSeq = FIELDS.column(
        NuclSeq(),
        nullable=False,
        description="Target sequence content"
    )
    ref_match_start: int = FIELDS.column(
        Integer(),
        nullable=False,
        description="Target sequence match start position",
    )
    ref_match_end: int = FIELDS.column(
        Integer(),
        nullable=False,
        description="Target sequence match end position",
    )
    ref_match_len: int = FIELDS.column(
        Integer(),
        nullable=False,
        description="Target sequence match length",
    )
    ### Query
    query_id: int = FIELDS.column(
        BigInteger(),
        nullable=True,
        description="Query sequence primary ID",
    )
    query_subseq: NuclSeq = FIELDS.column(
        NuclSeq(),
        nullable=False,
        description="Query sequence content"
    )
    query_match_start: int = FIELDS.column(
        Integer(),
        nullable=False,
        description="Query sequence match start position",
    )
    query_match_end: int = FIELDS.column(
        Integer(),
        nullable=False,
        description="Query sequence match end position",
    )
    query_match_len: int = FIELDS.column(
        Integer(),
        nullable=False,
        description="Query sequence match length",
    )
    ### BWA specific
    is_primary: bool = FIELDS.column(
        Boolean(),
        nullable=False,
        description="Is primary",
    )
    is_secondary: bool = FIELDS.column(
        Boolean(),
        nullable=False,
        description="Is secondary",
    )
    is_reverse: bool = FIELDS.column(
        Boolean(),
        nullable=False,
        description="Is reverse",
    )
    cigar: str = FIELDS.column(
        Text(),
        nullable=False,
        description="Cigar sequence alignment representation",
    )
    score: int = FIELDS.column(
        Integer(),
        nullable=False,
        description="Match score",
    )


    