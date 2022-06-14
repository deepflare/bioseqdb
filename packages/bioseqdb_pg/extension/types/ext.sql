CREATE TYPE nucl_seq;

CREATE FUNCTION nuclseq_in(cstring)
    RETURNS nucl_seq
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION nuclseq_out(nucl_seq)
    RETURNS CSTRING
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE nucl_seq (
    internallength = VARIABLE,
    storage = EXTENDED,
	alignment = double,
    input = nuclseq_in,
    output = nuclseq_out
);

-------

CREATE TYPE amb_aa_seq;

CREATE FUNCTION amb_aa_seq_in(cstring)
    RETURNS nucl_seq
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION amb_aa_seq_out(nucl_seq)
    RETURNS CSTRING
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE amb_aa_seq (
    internallength = VARIABLE,
    storage = EXTENDED,
	alignment = double,
    input = amb_aa_seq_in,
    output = amb_aa_seq_out
);

-------


CREATE TYPE aa_seq;

CREATE FUNCTION aa_seq_in(cstring)
    RETURNS aa_seq
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION aa_seq_out(aa_seq)
    RETURNS CSTRING
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE aa_seq (
    internallength = VARIABLE,
    storage = EXTENDED,
	alignment = double,
    input = aa_seq_in,
    output = aa_seq_out
);