#include "common.h"
#include "sequence.h"

// Lowercase nucleotides should not be allowed to be stored in the database. Their meaning in non-standardized, and some
// libraries can handle them poorly (for example, by replacing them with Ns). They should be handled before importing
// them into the database, in order to make the internals more robust and prevent accidental usage. A valid option when
// importing is replacing them with uppercase ones, as their most common use is for repeating but valid nucleotides.
POSTGRES_FUNCTION(nuclseq_in)
{
    std::string_view text = PG_GETARG_CSTRING(0);

    if (text.length() > INT32_MAX / 4)
        raise_pg_error(ERRCODE_INVALID_PARAMETER_VALUE, errmsg("provided sequence is too long"));

    for (char chr : text) {
        if (std::find(allowed_nucleotides.begin(), allowed_nucleotides.end(), chr) == allowed_nucleotides.end()) {
            raise_pg_error(ERRCODE_INVALID_TEXT_REPRESENTATION,
                    errmsg("invalid nucleotide in nuclseq_in: '%c'", chr));
        }
    }

    PG_RETURN_POINTER(nuclseq_from_text(text));
}
POSTGRES_FUNCTION_END()

POSTGRES_FUNCTION(nuclseq_out)
{
    auto nucls = reinterpret_cast<const NucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)));
    PG_RETURN_CSTRING(nucls->to_text_palloc());
}
POSTGRES_FUNCTION_END()

/////////

POSTGRES_FUNCTION(aa_seq_in)
{
    std::string_view text = PG_GETARG_CSTRING(0);

    if (text.length() > INT32_MAX / 4)
        raise_pg_error(ERRCODE_INVALID_PARAMETER_VALUE, errmsg("provided sequence is too long"));

    for (char chr : text) {
        if (std::find(allowed_nucleotides.begin(), allowed_nucleotides.end(), chr) == allowed_nucleotides.end()) {
            raise_pg_error(ERRCODE_INVALID_TEXT_REPRESENTATION,
                    errmsg("invalid nucleotide in nuclseq_in: '%c'", chr));
        }
    }

    PG_RETURN_POINTER(nuclseq_from_text(text));
}
POSTGRES_FUNCTION_END()

POSTGRES_FUNCTION(aa_seq_out)
{
    auto nucls = reinterpret_cast<const NucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)));
    PG_RETURN_CSTRING(nucls->to_text_palloc());
}
POSTGRES_FUNCTION_END()

/////

POSTGRES_FUNCTION(amb_aa_seq_in)
{
    std::string_view text = PG_GETARG_CSTRING(0);

    if (text.length() > INT32_MAX / 4)
        raise_pg_error(ERRCODE_INVALID_PARAMETER_VALUE, errmsg("provided sequence is too long"));

    for (char chr : text) {
        if (std::find(allowed_nucleotides.begin(), allowed_nucleotides.end(), chr) == allowed_nucleotides.end()) {
            raise_pg_error(ERRCODE_INVALID_TEXT_REPRESENTATION,
                    errmsg("invalid nucleotide in nuclseq_in: '%c'", chr));
        }
    }

    PG_RETURN_POINTER(nuclseq_from_text(text));
}
POSTGRES_FUNCTION_END()

POSTGRES_FUNCTION(amb_aa_seq_out)
{
    auto nucls = reinterpret_cast<const NucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)));
    PG_RETURN_CSTRING(nucls->to_text_palloc());
}
POSTGRES_FUNCTION_END()