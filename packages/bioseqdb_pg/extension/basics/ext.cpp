#include "common.h"
#include "sequence.h"

POSTGRES_FUNCTION(nuclseq_len)
{
    auto nucls = reinterpret_cast<const NucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)));
    PG_RETURN_UINT64(nucls->length());
}
POSTGRES_FUNCTION_END()

POSTGRES_FUNCTION(nuclseq_content)
{
    auto nucls = reinterpret_cast<NucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)));
    std::string_view needle = PG_GETARG_CSTRING(1);
    if (needle.length() != 1 || std::find(allowed_nucleotides.begin(), allowed_nucleotides.end(), needle[0]) == allowed_nucleotides.end()) {
        raise_pg_error(ERRCODE_INVALID_PARAMETER_VALUE,
                errmsg("invalid nucleotide in nuclseq_content: '%s'", needle.data()));
    }

    auto matches = static_cast<double>(nucls->occurences(needle[0]));
    PG_RETURN_FLOAT8(matches / nucls->length());
}
POSTGRES_FUNCTION_END()

POSTGRES_FUNCTION(nuclseq_complement)
{
    auto nucls = reinterpret_cast<const NucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)));
    PG_RETURN_POINTER(nucls->complement());
}
POSTGRES_FUNCTION_END()

POSTGRES_FUNCTION(nuclseq_reverse)
{
    auto nucls = reinterpret_cast<const NucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)));
    PG_RETURN_POINTER(nucls->reverse());
}
POSTGRES_FUNCTION_END()
