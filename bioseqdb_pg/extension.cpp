#include <algorithm>
#include <string_view>

extern "C" {
#include <postgres.h>
#include <fmgr.h>
#include <funcapi.h>
}

struct PgNucleotideSequence {
    char vl_len[4];
    char nucleotides[];

    std::string_view text() const {
        return {nucleotides, VARSIZE(this) - 4};
    }

    static PgNucleotideSequence* palloc(size_t len) {
        auto ptr = static_cast<PgNucleotideSequence*>(::palloc(4 + len));
        SET_VARSIZE(ptr, 4 + len);
        return ptr;
    }
};

extern "C" {

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(nuclseq_in);
Datum nuclseq_in(PG_FUNCTION_ARGS) {
    std::string_view text = PG_GETARG_CSTRING(0);
    for (char chr : text) {
        if (chr != 'A' && chr != 'C' && chr != 'G' && chr != 'T') {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION)),
                    errmsg("invalid nucleotide in nuclseq_in: '%c'", chr));
        }
    }

    auto nuclseq = PgNucleotideSequence::palloc(text.size());
    std::copy(text.begin(), text.end(), nuclseq->nucleotides);
    PG_RETURN_POINTER(nuclseq);
}

PG_FUNCTION_INFO_V1(nuclseq_out);
Datum nuclseq_out(PG_FUNCTION_ARGS) {
    auto nucls = reinterpret_cast<PgNucleotideSequence*>(PG_GETARG_POINTER(0))->text();
    auto text = static_cast<char*>(palloc(nucls.size() + 1));
    std::copy(nucls.begin(), nucls.end(), text);
    text[nucls.size()] = '\0';
    PG_RETURN_CSTRING(text);
}

PG_FUNCTION_INFO_V1(nuclseq_len);
Datum nuclseq_len(PG_FUNCTION_ARGS) {
    auto nucls = reinterpret_cast<PgNucleotideSequence*>(PG_GETARG_POINTER(0))->text();
    PG_RETURN_UINT64(nucls.size());
}

PG_FUNCTION_INFO_V1(nuclseq_content);
Datum nuclseq_content(PG_FUNCTION_ARGS) {
    auto nucls = reinterpret_cast<PgNucleotideSequence*>(PG_GETARG_POINTER(0))->text();
    std::string_view needle = PG_GETARG_CSTRING(1);
    if (needle != "A" && needle != "C" && needle != "G" && needle != "T") {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE)),
                errmsg("invalid nucleotide in nuclseq_content: '%s'", needle.data()));
    }

    auto matches = static_cast<double>(std::count(nucls.begin(), nucls.end(), needle[0]));
    PG_RETURN_FLOAT8(matches / nucls.size());
}

PG_FUNCTION_INFO_V1(nuclseq_complement);
Datum nuclseq_complement(PG_FUNCTION_ARGS) {
    auto nucls = reinterpret_cast<PgNucleotideSequence*>(PG_GETARG_POINTER(0))->text();
    auto complement = PgNucleotideSequence::palloc(nucls.size());
    for (size_t i=0; i<nucls.size(); ++i) {
        if (nucls[i] == 'A') {
            complement->nucleotides[i] = 'C';
        } else if (nucls[i] == 'C') {
            complement->nucleotides[i] = 'A';
        } else if (nucls[i] == 'T') {
            complement->nucleotides[i] = 'G';
        } else if (nucls[i] == 'G') {
            complement->nucleotides[i] = 'T';
        }
    }
    PG_RETURN_POINTER(complement);
}

PG_FUNCTION_INFO_V1(yoyo_v1);
Datum yoyo_v1(PG_FUNCTION_ARGS) {
    if (SRF_IS_FIRSTCALL()) {
        FuncCallContext* funcctx = SRF_FIRSTCALL_INIT();
        MemoryContext oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

        int num_tuples = PG_GETARG_INT32(0);
        if (num_tuples < 0) {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE)),
                    errmsg("number of rows cannot be negative"));
        }
        funcctx->max_calls = num_tuples;

        MemoryContextSwitchTo(oldcontext);
    }

    FuncCallContext* funcctx = SRF_PERCALL_SETUP();
    uint64_t result = funcctx->call_cntr;
    if (funcctx->call_cntr < funcctx->max_calls) {
        SRF_RETURN_NEXT(funcctx, UInt64GetDatum(result));
    } else {
        SRF_RETURN_DONE(funcctx);
    }
}

}
