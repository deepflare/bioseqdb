#include <algorithm>
#include <array>
#include <chrono>
#include <charconv>
#include <string>
#include <string_view>
#include <optional>
#include <stdint.h>
#include <cstdlib>
#include <mutex>
#include <queue>
#include <thread>
#include <condition_variable>
#include <algorithm>

extern "C" {
#include <postgres.h>
#include <fmgr.h>
#include <funcapi.h>
#include <miscadmin.h>
#include <executor/spi.h>
#include <catalog/pg_type.h>
}

#include "bwa.h"
#include "sequence.h"

#define raise_pg_error(code, msg) ereport(ERROR, (errcode(code)), msg);

namespace {

text *string_view_to_text(std::string_view s) {
    text *result = (text *) palloc(s.size() + VARHDRSZ);
    SET_VARSIZE(result, s.size() + VARHDRSZ);
    memcpy(VARDATA(result), s.data(), s.size());
    return result;
}

}

extern "C" {

PG_MODULE_MAGIC;

// Lowercase nucleotides should not be allowed to be stored in the database. Their meaning in non-standardized, and some
// libraries can handle them poorly (for example, by replacing them with Ns). They should be handled before importing
// them into the database, in order to make the internals more robust and prevent accidental usage. A valid option when
// importing is replacing them with uppercase ones, as their most common use is for repeating but valid nucleotides.
PG_FUNCTION_INFO_V1(nuclseq_in);
Datum nuclseq_in(PG_FUNCTION_ARGS) {
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

PG_FUNCTION_INFO_V1(nuclseq_out);
Datum nuclseq_out(PG_FUNCTION_ARGS) {
    auto nucls = reinterpret_cast<RawNucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)))->wrapped();
    PG_RETURN_CSTRING(nucls.to_text_malloc());
}

PG_FUNCTION_INFO_V1(nuclseq_len);
Datum nuclseq_len(PG_FUNCTION_ARGS) {
    auto nucls = reinterpret_cast<RawNucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)))->wrapped();
    PG_RETURN_UINT64(nucls.length());
}

PG_FUNCTION_INFO_V1(nuclseq_content);
Datum nuclseq_content(PG_FUNCTION_ARGS) {
    auto nucls = reinterpret_cast<RawNucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)))->wrapped();
    std::string_view needle = PG_GETARG_CSTRING(1);
    if (needle.length() != 1 || std::find(allowed_nucleotides.begin(), allowed_nucleotides.end(), needle[0]) == allowed_nucleotides.end()) {
        raise_pg_error(ERRCODE_INVALID_PARAMETER_VALUE,
                errmsg("invalid nucleotide in nuclseq_content: '%s'", needle.data()));
    }

    auto matches = static_cast<double>(nucls.occurences(needle[0]));
    PG_RETURN_FLOAT8(matches / nucls.length());
}

PG_FUNCTION_INFO_V1(nuclseq_complement);
Datum nuclseq_complement(PG_FUNCTION_ARGS) {
    auto nucls = reinterpret_cast<RawNucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)))->wrapped();
    PG_RETURN_POINTER(nucls.complement());
}

PG_FUNCTION_INFO_V1(nuclseq_reverse);
Datum nuclseq_reverse(PG_FUNCTION_ARGS) {
    auto nucls = reinterpret_cast<RawNucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)))->wrapped();
    PG_RETURN_POINTER(nucls.reverse());
}

}

namespace {

template<typename F>
Portal iterate_nuclseq_table(const char* sql, Oid nuclseq_oid, F f) {
    Portal portal = SPI_cursor_open_with_args(nullptr, sql, 0, nullptr, nullptr, nullptr, true, 0);
    long batch_size = 1;

    SPI_cursor_fetch(portal, true, batch_size);
    while (SPI_processed > 0 && SPI_tuptable != NULL) {
        int n = SPI_processed;
        SPITupleTable* tuptable = SPI_tuptable;
        TupleDesc tupdesc = tuptable->tupdesc;

        switch(SPI_gettypeid(tupdesc, 1)) {
            case INT2OID:
            case INT4OID:
            case INT8OID:
                break;
            default:
            raise_pg_error(ERRCODE_DATATYPE_MISMATCH, errmsg("expected column of integers"));
        }

        if (SPI_gettypeid(tupdesc, 2) != nuclseq_oid)
            raise_pg_error(ERRCODE_DATATYPE_MISMATCH, errmsg("expected column of nuclseqs"));

        for(int i = 0 ; i < n; i++) {
            HeapTuple tup = tuptable->vals[i];
            bool null_id = false, null_seq = false;

            Datum id = SPI_getbinval(tup, tupdesc, 1, &null_id);
            Datum nucls = SPI_getbinval(tup, tupdesc, 2, &null_seq);

            elog(NOTICE, "bwa search supplied with null id or sequence");

            if (!null_id && !null_seq)
                f(static_cast<int64_t>(id), reinterpret_cast<RawNucleotideSequence*>(PG_DETOAST_DATUM(nucls))->wrapped());
        }

        SPI_freetuptable(tuptable);
        SPI_cursor_fetch(portal, true, batch_size);
    }
    return portal;

}


BwaIndex bwa_index_from_query(const char* sql, Oid nuclseq_oid) {
    BwaIndex bwa;
    size_t count = 0;

    Portal portal = iterate_nuclseq_table(sql, nuclseq_oid, [&](auto id, auto nucls){
        bwa.add_ref_sequence(id, nucls);
        count++;
    });
    SPI_cursor_close(portal);
    bwa.options->max_occ = std::max<int>(500, count * 2);
    bwa.build();

    return bwa;
}

void assert_can_return_set(ReturnSetInfo* rsi) {
    if (rsi == NULL || !IsA(rsi, ReturnSetInfo)) {
        raise_pg_error(ERRCODE_FEATURE_NOT_SUPPORTED,
                errmsg("set-valued function called in context that cannot accept a set"));
    }
    if (!(rsi->allowedModes & SFRM_Materialize)) {
        raise_pg_error(ERRCODE_FEATURE_NOT_SUPPORTED,
                errmsg("materialize mode required, but it is not allowed in this context"));
    }
}

TupleDesc get_retval_tupledesc(const FunctionCallInfo& fcinfo) {
    TupleDesc tupledesc;

    switch (get_call_result_type(fcinfo, nullptr, &tupledesc)) {
        case TYPEFUNC_COMPOSITE:
            break;
        case TYPEFUNC_RECORD:
            raise_pg_error(ERRCODE_FEATURE_NOT_SUPPORTED,
                    errmsg("function returning record called in context that cannot accept type record"));
            break;
        default:
            raise_pg_error(ERRCODE_FEATURE_NOT_SUPPORTED,
                    errmsg("return type must be a row type"));
            break;
    }

    return tupledesc;
}

Tuplestorestate* create_tuplestore(ReturnSetInfo* rsi, TupleDesc& tupledesc) {
    MemoryContext per_query_ctx = rsi->econtext->ecxt_per_query_memory;
    MemoryContext old_ctx = MemoryContextSwitchTo(per_query_ctx);

    tupledesc = CreateTupleDescCopy(tupledesc);
    tupledesc = BlessTupleDesc(tupledesc);

    Tuplestorestate* tupstore = tuplestore_begin_heap(SFRM_Materialize_Random, false, work_mem);
    MemoryContextSwitchTo(old_ctx);

    return tupstore;
}

HeapTuple build_tuple_bwa(std::optional<int64_t> query_id, const BwaMatch& match, TupleDesc& tupledesc) {
    std::array<bool, 15> nulls;
    std::array<Datum, 15> values { {
        Int64GetDatum(match.ref_id),
        PointerGetDatum(nuclseq_from_text(match.ref_subseq)),
        Int32GetDatum(match.ref_match_begin),
        Int32GetDatum(match.ref_match_end),
        Int32GetDatum(match.ref_match_len),
        Int64GetDatum(query_id.value_or(0)),
        PointerGetDatum(nuclseq_from_text(match.query_subseq)),
        Int32GetDatum(match.query_match_begin),
        Int32GetDatum(match.query_match_end),
        Int32GetDatum(match.query_match_len),
        BoolGetDatum(match.is_primary),
        BoolGetDatum(match.is_secondary),
        BoolGetDatum(match.is_reverse),
        PointerGetDatum(string_view_to_text(match.cigar)),
        Int32GetDatum(match.score),
    } };

    nulls.fill(false);
    nulls[5] = !query_id.has_value();

    return heap_form_tuple(tupledesc, values.data(), nulls.data());
}

}

extern "C" {

PG_FUNCTION_INFO_V1(nuclseq_search_bwa);
Datum nuclseq_search_bwa(PG_FUNCTION_ARGS) {
    ReturnSetInfo* rsi = reinterpret_cast<ReturnSetInfo*>(fcinfo->resultinfo);
    assert_can_return_set(rsi);

    auto nucls = reinterpret_cast<RawNucleotideSequence*>(PG_DETOAST_DATUM(PG_GETARG_POINTER(0)))->wrapped();
    const char* reference_sql = PG_GETARG_CSTRING(1);

    if (int ret = SPI_connect(); ret < 0)
        elog(ERROR, "connectby: SPI_connect returned %d", ret);

    TupleDesc ret_tupdesc = get_retval_tupledesc(fcinfo);

    Oid nuclseq_oid = TupleDescAttr(ret_tupdesc, 1)->atttypid;
    BwaIndex bwa = bwa_index_from_query(reference_sql, nuclseq_oid);
    SPI_finish();

    Tuplestorestate* ret_tupstore = create_tuplestore(rsi, ret_tupdesc);
    AttInMetadata* attr_input_meta = TupleDescGetAttInMetadata(ret_tupdesc);

    std::vector<BwaMatch> aligns = bwa.align_sequence(nucls);

    for (BwaMatch& row : aligns) {
        HeapTuple tuple = build_tuple_bwa(std::nullopt, row, ret_tupdesc);
        tuplestore_puttuple(ret_tupstore, tuple);
        heap_freetuple(tuple);
    }

    rsi->returnMode = SFRM_Materialize;
    rsi->setResult = ret_tupstore;
    rsi->setDesc = ret_tupdesc;
    return (Datum) nullptr;
}

PG_FUNCTION_INFO_V1(nuclseq_multi_search_bwa);
Datum nuclseq_multi_search_bwa(PG_FUNCTION_ARGS) {
    ReturnSetInfo* rsi = reinterpret_cast<ReturnSetInfo*>(fcinfo->resultinfo);
    assert_can_return_set(rsi);

    const char* query_sql = PG_GETARG_CSTRING(0);
    const char* reference_sql = PG_GETARG_CSTRING(1);

    if (int ret = SPI_connect(); ret < 0)
        elog(ERROR, "connectby: SPI_connect returned %d", ret);

    TupleDesc ret_tupdesc = get_retval_tupledesc(fcinfo);
    Oid nuclseq_oid = TupleDescAttr(ret_tupdesc, 1)->atttypid;
    BwaIndex bwa = bwa_index_from_query(reference_sql, nuclseq_oid);
    Tuplestorestate* ret_tupstore = create_tuplestore(rsi, ret_tupdesc);
    AttInMetadata* attr_input_meta = TupleDescGetAttInMetadata(ret_tupdesc);


    std::mutex m_queries;
    std::queue<std::tuple<int64_t, char*, uint32_t>> queries;

    iterate_nuclseq_table(query_sql, nuclseq_oid, [&](auto id, auto nuclseq){
        queries.emplace(id, nuclseq.to_text_malloc(), nuclseq.length());
    });
    size_t queries_cnt = queries.size();

    unsigned threads_num = 6;
    std::vector<std::unique_ptr<std::thread>> threads;

    std::mutex m_matches;
    std::queue<std::pair<int64_t, std::vector<BwaMatch>>> matches;
    std::condition_variable cv;

    for(unsigned i = 0 ; i < threads_num ; i++) {
        threads.push_back(std::make_unique<std::thread>([&]() {
            m_queries.lock();
            while (!queries.empty()) {
                auto [id, seq, len] = queries.front();
                queries.pop();
                m_queries.unlock();

                std::vector<BwaMatch> aligns = bwa.align_string_view(std::string_view(seq, len));
                free(seq);

                m_matches.lock();
                matches.emplace(id, move(aligns));
                cv.notify_one();
                m_matches.unlock();

                m_queries.lock();
            }
            m_queries.unlock();
        }));
    }

    size_t matches_cnt = 0;
    while (matches_cnt < queries_cnt) {
        std::unique_lock lk(m_matches);
        cv.wait(lk, [&] { return !matches.empty(); });
        while(!matches.empty()) {
            const auto& [id, aligns] = matches.front();
            for (const BwaMatch& row : aligns) {
                HeapTuple tuple = build_tuple_bwa(id, row, ret_tupdesc);
                tuplestore_puttuple(ret_tupstore, tuple);
                heap_freetuple(tuple);
            }
            matches.pop();
            matches_cnt++;
        }

    }

    for(auto& t: threads)
        t->join();

    SPI_finish();

    rsi->returnMode = SFRM_Materialize;
    rsi->setResult = ret_tupstore;
    rsi->setDesc = ret_tupdesc;
    return (Datum) nullptr;
}

}
