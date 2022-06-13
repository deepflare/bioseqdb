#include <vector>
#include <string>
#include <exception>
#include <signal.h>
#include "../common/mmseq2lib.h"
#include "rpc/client.h"
#include "rpc/rpc_error.h"
#include "config.h"

extern "C"
{
#include <postgres.h>
#include <fmgr.h>
#include <executor/spi.h>    // SPI_*
#include <funcapi.h>         // get_call_result_type, AttInMetadata, TupleDescGetAttInMetadata, BuildTupleFromCStrings
#include <miscadmin.h>       // work_mem
#include <utils/array.h>     // deconstruct_array, ARR_ELEMTYPE, PG_GETARG_ARRAYTYPE_P
#include <utils/builtins.h>  // text_to_cstring
#include <utils/lsyscache.h> // get_typlenbyvalalign
}

#define TEST_ARGS_FOR_NULL(n, pred) \
for (uint32_t i = 0; i < n; i++) \
    if (PG_ARGISNULL(i) && pred) \
        PG_RETURN_NULL()

#define DECLARE_VECTORS_AND_TARGET_NAMES() \
common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{}); \
common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{}); \
common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{}); \
common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{}); \
std::optional<std::string> target_tblname = std::nullopt; \
std::optional<std::string> target_colname = std::nullopt; \
bool all_targets = false; \
uint32_t tCount

#define GET_OPTIONAL_PARAMS(n) \
std::string subst_matrix_name = text_to_cstring(PG_GETARG_TEXT_PP(n)); \
uint32_t kmer_length = PG_GETARG_INT32(n + 1); \
uint32_t kmer_gen_threshold = PG_GETARG_INT32(n + 2); \
uint32_t ungapped_alignment_score = PG_GETARG_INT32(n + 3); \
double eval_threshold = PG_GETARG_FLOAT8(n + 4); \
uint32_t gap_open_cost = PG_GETARG_INT32(n + 5); \
uint32_t gap_penalty_cost = PG_GETARG_INT32(n + 6); \
uint32_t thread_number = PG_GETARG_INT32(n + 7)

#define ADD_QUERIES_ONE(n) add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(n)), qIds, queries)

#define ADD_TARGETS_ONE(n) add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(n)), tIds, targets)

#define ADD_QUERIES_ARR(n) add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(n), qIds, queries)

#define ADD_TARGETS_ARR(n) add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(n), tIds, targets)

#define ADD_QUERIES_DB(n, a, b) \
std::string query_tblname = text_to_cstring(PG_GETARG_TEXT_PP(a)); \
std::string query_colname = text_to_cstring(PG_GETARG_TEXT_PP(b)); \
if (!check_valid_name(query_tblname)) \
    elog(ERROR, "%s", "Invalid query table name!"); \
if (!check_valid_name(query_colname)) \
    elog(ERROR, "%s", "Invalid query column name!"); \
if (PG_ARGISNULL(n)) \
    add_queries_all(query_tblname, query_colname, qIds, queries); \
else \
    add_queries_with_ids(query_tblname, query_colname, qIds, queries, PG_GETARG_ARRAYTYPE_P(n))

#define ADD_TARGETS_DB(n, a, b) \
target_tblname = text_to_cstring(PG_GETARG_TEXT_PP(a)); \
target_colname = text_to_cstring(PG_GETARG_TEXT_PP(b)); \
if (!check_valid_name(target_tblname.value())) \
    elog(ERROR, "%s", "Invalid target table name!"); \
if (!check_valid_name(target_colname.value())) \
    elog(ERROR, "%s", "Invalid target column name!"); \
tCount = count_targets(target_tblname.value()); \
if (PG_ARGISNULL(n)) \
{ \
    all_targets = true; \
} \
else \
{ \
    add_targets_with_ids(target_tblname.value(), target_colname.value(), tIds, PG_GETARG_ARRAYTYPE_P(n)); \
}

#define SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION() \
Tuplestorestate *tupstore; \
AttInMetadata *attinmeta; \
enable_materialize_mode(fcinfo, tupstore, attinmeta); \
seq_search_mmseqs_main(target_tblname, target_colname, \
                       qIds, tIds, queries, targets, tCount, \
                       kmer_length, \
                       subst_matrix_name, \
                       kmer_gen_threshold, \
                       ungapped_alignment_score, \
                       eval_threshold, \
                       gap_open_cost, \
                       gap_penalty_cost, \
                       thread_number, \
                       all_targets, \
                       seqType, \
                       enableAmbiguity, \
                       tupstore, \
                       attinmeta \
                       ); \
return (Datum)0

#define DEFINE_PG_FUNCTIONS_FOR_TYPE(prefix, seqType, ambiguity) \
PG_FUNCTION_INFO_V1(prefix##_search_one_to_one); \
PG_FUNCTION_INFO_V1(prefix##_search_one_to_arr); \
PG_FUNCTION_INFO_V1(prefix##_search_one_to_db); \
PG_FUNCTION_INFO_V1(prefix##_search_arr_to_one); \
PG_FUNCTION_INFO_V1(prefix##_search_arr_to_arr); \
PG_FUNCTION_INFO_V1(prefix##_search_arr_to_db); \
PG_FUNCTION_INFO_V1(prefix##_search_db_to_one); \
PG_FUNCTION_INFO_V1(prefix##_search_db_to_arr); \
PG_FUNCTION_INFO_V1(prefix##_search_db_to_db); \
Datum prefix##_search_one_to_one(PG_FUNCTION_ARGS) \
{ \
    return seq_search_mmseqs_one_to_one(fcinfo, seqType, ambiguity); \
} \
Datum prefix##_search_one_to_arr(PG_FUNCTION_ARGS) \
{ \
    return seq_search_mmseqs_one_to_arr(fcinfo, seqType, ambiguity); \
} \
Datum prefix##_search_one_to_db(PG_FUNCTION_ARGS) \
{ \
    return seq_search_mmseqs_one_to_db(fcinfo, seqType, ambiguity); \
} \
Datum prefix##_search_arr_to_one(PG_FUNCTION_ARGS) \
{ \
    return seq_search_mmseqs_arr_to_one(fcinfo, seqType, ambiguity); \
} \
Datum prefix##_search_arr_to_arr(PG_FUNCTION_ARGS) \
{ \
    return seq_search_mmseqs_arr_to_arr(fcinfo, seqType, ambiguity); \
} \
Datum prefix##_search_arr_to_db(PG_FUNCTION_ARGS) \
{ \
    return seq_search_mmseqs_arr_to_db(fcinfo, seqType, ambiguity); \
} \
Datum prefix##_search_db_to_one(PG_FUNCTION_ARGS) \
{ \
    return seq_search_mmseqs_db_to_one(fcinfo, seqType, ambiguity); \
} \
Datum prefix##_search_db_to_arr(PG_FUNCTION_ARGS) \
{ \
    return seq_search_mmseqs_db_to_arr(fcinfo, seqType, ambiguity); \
} \
Datum prefix##_search_db_to_db(PG_FUNCTION_ARGS) \
{ \
    return seq_search_mmseqs_db_to_db(fcinfo, seqType, ambiguity); \
}

namespace
{
    constexpr uint32_t OUT_TUPLE_ARITY = 18;

    uint32_t count_targets(const std::string &target_tblname)
    {
        const std::string countTargetsQuery =
            std::string("SELECT max(id) + 1 FROM ") +
            target_tblname +
            std::string(";");
        
        SPI_connect();
        SPI_exec(countTargetsQuery.data(), 0);

        const TupleDesc spi_tupdesc = SPI_tuptable->tupdesc;
        const HeapTuple spi_tuple = SPI_tuptable->vals[0];
        const std::string tCountStr = SPI_getvalue(spi_tuple, spi_tupdesc, 1);
        const uint32_t tCount = std::stol(tCountStr);

        SPI_finish();
        return tCount;
    }

    void add_one_sequence(const std::string &sequence,
                          const common::InputParams::Vec64Ptr &ids, const common::InputParams::VecStrPtr &sequences)
    {
        ids->push_back(1);
        sequences->push_back(std::make_shared<std::string>(sequence));
    }

    void add_array_of_sequences(ArrayType *const sequences_array_type,
                                const common::InputParams::Vec64Ptr &ids, const common::InputParams::VecStrPtr &sequences)
    {
        // Get the element type of the array of queries
        const Oid element_type = ARR_ELEMTYPE(sequences_array_type);
        int16 elem_type_width;
        bool elem_type_by_val;
        char elem_type_align;

        // Deconstruct the array
        Datum *sequences_array;
        int sequences_num;
        get_typlenbyvalalign(element_type, &elem_type_width, &elem_type_by_val, &elem_type_align);
        deconstruct_array(sequences_array_type, element_type, elem_type_width, elem_type_by_val, elem_type_align,
                          &sequences_array, NULL, &sequences_num);

        for (int i = 0; i < sequences_num; i++)
        {
            ids->push_back(i + 1);
            const std::string sequence = TextDatumGetCString(sequences_array[i]);
            sequences->push_back(std::make_shared<std::string>(sequence));
        }
    }

    void add_queries_all(const std::string &query_tblname, const std::string &query_colname,
                         const common::InputParams::Vec64Ptr &qIds, const common::InputParams::VecStrPtr &queries)
    {
        // Perform a query via SPI
        const std::string getQueriesQuery =
            std::string("SELECT id, ") +
            query_colname +
            std::string(" FROM ") +
            query_tblname +
            std::string(";");
        SPI_connect();
        SPI_exec(getQueriesQuery.data(), 0);

        // Move data from SPI_tuptable to std::vectors qIds and queries
        const TupleDesc spi_tupdesc = SPI_tuptable->tupdesc;
        const uint64_t processed = SPI_processed;
        for (uint32_t i = 0; i < processed; i++)
        {
            const HeapTuple spi_tuple = SPI_tuptable->vals[i];
            const std::string id_str = SPI_getvalue(spi_tuple, spi_tupdesc, 1);
            const std::string sequence = SPI_getvalue(spi_tuple, spi_tupdesc, 2);
            qIds->push_back(std::stol(id_str));
            queries->push_back(std::make_shared<std::string>(sequence));
        }

        SPI_finish();
    }

    void add_queries_with_ids(const std::string &query_tblname, const std::string &query_colname,
                              const common::InputParams::Vec64Ptr &qIds, const common::InputParams::VecStrPtr &queries,
                              ArrayType *const queries_array_type)
    {
        // Get the element type of the array of queries
        const Oid element_type = ARR_ELEMTYPE(queries_array_type);
        int16 elem_type_width;
        bool elem_type_by_val;
        char elem_type_align;

        // Deconstruct the array
        uint64_t *queries_array;
        int queries_num;
        get_typlenbyvalalign(element_type, &elem_type_width, &elem_type_by_val, &elem_type_align);
        deconstruct_array(queries_array_type, element_type, elem_type_width, elem_type_by_val, elem_type_align,
                          &queries_array, NULL, &queries_num);

        SPI_connect();
        for (int i = 0; i < queries_num; i++)
        {
            qIds->push_back(queries_array[i]);
            const std::string getOneQueryQuery =
                std::string("SELECT ") +
                query_colname +
                std::string(" FROM ") +
                query_tblname +
                std::string(" WHERE id = ") +
                std::to_string(queries_array[i]) +
                std::string(";");
            SPI_exec(getOneQueryQuery.data(), 0);
            const TupleDesc spi_tupdesc = SPI_tuptable->tupdesc;
            const HeapTuple spi_tuple = SPI_tuptable->vals[0];
            const std::string sequence = SPI_getvalue(spi_tuple, spi_tupdesc, 1);
            queries->push_back(std::make_shared<std::string>(sequence));
        }
        SPI_finish();

        // Free the allocated memory
        pfree(queries_array);
    }

    void add_targets_all(const std::string &target_tblname, const std::string &target_colname,
                         const common::InputParams::Vec64Ptr &tIds)
    {
        const std::string getTargetsQuery =
            std::string("SELECT id FROM ") +
            target_tblname +
            std::string(";");
        SPI_connect();
        SPI_exec(getTargetsQuery.data(), 0);

        const TupleDesc spi_tupdesc = SPI_tuptable->tupdesc;
        const uint64_t processed = SPI_processed;
        for (uint32_t i = 0; i < processed; i++)
        {
            const HeapTuple spi_tuple = SPI_tuptable->vals[i];
            const std::string id_str = SPI_getvalue(spi_tuple, spi_tupdesc, 1);
            tIds->push_back(std::stol(id_str));
        }
        SPI_finish();
    }

    void add_targets_with_ids(const std::string &target_tblname, const std::string &target_colname,
                              const common::InputParams::Vec64Ptr &tIds,
                              ArrayType *const targets_array_type)
    {
        const Oid element_type = ARR_ELEMTYPE(targets_array_type);
        int16 elem_type_width;
        bool elem_type_by_val;
        char elem_type_align;

        uint64_t *targets_array;
        int targets_num;
        get_typlenbyvalalign(element_type, &elem_type_width, &elem_type_by_val, &elem_type_align);
        deconstruct_array(targets_array_type, element_type, elem_type_width, elem_type_by_val, elem_type_align,
                          &targets_array, NULL, &targets_num);

        for (int i = 0; i < targets_num; i++)
        {
            tIds->push_back(targets_array[i]);
        }
    }

    bool check_valid_name(std::string s)
    {
        if (s.empty() || (s[0] != '_' && !isalpha(s[0])))
            return false;
        for (const char c : s)
        {
            if (c != '_' && !isalpha(c) && !isdigit(c))
                return false;
        }
        return true;
    }

    void sigint_handler(int signo)
    {
        if (signo == SIGINT)
            elog(ERROR, "%s", "Received SIGINT");
    }

    void seq_search_mmseqs_main(std::optional<std::string> target_tblname,
                                std::optional<std::string> target_colname,
                                common::InputParams::Vec64Ptr qIds,
                                common::InputParams::Vec64Ptr tIds,
                                common::InputParams::VecStrPtr queries,
                                common::InputParams::VecStrPtr targets,
                                const uint32_t tCount,
                                const uint32_t kmer_length,
                                const std::string &subst_matrix_name,
                                const uint32_t kmer_gen_threshold,
                                const uint32_t ungapped_alignment_score,
                                const double eval_threshold,
                                const uint32_t gap_open_cost,
                                const uint32_t gap_penalty_cost,
                                const uint32_t thread_number,
                                const bool all_targets,
                                char seqType,
                                bool enableAmbiguity,
                                Tuplestorestate *tupstore,
                                AttInMetadata *attinmeta)
    {
        // Add a sigint handler
        struct sigaction act;
        act.sa_handler = sigint_handler;
        if (sigaction(SIGINT, &act, NULL) == -1)
            elog(ERROR, "%s", "Couldn't set SIGINT handler");

        // Ad-hoc/local targets
        const bool local_targets = target_tblname == std::nullopt;
        uint32_t no_of_targets = 0;
        if (local_targets)
        {
            target_tblname = "";
            target_colname = "";
            no_of_targets = tIds->size();
        }
        else 
        {
            no_of_targets = tCount;
        }

        // Input params
        common::InputParams input_params(qIds->size(), no_of_targets, qIds, tIds, queries,
                                         all_targets, local_targets, targets,
                                         std::make_shared<std::string>(target_tblname.value()),
                                         std::make_shared<std::string>(target_colname.value()),
                                         std::make_shared<std::string>(subst_matrix_name),
                                         kmer_length,
                                         kmer_gen_threshold,
                                         ungapped_alignment_score,
                                         eval_threshold,
                                         gap_open_cost,
                                         gap_penalty_cost,
                                         thread_number,
                                         seqType,
                                         enableAmbiguity);

        // Client
        rpc::client c(MMSEQ_HOSTNAME, MMSEQ_PORT);
        c.set_timeout(MMSEQ_TIMEOUT);
        common::VecRes mmseq_result;
        try {
            c.call("mmseq2", input_params).get().convert(mmseq_result);
        }
        catch (rpc::rpc_error &e) {
            std::string error_msg;
            e.get_error().get().convert(error_msg);
            elog(ERROR, "%s", error_msg.data());
        }
        catch (const std::exception &e) {
            std::string error_msg = e.what();
            elog(ERROR, "%s", error_msg.data());
        }
        uint32_t n = mmseq_result.size();

        // Build the output table
        for (uint32_t i = 0; i < n; i++)
        {
            common::MmseqResult t = mmseq_result[i];
            char **values = (char **)palloc0(sizeof(char *) * OUT_TUPLE_ARITY);
            values[0] = std::to_string(t.getQueryId()).data();
            values[1] = std::to_string(t.getTargetId()).data();
            values[2] = std::to_string(t.getRawScore()).data();
            values[3] = std::to_string(t.getBitScore()).data();
            values[4] = std::to_string(t.getEValue()).data();
            values[5] = std::to_string(t.getQStart()).data();
            values[6] = std::to_string(t.getQEnd()).data();
            values[7] = std::to_string(t.getQLen()).data();
            values[8] = std::to_string(t.getTStart()).data();
            values[9] = std::to_string(t.getTEnd()).data();
            values[10] = std::to_string(t.getTLen()).data();

            std::string qAln = t.getQAln();
            std::string tAln = t.getTAln();
            std::string cigar = t.getCigar();
            values[11] = qAln.data();
            values[12] = tAln.data();
            values[13] = cigar.data();

            values[14] = std::to_string(t.getAlnLen()).data();
            values[15] = std::to_string(t.getMismatch()).data();
            values[16] = std::to_string(t.getGapOpen()).data();
            values[17] = std::to_string(t.getPident()).data();

            const HeapTuple tuple = BuildTupleFromCStrings(attinmeta, values);
            tuplestore_puttuple(tupstore, tuple);
            heap_freetuple(tuple);
            pfree(values);
        }
    }

    void enable_materialize_mode(FunctionCallInfo &fcinfo,
                                 Tuplestorestate * &tupstore,
                                 AttInMetadata * &attinmeta)
    {
        // Get rsinfo and per_query_ctx
        ReturnSetInfo *rsinfo = (ReturnSetInfo *)fcinfo->resultinfo;
        MemoryContext per_query_ctx = rsinfo->econtext->ecxt_per_query_memory;
        MemoryContext oldcontext;

        // Tuple descriptor
        TupleDesc tupdesc;
        get_call_result_type(fcinfo, NULL, &tupdesc);

        // Get the AttInMetadata (out parameter)
        attinmeta = TupleDescGetAttInMetadata(tupdesc);

        // Switch to per_query memory context
        oldcontext = MemoryContextSwitchTo(per_query_ctx);

        // Initialize the output table as empty
        tupdesc = CreateTupleDescCopy(tupdesc);
        tupstore = tuplestore_begin_heap(
            rsinfo->allowedModes & SFRM_Materialize_Random,
            false,
            work_mem);

        // Switch back to the old memory context
        MemoryContextSwitchTo(oldcontext);

        // Set return mode to materialize and make the function return the tupstore
        rsinfo->returnMode = SFRM_Materialize;
        rsinfo->setResult = tupstore;
        rsinfo->setDesc = tupdesc;
    }

    Datum seq_search_mmseqs_one_to_one(FunctionCallInfo &fcinfo, char seqType, bool enableAmbiguity)
    {
        TEST_ARGS_FOR_NULL(10, true);
        DECLARE_VECTORS_AND_TARGET_NAMES();
        GET_OPTIONAL_PARAMS(2);
        ADD_QUERIES_ONE(0);
        ADD_TARGETS_ONE(1);
        SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION();
    }

    Datum seq_search_mmseqs_one_to_arr(FunctionCallInfo &fcinfo, char seqType, bool enableAmbiguity)
    {
        TEST_ARGS_FOR_NULL(10, true);
        DECLARE_VECTORS_AND_TARGET_NAMES();
        GET_OPTIONAL_PARAMS(2);
        ADD_QUERIES_ONE(0);
        ADD_TARGETS_ARR(1);
        SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION();
    }

    Datum seq_search_mmseqs_one_to_db(FunctionCallInfo &fcinfo, char seqType, bool enableAmbiguity)
    {
        TEST_ARGS_FOR_NULL(12, i != 3);
        DECLARE_VECTORS_AND_TARGET_NAMES();
        GET_OPTIONAL_PARAMS(4);
        ADD_QUERIES_ONE(0);
        ADD_TARGETS_DB(3, 1, 2);
        SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION();
    }

    Datum seq_search_mmseqs_arr_to_one(FunctionCallInfo &fcinfo, char seqType, bool enableAmbiguity)
    {
        TEST_ARGS_FOR_NULL(10, true);
        DECLARE_VECTORS_AND_TARGET_NAMES();
        GET_OPTIONAL_PARAMS(2);
        ADD_QUERIES_ARR(0);
        ADD_TARGETS_ONE(1);
        SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION();
    }

    Datum seq_search_mmseqs_arr_to_arr(FunctionCallInfo &fcinfo, char seqType, bool enableAmbiguity)
    {
        TEST_ARGS_FOR_NULL(10, true);
        DECLARE_VECTORS_AND_TARGET_NAMES();
        GET_OPTIONAL_PARAMS(2);
        ADD_QUERIES_ARR(0);
        ADD_TARGETS_ARR(1);
        SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION();
    }

    Datum seq_search_mmseqs_arr_to_db(FunctionCallInfo &fcinfo, char seqType, bool enableAmbiguity)
    {
        TEST_ARGS_FOR_NULL(12, i != 3);
        DECLARE_VECTORS_AND_TARGET_NAMES();
        GET_OPTIONAL_PARAMS(4);
        ADD_QUERIES_ARR(0);
        ADD_TARGETS_DB(3, 1, 2);
        SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION();
    }

    Datum seq_search_mmseqs_db_to_one(FunctionCallInfo &fcinfo, char seqType, bool enableAmbiguity)
    {
        TEST_ARGS_FOR_NULL(12, i != 3);
        DECLARE_VECTORS_AND_TARGET_NAMES();
        GET_OPTIONAL_PARAMS(4);
        ADD_QUERIES_DB(3, 0, 1);
        ADD_TARGETS_ONE(2);
        SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION();
    }

    Datum seq_search_mmseqs_db_to_arr(FunctionCallInfo &fcinfo, char seqType, bool enableAmbiguity)
    {
        TEST_ARGS_FOR_NULL(12, i != 3);
        DECLARE_VECTORS_AND_TARGET_NAMES();
        GET_OPTIONAL_PARAMS(4);
        ADD_QUERIES_DB(3, 0, 1);
        ADD_TARGETS_ARR(2);
        SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION();
    }

    Datum seq_search_mmseqs_db_to_db(FunctionCallInfo &fcinfo, char seqType, bool enableAmbiguity)
    {
        TEST_ARGS_FOR_NULL(14, i != 4 && i != 5);
        DECLARE_VECTORS_AND_TARGET_NAMES();
        GET_OPTIONAL_PARAMS(6);
        ADD_QUERIES_DB(4, 0, 1);
        ADD_TARGETS_DB(5, 2, 3);
        SET_MATERIALIZE_AND_CALL_MAIN_FUNCTION();
    }
}

extern "C"
{
    PG_MODULE_MAGIC;
    DEFINE_PG_FUNCTIONS_FOR_TYPE(nucl, 'n', false)
    DEFINE_PG_FUNCTIONS_FOR_TYPE(aa, 'a', false)
    DEFINE_PG_FUNCTIONS_FOR_TYPE(amb_aa, 'a', true)
}