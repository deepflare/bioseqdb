#ifndef PSQL_INTERFACE_H
#define PSQL_INTERFACE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

    // For logging in Postgres.
    void log_from_cpp(const char *str);

    // For fetching targets (or queries if you'd like to).
    // Overwrites SPI_tuptable.
    const char *get_sequence(const char *table_name, const char *column_name, uint64_t id);

    // Fetches indexes for a given kmer into SPI_tuptable.
    // To access them from C++ use get_ith_index() but you have to do so
    // before calling get_indexes() or get_sequence() again.
    uint64_t get_indexes(const char *column_name, const char *kmer);

    // Fetches i-th index from SPI_tuptable (assuming SPI_tuptable contains indexes).
    void get_ith_index(int i, uint64_t *target_id, uint32_t *position);

#ifdef __cplusplus
}
#endif

#endif // PSQL_INTERFACE_H
