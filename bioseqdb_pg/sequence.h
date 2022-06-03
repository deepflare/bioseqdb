#pragma once

#include <cstdint>
#include <string_view> 

extern "C" {
#include <bwa/bwt.h>
#include <bwa/bwamem.h>
#include <postgres.h>
#include <catalog/pg_type.h>
}

static_assert(sizeof(bntamb1_t) == 16, "This should not happen");
static_assert(alignof(bntamb1_t) == 8, "This should not happen");

constexpr std::string_view allowed_nucleotides = "ACGTNWSMKRYBDHV";

struct RawPgNucleotideSequence;

struct PgNucleotideSequence {
    char* to_palloc_text() const;
    char* to_malloc_text() const;
    uint32_t occurences(char symbol) const;
    size_t length() const { return len; }
    RawPgNucleotideSequence* complement() const; 

    uint32_t holes_num;
    uint32_t len;
    uint32_t padded_len;
    bntamb1_t* holes;
    ubyte_t* pac;
};

struct RawPgNucleotideSequence {
    PgNucleotideSequence wrapped();

    char vl_len[4];
    uint32_t holes_num;
    uint32_t len;
    uint32_t padded_len;
    ubyte_t data[];
};

RawPgNucleotideSequence* nuclseq_from_text(std::string_view str);

static inline size_t pac_byte_size(size_t x) { return x / 4 + (x % 4 != 0 ? 1 : 0); }

static inline int32_t nuclcode_from_char(char chr) {
    return nst_nt4_table[chr];
}

static inline uint8_t pac_raw_get(const ubyte_t* pac, size_t index) {
    return pac[index >> 2] >> ((~index & 3) << 1) & 3;
}

static inline void pac_raw_set(ubyte_t* pac, size_t index, uint8_t value) {
    pac[index >> 2] |= value << ((~index & 3) << 1);
}


