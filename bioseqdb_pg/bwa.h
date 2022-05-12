#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

extern "C" {
#include <bwa/bwamem.h>
}

struct BwaSequence {
    std::string_view id;
    std::string_view seq;
};

struct BwaMatch {
    std::string ref_id;
    int64_t ref_match_begin;
    int64_t ref_match_end;
    int64_t ref_match_len;
    std::string query_subseq;
    int query_match_begin;
    int query_match_end;
    int query_match_len;
    bool is_primary;
    bool is_secondary;
    bool is_reverse;
    std::string cigar;
    int score;
};

class BwaIndex {
public:
    BwaIndex();

    ~BwaIndex();

    void build_index(const std::vector<BwaSequence>& v);

    std::vector<BwaMatch> align_sequence(std::string_view read_nucleotides) const;

private:
    // Store the options in memory
    mem_opt_t * memopt;

    // hold the full index structure
    bwaidx_t* idx;
};
