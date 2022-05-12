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
    explicit BwaIndex(const std::vector<BwaSequence>& ref_seqs);
    ~BwaIndex();

    std::vector<BwaMatch> align_sequence(std::string_view query) const;

private:
    bwaidx_t* index;
    mem_opt_t* options;
};
