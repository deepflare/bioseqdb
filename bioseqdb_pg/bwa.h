#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

extern "C" {
#include <bwa/bwt.h>
#include <bwa/bwamem.h>
}

#include "sequence.h"

struct BwaMatch {
    std::string_view ref_id;
    std::string ref_subseq;
    int64_t ref_match_begin;
    int64_t ref_match_end;
    int64_t ref_match_len;
    std::string_view query_subseq;
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
    explicit BwaIndex();
    ~BwaIndex();

    std::vector<BwaMatch> align_sequence(const PgNucleotideSequence& seq) const;
    void build();
    void add_ref_sequence(std::string_view id, const PgNucleotideSequence& seq);

    mem_opt_t* options;

private:
    std::vector<ubyte_t> pac_forward;
    std::vector<bntamb1_t> holes;
    std::vector<bntann1_t> annotations; 
    bwaidx_t* index;
};

