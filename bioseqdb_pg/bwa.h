#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>

extern "C" {
#include <bwa/bwamem.h>
#include <bwa/kseq.h>
int is_bwt(ubyte_t *T, int n);
KSEQ_DECLARE(gzFile)
}

#define MEM_F_SOFTCLIP  0x200

struct UnalignSequence {
    std::string_view name;
    std::string_view nucleotides;
};

struct AlignMatch {
    int reference_id;
    bool is_secondary;
};

struct BioseqdbBWA {
    BioseqdbBWA() {
        idx = nullptr;
        memopt = mem_opt_init();
        memopt->flag |= MEM_F_SOFTCLIP;
    }

    ~BioseqdbBWA() {
        if (idx)
            bwa_idx_destroy(idx);
        if (memopt)
            free(memopt);
    }

    void build_index(const std::vector<UnalignSequence>& v);

    std::vector<AlignMatch> align_sequence(std::string_view read_nucleotides, bool hardclip) const;

private:

    // Store the options in memory
    mem_opt_t * memopt;

    // hold the full index structure
    bwaidx_t* idx;

    // overwrite the bwa bwt_pac2pwt function
    bwt_t *seqlib_bwt_pac2bwt(const uint8_t *pac, int bwt_seq_lenr);

    // add an anns (chr annotation structure)
    bntann1_t* seqlib_add_to_anns(std::string_view name, std::string_view seq, bntann1_t * ann, size_t offset);

    // overwrite the bwa-mem add1 function, which takes a sequence and adds to pac
    uint8_t* seqlib_add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q);

    // make the pac structure (2-bit encoded packed sequence)
    uint8_t* seqlib_make_pac(const std::vector<UnalignSequence>& v, bool for_only);
};
