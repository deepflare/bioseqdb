#include <algorithm>
#include <cstdint>
#include <random>

#include "sequence.h"

inline namespace {

char complement_symbol(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        case 'W': return 'W';
        case 'S': return 'S';
        case 'M': return 'K';
        case 'K': return 'M';
        case 'R': return 'Y';
        case 'Y': return 'R';
        case 'B': return 'V';
        case 'D': return 'H';
        case 'H': return 'D';
        case 'V': return 'B';
    }

    return 'N';
}

template<typename F>
void for_each_block(const NucletideSequence& nucls, F f) {
    uint32_t p = 0;
    for(uint32_t i = 0 ; i < nucls.holes_num ; i++) {
        const auto& hole = nucls.holes[i];

        if(hole.offset > p)
            f(p, hole.offset);
        p = (uint32_t) hole.offset + hole.len;
    }

    if (p < nucls.len)
        f(p, nucls.len);
}

uint32_t calculate_num_of_holes(std::string_view str) {
    uint32_t count = 0;
    char prev_chr = 0;

    for(auto chr : str) {
        if (prev_chr != chr && nuclcode_from_char(chr) >= 4)
            count += 1;
        prev_chr = chr;
    }

    return count;
}

RawNucleotideSequence* alloc_raw_nucls(uint32_t holes_num, uint32_t len) {
    // Postgresql requires logicaly same values to have same bits, so we use zero alloc to fill paddings of holes.
    const auto size = 4 * sizeof(uint32_t) + holes_num * sizeof(bntamb1_t) + pac_byte_size(len);
    const auto ptr = static_cast<RawNucleotideSequence*>(palloc0(size));

    SET_VARSIZE(ptr, size);
    ptr->holes_num = holes_num;
    ptr->len = len;
    ptr->padded_len = pac_byte_size(len) * 4;

    return ptr;
}

void inplace_to_text(const NucletideSequence& nucls, char* text) {
    for(uint32_t i = 0 ; i < nucls.len; i++)
        text[i] = "ACGT"[pac_raw_get(nucls.pac, i)];

    for(bntamb1_t* hole = nucls.holes ; hole < nucls.holes + nucls.holes_num ; hole++)
        std::fill(text + hole->offset, text + hole->offset + hole->len, hole->amb);

    text[nucls.len] = '\0';
}

}

uint32_t NucletideSequence::occurences(char chr) const {
    const ubyte_t code = nuclcode_from_char(chr);
    uint32_t count = 0;

    if (code >= 4) {
        for(uint32_t i = 0 ; i < holes_num ; i++) {
            if (holes[i].amb == chr)
                count += holes[i].len;
        }
    } else {
        for_each_block(*this, [&](uint32_t p, uint32_t q) {
            for(uint32_t i = p; i < q ; i++) {
                if(pac_raw_get(pac, i) == code)
                    count++;
            }
        });
    }

    return count;
};

RawNucleotideSequence* NucletideSequence::complement() const {
    auto cptr = alloc_raw_nucls(holes_num, len);
    auto cnucls = cptr->wrapped();

    std::copy_n(holes, holes_num, cnucls.holes);
    for(uint32_t i = 0 ; i < holes_num ; i++) {
        const auto& hole = holes[i];
        cnucls.holes[i].amb = complement_symbol(hole.amb);
        for(uint32_t j = hole.offset ; j < hole.offset + hole.len ; j++)
            pac_raw_set(cnucls.pac, j, pac_raw_get(pac, j));
    }

    for_each_block(*this, [&](uint32_t p, uint32_t q) {
        for(uint32_t i = p; i < q ; i++) {
            pac_raw_set(cnucls.pac, i, 0b11 - pac_raw_get(pac, i));
        }
    });

    return cptr;
};

RawNucleotideSequence* NucletideSequence::reverse() const {
    auto rev_ptr = alloc_raw_nucls(holes_num, len);
    auto rev_nucls = rev_ptr->wrapped();
    std::minstd_rand rng(holes_num ^ len);

    for_each_block(*this, [&](uint32_t p, uint32_t q) {
        for(uint32_t i = p; i < q ; i++) {
            pac_raw_set(rev_nucls.pac, len - i - 1, pac_raw_get(pac, i));
        }
    });

    for(uint32_t i = 0  ; i < holes_num ; i++) {
        const auto& hole = holes[i];
        auto& rev_hole = rev_nucls.holes[holes_num - i - 1];
        rev_hole = hole;
        rev_hole.offset = len - hole.offset - 1;

        for(uint32_t j = rev_hole.offset ; j < rev_hole.offset + rev_hole.len ; j++)
            pac_raw_set(rev_nucls.pac, j, rng() & 0b11);
    }

    for(uint32_t i = len; i < padded_len ; i++)
        pac_raw_set(rev_nucls.pac, i, rng() & 0b11);

    return rev_ptr;
}

char* NucletideSequence::to_text_palloc() const {
    auto text = reinterpret_cast<char*>(palloc(len + 1));
    inplace_to_text(*this, text);
    return text;
}

char* NucletideSequence::to_text_malloc() const {
    auto text = reinterpret_cast<char*>(malloc(len + 1));
    inplace_to_text(*this, text);
    return text;
}

RawNucleotideSequence* nuclseq_from_text(std::string_view str) {
    uint32_t holes_num = calculate_num_of_holes(str);
    RawNucleotideSequence* ptr = alloc_raw_nucls(holes_num, str.size());
    NucletideSequence nucls = ptr->wrapped();

    // libbwa requires random values inside holes, but again we want them to be deterministic => lcg
    std::minstd_rand rng(holes_num ^ str.size());
    bntamb1_t* hole = nucls.holes;
    uint32_t holes_cnt = 0;
    char prev_chr = 0;

    for(uint32_t idx = 0 ; idx < str.size() ; idx++) {
        char chr = str[idx];
        ubyte_t code = nuclcode_from_char(chr);

        if (code >= 4) {
            if (prev_chr == chr) {
                hole->len++;
            } else {
                hole = nucls.holes + holes_cnt++;
                hole->amb = chr;
                hole->offset = idx;
                hole->len = 1;
            }
            pac_raw_set(nucls.pac, idx, rng() & 0b11);
        }
        else {
            pac_raw_set(nucls.pac, idx, code);
        }

        prev_chr = chr;
    }

    for(uint32_t i = ptr->len ; i < ptr->padded_len ; i++)
        pac_raw_set(nucls.pac, i, rng() & 0b11);

    return ptr;
}

NucletideSequence RawNucleotideSequence::wrapped() {
    uint32_t holes_size = sizeof(bntamb1_t) * holes_num;

    return NucletideSequence {
        .holes_num = holes_num,
        .len = len,
        .padded_len = padded_len,
        .holes = reinterpret_cast<bntamb1_t*>(data),
        .pac = reinterpret_cast<ubyte_t*>(data + holes_size),
    };
}
