#include "bwa.h"

#include <algorithm>
#include <cstring>
#include <stdexcept>

#include <htslib/htslib/sam.h>
extern "C" {
// Internal libbwa symbol, not exported through any of the headers.
int is_bwt(ubyte_t *T, int n);
}

inline namespace {

class PacVector {
public:
    uint8_t get(uint64_t index) const {
        return pac[index >> 2] >> ((~index & 3) << 1) & 3;
    }

    void set(uint64_t index, uint8_t value) {
        pac[index >> 2] |= value << ((~index & 3) << 1);
    }

    void push_back(uint8_t value) {
        if (n % 4 == 0)
            pac.push_back(0);
        set(n++, value);
    }

    uint8_t* data() {
        return pac.data();
    }

    size_t size() const {
        return n;
    }

    size_t byte_size() const {
        return pac.size();
    }

private:
    std::vector<uint8_t> pac;
    size_t n = 0;
};

struct CompressedReference {
    PacVector pac_bwa;
    PacVector pac_forward;
    bntseq_t* bns;
};

char* make_c_string(std::string_view str_view) {
    char* c_str = (char*) malloc(str_view.length() + 1);
    std::copy(str_view.begin(), str_view.end(), c_str);
    c_str[str_view.length()] = '\0';
    return c_str;
}

void compress_reference_seq(const BwaSequence& ref, bntseq_t* bns, PacVector& pac, std::vector<bntamb1_t>& holes) {
    bntann1_t& ann = bns->anns[bns->n_seqs];
    ann.name = make_c_string(ref.id);
    ann.anno = make_c_string("(null)");
    ann.gi = 0;
    ann.len = ref.seq.length();
    ann.offset = bns->n_seqs == 0 ? 0 : bns->anns[bns->n_seqs - 1].offset + bns->anns[bns->n_seqs - 1].len;
    ann.n_ambs = 0;

    bntamb1_t* current_hole = nullptr;
    char prev_char = 0;
    for (char chr : ref.seq) {
        // ACGT decode to 0123 respectively, N and unknown characters decode to 4, and minus decodes to 5. This encoding
        // lets you efficiently compute the complement.
        int c = nst_nt4_table[(int) chr];

        if (c >= 4) {
            if (prev_char == chr)
                ++current_hole->len;
            else {
                holes.push_back({
                    .offset = current_hole->offset,
                    .len = 1,
                    .amb = chr,
                });
                current_hole = &holes.back();
                ++ann.n_ambs;
                ++bns->n_holes;
            }
        }

        prev_char = chr;

        if (c >= 4)
            c = lrand48() & 3;
        pac.push_back(c);
        ++bns->l_pac;
    }

    ++bns->n_seqs;
}

CompressedReference compress_reference(const std::vector<BwaSequence>& refs) {
    bntseq_t * bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
    // BWA encodes holes (unknown nucleotides) as a random valid nucleotide in its 2-bit encoding. This can produce
    // nondeterministic results, so the RNG seed is fixed to 11 to prevent this.
    bns->seed = 11;
    bns->anns = (bntann1_t*) calloc(refs.size(), sizeof(bntann1_t));

    PacVector pac;

    std::vector<bntamb1_t> holes;
    holes.reserve(8);

    for (const BwaSequence& ref : refs)
        compress_reference_seq(ref, bns, pac, holes);

    bns->ambs = (bntamb1_t*) calloc(holes.size(), sizeof(bntamb1_t));
    std::copy(holes.begin(), holes.end(), bns->ambs);

    PacVector pac_bwa = pac;
    for (int64_t l = bns->l_pac - 1; l >= 0; --l)
        pac_bwa.push_back(3 - pac_bwa.get(l));

    return {
        .pac_bwa = pac_bwa,
        .pac_forward = pac,
        .bns = bns,
    };
}

// modified from bwa (heng li)
bwt_t *seqlib_bwt_pac2bwt(const uint8_t *pac, int bwt_seq_lenr)
{

    bwt_t *bwt;
    ubyte_t *buf;
    int i;
    //FILE *fp;

    // initialization
    bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
    bwt->seq_len = bwt_seq_lenr; //bwa_seq_len(fn_pac); //dummy
    bwt->bwt_size = (bwt->seq_len + 15) >> 4;
    //fp = xopen(fn_pac, "rb");

    // prepare sequence
    //pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
    //buf2 = (ubyte_t*)calloc(pac_size, 1);
    //err_fread_noeof(buf2, 1, pac_size, fp);
    //err_fclose(fp);
    memset(bwt->L2, 0, 5 * 4);
    buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
    for (i = 0; i < (int)bwt->seq_len; ++i) {
        buf[i] = pac[i>>2] >> ((3 - (i&3)) << 1) & 3;
        ++bwt->L2[1+buf[i]];
    }
    for (i = 2; i <= 4; ++i)
        bwt->L2[i] += bwt->L2[i-1];
    //free(buf2);

    // Burrows-Wheeler Transform
    bwt->primary = is_bwt(buf, bwt->seq_len);
    bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4);
    for (i = 0; i < (int)bwt->seq_len; ++i)
        bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
    free(buf);
    return bwt;
}

std::string cigar_compressed_to_string(const uint32_t *raw, int len) {
    std::string cigar;
    for (int i = 0; i < len; ++i) {
        cigar += std::to_string(bam_cigar_oplen(raw[i]));
        cigar += bam_cigar_opchr(raw[i]);
    }
    return cigar;
}

}

BwaIndex::BwaIndex(const std::vector<BwaSequence>& ref_seqs): index(nullptr), options(mem_opt_init()) {
    assert(!ref_seqs.empty());
    for (auto&& ref_seq : ref_seqs)
        assert(!(ref_seq.id.empty() || ref_seq.seq.empty()));

    CompressedReference ref_compressed = compress_reference(ref_seqs);

    bwt_t *bwt = seqlib_bwt_pac2bwt(ref_compressed.pac_bwa.data(), ref_compressed.pac_bwa.size());
    bwt_bwtupdate_core(bwt);
    bwt_cal_sa(bwt, 32);
    bwt_gen_cnt_table(bwt);

    index = (bwaidx_t*) calloc(1, sizeof(bwaidx_t));
    index->bwt = bwt;
    index->bns = ref_compressed.bns;
    index->pac = (uint8_t*) malloc(ref_compressed.pac_forward.size());
    std::copy_n(ref_compressed.pac_forward.data(), ref_compressed.pac_forward.byte_size(), index->pac);
}

BwaIndex::~BwaIndex() {
    if (index)
        bwa_idx_destroy(index);
    if (options)
        free(options);
}

std::vector<BwaMatch> BwaIndex::align_sequence(std::string_view query) const {
    mem_alnreg_v ar = mem_align1(options, index->bwt, index->bns, index->pac, query.length(), query.data()); // get all the hits (was c_str())

    // TODO: Revert the CIGAR string when the match is reversed?

    std::vector<BwaMatch> matches;
    for (mem_alnreg_t* alignment = ar.a; alignment != ar.a + ar.n; ++alignment) {
        mem_aln_t a = mem_reg2aln(options, index->bns, index->pac, query.length(), query.data(), alignment);
        matches.push_back({
            .ref_id = std::string_view(index->bns->anns[alignment->rid].name),
            .ref_match_begin = alignment->rb,
            .ref_match_end = alignment->re,
            .ref_match_len = alignment->re - alignment->rb,
            .query_subseq = query.substr(alignment->qb, alignment->qe - alignment->qb),
            .query_match_begin = alignment->qb,
            .query_match_end = alignment->qe,
            .query_match_len = alignment->qe - alignment->qb,
            .is_primary = (a.flag & BAM_FSECONDARY) == 0,
            .is_secondary = (a.flag & BAM_FSECONDARY) != 0,
            .is_reverse = a.is_rev != 0,
            .cigar = cigar_compressed_to_string(a.cigar, a.n_cigar),
            .score = a.score,
        });
        free(a.cigar);
    }
    free(ar.a);
    return matches;
}
