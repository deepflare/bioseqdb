#include "bwa.h"

#include <cstring>
#include <stdexcept>

#include <htslib/htslib/sam.h>
extern "C" {
// Internal libbwa symbol, not exported through any of the headers.
int is_bwt(ubyte_t *T, int n);
}

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

inline namespace {

struct CompressedReference {
    std::vector<uint8_t> pac_bwa;
    uint8_t* pac_forward;
    bntseq_t* bns;
};

char* make_c_string(std::string_view str_view) {
    char* c_str = (char*) malloc(str_view.length() + 1);
    std::copy(str_view.begin(), str_view.end(), c_str);
    c_str[str_view.length()] = '\0';
    return c_str;
}

uint8_t* seqlib_add1(const BwaSequence& seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q)
{
    bntann1_t *p;
    int lasts;
    if (bns->n_seqs == *m_seqs) {
        *m_seqs <<= 1;
        bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
    }
    p = bns->anns + bns->n_seqs;
    p->name = make_c_string(seq.id);
    p->anno = make_c_string("(null)");
    p->gi = 0;
    p->len = seq.seq.length();
    p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
    p->n_ambs = 0;
    for (size_t i = lasts = 0; i < seq.seq.length(); ++i) {
        int c = nst_nt4_table[(int)seq.seq[i]];
        if (c >= 4) { // N
            if (lasts == seq.seq[i]) { // contiguous N
                ++(*q)->len;
            } else {
                if (bns->n_holes == *m_holes) {
                    (*m_holes) <<= 1;
                    bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
                }
                *q = bns->ambs + bns->n_holes;
                (*q)->len = 1;
                (*q)->offset = p->offset + i;
                (*q)->amb = seq.seq[i];
                ++p->n_ambs;
                ++bns->n_holes;
            }
        }
        lasts = seq.seq[i];
        { // fill buffer
            if (c >= 4) c = lrand48()&3;
            if (bns->l_pac == *m_pac) { // double the pac size
                *m_pac <<= 1;
                pac = (uint8_t*)realloc(pac, *m_pac/4);
                memset(pac + bns->l_pac/4, 0, (*m_pac - bns->l_pac)/4);
            }
            _set_pac(pac, bns->l_pac, c);
            ++bns->l_pac;
        }
    }
    ++bns->n_seqs;

    return pac;
}

CompressedReference compress_reference(const std::vector<BwaSequence>& ref) {
    bntseq_t * bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
    uint8_t *pac = 0;
    int32_t m_seqs, m_holes;
    int64_t m_pac;
    bntamb1_t *q;

    bns->seed = 11; // fixed seed for random generator
    m_seqs = m_holes = 8; m_pac = 0x10000;
    bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
    bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
    pac = (uint8_t*) calloc(m_pac/4, 1);
    q = bns->ambs;

    for (size_t k = 0; k < ref.size(); ++k)
        pac = seqlib_add1(ref[k], bns, pac, &m_pac, &m_seqs, &m_holes, &q);

    std::vector<uint8_t> pac_bwa(pac, pac + m_pac);
    int64_t m_pac_bwa = (bns->l_pac * 2 + 3) / 4 * 4;
    pac_bwa.resize(m_pac_bwa / 4);
    std::fill_n(pac_bwa.begin() + (bns->l_pac + 3) / 4, (m_pac_bwa - (bns->l_pac + 3) / 4 * 4) / 4, 0);
    uint64_t bns_l_pac_bwa = bns->l_pac;
    for (int64_t l = bns->l_pac - 1; l >= 0; --l, ++bns_l_pac_bwa)
        _set_pac(pac_bwa, bns_l_pac_bwa, 3 - _get_pac(pac_bwa, l));

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

std::string_view slice_match(std::string_view nucleotides, const uint32_t *cigar_raw, int cigar_len) {
    size_t start = 0;
    size_t len = 0;
    // Skip first N block, then use bam_cigar_type to check whether further blocks use seq from the query to
    // calculate the length of the match in the queried sequence.
    for (int i = 0; i < cigar_len; ++i) {
        if (i == 0 && bam_cigar_op(cigar_raw[i]) == BAM_CREF_SKIP)
            start = bam_cigar_oplen(cigar_raw[i]);
        else if (bam_cigar_type(bam_cigar_op(cigar_raw[i])) & 1)
            len += bam_cigar_oplen(cigar_raw[i]);
    }
    return nucleotides.substr(start, len);
}

}

BwaIndex::BwaIndex() {
    idx = nullptr;
    memopt = mem_opt_init();
    memopt->flag |= MEM_F_SOFTCLIP;
}

BwaIndex::~BwaIndex() {
    // TODO: Free idx->bns, hitting a double free right now.
    if (idx)
        bwa_idx_destroy(idx);
    if (memopt)
        free(memopt);
}

void BwaIndex::build_index(const std::vector<BwaSequence>& ref_seqs) {
    assert(ref_seqs.size() > 0);
    for (auto&& ref_seq : ref_seqs)
        assert(!(ref_seq.id.empty() || ref_seq.seq.empty()));
    assert(idx == nullptr);

    // allocate memory for idx
    idx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));;

    CompressedReference ref_compressed = compress_reference(ref_seqs);

    size_t tlen = 0;
    for (auto&& ref_seq : ref_seqs)
        tlen += ref_seq.seq.length();

#ifdef DEBUG_BWATOOLS
    std::cerr << "ref seq length: " << tlen << std::endl;
#endif

    // make the bwt
    bwt_t *bwt;
    bwt = seqlib_bwt_pac2bwt(ref_compressed.pac_bwa.data(), tlen*2); // *2 for fwd and rev
    bwt_bwtupdate_core(bwt);

    // construct sa from bwt and occ. adds it to bwt struct
    bwt_cal_sa(bwt, 32);
    bwt_gen_cnt_table(bwt);

    // make the in-memory idx struct
    idx->bwt = bwt;
    idx->bns = ref_compressed.bns;
    idx->pac = ref_compressed.pac_forward;

}

std::vector<BwaMatch> BwaIndex::align_sequence(std::string_view read_nucleotides) const {
    assert(idx != nullptr);

    mem_alnreg_v ar = mem_align1(memopt, idx->bwt, idx->bns, idx->pac, read_nucleotides.length(), read_nucleotides.data()); // get all the hits (was c_str())

    // TOOD: Free memory.

    std::vector<BwaMatch> matches;
    for (mem_alnreg_t* alignment = ar.a; alignment != ar.a + ar.n; ++alignment) {
        mem_aln_t a = mem_reg2aln(memopt, idx->bns, idx->pac, read_nucleotides.length(), read_nucleotides.data(), alignment);
        matches.push_back({
            .ref_id = std::string(idx->bns->anns[alignment->rid].name),
            .ref_match_begin = alignment->rb,
            .ref_match_end = alignment->re,
            .ref_match_len = alignment->re - alignment->rb,
            .query_subseq = std::string(slice_match(read_nucleotides, a.cigar, a.n_cigar)),
            .query_match_begin = alignment->qb,
            .query_match_end = alignment->qe,
            .query_match_len = alignment->qe - alignment->qb,
            .is_primary = (a.flag & BAM_FSECONDARY) == 0,
            .is_secondary = (a.flag & BAM_FSECONDARY) != 0,
            .is_reverse = a.is_rev != 0,
            .cigar = cigar_compressed_to_string(a.cigar, a.n_cigar),
            .score = a.score,
        });
    }
    return matches;

//        // if alignment is reverse, set it
//        if (a[i].is_rev)
//            b.b->core.flag |= BAM_FREVERSE;
//
//        // allocate the cigar. Reverse if aligned to neg strand, since mem_aln_t stores
//        // cigars relative to referemce string oreiatnion, not forward alignment
//        memcpy(b.b->data + b.b->core.l_qname, (uint8_t*)a[i].cigar, a[i].n_cigar<<2);
//
//        // convert N to S or H
//        int new_val = hardclip ? BAM_CHARD_CLIP : BAM_CSOFT_CLIP;
//        uint32_t * cigr = bam_get_cigar(b.b);
//        for (int k = 0; k < b.b->core.n_cigar; ++k) {
//            if ( (cigr[k] & BAM_CIGAR_MASK) == BAM_CREF_SKIP) {
//                cigr[k] &= ~BAM_CIGAR_MASK;
//                cigr[k] |= new_val;
//            }
//        }
//
//        int slen = new_seq.length();
//        int j = 0;
//        if (a[i].is_rev) {
//            for (int i = slen-1; i >= 0; --i) {
//
//                // bad idea but works for now
//                // this is REV COMP things
//                uint8_t base = 15;
//                if (new_seq.at(i) == 'T')
//                    base = 1;
//                else if (new_seq.at(i) == 'G')
//                    base = 2;
//                else if (new_seq.at(i) == 'C')
//                    base = 4;
//                else if (new_seq.at(i) == 'A')
//                    base = 8;
//
//                m_bases[j >> 1] &= ~(0xF << ((~j & 1) << 2));   ///< zero out previous 4-bit base encoding
//                m_bases[j >> 1] |= base << ((~j & 1) << 2);  ///< insert new 4-bit base encoding
//                ++j;
//            }
//        } else {
//            for (int i = 0; i < slen; ++i) {
//                // bad idea but works for now
//                uint8_t base = 15;
//                if (new_seq.at(i) == 'A')
//                    base = 1;
//                else if (new_seq.at(i) == 'C')
//                    base = 2;
//                else if (new_seq.at(i) == 'G')
//                    base = 4;
//                else if (new_seq.at(i) == 'T')
//                    base = 8;
//
//                m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
//                m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding
//
//            }
//        }
//
}
