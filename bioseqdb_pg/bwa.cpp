#include "bwa.h"

#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <htslib/htslib/sam.h>

//#define DEBUG_BWATOOLS 1

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

inline namespace {

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
    // Skip first N block, then use bam_cigar_type to check whether further blocks use nucleotides from the query to
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

void BioseqdbBWA::build_index(const std::vector<UnalignSequence>& v) {

    if (!v.size())
        return;

    // check the integrity of the input data
    for (std::vector<UnalignSequence>::const_iterator i = v.begin(); i != v.end(); ++i)
        if (i->name.empty() || i->nucleotides.empty())
            throw std::invalid_argument("BWAWrapper::constructIndex - Reference sequences must have non-empty name and seq");

    if (idx) {
        std::cerr << "...clearing old index" << std::endl;
        bwa_idx_destroy(idx);
        idx = 0;
    }

    // allocate memory for idx
    idx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));;

    // construct the forward-only pac
    uint8_t* fwd_pac = seqlib_make_pac(v, true); //true->for_only

    // construct the forward-reverse pac ("packed" 2 bit sequence)
    uint8_t* pac = seqlib_make_pac(v, false); // don't write, becasue only used to make BWT

    size_t tlen = 0;
    for (std::vector<UnalignSequence>::const_iterator i = v.begin(); i != v.end(); ++i)
        tlen += i->nucleotides.length();

#ifdef DEBUG_BWATOOLS
    std::cerr << "ref seq length: " << tlen << std::endl;
#endif

    // make the bwt
    bwt_t *bwt;
    bwt = seqlib_bwt_pac2bwt(pac, tlen*2); // *2 for fwd and rev
    bwt_bwtupdate_core(bwt);
    free(pac); // done with fwd-rev pac

    // construct sa from bwt and occ. adds it to bwt struct
    bwt_cal_sa(bwt, 32);
    bwt_gen_cnt_table(bwt);

    // make the bns
    bntseq_t * bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
    bns->l_pac = tlen;
    bns->n_seqs = v.size();
    bns->seed = 11;
    bns->n_holes = 0;

    // make the anns
    bns->anns = (bntann1_t*)calloc(v.size(), sizeof(bntann1_t));
    size_t offset = 0;
    for (size_t k = 0; k < v.size(); ++k) {
        seqlib_add_to_anns(v[k].name, v[k].nucleotides, &bns->anns[k], offset);
        offset += v[k].nucleotides.length();
    }

    //ambs is "holes", like N bases
    bns->ambs = 0; //(bntamb1_t*)calloc(1, sizeof(bntamb1_t));

    // make the in-memory idx struct
    idx->bwt = bwt;
    idx->bns = bns;
    idx->pac = fwd_pac;

}

std::vector<AlignMatch> BioseqdbBWA::align_sequence(std::string_view read_nucleotides) const {

    // we haven't made an index, just return
    if (!idx)
        throw std::runtime_error("bwa index not built");

    mem_alnreg_v ar;
    ar = mem_align1(memopt, idx->bwt, idx->bns, idx->pac, read_nucleotides.length(), read_nucleotides.data()); // get all the hits (was c_str())

#ifdef DEBUG_BWATOOLS
    std::cerr << "num hits: " << ar.n << " " << name << std::endl;
      size_t secondary_count_debug = 0;
      for (size_t i = 0; i < ar.n; ++i) {
	if (ar.a[i].secondary >= 0)
	  secondary_count_debug++;
      }
      std::cerr << "secondary count " << secondary_count_debug << std::endl;
    }
#endif

//    double primary_score = 0;
//
//    int secondary_count = 0;

    //setup a vector to store the hits. Why remake a vector? Because
    //I want them to be in the same order for all runs. With BWA multithreading,
    //it may not be so in the initial ar.a vector.
//    std::vector<mem_aln_t> a;
    std::vector<AlignMatch> matches;

    //size_t num_secondary = 0;
    // loop through the hits
    for (size_t i = 0; i < ar.n; ++i) {
        mem_aln_t a = mem_reg2aln(memopt, idx->bns, idx->pac, read_nucleotides.length(), read_nucleotides.data(), &ar.a[i]);
        matches.push_back({
            .ref_id_index = ar.a[i].rid,
            .ref_match_begin = ar.a[i].rb,
            .ref_match_end = ar.a[i].re,
            .query_subseq = std::string(slice_match(read_nucleotides, a.cigar, a.n_cigar)),
            .query_match_begin = ar.a[i].qb,
            .query_match_end = ar.a[i].qe,
            .is_primary = (a.flag & BAM_FSECONDARY) == 0,
            .is_secondary = (a.flag & BAM_FSECONDARY) != 0,
            .cigar = cigar_compressed_to_string(a.cigar, a.n_cigar),
            .score = a.score,
        });

//        if (ar.a[i].secondary >= 0 && (keep_sec_with_frac_of_primary_score < 0 || keep_sec_with_frac_of_primary_score > 1))
//            continue; // skip secondary alignments

        // get forward-strand position and CIGAR
//        mem_aln_t a_aln;

//        a_aln = mem_reg2aln(memopt, idx->bns, idx->pac, seq.length(), seq.c_str(), &ar.a[i]);

//        a.push_back(a_aln);
    }

    return matches;

    // sort it
//    std::sort(a.begin(), a.end(), aln_sort);
//
//    for (size_t i = 0; i < a.size(); ++i) {
//
//        // if score not sufficient or past cap, continue
//        //bool sec_and_low_score =  ar.a[i].secondary >= 0 && (primary_score * keep_sec_with_frac_of_primary_score) > a.score;
//        //bool sec_and_cap_hit = ar.a[i].secondary >= 0 && (int)i > max_secondary;
//        bool sec_and_low_score = (a[i].flag&BAM_FSECONDARY) && (primary_score * keep_sec_with_frac_of_primary_score) > a[i].score;
//        bool sec_and_cap_hit =   (a[i].flag&BAM_FSECONDARY) && (int)i > max_secondary;
//        if (sec_and_low_score || sec_and_cap_hit) {
//            free(a[i].cigar);
//            continue;
//        } else if (!(a[i].flag&BAM_FSECONDARY)) {
//            primary_score = a[i].score;
//            //num_secondary = 0;
//        }
//
//        // instantiate the read
//        BamRecord b;
//        b.init();
//
//        b.b->core.tid = a[i].rid;
//        b.b->core.pos = a[i].pos;
//        b.b->core.qual = a[i].mapq;
//        b.b->core.flag = a[i].flag;
//        b.b->core.n_cigar = a[i].n_cigar;
//
//        // set dumy mate
//        b.b->core.mtid = -1;
//        b.b->core.mpos = -1;
//        b.b->core.isize = 0;
//
//        // if alignment is reverse, set it
//        if (a[i].is_rev)
//            b.b->core.flag |= BAM_FREVERSE;
//
//        std::string new_seq = seq;
//        // if hardclip, figure out what to clip
//        if (hardclip) {
//            size_t tstart = 0;
//            size_t len = 0;
//            for (int i = 0; i < a[i].n_cigar; ++i) {
//                if (i == 0 && bam_cigar_op(a[i].cigar[i]) == BAM_CREF_SKIP) // first N (e.g. 20N50M)
//                    tstart = bam_cigar_oplen(a[i].cigar[i]);
//                else if (bam_cigar_type(bam_cigar_op(a[i].cigar[i]))&1) // consumes query, but not N
//                    len += bam_cigar_oplen(a[i].cigar[i]);
//            }
//            assert(len > 0);
//            assert(tstart + len <= seq.length());
//            new_seq = seq.substr(tstart, len);
//        }
//
//        // allocate all the data
//        b.b->core.l_qname = name.length() + 1;
//        b.b->core.l_qseq = new_seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
//        b.b->l_data = b.b->core.l_qname + (a[i].n_cigar<<2) + ((b.b->core.l_qseq+1)>>1) + (b.b->core.l_qseq);
//        b.b.get()->data = (uint8_t*)malloc(b.b.get()->l_data);
//
//        // allocate the qname
//        memcpy(b.b->data, name.c_str(), name.length() + 1);
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
//        // allocate the sequence
//        uint8_t* m_bases = b.b->data + b.b->core.l_qname + (b.b->core.n_cigar<<2);
//
//        // TODO move this out of bigger loop
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
//        // allocate the quality to NULL
//        uint8_t* s = bam_get_qual(b.b);
//        s[0] = 0xff;
//
//        b.AddIntTag("NA", ar.n); // number of matches
//        b.AddIntTag("NM", a[i].NM);
//
//        if (a[i].XA)
//            b.AddZTag("XA", std::string(a[i].XA));
//
//        // add num sub opt
//        //b.AddIntTag("SB", ar.a[i].sub_n);
//        b.AddIntTag("AS", a[i].score);
//
//        // count num secondaries
//        if (b.SecondaryFlag())
//            ++secondary_count;
//
//        vec.push_back(b);

#ifdef DEBUG_BWATOOLS
        // print alignment
      printf("\t%c\t%s\t%ld\t%d\t", "+-"[a[i].is_rev], idx->bns->anns[a[i].rid].name, (long)a[i].pos, a[i].mapq);
      for (int k = 0; k < a[i].n_cigar; ++k) // print CIGAR
      	printf("%d%c", a[i].cigar[k]>>4, "MIDSH"[a[i].cigar[k]&0xf]);
      printf("\t%d\n", a[i].NM); // print edit distance
#endif

//        free(a[i].cigar); // don't forget to deallocate CIGAR
//    }
//
//    free (ar.a); // dealloc the hit list
//
//    // add the secondary counts
//    for (BamRecordVector::iterator i = vec.begin(); i != vec.end(); ++i)
//        i->AddIntTag("SQ", secondary_count);

}

// modified from bwa (heng li)
uint8_t* BioseqdbBWA::seqlib_add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q)
{
    bntann1_t *p;
    int lasts;
    if (bns->n_seqs == *m_seqs) {
        *m_seqs <<= 1;
        bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
    }
    p = bns->anns + bns->n_seqs;
    p->name = strdup((char*)seq->name.s);
    p->anno = seq->comment.l > 0? strdup((char*)seq->comment.s) : strdup("(null)");
    p->gi = 0; p->len = seq->seq.l;
    p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
    p->n_ambs = 0;
    for (size_t i = lasts = 0; i < seq->seq.l; ++i) {
        int c = nst_nt4_table[(int)seq->seq.s[i]];
        if (c >= 4) { // N
            if (lasts == seq->seq.s[i]) { // contiguous N
                ++(*q)->len;
            } else {
                if (bns->n_holes == *m_holes) {
                    (*m_holes) <<= 1;
                    bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
                }
                *q = bns->ambs + bns->n_holes;
                (*q)->len = 1;
                (*q)->offset = p->offset + i;
                (*q)->amb = seq->seq.s[i];
                ++p->n_ambs;
                ++bns->n_holes;
            }
        }
        lasts = seq->seq.s[i];
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

// modified from bwa (heng li)
uint8_t* BioseqdbBWA::seqlib_make_pac(const std::vector<UnalignSequence>& v, bool for_only)
{

    bntseq_t * bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
    uint8_t *pac = 0;
    int32_t m_seqs, m_holes;
    int64_t m_pac, l;
    bntamb1_t *q;

    bns->seed = 11; // fixed seed for random generator
    m_seqs = m_holes = 8; m_pac = 0x10000;
    bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
    bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
    pac = (uint8_t*) calloc(m_pac/4, 1);
    q = bns->ambs;

    // move through the unaligned sequences
    for (size_t k = 0; k < v.size(); ++k) {

        // make the ref name kstring
        kstring_t * name = (kstring_t*)malloc(1 * sizeof(kstring_t));
        name->l = v[k].name.length() + 1;
        name->m = v[k].name.length() + 3;
        name->s = (char*)calloc(name->m, sizeof(char));
        memcpy(name->s, v[k].name.data(), v[k].name.length());
        name->s[v[k].name.length()] = '\0';

        // make the sequence kstring
        kstring_t * t = (kstring_t*)malloc(sizeof(kstring_t));
        t->l = v[k].nucleotides.length();
        t->m = v[k].nucleotides.length() + 2;
        //t->s = (char*)calloc(v[k].Seq.length(), sizeof(char));
        t->s = (char*)malloc(t->m);
        memcpy(t->s, v[k].nucleotides.data(), v[k].nucleotides.length());

        // put into a kstring
        kseq_t *ks = (kseq_t*)calloc(1, sizeof(kseq_t));
        ks->seq = *t;
        ks->name = *name;

        // make the forward only pac
        pac = seqlib_add1(ks, bns, pac, &m_pac, &m_seqs, &m_holes, &q);

        // clear it out
        free(name->s);
        free(name);
        free(t->s);
        free(t);
        //free(ks->name.s);
        //free(ks->seq.s);
        //free(ks->f->buf);
        //free(
        free(ks);
        // NOTE free kstring_t?
        //kseq_destroy(s);
    }

    if (!for_only)
    {
        // add the reverse complemented sequence
        m_pac = (bns->l_pac * 2 + 3) / 4 * 4;
        pac = (uint8_t*)realloc(pac, m_pac/4);
        memset(pac + (bns->l_pac+3)/4, 0, (m_pac - (bns->l_pac+3)/4*4) / 4);
        for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac)
            _set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
    }

    bns_destroy(bns);

    return pac;
}

// modified from bwa (heng li)
bwt_t *BioseqdbBWA::seqlib_bwt_pac2bwt(const uint8_t *pac, int bwt_seq_lenr)
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

// modified from bwa (heng li)
bntann1_t* BioseqdbBWA::seqlib_add_to_anns(std::string_view name, std::string_view seq, bntann1_t* ann, size_t offset)
{

    ann->offset = offset;
    ann->name = (char*)malloc(name.length()+1); // +1 for \0
    strncpy(ann->name, name.data(), name.length());
    ann->name[name.length()] = '\0';
    ann->anno = (char*)malloc(7);
    strcpy(ann->anno, "(null)\0");
    ann->len = seq.length();
    ann->n_ambs = 0; // number of "holes"
    ann->gi = 0; // gi?
    ann->is_alt = 0;

    return ann;
}
