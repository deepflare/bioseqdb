#include <fmt/format.h>
#include <iostream>
#include <bioseq_lib/bioseq_lib.h>
#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
using namespace SeqLib;

void bioseq_lib_hello_world() {
  std::cout << "Hello world! :D";
  //elog(WARNING, "PAPIEŻ ZAWADIAKA W POSTGRESIE MORDO");
}

void bioseq_lib_align(std::string querySeq = "CAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAG", std::string filepath="hg19.fasta") {
    RefGenome ref;
    ref.LoadIndex(filepath);

    // get sequence at given locus
    std::string seq = ref.QueryRegion("1", 1000000,1001000);

    // Make an in-memory BWA-MEM index of region
    BWAWrapper bwa;
    UnalignedSequenceVector usv = {{"chr_reg1", seq}};
    bwa.ConstructIndex(usv);

    // align an example string with BWA-MEM
    BamRecordVector results;
    // hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
    bwa.AlignSequence(querySeq, "my_seq", results, false, 0.9, 10); 

    // print results to stdout
    for (auto& i : results)
        std::cout << i << std::endl;
}
