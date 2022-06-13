#include "mmseq2.h"
#include "rpc/server.h"

namespace {
    void runMMSeq2(std::vector<std::string> &&querySequences, std::vector<std::string> &&targetSequences,
                   std::string&& matrixName, uint32_t kMerLength,
                   int32_t kMerGenThreshold, int32_t ungappedAlignmentScore, double evalTreshold,
                   int32_t gapOpenCost, int32_t gapPenaltyCost, uint32_t threadNumber) {

        uint32_t qLen = querySequences.size();
        uint32_t tLen = targetSequences.size();

        common::InputParams::Vec64Ptr qIds = std::make_shared<std::vector<uint64_t>>(qLen, 0);
        for (uint32_t i = 0; i < qLen; i++) {
            (*qIds)[i] = (uint64_t)i;
        }
        common::InputParams::Vec64Ptr tIds = std::make_shared<std::vector<uint64_t>>(tLen, 0);
        for (uint32_t i = 0; i < tLen; i++) {
            (*tIds)[i] = (uint64_t)i;
        }
        common::InputParams::VecStrPtr queries = std::make_shared<std::vector<common::InputParams::StrPtr>>(qLen, nullptr);
        for (uint32_t i = 0; i < qLen; i++) {
            (*queries)[i] = std::make_shared<std::string>(querySequences[i]);
        }
        common::InputParams::VecStrPtr targets = std::make_shared<std::vector<common::InputParams::StrPtr>>(tLen, nullptr);
        for (uint32_t i = 0; i < tLen; i++) {
            (*targets)[i] = std::make_shared<std::string>(targetSequences[i]);
        }

        common::InputParams::StrPtr targetTableName = std::make_shared<std::string>("target");
        common::InputParams::StrPtr targetColumnName = std::make_shared<std::string>("seq");
        common::InputParams::StrPtr substitutionMatrixName = std::make_shared<std::string>(matrixName);

        bool allTargets = true; // sth for psql backend to know if there is sth query filtering
        bool localTargets = true; // local getters will be used

        common::InputParams inputParams(qLen, tLen, qIds, tIds, queries, allTargets, localTargets, targets,
                                        targetTableName, targetColumnName, substitutionMatrixName, kMerLength,
                                        kMerGenThreshold, ungappedAlignmentScore, evalTreshold,
                                        gapOpenCost, gapPenaltyCost, threadNumber, 'a', true);

        auto res = mmseq2::MMSeq2(inputParams);
        for (auto& el : res) {
            std::cout << "qId: " << el.getQueryId() << ", ";
            std::cout << "tId: " << el.getTargetId() << "\n";
            std::cout << "rawScore: " << el.getRawScore() << ", ";
            std::cout << "bitScore: " << el.getBitScore() << ", ";
            std::cout << "eValue: " << el.getEValue() << "\n";
            std::cout << "qStart: " << el.getQStart() << ", ";
            std::cout << "qEnd: " << el.getQEnd() << ", ";
            std::cout << "qLen: " << el.getQLen() << "\n";
            std::cout << "tStart: " << el.getTStart() << ", ";
            std::cout << "tEnd: " << el.getTEnd() << ", ";
            std::cout << "tLen: " << el.getTLen() << "\n";
            std::cout << "qAln: " << el.getQAln() << "\n";
            std::cout << "tAln: " << el.getTAln() << "\n";
            std::cout << "cigar: " << el.getCigar() << "\n";
            std::cout << "alnLen: " << el.getAlnLen() << "\n";
            std::cout << "mismatch: " << el.getMismatch() << "\n";
            std::cout << "gapOpen: " << el.getGapOpen() << "\n";
            std::cout << "pident: " << el.getPident() << "\n";
            std::cout << "\n";
        }
    }
}

int main()
{
   common::InputParams b;

   uint32_t port;
   std::cin >> port;

   rpc::server srv(port);

   srv.bind("mmseq2", [](common::InputParams inputParams)
            { return mmseq2::MMSeq2(inputParams); });

   srv.run();

//   runMMSeq2({"DDDDDDDDDCCGGGGGGGAA", "AAADDDDDDDCCGGGGGGGDD"},
//             {"DDDDDDDAAGGGGGGG"},
//             "blosum62", 7, 18, 0, 1, 4, 1, 2);
//
//    runMMSeq2({"MKFKLKTTSGTARRGEMTFSRPQGEFTVQTPAFMPVGTYGTVKGMTPEEVRATGAEILLGNTFHLWLRPGQEVMRKHGDLHDFMQWHRPILTDSGGFQVFSLGKLRKITEEGVKFQNPINGERIFLSPEKSMEIQYDLGSDIVMIFDECTPYPATFDYAKKSMEMSLRWAQRSRDRFDELGNKNALFGIVQGGTFEELRKVSAEGLVDIGFDGYAVGGLAVGEPKEEMHRILEFTTPLLPADKPRYLMGVGKPEDLVEGVRRGIDMFDCVMPTRNARNGHLFVTDGIVKIRNAKYRDDTSALDPHCDCYTCRHYTKSYLYHLDKCGEILGARLNTIHNLRYYQRLMEEIRQAIEEDRFDDFVVEFYTRMGKEVPPLQKP"},
//             {"MKFELDTTDGRARRGRLIFERGTVETPAFMPVGTYGTVKGMTPEEVRATGADILLGNTFHLWLRPGEEIMRKHGDLHDFMNWQRPILTDSGGFQVFSLGDIRKITEEGVHFRSPINGEKIFLDPEKSMQIQHALGSDVVMIFDECTPYPATEDEARKSMQMSLRWAKRSRDEFDRLENPNSLFGIIQGSVYEDLRDESLKGLVEIGFDGYAVGGLAVGEPKEDMHRILEHVCPQIPADKPRYLMGVGKPEDLVEGVRRGIDMFDCVMPTRNARNGHLFTSEGVIKIRNARHRDDTSPLDPKCDCYTCKNYSRAYLYHLDRCNEILGARLNTIHNLRYYQMLMEGLRGAIETGTLDAFVKDFYTSQGREVPELVD",
//              "MYKLLKKSGKARRGEFTTPHGVIQTPVFMNVGTLAAIKGAVSSMDLKEIGCQVELSNTYHLHLRPGDEVVKKMGGLHKFMNWDRPILTDSGGFQVFSLSKIRKIQEEGVYFNSHIDGRKIFMGPEESMRIQSNLASTIAMAFDECVENPAPREYVEKSVERTTRWLHRCKDEMNRLNSLPDTINNKQMLFGINQGGTYEDIRIEHAKTIAKMDLDGYAIGGLAVGESHEDMYRIIDAVVPHLPEDKPIYLMGVGIPSNILEAVDRGVDFFDCVLPARNGRHAHVFTKEGKINLLNAKFELDDRPIDEGCQCPACKHYTRSYIRHLFKAKEMLAMRLCVLHNLYFYNNLMEEIRDAIDGDYFKEYKERKLKEWGGRA",
//              "MYKLLKKSGKARRGEFTTPHGVIQTPVFMNVGTLAAIKGAVSSMDLKEIGCQVELSNTYHLHLRPGDEVVKKMGGLHKFMNWDRPILTDSGGFQVFSLSKIRKIQEEGVYFNSHIDGRKIFMGPEESMRIQSNLASTIAMAFDECVENPAPREYVEKSVERTTRWLHRCKDEMNRLNSLPDTINNKQMLFGINQGGTYEDIRIEHAKTIAKMDLDGYAIGGLAVGESHEDMYRIIDAVVPHLPEDKPIYLMGVGIPSNILEAVDRGVDFFDCVLPARNGRHAHVFTKEGKINLLNAKFELDDRPIDEGCQCPACKHYTRSYIRHLFKAKEMLAMRLCVLHNLYFYNNLMEEIRDAIDGDYFKEYKERKLKEWGGRA"},
//              "blosum62", 7, 18, 15, 0.001, 11, 1, 1);
//
//    runMMSeq2({"AAADDDDDDDCCGGGGGGGDD"},
//            {"DDDDDDDAAGGGGGGG", "DDDDDDDDDCCGGGGGGGAA"},
//            "blosum62", 7, 20, 0, 1, 4, 1, 1);
//
//    runMMSeq2({"AACCTTGG", "ACTGACTGACTG", "TACTCAT"},
//            {"TACGGTAGCTTACTGA", "CTAGCTTACGATGCAAG", "CTTACAGCATACAGCATCGAT"},
//            "blosum62", 5, 10, 10, 1, 4, 1, 5);

    return 0;
}
