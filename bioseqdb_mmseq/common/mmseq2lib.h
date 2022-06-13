#ifndef MMSEQ2_LIB_H
#define MMSEQ2_LIB_H

#include <memory>
#include <vector>
#include <string>
#include <stdint.h>
#include <iostream>
#include "rpc/server.h"

namespace common
{
    class InputParams
    {
    public:
        using Vec64Ptr = std::shared_ptr<std::vector<uint64_t>>;
        using StrPtr = std::shared_ptr<std::string>;
        using VecStrPtr = std::shared_ptr<std::vector<StrPtr>>;
        using InputParamsPtr = std::shared_ptr<InputParams>;
        using Str = std::string;

        InputParams() {}

        InputParams(uint32_t qLen, uint32_t tLen, common::InputParams::Vec64Ptr qIds,
                    common::InputParams::Vec64Ptr tIds, common::InputParams::VecStrPtr queries,
                    bool allTargets, bool localTargets, common::InputParams::VecStrPtr targets,
                    common::InputParams::StrPtr targetTableName,
                    common::InputParams::StrPtr targetColumnName,
                    const common::InputParams::StrPtr &substitutionMatrixName, uint32_t kMerLength,
                    int32_t kMerGenThreshold, int32_t ungappedAlignmentScore, double evalTreshold,
                    int32_t gapOpenCost, int32_t gapPenaltyCost, uint32_t threadNumber,
                    char sequenceType, bool enableAmbiguity) : qLen{qLen}, tLen{tLen}, qIds{std::move(qIds)}, tIds{std::move(tIds)},
                                                                                          queries{std::move(queries)}, allTargets{allTargets}, localTargets{localTargets}, targets{std::move(targets)},
                                                                                          targetTableName{std::move(targetTableName)},
                                                                                          targetColumnName{std::move(targetColumnName)},
                                                                                          substitutionMatrixName{substitutionMatrixName}, kMerLength{kMerLength},
                                                                                          kMerGenThreshold{kMerGenThreshold}, ungappedAlignmentScore{ungappedAlignmentScore},
                                                                                          evalTreshold{evalTreshold}, gapOpenCost{gapOpenCost}, gapPenaltyCost{gapPenaltyCost},
                                                                                          threadNumber{threadNumber}, sequenceType{sequenceType}, enableAmbiguity{enableAmbiguity}, executionInterrupted{false}
        {
            if (sequenceType == 'n') {
                substitutionMatrixId = 0;
                return;
            }

            if (substitutionMatrixName.get()->length() != 8 || substitutionMatrixName.get()->compare(0, 6, "blosum") != 0)
            {
                throw std::invalid_argument("Wrong substitution matrix name");
            }

            uint32_t blosumId = 10 * (substitutionMatrixName.get()->at(6) - '0') + (substitutionMatrixName.get()->at(7) - '0');

            switch (blosumId)
            {
            case 45:
                substitutionMatrixId = 0;
                break;
            case 50:
                substitutionMatrixId = 1;
                break;
            case 62:
                substitutionMatrixId = 2;
                break;
            case 80:
                substitutionMatrixId = 3;
                break;
            case 90:
                substitutionMatrixId = 4;
                break;
            default:
                throw std::invalid_argument("Wrong blosum number");
            }
        }

        [[nodiscard]] char getSequenceType() const {
            return sequenceType;
        }

        [[nodiscard]] bool getEnableAmbiguity() const {
            return enableAmbiguity;
        }

        [[nodiscard]] uint32_t getQLen() const
        {
            return qLen;
        }

        [[nodiscard]] uint32_t getTLen() const
        {
            return tLen;
        }

        Vec64Ptr getQIds()
        {
            return qIds;
        }

        [[nodiscard]] Vec64Ptr getTIds() const
        {
            return tIds;
        }

        [[nodiscard]] VecStrPtr getQueries() const
        {
            return queries;
        }

        [[nodiscard]] bool getAllTargets() const
        {
            return allTargets;
        }

        [[nodiscard]] bool getLocalTargets() const
        {
            return localTargets;
        }

        [[nodiscard]] VecStrPtr getTargetsPtr() const
        {
            return targets;
        }

        [[nodiscard]] StrPtr getTargetTableName() const
        {
            return targetTableName;
        }

        [[nodiscard]] StrPtr getTargetColumnName() const
        {
            return targetColumnName;
        }

        [[nodiscard]] uint32_t getKMerLength() const
        {
            return kMerLength;
        }

        [[nodiscard]] int32_t getKMerGenThreshold() const
        {
            return kMerGenThreshold;
        }

        [[nodiscard]] int32_t getUngappedAlignmentScore() const
        {
            return ungappedAlignmentScore;
        }

        [[nodiscard]] double getEvalThreshold() const
        {
            return evalTreshold;
        }

        [[nodiscard]] int32_t getGapOpenCost() const
        {
            return gapOpenCost;
        }

        [[nodiscard]] int32_t getGapPenaltyCost() const
        {
            return gapPenaltyCost;
        }

        [[nodiscard]] uint32_t getThreadNumber() const
        {
            return threadNumber;
        }

        [[nodiscard]] uint32_t getSubstitutionMatrixId() const
        {
            return substitutionMatrixId;
        }

        void setInterruption(std::exception& e) {
            executionInterrupted = true;
            exceptionRaisedMessage = e.what();
        }

        bool getExecutionInterrupted() const {
            return executionInterrupted;
        }

        std::string& getExceptionRaisedMessage() {
            return exceptionRaisedMessage;
        }

        MSGPACK_DEFINE_MAP(qLen, tLen, qIds, tIds,
                           queries, allTargets, localTargets, targets, targetTableName, targetColumnName, substitutionMatrixName,
                           kMerLength, kMerGenThreshold, ungappedAlignmentScore, evalTreshold, gapOpenCost,
                           gapPenaltyCost, threadNumber, substitutionMatrixId, sequenceType, enableAmbiguity,
                           executionInterrupted, exceptionRaisedMessage);

    private:
        uint32_t qLen;
        uint32_t tLen;
        Vec64Ptr qIds;
        Vec64Ptr tIds;

        VecStrPtr queries;
        bool allTargets;
        bool localTargets;
        VecStrPtr targets;

        StrPtr targetTableName;
        StrPtr targetColumnName;
        StrPtr substitutionMatrixName;

        uint32_t kMerLength;
        int32_t kMerGenThreshold;
        int32_t ungappedAlignmentScore;
        double evalTreshold;
        int32_t gapOpenCost;
        int32_t gapPenaltyCost;
        uint32_t threadNumber;
        uint32_t substitutionMatrixId;

        char sequenceType;
        bool enableAmbiguity;

        bool executionInterrupted;
        Str exceptionRaisedMessage;
    };

    class MmseqResult
    {
    private:
        uint64_t queryId, targetId;
        double rawScore = 0.0, bitScore = 0.0, eValue = 0.0;
        uint32_t qStart = 0, qEnd = 0, qLen = 0, tStart = 0, tEnd = 0, tLen = 0;
        std::string qAln, tAln, cigar;
        uint32_t alnLen = 0, mismatch = 0, gapOpen = 0;
        double pident = 0.0;

    public:
        MmseqResult() = default;
        MmseqResult(uint64_t queryId, uint64_t targetId) : queryId{queryId}, targetId{targetId} {}

        uint64_t getQueryId() const
        {
            return queryId;
        }

        uint64_t getTargetId() const
        {
            return targetId;
        }

        double getRawScore() const
        {
            return rawScore;
        }

        double getBitScore() const
        {
            return bitScore;
        }

        double getEValue() const
        {
            return eValue;
        }

        uint32_t getQStart() const
        {
            return qStart;
        }

        uint32_t getQEnd() const
        {
            return qEnd;
        }

        uint32_t getQLen() const
        {
            return qLen;
        }

        uint32_t getTStart() const
        {
            return tStart;
        }

        uint32_t getTEnd() const
        {
            return tEnd;
        }

        uint32_t getTLen() const
        {
            return tLen;
        }

        std::string getQAln() const
        {
            return qAln;
        }

        std::string getTAln() const
        {
            return tAln;
        }

        std::string getCigar() const
        {
            return cigar;
        }

        uint32_t getAlnLen() const
        {
            return alnLen;
        }

        uint32_t getMismatch() const
        {
            return mismatch;
        }

        uint32_t getGapOpen() const
        {
            return gapOpen;
        }

        double getPident() const
        {
            return pident;
        }

        void setTargetId(uint64_t targetId) {
            MmseqResult::targetId = targetId;
        }

        void setQLen(uint32_t qLen)
        {
            MmseqResult::qLen = qLen;
        }

        void setTLen(uint32_t tLen)
        {
            MmseqResult::tLen = tLen;
        }

        void setQEnd(uint32_t qEnd)
        {
            this->qEnd = qEnd;
        }

        void setTEnd(uint32_t tEnd)
        {
            this->tEnd = tEnd;
        }

        void incrMismatch()
        {
            this->mismatch++;
        }

        void incrGapOpen()
        {
            this->gapOpen++;
        }

        void setQStart(uint32_t qStart)
        {
            this->qStart = qStart;
        }

        void setRawScore(double rawScore)
        {
            MmseqResult::rawScore = rawScore;
        }

        void setBitScore(double bitScore)
        {
            MmseqResult::bitScore = bitScore;
        }

        void setEValue(double eValue)
        {
            MmseqResult::eValue = eValue;
        }

        void setTStart(uint32_t tStart)
        {
            MmseqResult::tStart = tStart;
        }

        void setQAln(const std::string &qAln)
        {
            MmseqResult::qAln = qAln;
        }

        void setTAln(const std::string &tAln)
        {
            MmseqResult::tAln = tAln;
        }

        void setCigar(const std::string &cigar)
        {
            MmseqResult::cigar = cigar;
        }

        void setAlnLen(uint32_t alnLen)
        {
            MmseqResult::alnLen = alnLen;
        }

        void setPident(double pident)
        {
            MmseqResult::pident = pident;
        }

        MSGPACK_DEFINE_MAP(queryId, targetId,
                           rawScore, bitScore, eValue, qStart, qEnd, qLen, tStart, tEnd, tLen,
                           qAln, tAln, cigar,
                           alnLen, mismatch, gapOpen,
                           pident);
    };

    using VecRes = std::vector<common::MmseqResult>;
    using VecResPtr = std::shared_ptr<common::VecRes>;

    using SimKMers = std::vector<std::string>; // {kmer1, kmer2, kmer3, ...}
    using SimKMersPtr = std::shared_ptr<SimKMers>;
    using SimKMersHits = std::vector<std::pair<std::string, std::pair<uint64_t, uint32_t>>>; // {<kmer1, <tId1, pos1>>, ...}
    using SimKMersHitsPtr = std::shared_ptr<SimKMersHits>;
}

#endif // MMSEQ2_LIB_H
