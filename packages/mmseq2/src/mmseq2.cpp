#include "mmseq2.h"

#include <utility>
#include <thread>
#include <cmath>

namespace {
    // from hsp to bit score, using K, lambda
    double evalBitScore(double rawScore) {
        static double K = 0.041, lambda = 0.267, logK = std::log(0.041), log2 = std::log(2);
        return (lambda * rawScore - logK) / log2;
    }

    // number of expected hits of similar score
    double evalEValue(double bitScore, uint32_t m, uint32_t n) {
        return (double)m * (double)n * (std::pow(2.0, -bitScore));
    }
}

void processQueries(const common::InputParams::InputParamsPtr &inputParams, mmseq2::GetterInterface getterInterface,
                    std::mutex *mtx, uint32_t *nextQuery, std::mutex *resMtx, const common::VecResPtr &resultPtr);

void processSingleQuery(const mmseq2::GetterInterfacePtr &getterInterfacePtr, uint64_t qId, common::InputParams::StrPtr queryStr,
                        common::InputParams::InputParamsPtr inputParams, std::mutex *resMtx, const common::VecResPtr &resultPtr);

common::VecRes mmseq2::MMSeq2(common::InputParams inputParams) {
    common::InputParams::InputParamsPtr inputParamsPtr = std::make_shared<common::InputParams>(inputParams);
    std::vector<std::thread> workers{};
    std::mutex mtx, resMtx;
    uint32_t nextQuery = 0;
    common::VecResPtr resultPtr = std::make_shared<common::VecRes>();

    bool allTargets = inputParamsPtr.get()->getAllTargets(), localTargets = inputParamsPtr.get()->getLocalTargets();
    mmseq2::GetterInterface getterInterface(allTargets, localTargets);
    getterInterface.getTargetsPtr() = (*inputParamsPtr).getTargetsPtr();

    if (localTargets) {
        uint32_t kMerLength = inputParamsPtr.get()->getKMerLength();

        auto targetsPtr = inputParamsPtr.get()->getTargetsPtr();
        for (uint32_t i = 0; i < targetsPtr.get()->size(); i++) {
            std::string target = *targetsPtr.get()->at(i);
            if (target.size() < kMerLength) {
                continue;
            }
            for (uint j = 0; j <= target.size() - kMerLength; j++) {
                std::string kMer = target.substr(j, kMerLength);
                (*getterInterface.getIndexesMapPtr())[kMer].push_back({i, j});
            }
        }
    }

    for (uint32_t i = 0; i < inputParamsPtr.get()->getThreadNumber(); ++i) {
        workers.emplace_back(std::thread(
                processQueries, inputParamsPtr, getterInterface, &mtx, &nextQuery, &resMtx, resultPtr));
    }

    for (std::thread &worker : workers) {
        worker.join();
    }

    if (inputParamsPtr.get()->getExecutionInterrupted()) {
        rpc::this_handler().respond_error(inputParamsPtr.get()->getExceptionRaisedMessage());
    }
    else {
        if (localTargets) {
            for (uint32_t i = 0; i < resultPtr.get()->size(); i++) {
                uint64_t localTargetId = (*resultPtr)[i].getTargetId();
                (*resultPtr)[i].setTargetId(inputParamsPtr.get()->getTIds().get()->at(localTargetId));
            }
        }

        return *resultPtr;
    }

    return common::VecRes{};
}

// Returns the smallest id of query that was not processed yet
uint32_t getNextQuery(uint32_t q_len, std::mutex *mtx, uint32_t *nextQuery) {
    uint32_t res = -1;

    mtx->lock();
    uint32_t tmpNextQuery = *nextQuery;

    if (tmpNextQuery < q_len) {
        res = tmpNextQuery;
        *nextQuery = tmpNextQuery + 1;
    }

    mtx->unlock();

    return res;
}

void processQueries(const common::InputParams::InputParamsPtr &inputParams, mmseq2::GetterInterface getterInterface,
                    std::mutex *mtx, uint32_t *nextQuery, std::mutex *resMtx, const common::VecResPtr &resultPtr) {

    mmseq2::GetterInterfacePtr getterInterfacePtr = std::make_shared<mmseq2::GetterInterface>(getterInterface);

    try {
        if (!getterInterfacePtr.get()->getLocalTargets()) {
            (*getterInterfacePtr).getDBconnPtr() = std::make_shared<DB::DBconn>(*inputParams->getTargetTableName(), *inputParams->getTargetColumnName(),
                                                    inputParams->getAllTargets(), inputParams->getTIds());
        }

        uint32_t tmpNextQuery;

        while ((tmpNextQuery = getNextQuery(inputParams.get()->getQLen(), mtx, nextQuery)) != -1) {
            processSingleQuery(getterInterfacePtr, inputParams.get()->getQIds()->at(tmpNextQuery), inputParams.get()->getQueries()->at(tmpNextQuery), inputParams, resMtx, resultPtr);
        }
    }
    catch(std::exception& e) {
        inputParams.get()->setInterruption(e);
    }

    try {
        if (!getterInterfacePtr.get()->getLocalTargets()) {
            getterInterfacePtr.get()->getDBconnPtr().get()->CloseConnection();
        }
    }
    catch(std::exception& e) {} // error handled by db
}

void processSingleQuery(const mmseq2::GetterInterfacePtr &getterInterfacePtr, uint64_t qId, common::InputParams::StrPtr queryStr,
                        common::InputParams::InputParamsPtr inputParams, std::mutex *resMtx, const common::VecResPtr &resultPtr) {

    mmseq2::Query query{qId, std::move(queryStr), std::move(inputParams)};

    query.findPrefilterKmerStageResults(getterInterfacePtr);

    query.executeAlignment(getterInterfacePtr, resMtx, resultPtr);
}

void mmseq2::Query::findPrefilterKmerStageResults(const mmseq2::GetterInterfacePtr &getterInterfacePtr) {

    if (this->sequence.get()->size() < this->kMerLength) {
        return;
    }

    int32_t SMaxSuf = 0;

    for (uint32_t i = 0; i < this->kMerLength - 1; ++i) {
        SMaxSuf += bioSequence->getPenaltyByChars(this->substitutionMatrixId, this->sequence.get()->at(i), this->sequence.get()->at(i));
    }

    for (uint32_t kMerPos = 0; kMerPos + this->kMerLength <= this->sequence.get()->length(); ++kMerPos) {
        SMaxSuf += bioSequence->getPenaltyByChars(this->substitutionMatrixId, this->sequence.get()->at(kMerPos + this->kMerLength - 1), this->sequence.get()->at(kMerPos + this->kMerLength - 1));

        std::string kMer = this->sequence.get()->substr(kMerPos, this->kMerLength);

        // prepare environment
        (*getterInterfacePtr).getSimKMersPtr() = std::make_shared<common::SimKMers>();
        (*getterInterfacePtr).getSimKMersPosPtr() = std::make_shared<std::vector<uint32_t>>();
        this->processSimilarKMers(getterInterfacePtr, kMerPos, kMer, SMaxSuf);

        // get hits from db
        AddHitsFromSimilarKmers(getterInterfacePtr);

        SMaxSuf -= bioSequence->getPenaltyByChars(this->substitutionMatrixId, this->sequence.get()->at(kMerPos), this->sequence.get()->at(kMerPos));
    }
}

void mmseq2::Query::processSimilarKMers(const mmseq2::GetterInterfacePtr &getterInterfacePtr, uint32_t kMerPos, std::string &kMer, int32_t SMaxSuf,
                                        int32_t Spref, uint32_t indx) {
    if (indx == 0 && evalBitScore(Spref + SMaxSuf) >= this->kMerGenThreshold) {
        (*getterInterfacePtr).addSimKmerPos(kMerPos);
        (*getterInterfacePtr).addSimKMer(kMer);
    }

    if (indx >= kMer.size()) {
        return;
    }

    const char currentAA = kMer[indx];

    const uint32_t currentAAId = bioSequence->charToId(currentAA);

    SMaxSuf -= bioSequence->getPenaltyByIds(this->substitutionMatrixId, currentAAId, currentAAId);

    for (uint32_t aaId = 0; aaId < bioSequence->getAlphabetSize(); ++aaId) {
        Spref += bioSequence->getPenaltyByIds(this->substitutionMatrixId, currentAAId, aaId);

        if (evalBitScore(Spref + SMaxSuf) >= this->kMerGenThreshold) {
            kMer[indx] = bioSequence->idToChar(aaId);

            if (currentAAId != aaId) {
                (*getterInterfacePtr).addSimKmerPos(kMerPos);
                (*getterInterfacePtr).addSimKMer(kMer);
            }

            processSimilarKMers(getterInterfacePtr, kMerPos, kMer, SMaxSuf, Spref, indx + 1);
        }

        Spref -= bioSequence->getPenaltyByIds(this->substitutionMatrixId, currentAAId, aaId);
    }

    kMer[indx] = currentAA;
}

void mmseq2::Query::AddHitsFromSimilarKmers(const mmseq2::GetterInterfacePtr &getterInterfacePtr) {

    // similar kMers in simKMers from getterInterface are in right order
    // and hits from simKMersHits <kmer, <tId, pos>> should be also checked in the same kMer order

    common::SimKMersHitsPtr simKMersHitsPtr = std::make_shared<common::SimKMersHits>();
    getterInterfacePtr.get()->getSimKMersHits(simKMersHitsPtr);

    uint32_t simHitsSize = simKMersHitsPtr.get()->size();
    if (simHitsSize == 0) {
        return;
    }

    // map : <kmer, vector<tId, pos>> helpful for check hits in right order
    std::unordered_map<std::string, std::vector<std::pair<uint64_t, uint32_t>>> hitsMap;
    for (const auto &hit : *simKMersHitsPtr) {
        hitsMap[hit.first].push_back(hit.second);
    }

    for (uint32_t i = 0; i < getterInterfacePtr.get()->getSimKMersPtr().get()->size(); i++) {
        const auto &kMer = getterInterfacePtr.get()->getSimKMersPtr().get()->at(i);
        auto it = hitsMap.find(kMer);

        if (it != hitsMap.end()) {
            auto hits = it->second;
            uint32_t kMerPos = getterInterfacePtr.get()->getSimKMersPosPtr().get()->at(i);

            // sort by pos
            std::sort(hits.begin(), hits.end(), [](auto &left, auto &right){
                return left.second < right.second;
            });

            for (auto &hit : hits) {
                uint64_t tId = hit.first;
                uint32_t pos = hit.second;
                int32_t diagonal = (int32_t)pos - (int32_t)kMerPos;

                if (diagonalPreVVisited[tId] && diagonalPrev[tId] == diagonal) {
                    addMatch(tId, diagonal);
                }
                diagonalPrev[tId] = diagonal;
                diagonalPreVVisited[tId] = true;
            }
        }
    }
}

double mmseq2::Query::ungappedAlignment(const StrPtr &querySequence, const StrPtr &targetSequence, int32_t diagonal) const {
    uint32_t querySequenceLength = querySequence.get()->size(), targetSequenceLength = targetSequence.get()->size();
    int32_t queryPosition = 0, queryLastPosition = (int32_t)querySequenceLength - 1;
    int32_t targetPosition = queryPosition + diagonal;

    if (targetPosition < 0) {
        queryPosition -= targetPosition;
        targetPosition = 0;
    }
    if (queryLastPosition + diagonal >= targetSequenceLength) {
        queryLastPosition = (int32_t)targetSequenceLength - diagonal - 1;
    }

    int32_t score = 0, maxScore = 0;

    while (queryPosition <= queryLastPosition) {
        score += bioSequence->getPenaltyByChars(getSubstitutionMatrixId(), querySequence.get()->at(queryPosition), targetSequence.get()->at(targetPosition));
        score = std::max(score, 0);
        maxScore = std::max(maxScore, score);

        queryPosition++;
        targetPosition++;
    }

    return evalBitScore(maxScore);
}

void mmseq2::Query::gappedAlignment(const StrPtr &querySequence, const StrPtr &targetSequence, common::MmseqResult &result) const {
    uint32_t qSeqLen = querySequence.get()->size(), tSeqLen = targetSequence.get()->size();
    int32_t costOp = this->gapOpenCost, costEx = this->costGapExtended;
    // E - gap in row, F - gap in column, H - best score
    std::vector<std::vector<int32_t> > E(qSeqLen, std::vector<int>(tSeqLen, -costOp));
    std::vector<std::vector<int32_t> > F(qSeqLen, std::vector<int>(tSeqLen, -costOp));
    std::vector<std::vector<int32_t> > H(qSeqLen, std::vector<int>(tSeqLen, 0));

    int32_t bestScore = 0;
    uint32_t qPos = -1, tPos = -1;

    // compute matrix
    for (uint32_t qInd = 0; qInd < qSeqLen; qInd++) {
        for (uint32_t tInd = 0; tInd < tSeqLen; tInd++) {

            if (tInd > 0) {
                E[qInd][tInd] = std::max(E[qInd][tInd - 1] - costEx, H[qInd][tInd - 1] - costOp);
                H[qInd][tInd] = std::max(H[qInd][tInd], E[qInd][tInd]);
            }
            if (qInd > 0) {
                F[qInd][tInd] = std::max(F[qInd - 1][tInd] - costEx, H[qInd - 1][tInd] - costOp);
                H[qInd][tInd] = std::max(H[qInd][tInd], F[qInd][tInd]);
            }

            int32_t match = bioSequence->getPenaltyByChars(this->substitutionMatrixId, querySequence.get()->at(qInd), targetSequence.get()->at(tInd));
            H[qInd][tInd] = std::max(H[qInd][tInd], match);
            if (qInd > 0 && tInd > 0) {
                H[qInd][tInd] = std::max(H[qInd][tInd], H[qInd - 1][tInd - 1] + match);
            }

            if (H[qInd][tInd] > bestScore) {
                bestScore = H[qInd][tInd];
                qPos = qInd;
                tPos = tInd;
            }
        }
    }

    result.setQEnd(qPos);
    result.setTEnd(tPos);
    std::string qAl, tAl, cigar;
    int lastAction = -1, identicalMatches = 0;
    auto qInd = (int32_t)qPos, tInd = (int32_t)tPos;

    // backtrace
    while (qInd >= 0 && tInd >= 0 && H[qInd][tInd] > 0) {
        int32_t match = bioSequence->getPenaltyByChars(this->substitutionMatrixId,
                                              querySequence.get()->at(qInd), targetSequence.get()->at(tInd));
        if ((qInd > 0 && tInd > 0 && H[qInd][tInd] == H[qInd - 1][tInd - 1] + match) || (H[qInd][tInd] == match)) {
            if (querySequence.get()->at(qInd) != targetSequence.get()->at(tInd)) {
                result.incrMismatch();
            } else {
                identicalMatches += 1;
            }
            cigar += 'M';
            qAl += querySequence.get()->at(qInd);
            tAl += targetSequence.get()->at(tInd);
            qInd--;
            tInd--;
            lastAction = 0;
        } else if (H[qInd][tInd] == E[qInd][tInd]) {
            result.incrGapOpen();
            while (tInd > 0 && E[qInd][tInd] == E[qInd][tInd - 1] - costEx) {
                cigar += 'D';
                qAl += ' ';
                tAl += targetSequence.get()->at(tInd);
                tInd--;
            }
            if (tInd > 0 && E[qInd][tInd] == H[qInd][tInd - 1] - costOp) {
                cigar += 'D';
                qAl += ' ';
                tAl += targetSequence.get()->at(tInd);
                tInd--;
            }
            lastAction = 1;
        } else if (H[qInd][tInd] == F[qInd][tInd]) {
            result.incrGapOpen();
            while (qInd > 0 && F[qInd][tInd] == F[qInd - 1][tInd] - costEx) {
                cigar += 'I';
                tAl += ' ';
                qAl += querySequence.get()->at(qInd);
                qInd--;
            }
            if (qInd > 0 && F[qInd][tInd] == H[qInd - 1][tInd] - costOp) {
                cigar += 'I';
                tAl += ' ';
                qAl += querySequence.get()->at(qInd);
                qInd--;
            }
            lastAction = 2;
        }
    }

    std::reverse(qAl.begin(), qAl.end());
    std::reverse(tAl.begin(), tAl.end());
    std::reverse(cigar.begin(), cigar.end());

    result.setQStart(qInd + (lastAction % 2 == 0 ? 1 : 0));
    result.setTStart(tInd + (lastAction < 2 ? 1 : 0));
    result.setQAln(qAl);
    result.setTAln(tAl);
    result.setAlnLen(qAl.size());

    result.setRawScore(bestScore);
    result.setBitScore(evalBitScore(result.getRawScore()));
    result.setEValue(evalEValue(result.getBitScore(), result.getQLen(), result.getTLen()));

    result.setPident((double)identicalMatches / (double)qAl.size());
    result.setCigar(cigar);
}

void mmseq2::Query::executeAlignment(const mmseq2::GetterInterfacePtr &getterInterfacePtr, std::mutex *resMtx, const common::VecResPtr &mmseqResult) {
    const PrefilterKmerStageResults &kmerStageResults = mmseq2::Query::getPrefilterKmerStageResults();
    for (uint32_t i = 0; i < kmerStageResults.getTargetsNumber(); i++) {
        const int32_t diagonal = kmerStageResults.getDiagonal((int)i);
        const uint64_t targetId = kmerStageResults.getTargetId((int)i);

        StrPtr querySequence = this->sequence;
        StrPtr targetSequence = getterInterfacePtr.get()->getTargetById(targetId);

        if (ungappedAlignment(querySequence, targetSequence, diagonal) >= this->ungappedAlignmentScore && filteredTargetIds.find(targetId) == filteredTargetIds.end()) {
            filteredTargetIds.insert(targetId);

            common::MmseqResult result(this->queryId, targetId);
            result.setQLen(querySequence->size());
            result.setTLen(targetSequence->size());
            gappedAlignment(querySequence, targetSequence, result);

            if (result.getEValue() <= this->evalTreshold) {
                resMtx->lock();
                mmseqResult->push_back(result);
                resMtx->unlock();
            }
        }
    }
}
