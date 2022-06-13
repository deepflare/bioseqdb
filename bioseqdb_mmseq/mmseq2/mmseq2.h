#ifndef BIOSEQDB_MMSEQ2_H
#define BIOSEQDB_MMSEQ2_H

#include <stdint.h>
#include <memory>
#include <exception>
#include <map>
#include <utility>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <mutex>
#include "../common/mmseq2lib.h"
#include "dbconn.h"
#include "rpc/this_handler.h"

namespace mmseq2
{
    class BioSequence
    {
    protected:
        uint32_t alphabetSize;
    public:
        class CharNotInAlphabet: public std::exception {
            virtual const char* what() const throw()
            {
                return "Invalid char: char does not belong to the chosen alphabet";
            }
        } charNotInAlphabet;

        class BioSequenceNotAssigned : public std::exception {
            virtual const char* what() const throw()
            {
                return "Invalid usage of BioSequence class";
            }
        } bioSequenceNotAssigned;

        virtual uint32_t charToId(char aa_id) const {
            throw bioSequenceNotAssigned;
        }

        virtual char idToChar(uint32_t aa_char) const {
            throw bioSequenceNotAssigned;
        }

        virtual int32_t getPenaltyByIds(uint32_t matrixId, uint32_t currentId, uint32_t replacementId) const {
            throw bioSequenceNotAssigned;
        }

        virtual int32_t getPenaltyByChars(uint32_t matrixId, char currentId, char replacementId) const {
            throw bioSequenceNotAssigned;
        }

        uint32_t getAlphabetSize() {
            return alphabetSize;
        }
    };

    class Nucleotide : public BioSequence
    {
    public:
        Nucleotide(bool isAmbiguityEnabled) {
            if (isAmbiguityEnabled) {
                alphabetSize = 5;
            }
            else {
                alphabetSize = 4;
            }
        }

        uint32_t charToId(char aa_id) const {
            for (uint32_t id = 0; id < alphabetSize; id++)
            {
                if (charId[id] == aa_id)
                {
                    return id;
                }
            }
            throw charNotInAlphabet;
        }

        char idToChar(uint32_t aa_char) const {
            if (aa_char >= alphabetSize)
            {
                throw charNotInAlphabet;
            }

            return charId[aa_char];
        }

        int32_t getPenaltyByIds(uint32_t matrixId, uint32_t currentId, uint32_t replacementId) const {
            return substitutionMatrix[matrixId][currentId][replacementId];
        }

        int32_t getPenaltyByChars(uint32_t matrixId, char currentId, char replacementId) const {
            return substitutionMatrix[matrixId][charToId(currentId)][charToId(replacementId)];
        }
    private:
        static constexpr char charId[5] = {'A', 'C', 'T', 'G', 'X'};

        static constexpr int32_t substitutionMatrix[1][5][5] = {
                {
                        {2, -3, -3, -3, -3},
                        {-3, 2, -3, -3, -3},
                        {-3, -3, 2, -3, -3},
                        {-3, -3, -3, 2, -3},
                        {-3, -3, -3, -3, -3}
                }
        };
    };

    class AminoAcid : public BioSequence
    {
    public:
        AminoAcid() = delete;

        AminoAcid(bool isAmbiguityEnabled)
        {
            if (isAmbiguityEnabled) {
                this->alphabetSize = 25;
            }
            else {
                this->alphabetSize = 20;
            }
        }

        uint32_t charToId(char aa_id) const
        {
            auto it = mapOfCharIds.find(aa_id);

            if (it == mapOfCharIds.end() || it->second >= alphabetSize) {
                throw charNotInAlphabet;
            }

            return it->second;
        }

        char idToChar(uint32_t aa_char) const
        {
            if (aa_char >= alphabetSize)
            {
                throw charNotInAlphabet;
            }
            return charId[aa_char];
        }

        int32_t getPenaltyByIds(uint32_t matrixId, uint32_t currentId, uint32_t replacementId) const
        {
            return substitutionMatrix[matrixId][currentId][replacementId];
        }

        int32_t getPenaltyByChars(uint32_t matrixId, char currentId, char replacementId) const 
        {
            return substitutionMatrix[matrixId][charToId(currentId)][charToId(replacementId)];
        }

    private:
        static constexpr char charId[25] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X', '*'};

        std::map<char, int> mapOfCharIds = {
                {'A', 0},
                {'R', 1},
                {'N', 2},
                {'D', 3},
                {'C', 4},
                {'Q', 5},
                {'E', 6},
                {'G', 7},
                {'H', 8},
                {'I', 9},
                {'L', 10},
                {'K', 11},
                {'M', 12},
                {'F', 13},
                {'P', 14},
                {'S', 15},
                {'T', 16},
                {'W', 17},
                {'Y', 18},
                {'V', 19},
                {'B', 20},
                {'J', 21},
                {'Z', 22},
                {'X', 23},
                {'*', 24}
        };

        static constexpr int32_t substitutionMatrix[5][25][25] =
        {
            {{5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -2, -2, 0, -1, -1, -1, -1, -5},
             {-2, 7, 0, -1, -3, 1, 0, -2, 0, -3, -2, 3, -1, -2, -2, -1, -1, -2, -1, -2, -1, -3, 1, -1, -5},
             {-1, 0, 6, 2, -2, 0, 0, 0, 1, -2, -3, 0, -2, -2, -2, 1, 0, -4, -2, -3, 5, -3, 0, -1, -5},
             {-2, -1, 2, 7, -3, 0, 2, -1, 0, -4, -3, 0, -3, -4, -1, 0, -1, -4, -2, -3, 6, -3, 1, -1, -5},
             {-1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -2, -3, -1, -5},
             {-1, 1, 0, 0, -3, 6, 2, -2, 1, -2, -2, 1, 0, -4, -1, 0, -1, -2, -1, -3, 0, -2, 4, -1, -5},
             {-1, 0, 0, 2, -3, 2, 6, -2, 0, -3, -2, 1, -2, -3, 0, 0, -1, -3, -2, -3, 1, -3, 5, -1, -5},
             {0, -2, 0, -1, -3, -2, -2, 7, -2, -4, -3, -2, -2, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1, -5},
             {-2, 0, 1, 0, -3, 1, 0, -2, 10, -3, -2, -1, 0, -2, -2, -1, -2, -3, 2, -3, 0, -2, 0, -1, -5},
             {-1, -3, -2, -4, -3, -2, -3, -4, -3, 5, 2, -3, 2, 0, -2, -2, -1, -2, 0, 3, -3, 4, -3, -1, -5},
             {-1, -2, -3, -3, -2, -2, -2, -3, -2, 2, 5, -3, 2, 1, -3, -3, -1, -2, 0, 1, -3, 4, -2, -1, -5},
             {-1, 3, 0, 0, -3, 1, 1, -2, -1, -3, -3, 5, -1, -3, -1, -1, -1, -2, -1, -2, 0, -3, 1, -1, -5},
             {-1, -1, -2, -3, -2, 0, -2, -2, 0, 2, 2, -1, 6, 0, -2, -2, -1, -2, 0, 1, -2, 2, -1, -1, -5},
             {-2, -2, -2, -4, -2, -4, -3, -3, -2, 0, 1, -3, 0, 8, -3, -2, -1, 1, 3, 0, -3, 1, -3, -1, -5},
             {-1, -2, -2, -1, -4, -1, 0, -2, -2, -2, -3, -1, -2, -3, 9, -1, -1, -3, -3, -3, -2, -3, -1, -1, -5},
             {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -3, -1, -2, -2, -1, 4, 2, -4, -2, -1, 0, -2, 0, -1, -5},
             {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, 2, 5, -3, -1, 0, 0, -1, -1, -1, -5},
             {-2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2, 1, -3, -4, -3, 15, 3, -3, -4, -2, -2, -1, -5},
             {-2, -1, -2, -2, -3, -1, -2, -3, 2, 0, 0, -1, 0, 3, -3, -2, -1, 3, 8, -1, -2, 0, -2, -1, -5},
             {0, -2, -3, -3, -1, -3, -3, -3, -3, 3, 1, -2, 1, 0, -3, -1, 0, -3, -1, 5, -3, 2, -3, -1, -5},
             {-1, -1, 5, 6, -2, 0, 1, -1, 0, -3, -3, 0, -2, -3, -2, 0, 0, -4, -2, -3, 5, -3, 1, -1, -5},
             {-1, -3, -3, -3, -2, -2, -3, -4, -2, 4, 4, -3, 2, 1, -3, -2, -1, -2, 0, 2, -3, 4, -2, -1, -5},
             {-1, 1, 0, 1, -3, 4, 5, -2, 0, -3, -2, 1, -1, -3, -1, 0, -1, -2, -2, -3, 1, -2, 5, -1, -5},
             {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5},
             {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1}},
            {{5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -2, -2, -1, -1, -5},
             {-2, 7, -1, -2, -4, 1, 0, -3, 0, -4, -3, 3, -2, -3, -3, -1, -1, -3, -1, -3, -1, -3, 0, -1, -5},
             {-1, -1, 7, 2, -2, 0, 0, 0, 1, -3, -4, 0, -2, -4, -2, 1, 0, -4, -2, -3, 5, -4, 0, -1, -5},
             {-2, -2, 2, 8, -4, 0, 2, -1, -1, -4, -4, -1, -4, -5, -1, 0, -1, -5, -3, -4, 6, -4, 1, -1, -5},
             {-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -2, -3, -1, -5},
             {-1, 1, 0, 0, -3, 7, 2, -2, 1, -3, -2, 2, 0, -4, -1, 0, -1, -1, -1, -3, 0, -3, 4, -1, -5},
             {-1, 0, 0, 2, -3, 2, 6, -3, 0, -4, -3, 1, -2, -3, -1, -1, -1, -3, -2, -3, 1, -3, 5, -1, -5},
             {0, -3, 0, -1, -3, -2, -3, 8, -2, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -4, -1, -4, -2, -1, -5},
             {-2, 0, 1, -1, -3, 1, 0, -2, 10, -4, -3, 0, -1, -1, -2, -1, -2, -3, 2, -4, 0, -3, 0, -1, -5},
             {-1, -4, -3, -4, -2, -3, -4, -4, -4, 5, 2, -3, 2, 0, -3, -3, -1, -3, -1, 4, -4, 4, -3, -1, -5},
             {-2, -3, -4, -4, -2, -2, -3, -4, -3, 2, 5, -3, 3, 1, -4, -3, -1, -2, -1, 1, -4, 4, -3, -1, -5},
             {-1, 3, 0, -1, -3, 2, 1, -2, 0, -3, -3, 6, -2, -4, -1, 0, -1, -3, -2, -3, 0, -3, 1, -1, -5},
             {-1, -2, -2, -4, -2, 0, -2, -3, -1, 2, 3, -2, 7, 0, -3, -2, -1, -1, 0, 1, -3, 2, -1, -1, -5},
             {-3, -3, -4, -5, -2, -4, -3, -4, -1, 0, 1, -4, 0, 8, -4, -3, -2, 1, 4, -1, -4, 1, -4, -1, -5},
             {-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -3, -1, -1, -5},
             {1, -1, 1, 0, -1, 0, -1, 0, -1, -3, -3, 0, -2, -3, -1, 5, 2, -4, -2, -2, 0, -3, 0, -1, -5},
             {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5, -3, -2, 0, 0, -1, -1, -1, -5},
             {-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1, 1, -4, -4, -3, 15, 2, -3, -5, -2, -2, -1, -5},
             {-2, -1, -2, -3, -3, -1, -2, -3, 2, -1, -1, -2, 0, 4, -3, -2, -2, 2, 8, -1, -3, -1, -2, -1, -5},
             {0, -3, -3, -4, -1, -3, -3, -4, -4, 4, 1, -3, 1, -1, -3, -2, 0, -3, -1, 5, -3, 2, -3, -1, -5},
             {-2, -1, 5, 6, -3, 0, 1, -1, 0, -4, -4, 0, -3, -4, -2, 0, 0, -5, -3, -3, 6, -4, 1, -1, -5},
             {-2, -3, -4, -4, -2, -3, -3, -4, -3, 4, 4, -3, 2, 1, -3, -3, -1, -2, -1, 2, -4, 4, -3, -1, -5},
             {-1, 0, 0, 1, -3, 4, 5, -2, 0, -3, -3, 1, -1, -4, -1, 0, -1, -2, -2, -3, 1, -3, 5, -1, -5},
             {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5},
             {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1}},
            {{4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1, -1, -4},
             {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2, 0, -1, -4},
             {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 4, -3, 0, -1, -4},
             {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, -3, 1, -1, -4},
             {0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4},
             {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, -2, 4, -1, -4},
             {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, -3, 4, -1, -4},
             {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1, -4},
             {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, -3, 0, -1, -4},
             {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, 3, -3, -1, -4},
             {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, 3, -3, -1, -4},
             {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, -3, 1, -1, -4},
             {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, 2, -1, -1, -4},
             {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, 0, -3, -1, -4},
             {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4},
             {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, -2, 0, -1, -4},
             {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, -1, -1, -4},
             {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -2, -2, -1, -4},
             {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -1, -2, -1, -4},
             {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, 2, -2, -1, -4},
             {-2, -1, 4, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, -3, 0, -1, -4},
             {-1, -2, -3, -3, -1, -2, -3, -4, -3, 3, 3, -3, 2, 0, -3, -2, -1, -2, -1, 2, -3, 3, -3, -1, -4},
             {-1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -2, -2, -2, 0, -3, 4, -1, -4},
             {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4},
             {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1}},
            {
                {5, -2, -2, -2, -1, -1, -1, 0, -2, -2, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -2, -2, -1, -1, -6},
                {-2, 6, -1, -2, -4, 1, -1, -3, 0, -3, -3, 2, -2, -4, -2, -1, -1, -4, -3, -3, -1, -3, 0, -1, -6},
                {-2, -1, 6, 1, -3, 0, -1, -1, 0, -4, -4, 0, -3, -4, -3, 0, 0, -4, -3, -4, 5, -4, 0, -1, -6},
                {-2, -2, 1, 6, -4, -1, 1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4, 5, -5, 1, -1, -6},
                {-1, -4, -3, -4, 9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -4, -2, -4, -1,
                 -6},
                {-1, 1, 0, -1, -4, 6, 2, -2, 1, -3, -3, 1, 0, -4, -2, 0, -1, -3, -2, -3, 0, -3, 4, -1, -6},
                {-1, -1, -1, 1, -5, 2, 6, -3, 0, -4, -4, 1, -2, -4, -2, 0, -1, -4, -3, -3, 1, -4, 5, -1, -6},
                {0, -3, -1, -2, -4, -2, -3, 6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -1, -5, -3, -1,
                 -6},
                {-2, 0, 0, -2, -4, 1, 0, -3, 8, -4, -3, -1, -2, -2, -3, -1, -2, -3, 2, -4, -1, -4, 0, -1, -6},
                {-2, -3, -4, -4, -2, -3, -4, -5, -4, 5, 1, -3, 1, -1, -4, -3, -1, -3, -2, 3, -4, 3, -4, -1, -6},
                {-2, -3, -4, -5, -2, -3, -4, -4, -3, 1, 4, -3, 2, 0, -3, -3, -2, -2, -2, 1, -4, 3, -3, -1, -6},
                {-1, 2, 0, -1, -4, 1, 1, -2, -1, -3, -3, 5, -2, -4, -1, -1, -1, -4, -3, -3, -1, -3, 1, -1, -6},
                {-1, -2, -3, -4, -2, 0, -2, -4, -2, 1, 2, -2, 6, 0, -3, -2, -1, -2, -2, 1, -3, 2, -1, -1, -6},
                {-3, -4, -4, -4, -3, -4, -4, -4, -2, -1, 0, -4, 0, 6, -4, -3, -2, 0, 3, -1, -4, 0, -4, -1, -6},
                {-1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4, 8, -1, -2, -5, -4, -3, -2, -4, -2, -1,
                 -6},
                {1, -1, 0, -1, -2, 0, 0, -1, -1, -3, -3, -1, -2, -3, -1, 5, 1, -4, -2, -2, 0, -3, 0, -1, -6},
                {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2, 1, 5, -4, -2, 0, -1, -1, -1, -1, -6},
                {-3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2, 0, -5, -4, -4, 11, 2, -3, -5, -3, -3, -1,
                 -6},
                {-2, -3, -3, -4, -3, -2, -3, -4, 2, -2, -2, -3, -2, 3, -4, -2, -2, 2, 7, -2, -3, -2, -3, -1,
                 -6},
                {0, -3, -4, -4, -1, -3, -3, -4, -4, 3, 1, -3, 1, -1, -3, -2, 0, -3, -2, 4, -4, 2, -3, -1, -6},
                {-2, -1, 5, 5, -4, 0, 1, -1, -1, -4, -4, -1, -3, -4, -2, 0, -1, -5, -3, -4, 5, -4, 0, -1, -6},
                {-2, -3, -4, -5, -2, -3, -4, -5, -4, 3, 3, -3, 2, 0, -4, -3, -1, -3, -2, 2, -4, 3, -3, -1, -6},
                {-1, 0, 0, 1, -4, 4, 5, -3, 0, -4, -3, 1, -1, -4, -2, 0, -1, -3, -3, -3, 0, -3, 5, -1, -6},
                {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                 -6},
                {-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,
                 1},
            },
            {{5, -2, -2, -3, -1, -1, -1, 0, -2, -2, -2, -1, -2, -3, -1, 1, 0, -4, -3, -1, -2, -2, -1, -1, -6},
             {-2, 6, -1, -3, -5, 1, -1, -3, 0, -4, -3, 2, -2, -4, -3, -1, -2, -4, -3, -3, -2, -3, 0, -1, -6},
             {-2, -1, 7, 1, -4, 0, -1, -1, 0, -4, -4, 0, -3, -4, -3, 0, 0, -5, -3, -4, 5, -4, -1, -1, -6},
             {-3, -3, 1, 7, -5, -1, 1, -2, -2, -5, -5, -1, -4, -5, -3, -1, -2, -6, -4, -5, 5, -5, 1, -1, -6},
             {-1, -5, -4, -5, 9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -4, -2, -4, -2, -5, -1, -6},
             {-1, 1, 0, -1, -4, 7, 2, -3, 1, -4, -3, 1, 0, -4, -2, -1, -1, -3, -3, -3, -1, -3, 5, -1, -6},
             {-1, -1, -1, 1, -6, 2, 6, -3, -1, -4, -4, 0, -3, -5, -2, -1, -1, -5, -4, -3, 1, -4, 5, -1, -6},
             {0, -3, -1, -2, -4, -3, -3, 6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4, -5, -5, -2, -5, -3, -1, -6},
             {-2, 0, 0, -2, -5, 1, -1, -3, 8, -4, -4, -1, -3, -2, -3, -2, -2, -3, 1, -4, -1, -4, 0, -1, -6},
             {-2, -4, -4, -5, -2, -4, -4, -5, -4, 5, 1, -4, 1, -1, -4, -3, -1, -4, -2, 3, -5, 3, -4, -1, -6},
             {-2, -3, -4, -5, -2, -3, -4, -5, -4, 1, 5, -3, 2, 0, -4, -3, -2, -3, -2, 0, -5, 4, -4, -1, -6},
             {-1, 2, 0, -1, -4, 1, 0, -2, -1, -4, -3, 6, -2, -4, -2, -1, -1, -5, -3, -3, -1, -3, 1, -1, -6},
             {-2, -2, -3, -4, -2, 0, -3, -4, -3, 1, 2, -2, 7, -1, -3, -2, -1, -2, -2, 0, -4, 2, -2, -1, -6},
             {-3, -4, -4, -5, -3, -4, -5, -5, -2, -1, 0, -4, -1, 7, -4, -3, -3, 0, 3, -2, -4, 0, -4, -1, -6},
             {-1, -3, -3, -3, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4, 8, -2, -2, -5, -4, -3, -3, -4, -2, -1, -6},
             {1, -1, 0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2, 5, 1, -4, -3, -2, 0, -3, -1, -1, -6},
             {0, -2, 0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -3, -2, 1, 6, -4, -2, -1, -1, -2, -1, -1, -6},
             {-4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2, 0, -5, -4, -4, 11, 2, -3, -6, -3, -4, -1, -6},
             {-3, -3, -3, -4, -4, -3, -4, -5, 1, -2, -2, -3, -2, 3, -4, -3, -2, 2, 8, -3, -4, -2, -3, -1, -6},
             {-1, -3, -4, -5, -2, -3, -3, -5, -4, 3, 0, -3, 0, -2, -3, -2, -1, -3, -3, 5, -4, 1, -3, -1, -6},
             {-2, -2, 5, 5, -4, -1, 1, -2, -1, -5, -5, -1, -4, -4, -3, 0, -1, -6, -4, -4, 5, -5, 0, -1, -6},
             {-2, -3, -4, -5, -2, -3, -4, -5, -4, 3, 4, -3, 2, 0, -4, -3, -2, -3, -2, 1, -5, 4, -4, -1, -6},
             {-1, 0, -1, 1, -5, 5, 5, -3, 0, -4, -4, 1, -2, -4, -2, -1, -1, -4, -3, -3, 0, -4, 5, -1, -6},
             {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6},
             {-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1}}
        };
    };

    class PrefilterKmerStageResults
    {
    private:
        uint64_t queryId;
        std::vector<uint64_t> targetIds;
        std::vector<int32_t> diagonals;

    public:
        explicit PrefilterKmerStageResults(uint64_t queryId) : queryId{queryId}, targetIds{std::vector<uint64_t>{}},
                                                               diagonals{std::vector<int32_t>{}} { }

        PrefilterKmerStageResults() = delete;

        void addDiagonal(uint64_t targetId, int32_t diagonal)
        {
            targetIds.push_back(targetId);
            diagonals.push_back(diagonal);
        }

        [[nodiscard]] int32_t getDiagonal(int index) const
        {
            return diagonals[index];
        }

        [[nodiscard]] uint64_t getTargetId(int index) const
        {
            return targetIds[index];
        }

        [[nodiscard]] uint32_t getTargetsNumber() const
        {
            return targetIds.size();
        }
    };

    common::VecRes MMSeq2(common::InputParams inputParams);

    class GetterInterface;
    using GetterInterfacePtr = std::shared_ptr<GetterInterface>;

    class Query
    {
    public:
        using StrPtr = common::InputParams::StrPtr;
        using VecStrPtr = common::InputParams::VecStrPtr;

        Query() = delete;

        Query(uint64_t queryId, StrPtr query, common::InputParams::InputParamsPtr inputParams) : queryId{queryId}, sequence{std::move(query)},
                                                                                                 kMerLength{inputParams.get()->getKMerLength()},
                                                                                                 substitutionMatrixId{inputParams.get()->getSubstitutionMatrixId()},
                                                                                                 kMerGenThreshold{inputParams.get()->getKMerGenThreshold()},
                                                                                                 ungappedAlignmentScore{inputParams.get()->getUngappedAlignmentScore()},
                                                                                                 targetColumnName{inputParams.get()->getTargetColumnName()},
                                                                                                 evalTreshold{inputParams.get()->getEvalThreshold()},
                                                                                                 gapOpenCost{inputParams.get()->getGapOpenCost()},
                                                                                                 costGapExtended{inputParams.get()->getGapPenaltyCost()},
                                                                                                 targetTableName{inputParams.get()->getTargetTableName()},
                                                                                                 prefilterKmerStageResults{queryId},
                                                                                                 diagonalPreVVisited{std::vector<bool>(inputParams.get()->getTLen(), false)},
                                                                                                 diagonalPrev{std::vector<int32_t>(inputParams.get()->getTLen(), 0)} {
            if (inputParams->getSequenceType() == 'n') {
                bioSequence = std::make_shared<Nucleotide>(Nucleotide(inputParams->getEnableAmbiguity()));
            }
            else if (inputParams->getSequenceType() == 'a') {
                bioSequence = std::make_shared<AminoAcid>(AminoAcid(inputParams->getEnableAmbiguity()));
            }
            else {
                throw std::invalid_argument("Invalid Sequence type. Sequence has to be either 'a' (Aminoacid) or 'n' (Nucleotide)");
            }
        }

        [[nodiscard]] const PrefilterKmerStageResults &getPrefilterKmerStageResults() const
        {
            return prefilterKmerStageResults;
        }

        void addMatch(uint64_t targetId, int32_t diagonal)
        {
            prefilterKmerStageResults.addDiagonal(targetId, diagonal);
        }

        void findPrefilterKmerStageResults(const mmseq2::GetterInterfacePtr &getterInterfacePtr);

        void executeAlignment(const mmseq2::GetterInterfacePtr &getterInterfacePtr, std::mutex *resMtx, const common::VecResPtr &mmseqResult);

        [[nodiscard]] uint32_t getSubstitutionMatrixId() const
        {
            return substitutionMatrixId;
        }

        [[nodiscard]] uint32_t getKMerLength() const
        {
            return kMerLength;
        }

    private:
        std::shared_ptr<BioSequence> bioSequence;

        uint64_t queryId;
        StrPtr sequence;
        StrPtr targetTableName;

        uint32_t kMerLength;
        uint32_t substitutionMatrixId;
        int32_t kMerGenThreshold;
        int32_t ungappedAlignmentScore;
        StrPtr targetColumnName;
        double evalTreshold;
        int32_t gapOpenCost;
        int32_t costGapExtended;

        PrefilterKmerStageResults prefilterKmerStageResults;

        std::vector<int32_t> diagonalPrev;

        std::vector<bool> diagonalPreVVisited;

        std::set<uint64_t> filteredTargetIds;

        void processSimilarKMers(const mmseq2::GetterInterfacePtr &getterInterfacePtr, uint32_t diagonalNumber, std::string &kMer, int32_t SMaxSuf,
                                 int32_t Spref = 0, uint32_t indx = 0);

        void AddHitsFromSimilarKmers(const mmseq2::GetterInterfacePtr &getterInterfacePtr);

        [[nodiscard]] double ungappedAlignment(const StrPtr &querySequence, const StrPtr &targetSequence, int32_t diagonal) const;

        void gappedAlignment(const StrPtr &querySequence, const StrPtr &targetSequence, common::MmseqResult &mmseqResult) const;
    };

    class GetterInterface
    {
    public:
        using IndexesMap = std::unordered_map<std::string, std::vector<std::pair<uint64_t, uint32_t>>>;
        using IndexesMapPtr = std::shared_ptr<IndexesMap>;
        using DBconnPtr = std::shared_ptr<DB::DBconn>;
        using Vec32Ptr = std::shared_ptr<std::vector<uint32_t>>;

        GetterInterface(bool allTs, bool localTs)
        {
            indexesMapPtr = std::make_shared<IndexesMap>();
            dbconnPtr = nullptr;
            allTargets = allTs;
            localTargets = localTs;
        }

        [[nodiscard]] DBconnPtr &getDBconnPtr()
        {
            return dbconnPtr;
        }

        [[nodiscard]] IndexesMapPtr &getIndexesMapPtr()
        {
            return indexesMapPtr;
        }

        [[nodiscard]] common::InputParams::VecStrPtr &getTargetsPtr()
        {
            return targetsPtr;
        }

        [[nodiscard]] bool getLocalTargets() const
        {
            return localTargets;
        }

        [[nodiscard]] bool getAllTargets() const
        {
            return allTargets;
        }

        [[nodiscard]] Vec32Ptr &getSimKMersPosPtr()
        {
            return simKMersPosPtr;
        }

        void addSimKmerPos(uint32_t kMerPos)
        {
            (*simKMersPosPtr).push_back(kMerPos);
        }

        [[nodiscard]] common::SimKMersPtr &getSimKMersPtr()
        {
            return simKMersPtr;
        }

        void addSimKMer(const std::string &kMer)
        {
            (*simKMersPtr).push_back(kMer);
        }

        void getSimKMersHits(common::SimKMersHitsPtr &simKMersHitsPtr)
        {
            if (localTargets)
            {
                for (const auto &kMer : *simKMersPtr)
                {
                    auto it = indexesMapPtr.get()->find(kMer);
                    if (it != indexesMapPtr.get()->end())
                    {
                        for (const auto hit : it->second)
                        {
                            (*simKMersHitsPtr).emplace_back(kMer, hit);
                        }
                    }
                }
            }
            else
            {
                dbconnPtr.get()->GetSimKMersHits(simKMersPtr, simKMersHitsPtr);
            }
        }

        Query::StrPtr getTargetById(uint64_t id) {
            if (localTargets)
            {
                if (targetsPtr.get()->size() <= id)
                {
                    throw std::invalid_argument("out of range in targetsSequences");
                }
                return targetsPtr.get()->at(id);
            }
            return dbconnPtr.get()->GetTargetById(id);
        }

    private:
        bool localTargets = false;
        bool allTargets = false;
        DBconnPtr dbconnPtr;
        IndexesMapPtr indexesMapPtr;
        common::InputParams::VecStrPtr targetsPtr;
        Vec32Ptr simKMersPosPtr;
        common::SimKMersPtr simKMersPtr;
    };
}

#endif // BIOSEQDB_MMSEQ2_H
