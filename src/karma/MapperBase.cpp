/*
 * Copyright (c) 2009 Regents of the University of Michigan
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "GenomeSequence.h"
#include "ColorSpace.h"
#include "MappingStats.h"
#include "MatchedReadBase.h"
#include "MapperBase.h"
#include "MapperSEColorSpace.h"
#include "ReadsProcessor.h"
#include "MathConstant.h"
#include "Performance.h"
#include "SmithWaterman.h"
#include "Error.h"
#include "Util.h"

#include "../bam/SamFlag.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

unsigned int MapperBase::colorSpaceSNP[][3]={{ 5,10,15}, { 4,11,14}, { 7, 8,13}, { 6, 9,12},
                                             { 1,11,14}, { 0,10,15}, { 3, 9,12}, { 2, 8,13},
                                             { 2, 7,13}, { 3, 6,12}, { 0, 5,15}, { 1, 4,14},
                                             { 3, 6, 9}, { 2, 7, 8}, { 1, 4,11}, { 0, 5,10}
};

#include "debug.h"

/**
 * If trimming happened, restore to original read and quality and adjust genomeMatchingPosition
 * @param cigarRoller   store cigarRoller results
 * @param indexer       change read and quality in the indexer determined by bestMatch
 * @param genomeMatchPosition   store adjusted genomeMatchPosition
 */
void MapperBase::restoreTrimming(CigarRoller&   cigarRoller,
                                 ReadIndexer&   indexer,
                                 genomeIndex_t& genomeMatchPosition)
{
    // if trimming never occured, just return
    if (!(indexer.leftTrimmedBases.size() || indexer.rightTrimmedBases.size()))
        return;

    // rewrite the cigar string to include the two softclips
    CigarRoller newCigar;

    // recover the untrimmed read and quality string:
    std::string originalRead;
    std::string originalQualities;

    if (indexer.isForward)
    {
        originalRead = indexer.leftTrimmedBases + indexer.read + indexer.rightTrimmedBases;
        originalQualities = indexer.leftTrimmedQualities + indexer.phredQuality + indexer.rightTrimmedQualities;
        newCigar.Add(CigarRoller::softClip, indexer.leftTrimmedBases.size());
        newCigar.Add(cigarRoller);
        newCigar.Add(CigarRoller::softClip, indexer.rightTrimmedBases.size());
        // conditionally adjust the match postion for the soft clip on the left
        if (genomeMatchPosition != INVALID_GENOME_INDEX)
            genomeMatchPosition -= indexer.leftTrimmedBases.size();

    }
    else
    {
        originalRead = indexer.rightTrimmedBases + indexer.read + indexer.leftTrimmedBases;
        originalQualities = indexer.rightTrimmedQualities + indexer.phredQuality + indexer.leftTrimmedQualities;
        newCigar.Add(CigarRoller::softClip, indexer.rightTrimmedBases.size());
        newCigar.Add(cigarRoller);
        newCigar.Add(CigarRoller::softClip, indexer.leftTrimmedBases.size());
        // conditionally adjust the match postion for the soft clip on the left
        if (genomeMatchPosition != INVALID_GENOME_INDEX)
            genomeMatchPosition -= indexer.rightTrimmedBases.size();
    }

    cigarRoller = newCigar;
    indexer.read = originalRead;
    indexer.phredQuality = originalQualities;
}



//
// get the CIGAR string for this match position.
//
// Returns a match position offset that indicates how to adjust the
// given match position.  This allows for variations due to indels
// that may occur ahead of where the index word position is.
//
// When the match being printed (usually MapperBase::bestMatch) was
// a gapped match (i.e. Smith Waterman was run to find it), then we
// have a messy problem.
//
// The difficulty is that to generate the CIGAR string, we have to have
// the array H, which gets discarded after computing getSumQSW().  The
// reason is we keep possibly thousands of copies of the MatchedReadBase
// in paired end code, and the array is large and is statically allocated,
// so we can't expect it to be kept around very long.
//
// In order to obtain the CIGAR string, we have to rerun the
// Smith Waterman algorithm so we can do the walkback from the lower
// right portion of the array H.
//
//
void  MapperBase::populateCigarRollerAndGenomeMatchPosition()
{
    // getBestMatch is virtual, remember both MapperSE and MapperPE will call this funciton
    MatchedReadBase &bestMatch = getBestMatch();    

    if (bestMatch.qualityIsValid() && bestMatch.gappedAlignment)
    {
        int matchPositionOffset = 0;

        //
        // if gapped, and if a valid quality, go ahead and get the
        // CIGAR roller data structure.  If quality is invalid, the code
        // below is not safe to call.
        //
        // This code fixes the issue that the index alignment will
        // fix a read at a certain location, but in the event of a
        // structural variation that occurs before the index word used
        // to assign gemomeMatchPosition for the read, we need to
        // adjust the real read location by the amount of that
        // structural variation (zero or more, actually).
        //
        // In the calls below, bestMatch.whichWord tells us which
        // index word was used to assign genomeMatchPosition, and
        // therefore where we need to adjust to.
        //
        // For a visualization of this, imagine this read:
        //
        // ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
        //        ^       ^
        //        |       this is an example location where the index word starts (base #15)
        //        +- we'll imagine there is a single base deletion here
        // On fast index alignment, the index code places the start
        // of the read one single base earlier than it should be.
        //
        // Only by examining the cigar string can we know what the
        // adjustment is to be made.
        //
        // It is true that the Smith Waterman algorithm could incorporate
        // an offset argument that it would update according to its detection
        // of indels, but this adds even more complexity to an already complex
        // algorithm.
        //
        if (bestMatch.isForward())
        {
            //
            // To get the cigar roller, the following call is
            // essentially re-running Smith Waterman again, exactly
            // as it was done during the fast index lookup phase.
            //
            forward.getCigarRoller(bestMatch.genomeMatchPosition,
                                   bestMatch.whichWord,
                                   matchPositionOffset,
                                   cigarRoller);
        }
        else
        {
            backward.getCigarRoller(bestMatch.genomeMatchPosition,
                                    bestMatch.whichWord,
                                    matchPositionOffset,
                                    cigarRoller);
        }
        bestMatch.genomeMatchPosition += matchPositionOffset;
    }
    else if (!localGappedAlignment)
    {
        //
        // ungapped, or no match case, so don't update match position:
        //
        // Set a simple CIGAR string, basically "%dM" where %d
        // is the read length.
        //
        cigarRoller.clear();
        cigarRoller.Add(CigarRoller::match, forward.read.size());
    }
    else
    {
        // local gapped alignment yields a correct CIGAR string
        // and match position, so nothing to do here.
    }

    //
    // now unwind quality clipping changes that we made in ReadInder::setReadAndQuality
    //
    // re-assemble the indexer read (sequence) and quality strings and
    // also adjust genomeMatchPostion and the cigarRoller...
    //

    if (!bestMatch.qualityIsValid())
    {
        // XXX this is a hack to fix an earlier error
        // to fix it, put an assert() here and find out
        // what is causing the bug.
        bestMatch.indexer = &forward;
    }

    ReadIndexer &indexer = bestMatch.isForward() ? forward : backward;
    restoreTrimming(cigarRoller, indexer, bestMatch.genomeMatchPosition);
}


MapperBase::MapperBase() :  forward(mapperOptions), backward(mapperOptions)
{
    for (int i = 0; i < MAX_Q2P; i++)
        sumQualityToProb[i] = pow(10.0, (-0.1 * i));
    rand = &globalRandom;
    gs = NULL;
    wordIndex = NULL;
    line = 0;
    numberReads = 0;

    // lame:
    forward.isForward = true;
    backward.isForward = false;

    localGappedAlignment = false;
    //isProperAligned = false;

    samMateFlag = 0;        // for paired end, mark which is mate 1 versus mate 2
}

MapperBase::~MapperBase()
{
    // we don't own gs or wordIndex, so we don't close them here.
}

std::string MapperBase::Integer2Word(wordInteger_t n, unsigned int wordsize)
{
    assert(wordsize == wordIndex->wordSize);

    std::string word("NNNNN00000NNNNN");

    if (n > THIRTY_BIT_MASK)
        error("Invalid integer for a word of size %u: %d\n", wordsize, n);

    for (unsigned int i = 0; i < wordsize; i++)
    {
        word[wordsize - 1 - i] = MapIntegerToBase(n & 3);
        n >>= 2;
    }

    return word;
}

int MapperBase::processReadAndQuality(Fastq& fq)
{
    std::string tag = fq.tag.c_str();
    std::string read = fq.read.c_str();
    std::string qual = fq.qual.c_str();
    return processReadAndQuality(tag, read, qual);
}


// 
// Given a fastq record, we will trim it from left and/or right, then call setReadAndQuality() function;
// store original data and quality for color space reads
// @param fragmentTag, readFragment, dataQuality: a Fastq read.
// @return int 0: if setReadAndQuality() succeed; or 1 if failed
//
int MapperBase::processReadAndQuality(std::string& fragmentTag, std::string& readFragment, std::string& dataQuality)
{
#if 0
    // debug_by_tag
    // for debug a certain read
    String tag = fragmentTag.c_str();
    String toFind = "1028:13165";
    if (tag.Find(toFind) >= 0)
    {
        printf("1");
    }
#endif

    this->fragmentTag=fragmentTag;
    //
    // left and right truncate if requested
    //
    // It is possible that after trimming, the read
    // will be too short to map, but it is the job of
    // ::setReadAndQuality to determine that.
    //
    if (mapperOptions.trimLeft)
    {
        readFragment.erase(0, mapperOptions.trimLeft);
        dataQuality.erase(0, mapperOptions.trimLeft);
    }
    if (mapperOptions.trimRight)
    {
        throw std::logic_error("trim right option not tested");
        readFragment.erase(readFragment.size() - mapperOptions.trimRight);
        dataQuality.erase(dataQuality.size() - mapperOptions.trimRight);
    }

    // keep checks simple for now:
    if (readFragment.size() != dataQuality.size())
    {
        return 2;
    }

    //
    // on output of colorspace reads, we need the original read
    // for output to the SAM file.
    //
    if (gs->isColorSpace())
    {
        // save the read with the primer base
        ((MapperSEColorSpace *)this)->originalCSRead=readFragment;
        ((MapperSEColorSpace *)this)->originalCSQual=dataQuality;
    }

    // let setReadAndQuality conditionally truncate the first primer base
    // for color space reads
    //
    // setReadAndQuality has to be called once per read before
    // anything else is done with the class
    //

    if (setReadAndQuality(readFragment.c_str(), readFragment.size(), dataQuality.c_str()))
        return 1;

    return 0;
}


bool MapperBase::setReadAndQuality(const char *r, int len, const char *q)
{

    originalRead = r;
    originalQuality = q;

    // for color space reads,
    // we will convert char* r to base pair space
    // notice, color space and base space are treated differently:
    // base space:  forwardRead="ACGTA" backwardRead="ACGTA"
    // color space: forwardRead="ACGTA" backwardRead="TGCAT"
    // since if color space read is 01230, then reverse complement is 03210,
    // after map 0->A, ..., 3->T, actually, we are mapping forward ACGTA, the reverse complement ATGCA,
    // so we set the forward read to ACGTA, and set backward read to ATGCA.


    if (gs->isColorSpace() && isalpha(r[0]))
    {
        // Clip the primer base - it is not useful in our mapping code.
        // also clip the first color as it represents transition information
        // between the primer base and the first real base.
        forward.setReadAndQuality(r + 2, len - 2, q + 2);
        backward.setReadAndQuality(r + 2, len - 2, q + 2);
    }
    else
    {
        // base space - no translations are necessary - setReadAndQuality takes
        // care of converting the read to the backwards strand
        forward.setReadAndQuality(r,len,q);
        backward.setReadAndQuality(r,len,q);
    }

    //
    // Here we check the return code - if true, it means
    // we have too few valid index words for this read.
    //
    if (forward.setIndexStrategy() || backward.setIndexStrategy())
    {
        // TODO(zhanxw)
        // will remove this code eventually,
        // as we should call resetMapper() at begining of using Mapper.
        // when failing, make sure we clear out old data.
        getBestMatch().constructorClear();
        getBestMatch().indexer = &forward;
        return true;
    }
    return false;
}


//
// XXX we need two distinct modes of operation.
//
// 1) recurse to evaluate a single position
// 2) recurse to evaluate a group of positions
//
// We could:
// 1) add a second argument, evalFunctionType2, and if non-zero,
// we call it instead.
// 2) we pass a void * and a boolean to determine which one it is
// 3) turn these all into macros, which makes debug a nightmare,
//    but does solve the problem.
//
#define CANDIDATE_LIMIT 500
#if defined(CANDIDATE_LIMIT)
static int totalCandidateCount;
#endif

// a wrapper function
// for calling "evalBaseSpaceRead()" for forward & backward ReadIndex
void MapperBase::evalBaseSpaceReads(
                                    evalSinglePositionFunctionType onePositionMethod,
                                    evalAllPositionsFunctionType allPositionsMethod
                                    )
{
    assert((onePositionMethod==NULL) ^(allPositionsMethod==NULL));
#if defined(CANDIDATE_LIMIT)
    totalCandidateCount = 0;
#endif
    if (evalBaseSpaceRead(onePositionMethod, allPositionsMethod, forward)) return;
    if (evalBaseSpaceRead(onePositionMethod, allPositionsMethod, backward)) return;
}

// a wrapper function
// for calling "evalAllCandidatesForWord()" function on
// 1. every word stored in ReadIndexer class
// 2. every word stored (replace "N" with a letter) in ReadIndexer class
// 3. (optional) mutated every word stored in ReadIndexer class
bool MapperBase::evalBaseSpaceRead(
                                   evalSinglePositionFunctionType onePositionMethod,
                                   evalAllPositionsFunctionType allPositionsMethod,
                                   ReadIndexer &indexer)
{
    for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
    {
        if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, 0)) return true;
        if (indexer.wordNLocations[whichWord] > 0)
        {
            wordInteger_t mask1 = 1 << indexer.wordNLocations[whichWord];
            wordInteger_t mask2 = 2 << indexer.wordNLocations[whichWord];
            wordInteger_t mask3 = 3 << indexer.wordNLocations[whichWord];

            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask1)) return true;

            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask2)) return true;

            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask3)) return true;
        }
    }

    if (indexer.checkEdits)
    {

        wordInteger_t mask1 = 1;
        wordInteger_t mask2 = 2;
        wordInteger_t mask3 = 3;
        for (unsigned int i=0; i<wordIndex->wordSize; mask1<<=2, mask2<<=2, mask3<<=2, i++)
        {
            for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
            {

                // reads may now have invalid index words, but we
                // won't use them because of course they are invalid.
                //

                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask1)) return true;

                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask2)) return true;

                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask3)) return true;
            }
        }
    }
    return false;   // keep going
}

// a wrapper function
// for calling "evalBaseSpaceRead()" for forward & backward ReadIndex
void MapperBase::evalColorSpaceReads(
                                     evalSinglePositionFunctionType onePositionMethod,
                                     evalAllPositionsFunctionType allPositionsMethod
                                     )
{
    assert((onePositionMethod!=NULL) || (allPositionsMethod!=NULL));
    if (evalColorSpaceRead(onePositionMethod, allPositionsMethod, forward)) return;
    if (evalColorSpaceRead(onePositionMethod, allPositionsMethod, backward)) return;
}

// a wrapper function
// for calling "evalAllCandidatesForWord()" function on
// 1. every word stored in ReadIndexer class
// 2. every word stored (replace "N" with a letter) in ReadIndexer class
// 3. (optional) mutated every word stored in ReadIndexer class
bool MapperBase::evalColorSpaceRead(
                                    evalSinglePositionFunctionType onePositionMethod,
                                    evalAllPositionsFunctionType allPositionsMethod,
                                    ReadIndexer &indexer)
{
    for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
    {
        if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, 0)) return true;
        if (indexer.wordNLocations[whichWord] > 0)
        {
            wordInteger_t mask1 = 1 << indexer.wordNLocations[whichWord];
            wordInteger_t mask2 = 2 << indexer.wordNLocations[whichWord];
            wordInteger_t mask3 = 3 << indexer.wordNLocations[whichWord];

            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask1)) return true;
            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask2)) return true;
            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask3)) return true;
        }
    }

    if (indexer.checkEdits)
    {
        wordInteger_t mask1 = 1;
        wordInteger_t mask2 = 2;
        wordInteger_t mask3 = 3;
        for (unsigned int i=0; i<wordIndex->wordSize; mask1<<=2, mask2<<=2, mask3<<=2, i++)
        {
            for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
            {

                // reads may now have invalid index words, but we
                // won't use them because of course they are invalid.
                //
                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask1)) return true;
                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask2)) return true;
                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask3)) return true;
            }
        }

        // check single site SNPs
        unsigned int shift = 0;
        for (unsigned int i=0; i< wordIndex->wordSize; i++, shift+=2)
        {
            for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
            {
                if (evalAllCandidatesForColorSpaceWord(onePositionMethod, allPositionsMethod, indexer, whichWord, shift, 0)) return true;
                if (evalAllCandidatesForColorSpaceWord(onePositionMethod, allPositionsMethod, indexer, whichWord, shift, 1)) return true;
                if (evalAllCandidatesForColorSpaceWord(onePositionMethod, allPositionsMethod, indexer, whichWord, shift, 2)) return true;
            }
        }
    }
    return false;   // keep going
}

// a wrapper function
// to call either "allPositionMethod" (a function pointer passed in as parameter)
// or calling "evalAllCandidatePositions" function using parameter "onePositionMethod"
bool MapperBase::evalAllCandidatesForWord(
                                          evalSinglePositionFunctionType onePositionMethod,
                                          evalAllPositionsFunctionType allPositionsMethod,
                                          ReadIndexer &indexer,
                                          unsigned int whichWord,
                                          wordInteger_t xorMask)
{
    int count = 0;
    genomeIndex_t *candidates = NULL;

    if (wordIndex->wordReachedCutoff(indexer.wordInts[whichWord]^xorMask))
    {
        if (whichWord==indexer.wordInts.size()-1)
        {
            count = wordHashLeft->findGenomeLocations(indexer.wordInts[whichWord]^xorMask, indexer.wordInts[whichWord-1], candidates);
        }
        else
        {
            count = wordHashRight->findGenomeLocations(indexer.wordInts[whichWord]^xorMask, indexer.wordInts[whichWord+1], candidates);
        }
    }
    else
    {
        count = wordIndex->hashindices[(indexer.wordInts[whichWord]^xorMask) + 1] - wordIndex->hashindices[indexer.wordInts[whichWord]^xorMask];
        candidates = &wordIndex->wordpositions[wordIndex->hashindices[indexer.wordInts[whichWord]^xorMask]];
    }

#if  0
    // check how many count and where the word matched positions
    if (count > 0)
    {
        std::cerr << "count = " << count << " candidates = ";
        for (int i = 0; i< count; i++)
            std::cerr<< candidates[i] << ",";
        std::cerr<<" at " << __FILE__ << ":" << __LINE__ << std::endl;
    }
#endif

    if (allPositionsMethod)
    {
        return (*allPositionsMethod)(
                                     this,
                                     indexer,
                                     count,
                                     candidates,
                                     whichWord
                                     );
    }
    else
    {
        return evalAllCandidatePositions(
                                         onePositionMethod,
                                         indexer,
                                         whichWord,
                                         count,
                                         candidates
                                         );
    }

    return false;
}

// a wrapper function
// to call either "allPositionMethod" (a function pointer passed in as parameter)
// or calling "evalAllCandidatePositions" function using parameter "onePositionMethod"
bool MapperBase::evalAllCandidatesForColorSpaceWord(
                                                    evalSinglePositionFunctionType onePositionMethod,
                                                    evalAllPositionsFunctionType allPositionsMethod,
                                                    ReadIndexer &indexer,
                                                    unsigned int whichWord,
                                                    wordInteger_t shiftLocation,
                                                    wordInteger_t mutationIndex)
{

    wordInteger_t mask = 15;
    mask<<= shiftLocation;
    unsigned int adjacentColorCodes = (indexer.wordInts[whichWord] & mask) >>shiftLocation;
    wordInteger_t word = (indexer.wordInts[whichWord] & ~mask) |
        (colorSpaceSNP[adjacentColorCodes][mutationIndex] << shiftLocation);
    word &= wordIndex-> wordsCountIndexMask; // cap "word" to its maximum legal value

    int count = 0;
    genomeIndex_t *candidates = NULL;

    if (wordIndex->wordReachedCutoff(word))
    {
        if (whichWord==indexer.wordInts.size()-1)
        {
            count = wordHashLeft->findGenomeLocations(word, indexer.wordInts[whichWord-1], candidates);
        }
        else
        {
            count = wordHashRight->findGenomeLocations(word, indexer.wordInts[whichWord+1], candidates);
        }
    }
    else
    {
        count = wordIndex->hashindices[(word) + 1] - wordIndex->hashindices[word];
        candidates = &wordIndex->wordpositions[wordIndex->hashindices[word]];
    }

    if (allPositionsMethod)
    {
        return (*allPositionsMethod)(
                                     this,
                                     indexer,
                                     count,
                                     candidates,
                                     whichWord
                                     );
    }
    else
    {
        return evalAllCandidatePositions(
                                         onePositionMethod,
                                         indexer,
                                         whichWord,
                                         count,
                                         candidates
                                         );
    }

    return false;
}

inline bool MapperBase::evalAllCandidatePositions(
                                                  evalSinglePositionFunctionType pmethod,
                                                  ReadIndexer &indexer,
                                                  int     whichWord,
                                                  int     candidateCount,
                                                  genomeIndex_t *candidates)
{
    for (int i = 0; i < candidateCount; i ++)
    {
#if defined(LIMIT_CANDIDATES)
        if (++totalCandidateCount>LIMIT_CANDIDATES) return true;
#endif
        // avoid underflow in next subtract:
        if (indexer.wordPositions[whichWord] > candidates[i]) continue;

        // here, candidate position can (rarely) be smaller than the read
        // position ... see Test.cpp for test case
        //
        genomeIndex_t genomeMatchPosition =  candidates[i] - indexer.wordPositions[whichWord];

        // if genomeMatchPosition goes beyond the whole length of the genome
        if ((genomeMatchPosition + indexer.read.size()) > gs->sequenceLength())
            continue;

        // for color space, since we chopped the first two base, we should check the beginning part of the sequence
        if (genomeMatchPosition < 2 && indexer.gs->isColorSpace())
            continue;

        // check if this position has already been evalulated;
        if (indexer.checkedPositions.Find(genomeMatchPosition) != -1)
            continue;

        indexer.checkedPositions.Add(genomeMatchPosition);
        if ((*pmethod)(this, indexer, genomeMatchPosition, whichWord)) return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////
// Debug code
//////////////////////////////////////////////////////////////////////

void MapperBase::debugPrint(MatchedReadBase &matchedRead)
{
    std::string dataQuality;

    for (vector<uint8_t>::iterator it = forward.binaryQuality.begin() ; it < forward.binaryQuality.end(); it++)
        dataQuality.push_back(*it);

    if (!matchedRead.isForward())
        reverse(dataQuality.begin(), dataQuality.end());

    gs->debugPrintReadValidation(
                                 matchedRead.indexer->read,
                                 dataQuality,
                                 matchedRead.isForward() ? 'F' : 'R',
                                 matchedRead.genomeMatchPosition,
                                 matchedRead.quality,
                                 matchedRead.mismatchCount,
                                 false);
}

//////////////////////////////////////////////////////////////////////
// Obsolete code
//////////////////////////////////////////////////////////////////////
#ifdef COMPILE_OBSOLETE_CODE
int MapperBase::Word2Integer(std::string & word, unsigned int index, int &countNbases)
{
    assert(word.size() + index >= wordIndex->wordSize && index >= 0);

    int crtbaseinteger;
    int wordinteger = 0;

    //
    // XXX this is classic indication of a need to subclass, but I don't
    // know how practical it is for karma.
    //
    if (gs->isColorSpace())
    {
        //
        // NB: the + 1 is because the first character of the
        // color space read is a base pair key which we ignore here.
        // This is not a clean solution ... find something better.
        //
        for (unsigned int i = index + 1; i < index + wordIndex->wordSize + 1; i++)
        {
            wordinteger <<= 2;
            if (!isdigit(word[i])) return INVALID_WORDINDEX;
            crtbaseinteger = word[i] & 0xf; // works for ASCII and EBCDIC! whoo!
            wordinteger |= crtbaseinteger;
        }
    }
    else
    {
        for (unsigned int i = index; i < index + wordIndex->wordSize; i++)
        {
            wordinteger <<= 2;
            crtbaseinteger = GenomeSequence::base2int[(int) word[i]];

            switch (crtbaseinteger)
            {
            case GenomeSequence::baseXIndex:
                return INVALID_WORDINDEX;   // policy issue - should a read with a single invalid base be tossed?
            case GenomeSequence::baseNIndex:
                crtbaseinteger = 0;         // knowing this is 0, we'll visit all edits of this base later
            default:
                wordinteger |= crtbaseinteger;
            }
        }
    }

    return wordinteger;
}

//
// Given the short read input file in FASTQ format, populate the read fragment,
// data quality, fragment title (header/tag/whatever), and update the line #
// and # of reads.  Consider refactoring - a little messy looking.
//
// Returns:
//   1->read too short (deprecated)
//   2->read and quality lengths differ
//   3->too few index words
//
// Aborts on mal-formed FASTQ file.  This is a bit harsh, but does at least
// prevent ugly termination later when trying to read, e.g. a binary file.
//
int MapperBase::getReadAndQuality(IFILE f)
{
#pragma warning obsolete method, use FastqReader class 

    std::string ignore;
    std::string dataQuality;
    std::string readFragment;

    while (!ifeof(f))
    {
        line ++;

        // read line for readFragment name;
        f >> fragmentTag;

        // looking for first non-empty line
        if (fragmentTag.size() > 0) break;
    }

    if (ifeof(f)) return EOF;

    numberReads ++;

    if (fragmentTag[0] != '@')
        error("Tag line must start with [@]. check line %lu\n", line);

    // read actual readFragment
    line ++;
    f >> readFragment;
    if (readFragment.size() == 0)
        error("unexpected empty short read DNA strand line %lu\n", line);

    line++;
    f >> ignore; // ignore - it's a repeat of first line

    line++;
    f >> dataQuality; // quality data
    if (dataQuality.size() == 0)
        error("unexpected empty DNA quality line %lu\n", line);


    return processReadAndQuality(fragmentTag, readFragment, dataQuality);
}


#endif 
