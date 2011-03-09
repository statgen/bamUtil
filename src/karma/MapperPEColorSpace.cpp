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

#include "MapperPE.h"
#include "MappingStats.h"
#include "ReadsProcessor.h"
#include "MathConstant.h"
#include "Performance.h"
#include "Util.h"
#include "MapperSEColorSpace.h"
#include "MapperPEColorSpace.h"

#include <algorithm>
#include <vector>

MapperPEColorSpace::MapperPEColorSpace()
{
    //
    // This set of three vectors allows us to merge sort the sets of
    // words by shuffling pointers around to the match candidates
    // rather than the who MatchedReadPE object.
    //
    // This turned out to be weak because matchCandidatesPointers
    // will be invalidate if matchCandidates gets reallocated (grown), and
    // likewise the interators in matchCandidatesIndex will be invalidated
    // if matchCandidatesPointers is reallocated (grown).
    //
    matchCandidates.reserve(240*5000);
    matchCandidatesPointers.reserve(240*5000);
    matchCandidatesIndex.reserve(240);

    //
    // when we count the number of possible matches for a read, store
    // them in these counters - this will give us an idea of the cost
    // of doing the merge.
    //
    forwardCount = 0;
    backwardCount = 0;

    mapperSE = new MapperSEColorSpace;
    assert(mapperSE != NULL);

}

MapperPEColorSpace::~MapperPEColorSpace()
{
    if (mapperSE)
        delete mapperSE;
}

//
// this helper function is because getting pointer to method to work the way I want
// it to is a big pain in the backside.
//
static bool evalTrampoline(
                           MapperBase *mapper,
                           ReadIndexer &indexer,
                           int candidateCount,
                           genomeIndex_t *candidates,
                           int whichWord)
{
    return ((MapperPEColorSpace*)mapper)->mapReads(indexer, whichWord, candidateCount, candidates);
}

//
// mapReads:
//     pairedReadB -> the read with higher match counts
//
// We're basically taking a pair of reads, 'this' is the one with the lower count as
// determined by the above method MapperPEColorSpace::getMatchCountWithMutations.
//
// The goal is to iterate over all mutations (and originals) for first/second forward/reverse
// words in both reads.  This is not trivial to do as an iterator, but rather a combination
// of steps, so needs to be unrolled for each case.  That is to say, we can't just
// construct a simple iterator and iterate over all of the cases - the coding has to
// be explicit.
//
// The end result of mapping using this method should be the same regardless of
// the order paired ends are called here with.  The order may be performance
// sensitive, and we believe will be better if 'this' has a lower word match count.
//
void MapperPEColorSpace::mapReads(MapperPE *pairedReadB)
{
    this->mappingMethod = PE_MODE;
    pairedReadB->mappingMethod = PE_MODE;

    // feh... do this better...
    // need this so eval code can search sorted match candidates
    pairedReadB->otherEndMapper = this;
    otherEndMapper = pairedReadB;

    matchCandidate.constructorClear();
    pairedReadB->matchCandidate.constructorClear();
    bestMatch.constructorClear();
    pairedReadB->bestMatch.constructorClear();

    bestMatch.indexer = &forward;
    pairedReadB->bestMatch.indexer = &pairedReadB->forward;

    if (pairedReadB->populateMatchCandidates(true))
        return;

    // ad-hoc setting up mismatch cutoff values for both mappers
    forward.mismatchCutoff =
        backward.mismatchCutoff =
        2.5 + 4 * forward.read.size() * (mapperOptions.expectedErrorRate + mapperOptions.expectedSNPRate);


    pairedReadB->forward.mismatchCutoff =
        pairedReadB->backward.mismatchCutoff =
        2.5 + 4 * pairedReadB->forward.read.size() * (mapperOptions.expectedErrorRate + mapperOptions.expectedSNPRate);

    evalColorSpaceReads(NULL, evalTrampoline);
    return;

}

static MatchedReadPE   compareHelper;

bool MapperPEColorSpace::mapReads(
                                  ReadIndexer         &indexer,
                                  int                 whichWord,
                                  int                 candidateCount,
                                  genomeIndex_t       *candidates
                                  )
{


    //    otherEndMapper->printMatchCandidates();

    matchCandidate.indexer = &indexer;

    for (int i=0; i<candidateCount; i++)
    {
        //
        // avoid underflow if the match position is smaller than the
        // index of the word (ie the match position is <~15 or so):
        //
        if (candidates[i] < indexer.wordPositions[whichWord]) continue;

        matchCandidate.genomeMatchPosition = candidates[i] - indexer.wordPositions[whichWord];
        matchCandidate.quality = MatchedReadBase::UNSET_QUALITY;

        //
        // for this read, see if we already checked this genome match position.
        // This happens often, for example when the map is perfect, and both
        // words match, so the same location will come up twice.  We only want
        // to evaluate it once.
        //
        // Note that if we find we've been here, we don't need to continue the loop,
        // since all locations have been checked already.
        //
        if (indexer.checkedPositions.Find(matchCandidate.genomeMatchPosition) != -1)
        {
            continue;
        }
        indexer.checkedPositions.Add(matchCandidate.genomeMatchPosition);

        //
        // for each of the possible word positions from the read with the fewer number
        // of matches, we need to first obtain the match position, then compute the
        // quality of the match at that match position.
        //
        // XXX LAZY evaluate this as well!
        //
        MatchCandidatesPointers_t::iterator low, high, it;
        if (matchCandidate.genomeMatchPosition > mapperOptions.genomePositionFilterWidth)
            compareHelper.genomeMatchPosition = matchCandidate.genomeMatchPosition - mapperOptions.genomePositionFilterWidth;
        else
            compareHelper.genomeMatchPosition = 0;

        low = std::lower_bound(otherEndMapper->matchCandidatesPointers.begin(), otherEndMapper->matchCandidatesPointers.end(), &compareHelper, MatchedReadPEComp);
        compareHelper.genomeMatchPosition = matchCandidate.genomeMatchPosition + mapperOptions.genomePositionFilterWidth;
        high = std::upper_bound(otherEndMapper->matchCandidatesPointers.begin(), otherEndMapper->matchCandidatesPointers.end(), &compareHelper, MatchedReadPEComp);

        //
        // Here finally we examine a set of match locations that are
        // within genomePositionFilterWidth base pairs of each other.
        //
#if 0
        // debugging aids
        int candidates=high - low;
        if (candidates>100)
        {
            printf("candidates=%d\n", candidates);
        }
#endif
        for (it=low ; it < high; it++)
        {
            MatchedReadPE &matchCandidateB = **it;

            //
            // Here we are finally evaluating a pair of reads to see if they should
            // be considered.
            //
            // this->bestMatch and otherEndMapper->bestMatch should get updated if it is
            // a good candidate.
            //
            // Even if it is not a good candidate, we need to increment match
            // contributors for both probes, and update their posterior probabilities.
            //

#if 0
            if (candidates>1000)
            {
                printf("%ld: gmp=%u, matchCandidateB.gmp=%u (%d apart) word=%d, .whichWord=%d, mutIndex=%d\n",
                       it - low,
                       matchCandidate.genomeMatchPosition,
                       matchCandidateB.genomeMatchPosition,
                       matchCandidate.genomeMatchPosition- matchCandidateB.genomeMatchPosition,
                       matchCandidateB.word,
                       matchCandidateB.whichWord,
                       matchCandidateB.wordMutationIndex
                       );
            }
#endif
            //
            // lazy compute quality of read A match candidate:
            //
            if (matchCandidate.quality == MatchedReadBase::UNSET_QUALITY)
            {
                matchCandidate.quality = indexer.getColorSpaceSumQ(
                                                                   matchCandidate.genomeMatchPosition,
                                                                   matchCandidate.mismatchCount,
                                                                   bestMatch.quality,
                                                                   whichWord
                                                                   );
                matchCandidate.whichWord = whichWord;
                //
                // if we are too many edits away, just stop processing
                //
                if ((matchCandidate.mismatchCount > forward.mismatchCutoff))
                {
                    break;
                }

                //
                // if we have an invalid quality, we're done in this loop.
                //
                if (!matchCandidate.qualityIsValid()) break;

                // XXX as written, this is still accumlating some VERY small values... e.g. 10e-60
                if (matchCandidate.quality < MAX_Q2P)
                    bestMatch.cumulativePosteriorProbabilities += sumQualityToProb[matchCandidate.quality];

                // numMatchContributors should be recorded only when compute new sumQ. - Xiaowei
                // keep track of contributors
                bestMatch.numMatchContributors ++;
            }

            //
            // Here we are looking at (usually one, but can be
            // more than one) candidate for a match against mate A.
            //
            // Since not all read B match candidates will be referenced,
            // we lazy evaluate them here - this is known to be a slow
            // operation, so it is worth avoiding if we can.
            //
            // Note here we use a number of values out of the candidate
            // from read B, such as match position, direction, etc - these
            // are the arguments to the getQvalue method, and these are why
            // we keep all these values around rather than just the genome
            // match position.
            //
            if (matchCandidateB.quality == MatchedReadBase::UNSET_QUALITY)
            {
                matchCandidateB.quality = matchCandidateB.indexer->getColorSpaceSumQ(
                                                                                     matchCandidateB.genomeMatchPosition,
                                                                                     matchCandidateB.mismatchCount,
                                                                                     MAX_Q2P,        // XXX FIX THIS -
                                                                                     matchCandidateB.whichWord
                                                                                     );
                if (matchCandidateB.mismatchCount> forward.mismatchCutoff)
                {
                    continue;
                }

            }

            //
            // our termination condition for mateB is a little different than for
            // mate A - since we're iterating over mate B possibilities, we want
            // to continue on to the next one via continue.
            //
            if (!matchCandidateB.qualityIsValid() || matchCandidateB.quality >= MAX_Q2P) continue;

            // keep track of contributors
            otherEndMapper->bestMatch.numMatchContributors ++;
            otherEndMapper->bestMatch.cumulativePosteriorProbabilities += sumQualityToProb[matchCandidateB.quality];
#if 0
            printf("quality = %d, cumulativePosteriorProbabilities=%f\n", matchCandidateB.quality, otherEndMapper->bestMatch.cumulativePosteriorProbabilities);
#endif
#if 0
            otherEndMapper->debugPrint(matchCandidateB);
#endif

            // At this point, we have a match we wish to consider.  For each
            // mate, we've accumulated their contributors and posterior probabilties,
            // so all that is left is to appropriately filter and update the bestMatches
            // for both reads.


#if 0
            // debug
            //
            printf("A: %c%u W:%d WI:%d Q: %d", direction, matchCandidate.genomeMatchPosition, whichWord, wordMutationIndex, matchCandidate.quality);
            printf("\tB: %c%u W:%d WI:%d Q: %d\n", matchCandidateB.direction, matchCandidateB.genomeMatchPosition, matchCandidateB.whichWord, matchCandidateB.wordMutationIndex, matchCandidateB.quality);
#endif

            if (updateBestMatch(matchCandidateB) == false)
                continue;
        } // iterate on matchCandidateB
    } // iterate on (this->matchCandidate)
    // safe guard bestMatch & otherEndMapper->bestMatch indexer pointer NOT to bu
    if (bestMatch.quality == MatchedReadBase::UNSET_QUALITY &&
        otherEndMapper->bestMatch.quality == MatchedReadBase::UNSET_QUALITY)
    {
        bestMatch.indexer=& forward;
        otherEndMapper->bestMatch.indexer= & (otherEndMapper->forward);
    }
    return false;
}
bool MapperPEColorSpace::tryLocalAlign(MapperBase* anchor)
{
    genomeIndex_t localAlignmentSize;
    int64_t localAlignmentPosition;
    genomeIndex_t   optimalLocalAlignment;

    // here mapperSE serves as localMapper
    if (mapperSE->processReadAndQuality(fragmentTag, forward.read, forward.phredQuality) != 0)
    {
        // can't continue -- really the caller probably should
        // have already failed by now anyway.
        return false;
    }
    assert(mapperSE->processReadAndQuality(fragmentTag, forward.read, forward.phredQuality)==0);

    localAlignmentSize = mapperOptions.genomePositionFilterWidth +
        anchor->forward.readLength +
        mapperSE->forward.readLength;
    localAlignmentPosition = anchor->getBestMatch().genomeMatchPosition - localAlignmentSize/2;

    //
    // the following two checks are done using
    // 64 bit arithmetic...
    //
    // check for reads mapped before localAlignmentSize/2:
    //
    if (localAlignmentPosition < 0) localAlignmentPosition=0;

    //
    // check for reads near the end of the reference:
    //
    // XXX we need similar checks for begin/end of a chromosome,
    // since we don't want to consider alignments that cross
    // chromosome boundaries.  Alternatively, we can do cleanup
    // at the very end before printing the alignment by forcing
    // a clip of some kind for the portion of a read that crosses
    // a chromosome boundary.
    //
    if (localAlignmentPosition + localAlignmentSize >
        gs->sequenceLength())
        localAlignmentPosition=gs->sequenceLength() - localAlignmentSize;

    // TODO: uninitialized temporary variables, need to consider more!!!!
    int mismatchCount, quality, gapOpenCount, gapCloseCount, gapExtensionCount;
    if (!anchor->getBestMatch().indexer->isForward)
    {
        optimalLocalAlignment = this->mapperSE->forward.localAlignment(
                                                                       (genomeIndex_t) localAlignmentPosition,
                                                                       localAlignmentSize,
                                                                       mismatchCount,
                                                                       quality,
                                                                       gapOpenCount,
                                                                       gapCloseCount,
                                                                       gapExtensionCount,
                                                                       this->mapperSE->cigarRoller
                                                                       );
    }
    else
    {
        optimalLocalAlignment = this->mapperSE->backward.localAlignment(
                                                                        (genomeIndex_t) localAlignmentPosition,
                                                                        localAlignmentSize,
                                                                        mismatchCount,
                                                                        quality,
                                                                        gapOpenCount,
                                                                        gapCloseCount,
                                                                        gapExtensionCount,
                                                                        this->mapperSE->cigarRoller
                                                                        );
    };
    if (mapperSE->cigarRoller.size() > 5 ||
        mapperSE->cigarRoller.getExpectedQueryBaseCount() != (int) mapperSE->forward.read.size())
    {
        return false;
    }
    else
    {
        this->mapperSE->bestMatch.genomeMatchPosition = optimalLocalAlignment;
        this->mapperSE->bestMatch.numBest = 1;
        this->mapperSE->bestMatch.numMatchContributors = 1;
        this->mapperSE->bestMatch.mismatchCount = mismatchCount;
        this->mapperSE->bestMatch.quality = quality;
        if (!anchor->getBestMatch().indexer->isForward)
            this->mapperSE->bestMatch.indexer = &(this->mapperSE->forward);
        else
            this->mapperSE->bestMatch.indexer = &(this->mapperSE->backward);
        this->mapperSE->localGappedAlignment = true;   // this determines where we get the CIGAR string from
        return true;
    }
}


int MapperPEColorSpace::test(int testNum,
                             MapperPE    *otherMapper,
                             const char *read1, const char *qual1, char direction1, int chr1, genomeIndex_t index1, int misMatches1,
                             const char *read2, const char *qual2, char direction2, int chr2, genomeIndex_t index2, int misMatches2
                             )
{
    std::string r1 = read1;    // read
    std::string q1 = qual1;    // quality
    std::string t1 = "";      // tag

    std::string r2 = read2;    // read
    std::string q2 = qual2;    // quality
    std::string t2 = "";      // tag

    int rc = 0;

    init(r1, q1, t1);

    otherMapper->init(r2, q2, t2);

    mapReads(otherMapper);

    check(rc, testNum, "direction1", direction1, getBestMatch().isForward() ? 'F' : 'R');
    check(rc, testNum, "mismatch1", misMatches1, getBestMatch().mismatchCount);
    //    check(rc, testNum, "quality1", quality1, getBestMatch().quality);

    int c = gs->getChromosome(this->getBestMatch().genomeMatchPosition);
    genomeIndex_t i = getBestMatch().genomeMatchPosition - gs->getChromosomeStart(c) + 1;

    check(rc, testNum, "chromosome1", chr1-1, c);
    check(rc, testNum, "index1", index1, i);


    check(rc, testNum, "direction2", direction2, otherMapper->getBestMatch().isForward() ? 'F' : 'R');
    check(rc, testNum, "mismatch2", misMatches2, otherMapper->getBestMatch().mismatchCount);
    //    check(rc, testNum, "quality1", quality1, otherMapper->getBestMatch().quality);

    c = gs->getChromosome(otherMapper->getBestMatch().genomeMatchPosition);
    i = otherMapper->getBestMatch().genomeMatchPosition - gs->getChromosomeStart(c) + 1;

    check(rc, testNum, "chromosome2", chr2-1, c);
    check(rc, testNum, "index2", index2, i);

    return rc;

}
