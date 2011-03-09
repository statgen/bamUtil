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
#include "MapperSEBaseSpace.h"
#include "MapperPEBaseSpace.h"

#include <algorithm>
#include <vector>

MapperPEBaseSpace::MapperPEBaseSpace()
{
    mapperSE = new MapperSEBaseSpace;
    assert(mapperSE != NULL);
}

MapperPEBaseSpace::~MapperPEBaseSpace()
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
    return ((MapperPEBaseSpace*)mapper)->mapReads(indexer, whichWord, candidateCount, candidates);
}


//
// mapReads:
//     pairedReadB -> the read with higher match counts
//
// We're basically taking a pair of reads, 'this' is the one with the lower count as
// determined by the above method MapperPEBaseSpace::getMatchCountWithMutations.
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
void MapperPEBaseSpace::mapReads(MapperPE *pairedReadB)
{
    // XXX messy here and in single end, and worse, we handle it quite
    // differently in top level logic, but it should accomplish the same
    // thing...

    //
    // first --- if we aren't forcing smith waterman, attempt
    // a non gapped alignment
    //
    if (!mapperOptions.forceSmithWaterman)
        mapReads2(pairedReadB);

    //
    // second --- if we either are forcing smith waterman or we
    // failed to get a match and we are allowing it, get set to run
    // it.
    //
    if (mapperOptions.forceSmithWaterman || (mapperOptions.allowSmithWaterman && !bestMatch.qualityIsValid() && !pairedReadB->bestMatch.qualityIsValid()))
    {

        if (!mapperOptions.forceSmithWaterman)
        {
            //
            // clear our the old match information, we're starting over.
            //
            this->resetMapper();
            pairedReadB->resetMapper();
        }

        //
        // NB: yuuuuck - however, it gets the job done.
        // The underlying problem with these vile flags is
        // that if we don't do it this way, we have to pass
        // gapped flags down a ton of calls, and it gets tedious
        // and messy.
        //
        forward.useGapped = true;
        backward.useGapped = true;
        pairedReadB->forward.useGapped = true;
        pairedReadB->backward.useGapped = true;

        this->mappingMethod = PE_MODE;
        pairedReadB->mappingMethod = PE_MODE;
        mapReads2(pairedReadB); // now this does gapped alignment

        //
        // restore default non gapped behavior of the
        // indexers for both ends
        //
        forward.useGapped = false;
        backward.useGapped = false;
        pairedReadB->forward.useGapped = false;
        pairedReadB->backward.useGapped = false;

    }
}

void MapperPEBaseSpace::mapReads2(MapperPE *pairedReadB)
{
    // feh... do this better...
    // need this so eval code can search sorted match candidates
    pairedReadB->otherEndMapper = this;
    otherEndMapper = pairedReadB;


#if 0
    matchCandidate.constructorClear();
    pairedReadB->matchCandidate.constructorClear();     // XXX should be unused here
    bestMatch.constructorClear();
    pairedReadB->bestMatch.constructorClear();

    // assure that there is a default indexer hanging around in
    // case there are no maps:
    bestMatch.indexer = &forward;
    pairedReadB->bestMatch.indexer = &pairedReadB->forward;
#endif

    // Chain reaction here:
    // populateMatchCandidates()  in MapperPE                                                         =>
    // evalBaseSpaceReads(NULL, evalTrampoline)                                                       =>
    // evalBaseSpaceRead(NULL, evalTrampoline, forward/backward) in MapperBase                        =>
    // evalAllCandidatesForWord(NULL, evalTrampoline, indexer, whichWord, masks) in MapperBase        =>
    // ((MapperPE*)mapper)->populateMatchCandidates(indexer, whichWord, candidateCount, candidates);  by evalTrampoline
    //
    // this actually no longer fails, but check here anyway
    if (pairedReadB->populateMatchCandidates(false))
    {
        return;
    }

    forward.mismatchCutoff =
        backward.mismatchCutoff =
        2.5 + 4 * forward.read.size() * (mapperOptions.expectedErrorRate + mapperOptions.expectedSNPRate);


    pairedReadB->forward.mismatchCutoff =
        pairedReadB->backward.mismatchCutoff =
        2.5 + 4 * pairedReadB->forward.read.size() * (mapperOptions.expectedErrorRate + mapperOptions.expectedSNPRate);


    //  Another chain reaction:
    // evalBaseSpaceReads(NULL, evalTrampoline)                                                               =>
    // evalBaseSpaceRead(NULL, evalTrampoline, forward/backward) in MapperBase                                =>
    // evalAllCandidatesForWord(NULL, evalTrampoline, indexer, whichWord, masks) in Mapper Base               =>
    //     ((MapperPEBaseSpace*)mapper)->mapReads(indexer, whichWord, candidateCount, candidates) from evalTrampoline
    evalBaseSpaceReads(NULL, evalTrampoline);
}

//
// we don't want to keep calling the constructor for this,
// so make it static here (used only in ::mapReads below).
//
static MatchedReadPE   compareHelper;

//
// examine a candidate in read A (*this), checking against the
// map of candidates for read B.
//
// The set of matches for read A is not stored, only the
// running best match.
//
// The set of match candidates for read B is stored in matchCandidatesPointers,
// an ordered vector of pointers to the possible match candidates for read B.
// This gives us fast lookup for those possible positions on B that match our
// location.
//
// Normally, we don't expect there to be multiple matches for a given genome
// position, but because we are checking against a small range of positions
// on probe B (+/- genomePositionFilterWidth), there is a chance that there
// will be more than one match, therefore we use the STL ordered vector
// methods lower_bound and upper_bound to get all possible matches - again
// we normally expect zero or one, but there could be more.
//
// The read B quality score is not computed until we actually have a pair to
// evaluate.
//
// Performance tests to date (1/26/2009) show that evaluating list B is by
// far the slowest operation, so for the time being, no optimization has
// been done here.
//
bool MapperPEBaseSpace::mapReads(
                                 ReadIndexer         &indexer,
                                 int                 whichWord,
                                 int                 candidateCount,
                                 genomeIndex_t       *candidates
                                 )
{
    //    otherEndMapper->printMatchCandidates();
    //
    localGappedAlignment = false;
    otherEndMapper->localGappedAlignment = false;

    isProperAligned = true;
    otherEndMapper->isProperAligned = true;

    matchCandidate.indexer = &indexer;

    for (int i=0; i<candidateCount; i++)
    {
        // avoid underflow (genomeMatchPosition is unsigned,
        // so the result could be a very large # otherwise):
        if (candidates[i] < indexer.wordPositions[whichWord]) continue;
        matchCandidate.genomeMatchPosition = candidates[i] - indexer.wordPositions[whichWord];
        matchCandidate.quality = MatchedReadBase::UNSET_QUALITY;

        debugCheckPoint(matchCandidate.genomeMatchPosition, "::mapReads, top of loop for anchor");

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
            debugCheckPoint(matchCandidate.genomeMatchPosition, "::mapReads, location already checked");
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
        // NB: std::equal_range() returns a pair with the upper and lower
        // bound of interest to us.  It is about 50-60% the speed of just
        // calling them directly.  It might be due to the cost of implementation
        // or even just the additional overhead of creating the pair.
        //
        // I also attempted a shortcut for the case where there is only a single
        // entry in the table to be searched, but this was about 80% of the
        // simple minded code.
        //
        // In this instance, std::distance() is slower than (end-begin) by a wide margin.
        //
        //
        MatchCandidatesPointers_t::iterator low, high, it;
#if 0
        // old code always calls lower_bound, upper_bound - no clear performance
        // gain, but I need to check this change in.  We can revert later if need be.
        compareHelper.genomeMatchPosition = matchCandidate.genomeMatchPosition - mapperOptions.genomePositionFilterWidth;
        low = std::lower_bound(otherEndMapper->matchCandidatesPointers.begin(), otherEndMapper->matchCandidatesPointers.end(), &compareHelper, MatchedReadPEComp);
        compareHelper.genomeMatchPosition = matchCandidate.genomeMatchPosition + mapperOptions.genomePositionFilterWidth;
        high = std::upper_bound(otherEndMapper->matchCandidatesPointers.begin(), otherEndMapper->matchCandidatesPointers.end(), &compareHelper, MatchedReadPEComp);

#else
        // I don't normally try to hand optimize, but here is a case that
        // I'm not confident the compiler is doing good common subexpression
        // analysis on.
        //
        // XXX
        // XXX
        // XXX this code functions differently from the above code - FIND OUT WHY.
        // XXX
        // XXX
        //
        // When I test with paired end phiX reads, I get better and more reads mapped
        // using this little snippet.  This should NOT happen.
        //
        MatchCandidatesPointers_t::iterator begin = otherEndMapper->matchCandidatesPointers.begin();
        MatchCandidatesPointers_t::iterator end = otherEndMapper->matchCandidatesPointers.end();
        if ((end - begin) < 2)
        {
            low = begin;
            high = end;
        }
        else
        {
            compareHelper.genomeMatchPosition = matchCandidate.genomeMatchPosition - mapperOptions.genomePositionFilterWidth;
            low = std::lower_bound(begin, end, &compareHelper, MatchedReadPEComp);
            compareHelper.genomeMatchPosition = matchCandidate.genomeMatchPosition + mapperOptions.genomePositionFilterWidth;
            high = std::upper_bound(begin, end, &compareHelper, MatchedReadPEComp);
        }
#endif

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

            debugCheckPoint(matchCandidateB.genomeMatchPosition, "::mapReads, checking close by reads");

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
                matchCandidate.quality = indexer.getSumQ(
                                                         matchCandidate.genomeMatchPosition,
                                                         matchCandidate.mismatchCount,
                                                         bestMatch.quality,
                                                         whichWord
                                                         );
                matchCandidate.whichWord = whichWord;
#if 0
                //
                // if we are too many edits away, just stop processing
                //
                // XXX WRONG, MOST LIKELY:
                //
                if ((matchCandidate.mismatchCount > indexer.mismatchCutoff))
                {
                    debugCheckPoint(matchCandidate.genomeMatchPosition, "::mapReads, anchor has too many mismatches");
                    break;
                }
#endif

                //
                // if we have an invalid quality, we're done in this loop.
                //
                // likewise, if we have too low a quality, then there is no point
                // in continuing (check this assumption!)
                //
                // Tested against 100 base reads with no mutations, and this worked
                // great. Beware that terminating this search too soon could badly affect
                // some technologies like LS454.
                //
                // Gives a 20% speedup!
                //
                if (!matchCandidate.qualityIsValid() || matchCandidate.quality >= MAX_Q2P)
                {
                    debugCheckPoint(matchCandidate.genomeMatchPosition, "::mapReads, anchor too low quality");
                    break;
                }


                // XXX as written, this is still accumlating some VERY small values... e.g. 10e-60
                if (matchCandidate.quality < MAX_Q2P)
                    bestMatch.cumulativePosteriorProbabilities += sumQualityToProb[matchCandidate.quality];
            }
            // keep track of contributors
            bestMatch.numMatchContributors ++;

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
                matchCandidateB.quality = matchCandidateB.indexer->getSumQ(
                                                                           matchCandidateB.genomeMatchPosition,
                                                                           matchCandidateB.mismatchCount,
                                                                           MAX_Q2P,        // XXX FIX THIS -
                                                                           matchCandidateB.whichWord
                                                                           );

                if (!matchCandidateB.qualityIsValid() || matchCandidateB.quality >= MAX_Q2P)
                {
                    debugCheckPoint(matchCandidate.genomeMatchPosition, "::mapReads, matchCandidateB too low quality");
                    continue;
                }
#if 0
                //
                // XXX probably not what we want:
                //
                if (matchCandidateB.mismatchCount> matchCandidateB.indexer->mismatchCutoff*3)
                {
                    debugCheckPoint(matchCandidate.genomeMatchPosition, "::mapReads, matchCandidateB has too many mismatches");
                    continue;
                }
#endif

            }

#if 0
            //
            // debug code to show all possible matches... leave it here, please.
            //
            std::string s;
            gs->getChromosomeAndIndex(s, matchCandidateB.genomeMatchPosition);
            std::cout << "eval position " << matchCandidateB.genomeMatchPosition << " " << s << " this position quality = " << matchCandidateB.quality << std::endl;
#endif


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

    // safe guard bestMatch & otherEndMapper->bestMatch indexer pointer NOT to be null
    if (bestMatch.quality == MatchedReadBase::UNSET_QUALITY &&
        otherEndMapper->bestMatch.quality == MatchedReadBase::UNSET_QUALITY)
    {
        bestMatch.indexer=& forward;
        otherEndMapper->bestMatch.indexer= & (otherEndMapper->forward);
        isProperAligned = false;
        otherEndMapper->isProperAligned = false;
    }
    return false;
}

//
// XXX TODO: refactor this code.
//
// Attempt to do a local alignment.
//
// The goal appears to be to try local alignment so that the caller may
// determine whether or not it wishes to use the result.
//
bool MapperPEBaseSpace::tryLocalAlign(MapperBase* anchor)
{
    genomeIndex_t localAlignmentSize;
    int64_t localAlignmentPosition;
    genomeIndex_t   optimalLocalAlignment;

    // here mapperSE serves as localMapper
    if (mapperSE->processReadAndQuality(fragmentTag, originalRead, originalQuality) != 0)
    {
        // can't continue -- really the caller probably should
        // have already failed by now anyway.
        return false;
    }

    //
    // attempt to find a window of reasonable size
    //
    // mapperOptions.genomePositionFilterWidth is used, but is misnamed.
    //
    // the width should determine how far on either side of the anchor
    // postion we look, not just the overall width.  Otherwise, it is
    // misleading to most users (this one, for example).
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

    // TODO(zhanxw): uninitialized temporary variables, need to consider more!!!!
    int mismatchCount = 0, quality = 0, gapOpenCount = 0, gapCloseCount = 0, gapExtensionCount = 0;
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
        // fails obvious QC check here, so return failure
        return false;
    }
    else
    {
        //
        // XXX fill in gapOpenCount and mismatchCount based on cigar string?
        //
        // XXX technically, there are two cases this code should get called
        // 1) we get a good single end match, and guess that the other might
        // be nearby.
        // 2) we get two matches, one has high mismatch, and a realignment
        // picks up a local gapped alignment - in this case, we can re-use
        // the numBest and numMatchContributors values.
        //
        this->mapperSE->bestMatch.genomeMatchPosition = optimalLocalAlignment;
        this->mapperSE->bestMatch.numBest = 1;
        this->mapperSE->bestMatch.numMatchContributors = 1;
        this->mapperSE->bestMatch.mismatchCount = mismatchCount;
        this->mapperSE->bestMatch.quality = anchor->getBestMatch().quality; // quality should be at most anchor's quality
        if (!anchor->getBestMatch().indexer->isForward)
            this->mapperSE->bestMatch.indexer = &(this->mapperSE->forward);
        else
            this->mapperSE->bestMatch.indexer = &(this->mapperSE->backward);
        this->mapperSE->localGappedAlignment = true;   // this determines where we get the CIGAR string from
        return true;
    }
}

int MapperPEBaseSpace::test(int testNum,
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

    int c = gs->getChromosome(getBestMatch().genomeMatchPosition);
    genomeIndex_t i = getBestMatch().genomeMatchPosition - gs->getChromosomeStart(c) + 1;

    check(rc, testNum, "chromosome1", chr1, c + 1);
    check(rc, testNum, "index1", index1, i);


    check(rc, testNum, "direction2", direction2, otherMapper->getBestMatch().isForward() ? 'F' : 'R');
    check(rc, testNum, "mismatch2", misMatches2, otherMapper->getBestMatch().mismatchCount);
    //    check(rc, testNum, "quality1", quality1, otherMapper->getBestMatch().quality);

    c = gs->getChromosome(otherMapper->getBestMatch().genomeMatchPosition);
    i = otherMapper->getBestMatch().genomeMatchPosition - gs->getChromosomeStart(c) + 1;

    check(rc, testNum, "chromosome2", chr2, c + 1);
    check(rc, testNum, "index2", index2, i);

    return rc;

}
