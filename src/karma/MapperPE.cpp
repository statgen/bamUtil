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
#include "MatchedReadPE.h"
#include "MappingStats.h"
#include "ReadsProcessor.h"
#include "MathConstant.h"
#include "Performance.h"
#include "Util.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

MapperPE::MapperPE()
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
    // The original estimate of 240 * 5000 was reasonable, but with
    // color space matching, "around" 4X more are necessary.
    //
    // XXX we need to paramaterize this, and we need to defend better
    // against the capacity changing (reallocation breaks the pointers
    // from one vector into the other).
    //

#pragma message "reserve() does not work..."

#if 1
    matchCandidates.reserve(240*5000 * 4);
    matchCandidatesPointers.reserve(240*5000 * 4);
    matchCandidatesIndex.reserve(240 * 4);
#else
    matchCandidates.reserve(240*500 * 4);
    matchCandidatesPointers.reserve(240*500 * 4);
    matchCandidatesIndex.reserve(240 * 4);

#endif

    //
    // when we count the number of possible matches for a read, store
    // them in these counters - this will give us an idea of the cost
    // of doing the merge.
    //
    forwardCount = 0;
    backwardCount = 0;

    mappingMethod = PE_MODE;
}

void  MapperPE::init(std::string & readFragment, std::string &dataQuality, std::string &fragmentTag)
{
    setReadAndQuality(readFragment.c_str(), readFragment.size(), dataQuality.c_str());
    this->fragmentTag = fragmentTag;
    forwardCount = 0;
    backwardCount = 0;
}

void MapperPE::printMatchCandidates()
{
    MatchCandidates_t::iterator it;

    //
    //
    for (it= matchCandidates.begin() ; it !=  matchCandidates.end(); it++)
    {
        printf("%u\n", it->genomeMatchPosition);
    }
}

//
// This is here to call from gdb to verify that the resulting list is
// actually in order.
//
void MapperPE::testMatchCandidates()
{
    MapperPE::testMatchCandidates(matchCandidatesIndex.begin());

    fprintf(stderr, "matchCandidatePointers appears to be in order.\n");
    fflush(stderr);
    return;
}

void MapperPE::testMatchCandidates(MatchCandidatesIndex_t::iterator i1)
{
    MatchCandidatesPointers_t::iterator it;

    fprintf(stderr, "Checking matchCandidatePointers (%d elements):\n",
            (int)((*i1).second - (*i1).first));
    fflush(stderr);
    for (it=(*i1).first; it < (*i1).second - 1; it++)
    {
        if (((*it)->genomeMatchPosition > (*(it+1))->genomeMatchPosition))
        {
            //
            fprintf(stderr, "matchCandidatePointers [%d] is out of order:\n",
                    (int)(it - matchCandidatesPointers.begin()));
            fprintf(stderr, "it.direction = '%c', genome_match_pos = %u\n",
                    (*it)->isForward() ? 'F' : 'R', (*it)->genomeMatchPosition);
            fprintf(stderr, "it+1.direction = '%c', genome_match_pos = %u\n",
                    (*(it+1))->isForward() ? 'F' : 'R', (*(it+1))->genomeMatchPosition);
            //            printMatchCandidates(i1);
            fflush(stderr);
            throw std::logic_error("critical validation error in testMatchCandidates");
            return;
        }
    }
}

void MapperPE::printMatchCandidates(MatchCandidatesIndex_t::iterator i1)
{
    MatchCandidatesPointers_t::iterator it;

    FILE *out;

    out = fopen("genome_out.txt","a");

    fprintf(out, "Printing genomeMatchPositions (%d elements) (including one extra):\n",
            (int)((*i1).second - (*i1).first));
    fflush(out);
    for (it=(*i1).first; it < (*i1).second; it++)
    {
        fprintf(out,"%d: %c%u\t", (int)(it - (matchCandidatesPointers.begin())),
                (*it)->isForward() ? 'F' : 'R',
                (*it)->genomeMatchPosition);
    }
    fprintf(out,"\n");
    fflush(out);
    fclose(out);
}

// #define PARANOID_CHECK_ORDER

//
// This merges from i1 to i2 INCLUSIVE
//
// The expected result is that the range from i1.first to i2.second
// is merged (i1.first inclusive, i2.second EXCLUSIVE).
//
MapperPE::MatchCandidatesIndex_t::iterator
MapperPE::mergeSortedMatchCandidates(
                                     MatchCandidatesIndex_t::iterator i1,
                                     MatchCandidatesIndex_t::iterator i2)
{
    // STL algorithm ordered list merge references:
    //
    // see wikipedia article - std::merge and std::inplace_merge
    //
    // http://www.cplusplus.com/reference/algorithm/inplace_merge.html
    //
    // http://www.cplusplus.com/reference/algorithm/merge.html
    //
    //
    // if the slots are the same, don't merge
    //
    if (i1 == i2)
    {
        return i1;
    }
    //
    // if the regions are adjacent, merge them:
    //
    if ((*i1).second==(*i2).first)
    {
        //
        // This is currently the slowest part of paired read matching in karma -
        // this loop allocates intermediate results, then copies the final
        // result to the inplace destination.  It isn't clear that this is
        // what we want.
        //
        // For one thing, it entails a mostly redundant final copy that we don't
        // really benefit from.  For another, the allocation/deallocations are
        // redundant in a sense since we recursively merge the same data structure
        // over and over for each read with known upper resource bounds.
        //
        // One option is to declare an allocator for the list - if it uses
        // alloca(), there won't be any linked list operations, etc, for
        // allocating/deallocating, since it is on the local stack.
        //
        // However, here again, the issue seems to be the cost of iterating
        // and merging the elements themselves, rather than the cost of
        // alloc/dealloc.
        //
#if defined(PARANOID_CHECK_ORDER)
        testMatchCandidates(i1);
        testMatchCandidates(i2);
        //        printMatchCandidates(i1);
        //        printMatchCandidates(i2);
        assert((*i1).second == (*i2).first);
#endif

        std::inplace_merge((*i1).first, (*i1).second, (*i2).second, MatchedReadPEComp);
        (*i1).second = (*i2).second;
#if defined(PARANOID_CHECK_ORDER)
        testMatchCandidates(i1);
        printMatchCandidates(i1);
#endif
        return i1;
    }

    //
    // here we use recursion to simplify the task.
    // some optimizations may be possible, but get it right first.
    //
    MatchCandidatesIndex_t::iterator middle, final;
    if (i2<i1)
    {
        throw std::logic_error("MapperPE::mergeSortedMatchCandidates");
    }

    // We have a set of ranges to merge, starting
    // with i1 and ending with i2 inclusive.
    //
    // We need to divide them into two non-overlapping
    // sets, then merge those.
    //
    // If we simply divide the regions, we can get overlaps
    // because for example, the middle can wind up being == i1.
    //
    middle = i1 + (i2 - i1) / 2;
#if defined(PARANOID_CHECK_ORDER)
    testMatchCandidates(i1);
    testMatchCandidates(i2);
#endif

    mergeSortedMatchCandidates(i1, middle-1);
#if defined(PARANOID_CHECK_ORDER)
    testMatchCandidates(i1);
    assert((*i1).second==(*(middle-1)).second);
#endif

    mergeSortedMatchCandidates(middle, i2);
#if defined(PARANOID_CHECK_ORDER)
    testMatchCandidates(middle);
    assert((*middle).second==(*i2).second);

    // XXX this breaks: but end of first should be adjacent to start
    // of second merged result
    assert((*i1).second == (*middle).first);
#endif

    mergeSortedMatchCandidates(i1, middle);

#if defined(PARANOID_CHECK_ORDER)
    testMatchCandidates(i1);
#endif
    return final;
}

//
// This method is intended to return a count of approximately
// how many match positions this read has.  However, it is flawed
// because it ignores high repeat count words and the use of the
// secondary 30-mer hash.
//
// It doesn't matter too much yet because it is just used as
// a heuristic to guess which end of the pair to use for iteration.
// It should give the same result either way.
//
void MapperPE::getMatchCountWithMutations()
{
    forwardCount = 0;
    backwardCount = 0;
    //
    for (unsigned int i=0; i<forward.wordInts.size(); i++)
    {
        forwardCount += wordIndex->getMatchCountWithMutations(forward.wordInts[i]);
    }
    for (unsigned int i=0; i<backward.wordInts.size(); i++)
    {
        backwardCount += wordIndex->getMatchCountWithMutations(backward.wordInts[i]);
    }
    return;
}

void MapperPE::populateCigarRollerAndGenomeMatchPosition()
{
    //
    // overload the base because we have more work to do
    // for PE mapper...
    //
    if (mapperSE)
    {
        // XXX oh the horror
        mapperSE->populateCigarRollerAndGenomeMatchPosition();
    }
    MapperBase::populateCigarRollerAndGenomeMatchPosition();
}

//
// print the two best matches
//
void MapperPE::printBestReads(std::ostream &file, MapperPE *longerMapper)
{
    // obtain bestMatch for the pair
    MatchedReadBase* pBestMatch;
    MatchedReadBase* pOtherBestMatch;

    MapperBase* pThisMapper;
    MapperBase* pOtherMapper;
    if (this->mappingMethod == PE_MODE)
    {
        pThisMapper = this;
        pBestMatch = &bestMatch;
    }
    else
    {
        pThisMapper = this->mapperSE;
        pBestMatch = &(this->mapperSE->bestMatch);
    }

    if (longerMapper->mappingMethod == PE_MODE)
    {
        pOtherMapper = longerMapper;
        pOtherBestMatch = &(longerMapper->bestMatch);
    }
    else
    {
        pOtherMapper = longerMapper->mapperSE;
        pOtherBestMatch = &(longerMapper->mapperSE->bestMatch);
    }

    // calculate cigar string => TODO: may just need to set e.g. 119M in most cases.
    pThisMapper->populateCigarRollerAndGenomeMatchPosition();
    pOtherMapper->populateCigarRollerAndGenomeMatchPosition();

    // print the matched pairs
    pBestMatch->print(file,
                      pOtherBestMatch,   // pass it the mate
                      pThisMapper->fragmentTag,
                      mapperOptions.showReferenceBases,
                      pThisMapper->cigarRoller,
                      pThisMapper->isProperAligned,
                      pThisMapper->samMateFlag,
                      mapperOptions.readGroupID,
                      pThisMapper->alignmentPathTag);

    pOtherBestMatch->print(file,
                           pBestMatch,                 // pass it the mate
                           pOtherMapper->fragmentTag,
                           mapperOptions.showReferenceBases,
                           pOtherMapper->cigarRoller,
                           pOtherMapper->isProperAligned,
                           pOtherMapper->samMateFlag,
                           mapperOptions.readGroupID,
                           pOtherMapper->alignmentPathTag
                           );
}

//
// print the two best matches
//
void MapperPE::printCSBestReads(std::ostream &file,
                                GenomeSequence* gs,
                                GenomeSequence* csgs,
                                MapperPE *longerMapper)
{
    // obtain bestMatch for the pair
    MatchedReadBase* pBestMatch;
    MatchedReadBase* pOtherBestMatch;

    MapperBase* pThisMapper;
    MapperBase* pOtherMapper;
    if (this->mappingMethod == PE_MODE)
    {
        pThisMapper = this;
        pBestMatch = &bestMatch;
    }
    else
    {
        pThisMapper = this->mapperSE;
        pBestMatch = &(this->mapperSE->bestMatch);
    }

    if (longerMapper->mappingMethod == PE_MODE)
    {
        pOtherMapper = longerMapper;
        pOtherBestMatch = &(longerMapper->bestMatch);
    }
    else
    {
        pOtherMapper = longerMapper->mapperSE;
        pOtherBestMatch = &(longerMapper->mapperSE->bestMatch);
    }

    // calculate cigar string => TODO: may just need to set e.g. 119M in most cases.
    pThisMapper->populateCigarRollerAndGenomeMatchPosition();
    pOtherMapper->populateCigarRollerAndGenomeMatchPosition();

    // print the matched pairs
    pBestMatch->printColorSpace(file,
                                gs, csgs,
                                pOtherBestMatch,   // pass it the mate
                                ((MapperBase*)this)->originalCSRead,
                                ((MapperBase*)this)->originalCSQual,
                                pThisMapper->fragmentTag,
                                mapperOptions.showReferenceBases,
                                pThisMapper->cigarRoller,
                                pThisMapper->isProperAligned,
                                pThisMapper->samMateFlag,
                                mapperOptions.readGroupID,
                                pThisMapper->alignmentPathTag);

    pOtherBestMatch->printColorSpace(file,
                                     gs, csgs,
                                     pBestMatch,                 // pass it the mate
                                     ((MapperBase*) longerMapper)->originalCSRead,
                                     ((MapperBase*) longerMapper)->originalCSQual,
                                     pOtherMapper->fragmentTag,
                                     mapperOptions.showReferenceBases,
                                     pOtherMapper->cigarRoller,
                                     pOtherMapper->isProperAligned,
                                     pOtherMapper->samMateFlag,
                                     mapperOptions.readGroupID,
                                     pOtherMapper->alignmentPathTag
                                     );
}

void MapperPE::debugDump(std::ostream &out)
{
    out << "forwardCount = " << forwardCount << std::endl;
    out << "backwardCount = " << backwardCount << std::endl;
}

MatchedReadBase &MapperPE::getBestMatch()
{
    assert(mapperSE) ;
    if ( mappingMethod == PE_MODE)
        return bestMatch;
    else
        return mapperSE->bestMatch;
}

bool MapperPE::updateBestMatch(MatchedReadPE& matchCandidateB)
{
    //
    // Score the paired read - this is simple minded, but gets us something -
    // just add them together.
    //
    // If either mate is invalid, we're done.
    //
    int pairQuality;
    if (matchCandidate.qualityIsValid() && matchCandidateB.qualityIsValid())
        pairQuality = matchCandidate.quality + matchCandidateB.quality;
    else
    {
        pairQuality = MatchedReadBase::UNSET_QUALITY;
        return false;
    }

    // if the new pair map quality can affect the posterior probabilties, add it in
    if (pairQuality < MAX_Q2P)
    {
        bestMatch.pairCumulativePosteriorProbabilities += sumQualityToProb[pairQuality];
        otherEndMapper->bestMatch.pairCumulativePosteriorProbabilities += sumQualityToProb[pairQuality];
    }

    // if we have matches of equal quality, we want to keep
    // track of how many there are, and randomly assign one.
    if (bestMatch.qualityIsUnset() || pairQuality < bestMatch.pairQuality)
    {
        bestMatch.numBest = 1;
        otherEndMapper->bestMatch.numBest = 1;
        // XXX seems unclean, but set both mates pair
        // quality to the same value:
        matchCandidate.pairQuality = pairQuality;
        matchCandidateB.pairQuality = pairQuality;

        bestMatch.updateMatch(matchCandidate);
        otherEndMapper->bestMatch.updateMatch(matchCandidateB);

    }
    else if (pairQuality == bestMatch.pairQuality)
    {
#if 0
        printf("TagA: %s A: %c%u W:%d WI:%d Q: %d", fragmentTag.c_str(), direction, matchCandidate.genomeMatchPosition, whichWord, wordMutationIndex, matchCandidate.quality);
        printf("\tB: %c%u W:%d WI:%d Q: %d\n", matchCandidateB.direction, matchCandidateB.genomeMatchPosition, matchCandidateB.whichWord, matchCandidateB.wordMutationIndex, matchCandidateB.quality);
#endif
        bestMatch.numBest ++;
        otherEndMapper->bestMatch.numBest ++;
        double random = rand->Uniform();
        if (random < (1.0/ (double) bestMatch.numBest))
        {
            // XXX unclean:
            matchCandidate.pairQuality = pairQuality;
            matchCandidateB.pairQuality = pairQuality;

            bestMatch.updateMatch(matchCandidate);
            otherEndMapper->bestMatch.updateMatch(matchCandidateB);
        }
    }
    else
    {
        return false;   // XXX ugly coding, but this is what we want to do here
    }

    //
    // if we did an update in the bestMatch, we now want to see if we should
    // or can terminate early due to reaching Q_CUTOFF on both mates.
    //
    // XXX verify that this is right -- we changed the single
    // end code to compute Q_CUTOFF for every read...
    //
    // XXX also confirm when we get to this code.
    //
    // XXX I think we want to replace the code below with a single
    // check on pairQuality.
    //
    if (bestMatch.cumulativePosteriorProbabilities > Q_CUTOFF &&
        otherEndMapper->bestMatch.cumulativePosteriorProbabilities > Q_CUTOFF)
    {

        //
        // if BOTH mates are above the cutoff, and both mates have no
        // candidates, then we can terminate
        //
        // XXX this seems unlikely to be right.
        // I would think that what we really want to
        // terminate on is the pairQuality.
        //
        if (bestMatch.quality == MatchedReadBase::UNSET_QUALITY &&
            otherEndMapper->bestMatch.quality == MatchedReadBase::UNSET_QUALITY)
        {

            bestMatch.quality = otherEndMapper->bestMatch.quality = MatchedReadBase::EARLYSTOP_QUALITY;
            bestMatch.pairQuality = otherEndMapper->bestMatch.pairQuality = MatchedReadBase::EARLYSTOP_QUALITY;
            isProperAligned = otherEndMapper->isProperAligned = false;
            return true;
        }

        //
        // if BOTH mates are above the cutoff, and were also both above the
        // cutoff before we added in the new values, then we can terminate.
        //
        // XXX not sure it makes sense to track individual posterior
        // qualities.  The only reason it could help us is if we
        // could use one mate later to remap the other (we will do
        // this, but we have to ignore the scoring using this algorithm,
        // as it no longer applies).
        //
        if (
            (bestMatch.cumulativePosteriorProbabilities - sumQualityToProb[bestMatch.quality]) > Q_CUTOFF &&
            (otherEndMapper->bestMatch.cumulativePosteriorProbabilities - sumQualityToProb[otherEndMapper->bestMatch.quality]) > Q_CUTOFF)
        {

            bestMatch.quality = otherEndMapper->bestMatch.quality = MatchedReadBase::EARLYSTOP_QUALITY;
            bestMatch.pairQuality = otherEndMapper->bestMatch.pairQuality = MatchedReadBase::EARLYSTOP_QUALITY;
            isProperAligned = otherEndMapper->isProperAligned = false;

            return true;
        }


    } // check Q_CUTOFF

    return false;
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
    return ((MapperPE*)mapper)->populateMatchCandidates(indexer, whichWord, candidateCount, candidates);
}

//
// This populates several data structures optimized for rapid sorting of
// genome match candidate positions by their genome index.
//
// The Word Index contains genome indeces for words in ascending genome
// index order, so we choose to approach this problem as the final
// stage of a merge sort.
//
// The elements being sorted are matchCandidate objects.
//
// Because their size is significantly larger than one integer, we'll
// do the merge sort via pointers to them, rather than moving the
// elements themselves.
//
// To keep track of the sets of ordered words, we'll create a vector
// of pairs of indeces into the array of pointers.
//
// When doing the merge, we'll work through matchCandidatesIndex
// until it has one single pair of indeces - the first and last elements
// of the vector matchCandidatesPointers.
//
// At that point, matchCandidatesPointers will list the match candidates
// in genome sequence order so that we can do binary searches on that
// list for matching candidates on the other read.
//
bool MapperPE::populateMatchCandidates(bool isColorSpace)
{

    matchCandidates.clear();
    matchCandidatesPointers.clear();
    matchCandidatesIndex.clear();

    if (!isColorSpace)
        evalBaseSpaceReads(NULL, evalTrampoline);
    else
        evalColorSpaceReads(NULL, evalTrampoline);

    MatchCandidatesIndex_t::iterator finalIndex;
    //
    // Merge from start to end inclusive.  Check for an empty
    // list before we try to merge it!
    //
    if (matchCandidatesIndex.begin() != matchCandidatesIndex.end())
    {
        MatchCandidatesIndex_t::iterator begin = matchCandidatesIndex.begin();
        MatchCandidatesIndex_t::iterator end = matchCandidatesIndex.end() - 1;

#if 0
        printMatchCandidates(begin);
#endif

        finalIndex = mergeSortedMatchCandidates(begin, end);
    }
    return false;   // we don't fail anymore (used to fail on repeats)
}

//
// Build the list of match candidates for a read.
// Each candidate is defined by the arguments necessary
// to call getQValue(), so that we can
// later evaluate how good the match pair is.
//
bool MapperPE::populateMatchCandidates(
                                       ReadIndexer         &indexer,
                                       int                 whichWord,
                                       int                 candidateCount,
                                       genomeIndex_t       *candidates
                                       )
{
    std::pair<MatchCandidatesPointers_t::iterator, MatchCandidatesPointers_t::iterator> indexValue;
    MatchedReadPE matchCandidate;

    //
    // indexValue is a std::pair that contains two iterators
    // pointing to the first element of the range of sorted
    // genome postions, and the second points one beyond the
    // end.
    //
    indexValue.first = matchCandidatesPointers.end();
    indexValue.second = matchCandidatesPointers.end();

    matchCandidate.whichWord = whichWord;
    matchCandidate.quality = MatchedReadBase::UNSET_QUALITY;
    matchCandidate.indexer = &indexer;

    for (int i=0; i<candidateCount; i++)
    {
        matchCandidate.genomeMatchPosition = candidates[i] - indexer.wordPositions[whichWord];

        // genomeMatchPosition will be sometimes set to < 0 (== very large unsigned)
        // when candidates[i] is near zero, so check and ignore here:
        if (matchCandidate.genomeMatchPosition > gs->sequenceLength()) continue;

        //
        // for this read, see if we already checked this genome match position.
        // This happens often, for example when the map is perfect, and both
        // words match, so the same location will come up twice.  We only want
        // to evaluate it once.
        //
        // Note that if we find we've been here, we don't need to continue the loop,
        // since all locations have been checked already.
        //
        // This same logic is duplicated in MapperPEBaseSpace::mapPairedReads
        //

        if (indexer.checkedPositions.Find(matchCandidate.genomeMatchPosition) != -1) continue;
        indexer.checkedPositions.Add(matchCandidate.genomeMatchPosition);

        matchCandidates.push_back(matchCandidate);
        matchCandidatesPointers.push_back(&matchCandidates.back());
        indexValue.second++;
    }
    if (indexValue.first != indexValue.second) matchCandidatesIndex.push_back(indexValue);
    return false;
}

void MapperPE::checkHighMismatchMapping(MapperPE* otherMapper)
{
    bool properFlag = true;
    properFlag &= !this->clearHighMismatchMapping();
    properFlag &= !otherMapper->clearHighMismatchMapping();
    if (!properFlag)
    {
        this->isProperAligned
            = this->mapperSE->isProperAligned
            = otherMapper->isProperAligned
            = otherMapper->mapperSE->isProperAligned
            = false;
    }
    else
    {
        this->isProperAligned
            = this->mapperSE->isProperAligned
            = otherMapper->isProperAligned
            = otherMapper->mapperSE->isProperAligned
            = true;
    }
}

/**
 * If both mapper have results (no matter in PE_MODE or LOCAL_MODE), we will reset it when there are too many mismatches
 * NOTE: we do NOT count mismatches for 'N', '5' bases or bases with poorer quality than '#'
 */
bool MapperPE::clearHighMismatchMapping()
{
    int mismatchCount;
    if (bestMatch.qualityIsValid())
    {
        //
        // here we allow for a different mismatch count where we ignore
        // 'N' or '.' characters.
        //
        if (mappingMethod == LOCAL_MODE)
        {
            mismatchCount = mapperSE->bestMatch.mismatchCount;
            if (mismatchCount > forward.mismatchCutoff)
            {
                clearBestMatch();
                return true;
            }
        }
        else if (mappingMethod == PE_MODE)
        {
            if (!this->gs->isColorSpace())
                mismatchCount =
                    bestMatch.indexer -> getMismatchCount(bestMatch.genomeMatchPosition, 'N', '#');
            else
                mismatchCount =
                    bestMatch.indexer -> getMismatchCount(bestMatch.genomeMatchPosition, '5', '#');
            if (mismatchCount > forward.mismatchCutoff)
            {
                clearBestMatch();
                return true;
            }
        }
        else   // SE_MODE
        {
            mismatchCount =
                mapperSE->bestMatch.indexer -> getMismatchCount(bestMatch.genomeMatchPosition, 'N');
            if (mismatchCount > mapperSE->forward.mismatchCutoff)
            {
                clearBestMatch();
                return true;
            }
        }
    }
    return false;
}

void MapperPE::remapSingle(void)
{
    this->mapperSE->resetMapper();
    this->mapperSE->localGappedAlignment = false;
    if (forward.read.size()!=0 && mapperSE->processReadAndQuality(fragmentTag, originalRead, originalQuality)==0)
    {
        mapperSE->MapSingleRead();
    }

    // after running a single end alignment, we need to track
    // which direction it aligned with and mark that correctly
    // in our bestMatch
    if (mapperSE->getBestMatch().indexer == &mapperSE->forward)
    {
        this->bestMatch.indexer = & (this->forward);
    }
    else
    {
        this->bestMatch.indexer = & (this->backward);
    }
}

void MapperPE::setMappingMethodToPE(void)
{
    this->mappingMethod = PE_MODE;
}

void MapperPE::setMappingMethodToSE(void)
{
    this->mappingMethod = SE_MODE;
}

void MapperPE::setMappingMethodToLocal(void)
{
    this->mappingMethod = LOCAL_MODE;
}

//////////////////////////////////////////////////////////////////////
// OBSOLETE CODE
//////////////////////////////////////////////////////////////////////

#ifdef COMPILE_OBSOLETE_CODE
MatchedReadBase* MapperPE::chooseBestMatch()
{
    if (mappingMethod == PE_MODE)
    {
        return &(this->bestMatch);
    }
    else
    {
        // for mm == SE_MODE or mm == LOCAL_MODE, the mapping
        // result are both stored in mapperSE class, so just return it.
        return &(this->mapperSE->bestMatch);
    }
}
//
// XXX TODO: refactor this code.
//
// Attempt to do a local alignment.
//
// The goal appears to be to try local alignment so that the caller may
// determine whether or not it wishes to use the result.
//
bool MapperPE::tryLocalAlign(MapperPE* anchor)
{
    int64_t localAlignmentSize;
    int64_t localAlignmentPosition;
    genomeIndex_t   optimalLocalAlignment;

    // forward.read is the quality trimmed portion, and could be of zero length,
    // so check here:
    if (forward.read.size()==0) return false;

    this->mapperSE->resetMapper();
    this->mapperSE->localGappedAlignment = false; //may delete it!

    // here mapperSE serves as localMapper
    // XXX refactor this

    if (mapperSE->processReadAndQuality(fragmentTag, originalRead, originalQuality))
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
    //

    localAlignmentSize = mapperOptions.genomePositionFilterWidth +
        anchor->forward.readLength +
        mapperSE->forward.readLength;

    localAlignmentPosition = (int64_t) anchor->mapperSE->bestMatch.genomeMatchPosition - (int64_t) localAlignmentSize;

    localAlignmentSize *= 2;

    //
    // the following two checks are done using
    // 64 bit arithmetic...
    //
    // check for reads mapped before genome position localAlignmentSize/2:
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

    // XXX implement these or get rid of them:
    int mismatchCount = 0, quality = 0, gapOpenCount = 0, gapCloseCount = 0, gapExtensionCount = 0;

    // TODO(zhanxw):
    // XXX fragile code that is Illumina specific - this won't be correct
    // for colorspace.
    //
    if (anchor->mapperSE->bestMatch.isForward())
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
    }
    else
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
        this->mapperSE->bestMatch.quality = quality;
        if (!anchor->mapperSE->bestMatch.isForward())
            this->mapperSE->bestMatch.indexer = &(this->mapperSE->forward);
        else
            this->mapperSE->bestMatch.indexer = &(this->mapperSE->backward);
        this->mapperSE->localGappedAlignment = true;   // this determines where we get the CIGAR string from
        return true;
    }
}

#endif
