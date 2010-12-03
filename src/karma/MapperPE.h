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

#ifndef _MAPPERPE_H
#define _MAPPERPE_H

#include "GenomeSequence.h"
#include "MapperBase.h"
#include "MemoryMapArray.h"
#include "WordIndex.h"
#include "MapperSE.h"

#include <iostream>
#include <vector>

//
// these are actually match candidates, and serve three
// purposes:
// 1) running bestMatch
// 2) working currentMatch
// 3) part of a list of possible matches, which we
//    evaluate the Q score for in a lazy fashion.
//    In this instance, we need all args to pass to wordIndex::getQvalue
//
class MatchedReadPE : public MatchedReadBase
{
public:
    int             pairQuality;
    inline bool pairQualityIsValid() const
    {
        return pairQuality >=0;
    }
    double          pairCumulativePosteriorProbabilities;

    MatchedReadPE()
    {
        constructorClear();
    }
    void constructorClear();
    void printOptionalTags(std::ostream &, bool isProperAligned, const std::string &sampleGroupID, const std::string &alignmentPathTag);
    void updateMatch(MatchedReadPE &);
    double getQualityScore();
};

//
// Mapper algorithms for Paired End reads (PE)
//
class MapperPE : public MapperBase
{
public:
    typedef vector<MatchedReadPE>  matchCandidates_t;
    typedef vector<MatchedReadPE *>  matchCandidatesPointers_t;
    typedef vector<std::pair<
    matchCandidatesPointers_t::iterator,
    matchCandidatesPointers_t::iterator> >  matchCandidatesIndex_t;

    MapperPE    *otherEndMapper;
    MapperPE();

    void init(std::string & readFragment, std::string &dataQuality, std::string &fragmentTag);

    void populateCigarRollerAndGenomeMatchPosition();

    enum MappingMethod {PE_MODE, LOCAL_MODE, SE_MODE} mappingMethod;
    MapperSE* mapperSE;

    //
    // after all permutations of the words have been created
    // (forwardfirst, forwardsecond, backward...), count the
    // total number of matches for each direction.
    //
    void getMatchCountWithMutations();

    //
    // individual forward and backward
    // counts for getting a total mutation match count.
    //
    unsigned int forwardCount;
    unsigned int backwardCount;

    virtual void mapReads(MapperPE *)=0;
protected:
#if 0
    virtual bool mapPairedReads(
        ReadIndexer         &indexer,
        int                 whichWord,
        int                 candidateCount,
        genomeIndex_t       *candidates
    )=0;
#endif
public:
    MatchedReadPE   bestMatch;
    MatchedReadPE   matchCandidate;
    MatchedReadBase &getBestMatch();
    bool updateBestMatch(MatchedReadPE& matchCandidateB);
    //

    void debugDump(std::ostream &);
    void printBestReads(std::ostream &file, MapperPE *longerMapper);
    void printCSBestReads(std::ostream &file,
                          GenomeSequence* gs,
                          GenomeSequence* csgs,
                          MapperPE *longerMapper);

    // members and methods for managing and rapidly ordering a large
    // number of match candidates for a read.
    //
    // matchCandidates is a vector that is grown by sets of additional
    // match candidates obtained from either the pre-ordered word pointers
    // table (Word.h) or pre-ordered hashed pointers (WordHash.h).  These
    // are candidate positions with additional information necessary for
    // keeping track of mapped locations later (e.g. during printing).
    //
    matchCandidates_t           matchCandidates;
    //
    // matchCandidatesPointers is a vector containing pointers into
    // the above matchCandidates vector, and is used to make merge
    // sorting more efficient.
    //
    matchCandidatesPointers_t   matchCandidatesPointers;
    //
    // matchCandidatesIndex is a vector (initially un-ordered) of
    // pointers to begin and end subsets (which themselves are sorted),
    // and is used to do merging into progressively larger sets until
    // at the end, there is a top level iterator that points to the
    // entire ordered set.
    //
    matchCandidatesIndex_t      matchCandidatesIndex;


    bool populateMatchCandidates(bool isColorSpace);
    bool populateMatchCandidates(
        ReadIndexer &indexer,
        int whichWord,
        int candidateCount,
        genomeIndex_t *candidates
    );
    matchCandidatesIndex_t::iterator mergeSortedMatchCandidates(matchCandidatesIndex_t::iterator, matchCandidatesIndex_t::iterator);

    void testMatchCandidates();
    void testMatchCandidates(matchCandidatesIndex_t::iterator);
    void printMatchCandidates(matchCandidatesIndex_t::iterator);

    void printMatchCandidates();

    void clearBestMatch()
    {
        bestMatch.constructorClear();
        bestMatch.indexer = &forward;   // needs to point to something sane
        bestMatch.quality = MatchedReadBase::UNSET_QUALITY;
        bestMatch.pairQuality = MatchedReadBase::UNSET_QUALITY;
        bestMatch.genomeMatchPosition = INVALID_GENOME_INDEX;

        mapperSE->bestMatch.constructorClear();
        mapperSE->bestMatch.indexer = &forward;   // needs to point to something sane
        mapperSE->bestMatch.quality = MatchedReadBase::UNSET_QUALITY;
        mapperSE->bestMatch.genomeMatchPosition = INVALID_GENOME_INDEX;
    }
    void checkHighMismatchMapping(MapperPE* otherMapper);
    bool clearHighMismatchMapping();
    void remapSingle(void);
    virtual bool tryLocalAlign(MapperPE* anchor) = 0;
    void setMappingMethodToPE(void);
    void setMappingMethodToSE(void);
    void setMappingMethodToLocal(void);

    MatchedReadBase* chooseBestMatch();

};

//
// helper for inplace_merge algorithm used below
// in mergeSortedMatchCandidates:
//
// This is used because we don't want to merge the matchCandidatesPointers
// based on the pointer value, but rather the genome positiobn of the
// values they point to.
//
// The goal is to apply inplace_merge on pointers rather than the
// match candidtate structures, which are much larger than a single pointer.
//
inline static bool MatchedReadPEComp(MatchedReadPE *v1, MatchedReadPE *v2)
{
    return v1->genomeMatchPosition < v2->genomeMatchPosition;
}


#endif
