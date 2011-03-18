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
#include "MemoryMapArray.h"
#include "MapperBase.h"
#include "MapperSE.h"
#include "WordIndex.h"
#include "MatchedReadPE.h"

#include <iostream>
#include <vector>


//
// Mapper algorithms for Paired End reads (PE)
//
class MapperPE : public MapperBase
{
 public:
    // Define types
    typedef vector<MatchedReadPE>  MatchCandidates_t;
    typedef vector<MatchedReadPE *>  MatchCandidatesPointers_t;
    typedef MatchCandidatesPointers_t::iterator MatchCandidatesPointersIter_t;
    typedef vector< std::pair 
        <MatchCandidatesPointersIter_t, MatchCandidatesPointersIter_t> > MatchCandidatesIndex_t;
 public:

    MapperPE    *otherEndMapper;
    MapperPE();

    // reset method, it should be called before any alignment actually happen
    // note, this just reset one Mapper, you need to do the same for the other Mapper.
    void resetMapper() {
        this->clearBestMatch();
        this->matchCandidate.constructorClear();
        this->forward.checkedPositions.Clear();
        this->backward.checkedPositions.Clear();
    };

    void init(std::string & readFragment, std::string &dataQuality, std::string &fragmentTag);

    // calculate cigar for this PE mapper and (if SE mapper exists) the SE mapper.
    void populateCigarRollerAndGenomeMatchPosition();
    
    // for local align to single end alignment.
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

    MatchedReadBase &getBestMatch();
    bool updateBestMatch(MatchedReadPE& matchCandidateB);

    void printBestReads(std::ostream &file, MapperPE *longerMapper);
    void printCSBestReads(std::ostream &file,
                          GenomeSequence* gs,
                          GenomeSequence* csgs,
                          MapperPE *longerMapper);

    // reset bestMatch and assign its indexer to a valid ReadIndexer instance
    void clearBestMatch()
    {
        this->bestMatch.constructorClear();
        this->bestMatch.indexer = &forward;   // needs to point to something sane
        
        mapperSE->bestMatch.constructorClear();
        mapperSE->bestMatch.indexer = &forward;   // needs to point to something sane
    }

    // calling clearHighMismatchMapping() and set proper flag for every aligner
    void checkHighMismatchMapping(MapperPE* otherMapper);
    //clear bestMatch if the number of mismatches are higher than mismatch cutoff
    bool clearHighMismatchMapping();
    // try single end alignment for the read
    void remapSingle(void);
    virtual bool tryLocalAlign(MapperBase* anchor) = 0;
    // setter functions
    void setMappingMethodToPE(void);
    void setMappingMethodToSE(void);
    void setMappingMethodToLocal(void);

    enum MappingMethod {PE_MODE, LOCAL_MODE, SE_MODE} mappingMethod;
    MatchedReadPE   bestMatch;

    bool populateMatchCandidates(bool isColorSpace);
    bool populateMatchCandidates(
                                 ReadIndexer &indexer,
                                 int whichWord,
                                 int candidateCount,
                                 genomeIndex_t *candidates
                                 );
    MatchedReadPE   matchCandidate;
 protected:

    void debugDump(std::ostream &);

    // members and methods for managing and rapidly ordering a large
    // number of match candidates for a read.
    //
    // matchCandidates is a vector that is grown by sets of additional
    // match candidates obtained from either the pre-ordered word pointers
    // table (Word.h) or pre-ordered hashed pointers (WordHash.h).  These
    // are candidate positions with additional information necessary for
    // keeping track of mapped locations later (e.g. during printing).
    //
    MatchCandidates_t           matchCandidates;
 public:
    //
    // matchCandidatesPointers is a vector containing pointers into
    // the above matchCandidates vector, and is used to make merge
    // sorting more efficient.
    //
    MatchCandidatesPointers_t   matchCandidatesPointers;
 protected:
    //
    // matchCandidatesIndex is a vector (initially un-ordered) of
    // pointers to begin and end subsets (which themselves are sorted),
    // and is used to do merging into progressively larger sets until
    // at the end, there is a top level iterator that points to the
    // entire ordered set.
    //
    MatchCandidatesIndex_t      matchCandidatesIndex;

    // 
    MatchCandidatesIndex_t::iterator mergeSortedMatchCandidates(MatchCandidatesIndex_t::iterator, 
                                                                MatchCandidatesIndex_t::iterator);

    void testMatchCandidates();
    void testMatchCandidates(MatchCandidatesIndex_t::iterator);
    void printMatchCandidates(MatchCandidatesIndex_t::iterator);
    void printMatchCandidates();

#ifdef COMPILE_OBSOLETE_CODE
    // according to this->mappingMethod, return MatchedRead aligned in SE mode or PE mode
    MatchedReadBase* chooseBestMatch();
 protected:
    virtual bool mapPairedReads(
                                ReadIndexer         &indexer,
                                int                 whichWord,
                                int                 candidateCount,
                                genomeIndex_t       *candidates
                                )=0;
#endif

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
