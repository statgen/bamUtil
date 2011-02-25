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

#ifndef _MAPPERSE_H
#define _MAPPERSE_H

#include "GenomeSequence.h"
#include "MatchedReadSE.h"
#include "MapperBase.h"
#include "MemoryMapArray.h"
#include "WordIndex.h"

#include <iostream>
#include <vector>

//
// Mapper algorithms for Single End reads (SE)
//
class MapperSE : public MapperBase
{
public:
    MatchedReadSE bestMatch;

public:
    MapperSE();
    ~MapperSE();

    virtual void MapSingleRead()=0;
    // the following two are not claimed to be pure virtual,
    // since MapperSEColorSpace class has not ready to implement these 2 functions yet.
    virtual void MapSingleReadGapped()=0;
    virtual void MapSingleReadUnGapped()=0;

    void setBestMatch(MatchedReadSE& m)
    {
        this->bestMatch=m;
    };
    void clearBestMatch();

    int evaluateAllCandidates(
        ReadIndexer &indexer,
        int     whichWord,
        int     candidateCount,
        genomeIndex_t *candidates
    );

    MatchedReadBase &getBestMatch();

    //
    // for every single end match candidate, call
    // this method, which will compare the quality
    // of the overall best match with the current
    // match, and if necessary, update the best match.
    //
    // whichWord has to be passed through here because
    // this information is required when at the end
    // we wish to print the CIGAR string, since this
    // basically requires re-running the same code
    // getSumQ does (plus more to keep track of various
    // changes in the match).
    //
    bool updateBestMatch(
        ReadIndexer& indexer,
        int quality,
        int mismatchCount,
        int whichWord,
        genomeIndex_t genomeMatchPosition);
};

#endif
