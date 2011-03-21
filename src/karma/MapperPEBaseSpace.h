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

#ifndef _MAPPERPEBASESPACE_H_
#define _MAPPERPEBASESPACE_H_

#include "GenomeSequence.h"
#include "MapperBase.h"
#include "MemoryMapArray.h"
#include "WordIndex.h"

#include <iostream>
#include <vector>

//
// Mapper algorithms for Paired End reads (PE) in Base Space
//
class MapperPEBaseSpace : public MapperPE
{
 public:
    MapperPEBaseSpace();
    ~MapperPEBaseSpace();

    void mapReads(MapperPE *);
    void mapReads2(MapperPE *);
    bool mapReads(
                  ReadIndexer         &indexer,
                  int                 whichWord,
                  int                 candidateCount,
                  genomeIndex_t       *candidates
                  );
    bool tryLocalAlign(MapperBase* anchor);

 protected:
    void testMatchCandidates();
    void testMatchCandidates(MatchCandidatesIndex_t::iterator);
    void printMatchCandidates(MatchCandidatesIndex_t::iterator);

    void printMatchCandidates();

  private:
    MatchedReadPE   compareHelper;


 public:
    int test(
             int testNum,
             MapperPE    *otherMapper,
             const char *read1, const char *qual1, char direction1, int chr1, genomeIndex_t index1, int misMatches1,
             const char *read2, const char *qual2, char direction2, int chr2, genomeIndex_t index2, int misMatches2
             );

};


#endif /* _MAPPERPEBASESPACE_H_ */
