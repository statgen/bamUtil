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

#ifndef __READS_PROCESSOR_H__
#define __READS_PROCESSOR_H__

#include "Random.h"
#include "WordIndex.h"
#include "GenomeSequence.h"
#include "SamHeader.h"
#include "UserOptions.h"        // for MapperUserOptions

#include <vector>

// hack:
//#include "ReadIndexer.h"
#include "MapperBase.h"
#include "MapperSE.h"
#include "MapperSEBaseSpace.h"
#include "MapperPE.h"
#include "MapperPEBaseSpace.h"
#include "MappingStats.h"

using std::vector;

class ReadsProcessor
{
    GenomeSequence *gs;
    GenomeSequence *csgs;
    WordIndex *wi;
    WordHash  *wordHashLeft;
    WordHash  *wordHashRight;

    std::string outputBaseFilename;
    bool isColorSpace;
public:
    ReadsProcessor();
    ~ReadsProcessor();

    MapperUserOptions    mapperOptions;    // gets passed onto the mappers it uses
    uint64_t maxBases;
    uint64_t maxReads;
    uint64_t maxTotalReads;

    void setGenomeSequence(GenomeSequence *g)
    {
        gs = g;
    }
    void setColorSpaceGenomeSequence(GenomeSequence *g)
    {
        csgs = g;
    }
    void setWordIndex(WordIndex *w)
    {
        wi = w;
    }
    void setWordHashLeft(WordHash *w)
    {
        wordHashLeft = w;
    }
    void setWordHashRight(WordHash *w)
    {
        wordHashRight = w;
    }
    void setOutputBaseFilename(const char *s)
    {
        outputBaseFilename = s;
    }
    void setColorSpace(bool t)
    {
        isColorSpace = t;
    }

private:
    ///
    /// called by createSEMapper and createPEMapper to
    /// do common asserts on necessary word index, and
    /// left and right hashes.
    //
    void verifyHashesExist();
public:
    ///
    /// using the word index, left and righ hashes, etc,
    /// construct and return a single end mapper:
    ///
    /// @return a usable pointer to MapperSE class
    ///
    MapperSE* createSEMapper();
    //
    /// @return a usable pointer to MapperPE class
    ///
    MapperPE* createPEMapper();
    void MapPEReadsFromFiles(
        std::string filename1,
        std::string filename2,
        std::string of
    );

    void MapSEReadsFromFile(
        std::string filename,
        std::string of
    );

    void CalibratePairedReadsFiles(
        std::string filename1,
        std::string filename2
    );

    SamHeader header;

};

#endif
