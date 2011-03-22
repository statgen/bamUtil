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

#ifndef KARMA_TEST_H
#define KARMA_TEST_H

#include "GenomeSequence.h"
#include "WordHash.h"
#include "WordIndex.h"

class Test
{
    GenomeSequence  *gs;
    WordIndex       *wordIndex;
    WordHash        *wordHashLeft;
    WordHash        *wordHashRight;
public:
    void setGenomeSequence(GenomeSequence *g)
    {
        gs = g;
    }
    void setWordIndex(WordIndex *w)
    {
        wordIndex = w;
    }
    void setWordHashLeft(WordHash *w)
    {
        wordHashLeft = w;
    }
    void setWordHashRight(WordHash *w)
    {
        wordHashRight = w;
    }

    int test();

    int testSEBaseSpaceReads();
    int testSEColorSpaceReads();
    int testPEBaseSpaceReads();
    int testRemapReference(
        std::string &outputFile,
        std::string &chromosome,
        int readBaseCount,
        int skipBases
    );
};

#endif
