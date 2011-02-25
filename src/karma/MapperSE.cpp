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

#include "MapperSE.h"
#include "MappingStats.h"
#include "ReadsProcessor.h"
#include "Error.h"
#include "MathConstant.h"
#include "Performance.h"
#include "Util.h"

#include <algorithm>
#include <iostream>
#include <vector>



MatchedReadBase &MapperSE::getBestMatch()
{
    return bestMatch;
}

MapperSE::MapperSE()
{
}

MapperSE::~MapperSE()
{
}

// this function has 4 steps
// 1. adjust for invalid quality values ( <0, or cap to MAX_Q2P)
// 2. check mismatchCutoff
// 3. update bestMatch (randomly if tie happens)
// 4. adjust for huge posterior probability.
// XXX UPDATE whichWord as well, so getCigarStringOffset() can work
bool MapperSE::updateBestMatch(
    ReadIndexer& indexer,
    int quality,
    int mismatchCount,
    int whichWord,
    genomeIndex_t genomeMatchPosition)
{
#if 0
    std::cout<<"genomeMatchPosition="<<genomeMatchPosition<<" quality="<<quality<<std::endl;
#endif
    //
    // getQvalue returns -1 when it decides the read is no
    // longer worth computing the quality score for.  Here we
    // check for that, then drop the word if that is the case.
    //
    if (quality < 0) return false;

    if (mismatchCount > indexer.mismatchCutoff) return false;

    // topcode candidate quality:
    if (quality >= MAX_Q2P)
    {
        quality = MAX_Q2P - 1;
    }

    bestMatch.numMatchContributors ++;

    if (bestMatch.qualityIsUnset() || quality < bestMatch.quality)
    {
        bestMatch.numBest = 1;
        bestMatch.indexer = &indexer;
        bestMatch.whichWord = whichWord;
        bestMatch.genomeMatchPosition = genomeMatchPosition;
        bestMatch.quality = quality;
        bestMatch.mismatchCount = mismatchCount;
    }
    else if (quality == bestMatch.quality)
    {
        bestMatch.numBest ++;

        //
        // need to paramaterize this, but the basic idea should be sound...
        // After some number of duplicate positions all at top quality,
        // there really isn't any more one can do...
        //
        if (bestMatch.numBest > 5 && quality == 0)
        {
            bestMatch.genomeMatchPosition = 0;
            bestMatch.quality = MatchedReadBase::EARLYSTOP_DUPLICATES;
            return true;
        }
        // pick existing one w/ prob (bestMatch.numBest - 1) / bestMatch.numBest;
        // pick the crt one w/ prob 1 / bestMatch.numBest;
        double random = rand->Uniform();
        if (random < (1.0 / (double) bestMatch.numBest))
        {
            bestMatch.indexer = &indexer;
            bestMatch.genomeMatchPosition = genomeMatchPosition;
            bestMatch.quality = quality;
            bestMatch.mismatchCount = mismatchCount;
        }
    }

    // check early stop from reaching posterior probability.
    if (
        bestMatch.cumulativePosteriorProbabilities > cumulativePosteriorProbabilitiesCutoff &&
        sumQualityToProb[quality] > cumulativePosteriorProbabilitiesCutoff
    )
    {
#if 0
        // debug:
        double qs = bestMatch.getQualityScore(bestMatch.quality, bestMatch.cumulativePosteriorProbabilities);
        std::cout << "just rejected a read wtih qs = " << qs << std::endl;
#endif
        bestMatch.genomeMatchPosition = 0;
        bestMatch.quality = MatchedReadBase::EARLYSTOP_QUALITY;
        return true;
    }

    bestMatch.cumulativePosteriorProbabilities += sumQualityToProb[quality];

    return false;
}

void MapperSE::clearBestMatch()
{
    // re-init best match - we're starting over...
    bestMatch.constructorClear();
    // normally, mapper should only point to a map if we found one, but
    // some information still needs to be printed even if no match
    // positions were found.
    bestMatch.indexer = &forward;
    // only update the bestMatch if it has a better Q than this:
    bestMatch.quality = MatchedReadBase::UNSET_QUALITY;
    bestMatch.genomeMatchPosition = INVALID_GENOME_INDEX;
}
