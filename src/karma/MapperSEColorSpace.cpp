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
#include "MapperSEColorSpace.h"
#include "MappingStats.h"
#include "ReadsProcessor.h"
#include "MathConstant.h"
#include "Performance.h"
#include "Util.h"

#include <algorithm>
#include <iostream>
#include <vector>

/////////////////////////////////////////////////////////////////////
// Below I copied all functions from MapperSEBaseSpace class
// and rename it to MapperSEColorSpace
// meanwhile, replace getSumQ() and getSumQOrig() with suitable
// getColorSpaceSumQ() functions
/////////////////////////////////////////////////////////////////////
MapperSEColorSpace::MapperSEColorSpace()
{
}

//
// Wrapper routine for mapping a read.
//
// It's a wrapper because apart from ReadIndex::getSumQ calling
// ReadIndex::getSumQSW when we want gapped alignment, it is currently otherwise
// identical.
//
// However, we generally want to minimize calls to ReadIndex::getSumQSW(),
// we allow an option for calling it only if the ungapped mapper has
// failed in some horrible fashion.
//
void MapperSEColorSpace::MapSingleRead()
{
    if (mapperOptions.forceSmithWaterman)
    {
        //
        // fortune: you are a glutton for severe punishment,
        // you enjoy pain, you have a very large compute cluster
        // with nearly unlimited time and electrical power.
        //
        MapSingleReadGapped();
    }
    else
    {

        MapSingleReadUnGapped();

        //
        // If user or the technology (e.g. 454) suggests we can allow
        // gapped alignment, here we'll first check if ungapped alignment
        // failed badly enough that we wish to run it.
        //
        // XXX we'll probably want an option for a minimum quality threshhold,
        // provided that we trust map qualities enough
        //
        if (
            mapperOptions.allowSmithWaterman &&
            (bestMatch.quality == MatchedReadBase::REPEAT_QUALITY ||
             bestMatch.quality == MatchedReadBase::REPEAT_QUALITY)
            )
        {
            MapSingleReadGapped();
        }
    }
}

//
// This is horrible, but there is no way I'm going to duplicate all the
// iteration code in MapperSE::MapSingleReadUnGapped just so we can
// allow gapped sumQ.
//
// XXX figure out a way to replace this vile code.  Please!
//
void MapperSEColorSpace::MapSingleReadGapped()
{

    this->resetMapper();
    forward.useGapped = true;
    backward.useGapped = true;
    MapSingleReadUnGapped();
    bestMatch.gappedAlignment = true;
    forward.useGapped = false;
    backward.useGapped = false;
}

bool evalTrampolineEvalSEColorSpace(
                           MapperBase *mapper,
                           ReadIndexer &indexer,
                           genomeIndex_t genomeMatchPosition,
                           int whichWord)
{
    return ((MapperSEColorSpace*)mapper)->evalSEColorSpace(indexer, genomeMatchPosition, whichWord);
}

void MapperSEColorSpace::MapSingleReadUnGapped()
{
    this->resetMapper();
    cumulativePosteriorProbabilitiesCutoff = 100000;

    clearBestMatch();

    forward.setMismatchCutoff();
    backward.setMismatchCutoff();

    //evalBaseSpaceReads(evalTrampoline, NULL);
    evalColorSpaceReads(evalTrampolineEvalSEColorSpace, NULL);

    return;

}

inline bool MapperSEColorSpace::evalSEColorSpace(
                                                 ReadIndexer &indexer,
                                                 genomeIndex_t genomeMatchPosition,
                                                 unsigned int whichWord)
{
    int quality;
    int mismatchCount=-1; // in getSumQ(), this value should be calculated to be non-negative

    //
    // whichWord is used only by the Smith Waterman gapped
    // code, which is only called when indexer.useGapped is true.
    //
    // This is ugly, but gets the job done in the cleanest
    // and most expediate manner at the moment.
    //
    // Note, karma spends the most time by far inside this call,
    // even for exact match only, non-gapped alignments.
    //
    quality = indexer.getColorSpaceSumQ(
                                        genomeMatchPosition,
                                        mismatchCount,
                                        bestMatch.quality,
                                        whichWord
                                        );

    updateBestMatch(indexer, quality, mismatchCount, whichWord, genomeMatchPosition);

    return false;
}

int MapperSEColorSpace::test(int testNum, const char *read, const char *qual, char direction, int chr, genomeIndex_t index, int misMatches, int quality)
{
    setReadAndQuality(read, strlen(read), qual);
    fragmentTag = "";

    MapSingleRead();

    int rc = 0;

    check(rc, testNum, "direction", direction, getBestMatch().isForward() ? 'F' : 'R');
    check(rc, testNum, "mismatch", misMatches, getBestMatch().mismatchCount);
    check(rc, testNum, "quality", quality, getBestMatch().quality);

    int c = gs->getChromosome(this->getBestMatch().genomeMatchPosition);
    genomeIndex_t i = getBestMatch().genomeMatchPosition - gs->getChromosomeStart(c) + 1;
    check(rc, testNum, "chromosome", chr-1, c);
    check(rc, testNum, "index", index, i);

#if 0
    std::cout<< "read: "<<read << std::endl;
    std::cout<< "ref.: ";
    for (unsigned int i=0; i<  strlen(read); i++)
    {
        std::cout<< char(GenomeSequence::base2int[int((*gs)[getBestMatch().genomeMatchPosition+i])]+'0');
    }
    std::cout<< std::endl;
    std::cout<<"genomeMatchPosition: "<<getBestMatch().genomeMatchPosition<<std::endl;
    std::cout<<"chr#:pos "
             << gs->getChromosome(getBestMatch().genomeMatchPosition)
             << "\t"
             << getBestMatch().genomeMatchPosition - gs->getChromosomeStart(c) + 1
             << std::endl;
    std::cout<< std::endl;

    debugPrint(getBestMatch());
#endif
    return rc;
}


