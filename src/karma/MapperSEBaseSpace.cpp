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

#include "Generic.h"
#include "MapperSE.h"
#include "MapperSEBaseSpace.h"
#include "MappingStats.h"
#include "ReadsProcessor.h"
#include "Error.h"
#include "MathConstant.h"
#include "Performance.h"
#include "Util.h"

#include <algorithm>
#include <iostream>
#include <vector>


MapperSEBaseSpace::MapperSEBaseSpace()
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
void MapperSEBaseSpace::MapSingleRead()
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
             bestMatch.quality == MatchedReadBase::UNSET_QUALITY)
        )
        {
            //
            // bleh... gotta clean up some stuff in this mapper
            // so that we can re-run against the same read:
            //
            bestMatch.constructorClear();
            forward.checkedPositions.Clear();
            backward.checkedPositions.Clear();

            // now rerun, but do gapped alignment (slow)
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
void MapperSEBaseSpace::MapSingleReadGapped()
{
    forward.useGapped = true;
    backward.useGapped = true;
    MapSingleReadUnGapped();    // ugly -- now MapSingleReadUnGapped acts like a gapped mapper
    bestMatch.gappedAlignment = true;
    forward.useGapped = false;
    backward.useGapped = false;
}

//
// this helper function is because getting pointer to method to work the way I want
// it to is a big pain in the backside.
//
static bool evalTrampoline(
    MapperBase *mapper,
    ReadIndexer &indexer,
    genomeIndex_t genomeMatchPosition,
    int whichWord)
{
    return ((MapperSEBaseSpace*)mapper)->evalSEBaseSpace(indexer, genomeMatchPosition, whichWord);
}

void MapperSEBaseSpace::MapSingleReadUnGapped()
{
#if 0
    if (mapperOptions.minimumMapQuality>0)
//        cumulativePosteriorProbabilitiesCutoff = 1.0 - pow(10.0, -mapperOptions.minimumMapQuality/10.0);
        cumulativePosteriorProbabilitiesCutoff = 1.0 / (1 - pow(10.0, -mapperOptions.minimumMapQuality/10.0));
    else
#endif
        cumulativePosteriorProbabilitiesCutoff = 100000;

    clearBestMatch();

    //
    forward.setMismatchCutoff();
    backward.setMismatchCutoff();

    // a chain reaction is going to happen:
    // evalBaseSpaceReads(evalTrampoline, NULL)                                =>
    // evalBaseSpaceRead(evalTrampoline, NULL, forward/backward)               =>
    // evalAllCandidatesForWord(evalTrampoline, NULL, indxer, whichWord, 0)    =>
    // (MapperSEBaseSpace*)mapper->evalSEBaseSpace()   (from evalTrampoline)
    //
    evalBaseSpaceReads(evalTrampoline, NULL);

}

inline
bool MapperSEBaseSpace::evalSEBaseSpace(
    ReadIndexer &indexer,
    genomeIndex_t genomeMatchPosition,
    int whichWord)
{

    int quality;
    int mismatchCount;

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
    quality = indexer.getSumQ(
                  genomeMatchPosition,
                  mismatchCount,
                  bestMatch.quality,
                  whichWord
              );

    //
    // use a base class method to update the match --- shared
    // between paired and and single ended code... I think...
    // XXX incomplete - we need to pass in whichWord to keep track
    // of it for printing CIGAR strings when we're all done finding
    // the best match.
    //
    updateBestMatch(indexer, quality, mismatchCount, whichWord, genomeMatchPosition);

#if 0
    //
    // debug code to show all possible matches... leave it here, please.
    //
    std::string s;
    gs->getChromosomeAndIndex(s, genomeMatchPosition);
    std::cout << "eval position " << genomeMatchPosition << " " << s << " this position quality = " << quality << " bestmatch quality = " << bestMatch.quality << std::endl;
#endif

    return false;
}

//
// Test the read, and compare against the expected results.
// Returns 1 if there was a problem, 0 otherwise.
//
int MapperSEBaseSpace::test(int testNum, const char *read, const char *qual, char direction, int chr, genomeIndex_t index, int misMatches, int quality, const char *expectedCigarString)
{
    const char *cigarStringPtr;
    std::string cigarString;
    std::string expectedCigarStringLocal = expectedCigarString;
    setReadAndQuality(read, strlen(read), qual);
    fragmentTag = "";

    debugPosition = gs->getGenomePosition(chr-1, index);

    MapSingleRead();

    //
    // to deal with gapped alignment, we finally get the
    // CIGAR string, which in addition gets us an offset
    // adjustment for the best match index position.
    //
    populateCigarRollerAndGenomeMatchPosition();

    cigarStringPtr = cigarRoller.getString();  // gotta delete the result
    cigarString = cigarStringPtr;
    delete cigarStringPtr;

    int rc = 0;

    check(rc, testNum, "direction", direction, bestMatch.isForward() ? 'F' : 'R');
    check(rc, testNum, "mismatch", misMatches, bestMatch.mismatchCount);
    check(rc, testNum, "quality", quality, bestMatch.quality);

    int c = gs->getChromosome(bestMatch.genomeMatchPosition);
    genomeIndex_t i = bestMatch.genomeMatchPosition - gs->getChromosomeStart(c) + 1;
    check(rc, testNum, "chromosome", chr-1, c);
    check(rc, testNum, "index", index, i);
    check(rc, testNum, "cigarString", expectedCigarStringLocal, cigarString);

    if (rc)
    {
        //
        // XXX out of date comment here - offset is no longer
        // computed - we could add this back in, but not just
        // for this one debug statement:
        //
        // when offset is < 0, that means there are deletions
        // in the read, which means we want to print more of
        // the reference to see what the read really corresponds
        // to.
        //
        std::string referenceString;
        gs->getString(
            referenceString,
            bestMatch.genomeMatchPosition,
            strlen(read) /* - offset */);
        std::cout << "Want CIGAR:" << expectedCigarStringLocal << std::endl;
        std::cout << "Got CIGAR: " << cigarString << std::endl;
        std::cout << "Cigar Operations:" << std::endl;
        std::cout << cigarRoller;
//        std::cout << "Offset:    " << offset << std::endl;
        std::cout << "Read:      " << read << std::endl;
        std::cout << "Reference: " << referenceString << std::endl;

    }
    return rc;
}
