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

#ifndef _MAPPERBASE_H
#define _MAPPERBASE_H

#include "GenomeSequence.h"
#include "MemoryMapArray.h"
#include "MapperUserOption.h"
#include "Random.h"
#include "ReadIndexer.h"
#include "SmithWaterman.h"  // for CigarRoller
#include "WordIndex.h"
#include "MatchedReadBase.h"

#include <iostream>
#include <vector>
#include <algorithm>



// XXX this needs to go away.
#define    Q_CUTOFF (1.0 / 999.0)

class MapperBase
{

protected:
    WordIndex   *wordIndex;
    WordHash    *wordHashRight;
    WordHash    *wordHashLeft;

public:
    static const int MAX_Q2P = 1000;    // size of sumQualityToProb array

    double cumulativePosteriorProbabilitiesCutoff;

    // array of pre-calculated pow(10, - 0.1 * sumQuality) for sumQuality in [0, 999)
    double sumQualityToProb[MAX_Q2P];

    uint16_t    samMateFlag;
    std::string alignmentPathTag; // e.g. Set its value to "ABC" will result a SAM optional tag: XA:Z:alignmentPathTag

protected:
    GenomeSequence  *gs;
    Random *rand;
    MapperUserOption    mapperOptions;

// The following table is constructed to help detecting color space SNPs
// e.g.
// for index 0, we know its binary 0000, and possible SNPs are 0101, 1010, 1111, so we put 5, 10, 15
// for index 1, we know its binary 0001, and possible SNPs are 0100, 1011, 1110, so we put 4, 11, 14
    static unsigned int colorSpaceSNP[][3];

public:
    MapperBase();
    virtual ~MapperBase();

    int Word2Integer(std::string & word, unsigned int index, int &countNbases);
    std::string Integer2Word(wordInteger_t n, unsigned int wordsize);

    std::string fragmentTag;
    std::string originalRead;
    std::string originalQuality;
    int line;
    int numberReads;

    int getReadAndQuality(IFILE f);
    int processReadAndQuality(std::string& fragmentTag, std::string& readFragment, std::string& dataQuality);  // 0->success
    bool setReadAndQuality(const char *read, int len, const char *quality); // true -> failed

private:
    void setGenomeSequence(GenomeSequence *g)
    {
        gs = g;
        forward.gs = gs;
        backward.gs = gs;
    }
    void setWordIndex(WordIndex *w)
    {
        wordIndex = w;
        forward.wordIndex = w;
        backward.wordIndex = w;
    }

    void setWordHashLeft(WordHash *w)
    {
        wordHashLeft = w;
    }
    void setWordHashRight(WordHash *w)
    {
        wordHashRight = w;
    }
    void setMapperOptions(MapperUserOption &m)
    {
        //
        // NB forward::mapperOptions and backward::mapperOptions
        // are references to this member:
        //
        mapperOptions = m;      // copy values into here
    }

public:
    void initMapper(GenomeSequence *g, WordIndex *w, WordHash *wl, WordHash * wr, MapperUserOption &m)
    {
        this->setGenomeSequence(g);
        this->setWordIndex(w);
        this->setWordHashLeft(wl);
        this->setWordHashRight(wr);
        this->setMapperOptions(m);
    }

    ReadIndexer     forward;
    ReadIndexer     backward;

    virtual MatchedReadBase &getBestMatch() = 0;

    void    debugPrint(MatchedReadBase &);

    int getReadLength()
    {
        return forward.readLength;
    }

    void restoreTrimming(CigarRoller&   cigarRoller,
                         ReadIndexer&   indexer,
                         genomeIndex_t& genomeMatchPosition);

    //
    // Cigar string handling is tricky.  I'm uncomfortable
    // putting more complex data in MatchedReadBase because
    // it is a highly repeated and highly used structure in
    // the mapping process, and I doubt I can currently correctly
    // predict the cost of doing so.
    //
    // Because they are intended to be fast, SE and PE aligners
    // don't bother keeping track of the CigarRoller object that
    // is generated via the
    //
    //
    // if the bestMatch was gapped, we have lots
    // of work to do... see MapperBase.cpp.  This
    // method is here in MapperBase because getting
    // the gapped alignment Cigar string depends on
    // knowing a lot more about the read - all the
    // stuff in forward and backward, for example.
    //
    // virtual because PE mapper has additional perversions
    // to do because it now contains a nested SE mapper.
    //
    virtual void populateCigarRollerAndGenomeMatchPosition();

    //
    // When the cigarRoller object is valid (because of doing
    // local gapped alignment on this mate), then we need
    // to set localGappedAlignment to true.  This affects
    // operation of populateCigarRollerAndGenomeMatchPosition()
    // above.
    //
    bool        localGappedAlignment;
    bool        isProperAligned;
    CigarRoller cigarRoller;    // this is what printReads() uses to output the cigar string

    std::string originalCSRead;
    std::string originalCSQual;

    typedef bool (*evalSinglePositionFunctionType)(MapperBase *, ReadIndexer &, genomeIndex_t, int whichWord);
    typedef bool (*evalAllPositionsFunctionType)(MapperBase *, ReadIndexer &, int candidateCount, genomeIndex_t *candidates, int whichWord);

    void evalBaseSpaceReads(evalSinglePositionFunctionType, evalAllPositionsFunctionType);
    void evalColorSpaceReads(evalSinglePositionFunctionType, evalAllPositionsFunctionType);

private:
    bool evalBaseSpaceRead(
        evalSinglePositionFunctionType,
        evalAllPositionsFunctionType,
        ReadIndexer &indexer);
    bool evalColorSpaceRead(
        evalSinglePositionFunctionType,
        evalAllPositionsFunctionType,
        ReadIndexer &indexer);

    bool evalAllCandidatesForWord(
        evalSinglePositionFunctionType,
        evalAllPositionsFunctionType,
        ReadIndexer &indexer,
        unsigned int whichWord,
        wordInteger_t xorMask);

    bool evalAllCandidatesForColorSpaceWord(
        evalSinglePositionFunctionType,
        evalAllPositionsFunctionType,
        ReadIndexer &indexer,
        unsigned int whichWord,
        wordInteger_t shiftLocation,
        wordInteger_t mutationIndex);

    bool evalAllCandidatePositions(
        evalSinglePositionFunctionType,
        ReadIndexer &indexer,
        int     whichWord,
        int     candidateCount,
        genomeIndex_t *candidates
    );
public:
    //
    // Debug aid.  We'll sprinkle the mapper with some checkpoints
    // of interest - specifically for each condition that arises of
    // note, and conditionally print debug info.
    //
    // Performance impact: if enableMapperDebug==false and is a constant,
    // the inline code will be omitted.
    //
    static const bool enableMapperDebug = false;
    static const genomeIndex_t debugPositionWindow = 5000;
    genomeIndex_t   debugPosition;
    void debugCheckPoint(genomeIndex_t pos, const char *str)
    {
        if (enableMapperDebug)
        {
            if ((debugPosition - debugPositionWindow < pos) &&
                    (pos < debugPosition + debugPositionWindow))
            {
                std::string location;

                gs->getChromosomeAndIndex(location, pos);

                std::cout << str << ": genome index " << pos << " (" << location << ")" << std::endl;
            }
        }
    }
};


#endif
