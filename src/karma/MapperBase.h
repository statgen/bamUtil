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
#include "UserOptions.h"
#include "Random.h"
#include "ReadIndexer.h"
#include "SmithWaterman.h"  // for CigarRoller
#include "WordIndex.h"

#include <iostream>
#include <vector>
#include <algorithm>

// anything outside these values represents an invalid base
// base codes: 0-> A,    1-> C,     2-> G,      3-> T
// colorspace: 0-> blue, 1-> green, 2-> oragne, 3->red
//
// highest 2 bits: a base space nucleotide
// lowest 2 bits : a color space code
// Surprisingly, the following table is the same as const char fromBase2CS
static const char fromCS2Base[] =
{
    /* 0000 */ 0,   // A->A
    /* 0001 */ 1,   // A->C
    /* 0010 */ 2,   // A->G
    /* 0011 */ 3,   // A->T
    /* 0100 */ 1,   // C->C
    /* 0101 */ 0,   // C->A
    /* 0110 */ 3,   // C->T
    /* 0111 */ 2,   // C->G
    /* 1000 */ 2,   // G->G
    /* 1001 */ 3,   // G->T
    /* 1010 */ 0,   // G->A
    /* 1011 */ 1,   // G->C
    /* 1100 */ 3,   // T->T
    /* 1101 */ 2,   // T->G
    /* 1110 */ 1,   // T->C
    /* 1111 */ 0,   // T->A
};

//////////////////////////////////////////////////////////////////////
// functor class for get the compliment base
//////////////////////////////////////////////////////////////////////
struct BaseComplementor
{
    char operator()(char& c)
    {
        return GenomeSequence::base2complement[(int)c];
    };
    char operator()(char c)
    {
        return GenomeSequence::base2complement[(int)c];
    };
};

class MapperBase;
//
// When we get a match, mutated or not, keep records here
//
class MatchedReadBase
{
public:
    MatchedReadBase()
    {
        constructorClear();
    }
    virtual ~MatchedReadBase();
    void constructorClear();

    static const int UNSET_QUALITY = -1;         // default value
    static const int INVALID_DATA = -2;          // no used yet
    static const int EARLYSTOP_DUPLICATES = -3;  // early stop due to reaching max number of best matches, with all bases matched
    static const int EARLYSTOP_QUALITY = -4;     // early stop due to reaching max posterior quality
    static const int REPEAT_QUALITY = -5;        // not used yet

    // ---------------------------------
    //  Member variables
    // ---------------------------------
    //
    // whichWord is the index word
    // that was used to obtain this match.  When doing Smith Waterman
    // matching, it tells us where to anchor the forward and
    // reverse SW searches.
    //
    // When gappedAlignment is true, we got this match via using Smith
    // Waterman, which means the ::print method needs to correctly
    // print the corresponding cigar string.
    //
    int     whichWord;
    bool    gappedAlignment;    // Smith Waterman was used to align this result

    genomeIndex_t genomeMatchPosition;
    int mismatchCount;
    int quality;
    double cumulativePosteriorProbabilities;
    int numBest;
    int numMatchContributors;

    ReadIndexer *indexer;

    // if "other"  is given, include its quality as
    // a contributor
    double getQualityScore(double sumQ, double posteriorProbability);
    virtual double getQualityScore();

    bool isForward() const
    {
        return indexer?indexer->isForward:true;
    }
    virtual void updateMatch(MatchedReadBase &);


    // ---------------------------------
    //  SAM output functions:
    // ---------------------------------
    static void printHeader(std::ostream &);

    void calculateNAMEandPOS(const char* &chromosomeName, int& chromosomePosition,
                             const char* &mateChromosomeName, int& mateChromosomePosition,
                             const MatchedReadBase* mate);
    bool getBSSeq2print(std::string&    sequence2print,
                        CigarRoller&    cigarRoller,
                        bool            showReferenceBases = true);
    bool getCSSeqAndQual2print(std::string& sequence2print,
                               std::string& quality2print,
                               std::string&    cs_read_fragment,
                               std::string&    cs_data_quality,
                               GenomeSequence* gs,
                               GenomeSequence* csgs,
                               CigarRoller& cigarRoller,
                               bool         showReferenceBases = true);
    ///
    /// Write a single read to the stream.  It may
    /// have been sequenced as a pair, in which case
    /// mate points to the mate matchedReadBase.
    ///
    /// outputStream -> ostream that we send output to.
    ///
    /// mate -> if sequenced as a mate, this is the mate.
    /// Having mate!=NULL does not guarantee the mate was
    /// pair aligned (see isProperAligned).
    ///
    /// fragmentTag -> the tag for this read.  Since
    /// it is not included in the indexer class, we have
    /// to include it here.
    ///
    /// cigarRoller -> Vector of CIGAR operations indicating
    /// insert/delete/match/mismatch for the aligned read.
    /// This is populated by either the Smith Waterman aligner
    /// or by some other local alignment functions.
    ///
    /// isProperAligned -> indicates when a mate exists
    /// and is within the max insert distance specified
    /// by the command line arguments.  This affects the
    /// SAM flag bit 0x0002 and ISIZE fields.
    ///
    void print(
        std::ostream &outputStream,
        MatchedReadBase *mate,
        std::string &fragmentTag,
        bool showReferenceBases,
        CigarRoller  &cigarRoller,
        bool    isProperAligned,
        uint16_t    samMateFlag,
        const std::string &sampleGroupID,
        const std::string &alignmentPathTag
    );

    ///
    /// It is a similar procedure to print(), however, we implicitly translate
    /// color space reads to base space.
    ///
    void printColorSpace(
        std::ostream &file,
        GenomeSequence* gs,
        GenomeSequence* csgs,
        MatchedReadBase *mate,
        std::string &cs_read_fragment,
        std::string &cs_data_quality,
        std::string &fragmentTag,
        bool showReferenceBases,
        CigarRoller  &cigarRoller,
        bool    isProperAligned,
        uint16_t    samMateFlag,
        const std::string &sampleGroupID,
        const std::string &alignmentPathTag
    );

    virtual void printOptionalTags(std::ostream &, bool isProperAligned, const std::string &readGroupID, const std::string &alignmentPathTag);

    inline bool qualityIsValid() const
    {
        return quality >= 0;
    }

    inline bool qualityIsFailed() const
    {
        return quality < UNSET_QUALITY;
    }

    inline bool qualityIsUnset() const
    {
        return quality == UNSET_QUALITY;
    }

#if 0
    // not ever used funcitons
    // Xiaowei use #if to comment them out.
    virtual const char *getMateReferenceSequence();
    virtual genomeIndex_t getMatePosition();
    virtual int getISize();     // insertion size -> 0 for single, delta for paired
    const char *getSequence(std::string &sequence, genomeIndex_t genomeMatchPosition);
    virtual bool hasMate();
    virtual MatchedReadBase & getMate();
#endif
    // MapperBase has debugPrint also, but it does more validation.  This
    // method just prints the data to cout:
    void debugPrint();

    void fixBaseRange(int start, int end, /// inclusive boundaries
                      std::string& read_fragment, std::string& data_quality,
                      const std::string& cs_read_fragment, const std::string& cs_data_quality,
                      GenomeSequence* gs, GenomeSequence* csgs,
                      genomeIndex_t genomeMatchPosition,
                      bool isForwardStrand);

    void calibrateSequence(std::string& read_fragment, std::string& data_quality,
                           const std::string& cs_read_fragment, const std::string& cs_data_quality,
                           GenomeSequence* gs, GenomeSequence* csgs,
                           genomeIndex_t genomeMatchPosition,
                           bool isForwardStrand);

    // color space related variables and functions
    ///
    ///Convert 1 base and 1 color to 1 base (decode color space)
    ///
    void BaseAndColorToBase(const char& base, const char& color, char& secondBase)
    {
        int c = color - '0';
        BaseAndColorToBase(base, c, secondBase);
        return;
    };

    ///
    ///Convert 1 base and 1 color to 1 base (decode color space)
    ///
    void BaseAndColorToBase(const char& base, const int& color, char& secondBase)
    {
        int hi = GenomeSequence::base2int[(int)(base)];
        int lo = color;
        if (hi>= 4 || lo >= 4 || lo < 0)
        {
            secondBase = 'N';
            return;
        }
        secondBase= GenomeSequence::int2base[(int)(fromCS2Base[hi<<2 | lo])];
        return;

    };

    ///
    ///Convert 1 base and 1 base to the color of the transiation
    ///
    void BaseAndBaseToColor(const char& base1, const char& base2, int& color)
    {
        int hi = GenomeSequence::base2int[(int)(base1)];
        int lo = GenomeSequence::base2int[(int)(base2)];
        if (hi >= 4 || lo >= 4)
        {
            color = 5;
            return;
        }
        color = fromCS2Base[hi<<2 | lo];
        return;
    };

    ///
    ///Convert 1 base and 1 base to the color of the transiation
    ///
    void BaseAndBaseToColor(const char& base1, const char& base2, char& color)
    {
        int c;
        BaseAndBaseToColor(base1, base2, c);
        color='0'+c;
        return;
    };

    /*
     * convert a color space read to base space
     * Warning: if the input read contains '5', the unknown base, then the output cannot be predicted.
     * @param keepPrimer true: the leading primer the colr space will be stored as the first base of the oupt
     * @return translated read in base space
     */
    std::string convertCSRead(const std::string& csread, bool keepPrimer=false)
    {
        int cslength = csread.size();
        int bslength = (keepPrimer ? cslength : (cslength - 1));
        std::string bsread;
        bsread.resize(bslength);

        char lastBase = csread[0];
        bsread[0] = csread[0];
        for (int i= (keepPrimer ? 1: 0), j = 1;
                i < bslength; i++, j++)
        {
            BaseAndColorToBase(lastBase, csread[j], bsread[i]);
            lastBase = bsread[i];
        }
        return bsread;
    };

    std::string convertCSQuality(const std::string& csqual)
    {
        std::string bsqual = csqual.substr(1);
        return bsqual;
    };

    /*
     * get reverse of a input read (in base space or in color space)
     * @param read input read
     * @return reverse of the input read
     */
    std::string getReverse(const std::string& read)
    {
        std::string ret = read;
        std::reverse(ret.begin(), ret.end());
        return ret;
    };

    /*
     * get compement of a input read (in base space or in color space)
     * @param read input read
     * @return complement of the input read
     */
    std::string getComplement(const std::string& read)
    {
        std::string ret = read;
        // here we can user functor+transform() to speed things up
        static BaseComplementor baseComplementor;
        transform(read.begin(), read.end(), ret.begin(), baseComplementor);
        return ret;
    };

    /*
     * get reverse complement of a input read (in base space or in color space)
     * @param read input read
     * @return reverse complement of the input read
     */
    std::string getReverseComplement(const std::string& read)
    {
        return getReverse(getComplement(read));
    };
};

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
    std::string alignmentPathTag;
protected:

    GenomeSequence  *gs;
    Random *rand;
    MapperUserOptions    mapperOptions;

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
    void setMapperOptions(MapperUserOptions &m)
    {
        //
        // NB forward::mapperOptions and backward::mapperOptions
        // are references to this member:
        //
        mapperOptions = m;      // copy values into here
    }

public:
    void initMapper(GenomeSequence *g, WordIndex *w, WordHash *wl, WordHash * wr, MapperUserOptions &m)
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
