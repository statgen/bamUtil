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

#ifndef _READ_INDEXER
#define _READ_INDEXER

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

#include <string>
#include <vector>
#include "GenomeSequence.h"
#include "MapperUserOption.h"
#include "SmithWaterman.h"
#include "WordIndex.h"  // for wordInteger_t
//#include "MapperBase.h" // for using MatchedReadBase::UNSET_QUALITY

//
// for a given read, create and contain the various
// word index information, including:
//
//     word position (for faster indexing into genome match position)
//     word int (to hold the actual index itself)
//     count of N bases for each word
//     location of first N base for each word
//
class ReadIndexer
{
    int Word2Integer(const char *read, int offset, int &wordNLocation);
public:
    //
    // mismatchCutoff is calculated each read based on the
    // length of the read, and the expected SNP and machine read
    // error rates.  It is somewhat ad-hoc, but far better than
    // setting a fixed value.
    //
    int mismatchCutoff;

    GenomeSequence  *gs;
    WordIndex       *wordIndex;
    bool    isForward;
    int readLength;     // this helps us compute the reverse wordPositions
    ReadIndexer(MapperUserOption &m);

    //
    // Q: What is this enum structure for?
    //
    enum MatchStrategy {unable, noEdits, edits, secondaryHash} matchStrategy;
    bool setIndexStrategy();    // true -> too few valid index words
    void setReadAndQuality(const char *read, int readLength, const char *quality);
    void resize(uint32_t newWordCount, int offset);

    std::string     read;
    std::string     phredQuality;
    PackedRead      packedRead;
    PackedRead      packedReadWithPad;
    vector<uint8_t> binaryQuality;
    // set per read according to the qualityTrim argument:
    std::string     leftTrimmedBases;
    std::string     rightTrimmedBases;
    std::string     leftTrimmedQualities;
    std::string     rightTrimmedQualities;

    struct Word
    {
        wordInteger_t   wordInt;
        uint32_t        position;
        int             NLocation;
        uint32_t        sortKey;    // normally just the count of genome locations for wordInt
    };

    std::vector<Word> words;

    std::vector<wordInteger_t> wordInts;     // for each index word location, this is the word itself
    std::vector<uint32_t> wordPositions;     // precomputed offset into the read of the word
    std::vector<int> wordNLocations;         // for each index, where is the N if there is one?

    IntHash         checkedPositions;

    // I comment this out, since Paul have new getSumQ function
    // int getSumQ(genomeIndex_t, int &, int);
    int getColorSpaceSumQ(genomeIndex_t, int &, int, int);
    int getColorSpaceSumQOrig(genomeIndex_t, int &, int, int);
    bool checkColorSpaceSNP(const char& reference1,
                            const char& reference2,
                            const char& read1,
                            const char& read2);

    int getSumQ(genomeIndex_t, int &, int, int);
    int getSumQSW(genomeIndex_t, int &misMatches, int bestMatchSumQ, int whichWord);
    int getCigarRoller(genomeIndex_t, int whichWord, int &matchPositionOffset, CigarRoller &cigarRollger);
    int getSumQOrig(genomeIndex_t, int &, int, int);
    int getMismatchCount(genomeIndex_t location, char exclude='\0', char phredQualityThreshold='#')
    {
        int mismatchCount = 0;
        for (uint32_t i=0; i<read.size(); i++)
            if (read[i]!=exclude && read[i]!=(*gs)[location + i] && phredQuality[i] > phredQualityThreshold)
                mismatchCount ++ ;
        return mismatchCount;
    };
    MapperUserOption   &mapperOptions;

    bool checkEdits;
    bool useGapped;     // true->use gapped checking in ::getSumQ()

    void dump();        // debug dump to console for use in gdb - e.g. "call forward.dump()"

    //
    // this method sets a fairly critical value, which is the number of
    // allowed mismatches a read may have before we consider it unmappable
    // at a given location.
    //
    // Historically, this was a constant, but that is awkward when we're trying
    // to map files from different platforms.
    //
    // So, here we by default compute the mismatch cutoff ourselves based on
    // platform specific values, but allow the user to override that value.
    //
    void setMismatchCutoff()
    {
        if (mapperOptions.mismatchCutoff<0)
        {
            //
            // The calculation here is fairly abitrary.  I chose it empirically
            // based on watching what workeed well for Illumina map rates
            // based on the default SNP rate of .001 and read error rate of .01.
            //
            // I have no doubt that it could be improved, and actually a good
            // choice here will allow mapping performance to increase without
            // sacrificing sensitivity.
            //
            mismatchCutoff =
                2.5 +
                8 * read.size() *
                (mapperOptions.expectedErrorRate + mapperOptions.expectedSNPRate);
        }
        else
        {
            mismatchCutoff = mapperOptions.mismatchCutoff;
        }

    }

    //
    // run a local alignment algorithm on Base Space sequences.
    //
    // update: localAlignmentWindowPostion
    // input: localAlignmentWindowSize
    // output: mismatchCount, gapOpenCount, gapCloseCount, gapExtensionCount, cigarRoller
    //
    // returns the lowest cost postion of the read within the given reference window.
    //
    //
    genomeIndex_t localAlignment(
        genomeIndex_t localAlignmentWindowPostion,
        genomeIndex_t localAlignmentWindowSize,
        int &mismatchCount,
        int &quality,
        int &gapOpenCount,
        int &gapCloseCount,
        int &gapExtensionCount,
        CigarRoller &cigarRoller
    );


    //////////////////////////////////////////////////
    // Debug code
    void dumpWords(const std::vector<ReadIndexer::Word>& w) 
    {
        for (unsigned int i = 0; i < w.size(); i++) 
        {
            std::cout << "i=" << i << " ";
            std::cout << "sortKey=" << w[i].sortKey << " ";
            std::cout << std::endl;
        }
    }
    //////////////////////////////////////////////////
};

#if 1
//
// This is a fairly finicky test.
//
// The desire is to stop computing sumQ when we know it can't contribute
// to the final result of the best match.  The QUALITY_MARGIN macro helps
// determine this, but it is totally arbitrary.
//
// It turns out that we want to disregard the mismatchCutoff when effectively
// there is no existing good match.  An example of this is in Test.cpp where
// the mismatches are very low quality scores (sum=66), but the total mismatch
// count for the 108 base read is 12.  Most reasonable fixed cutoff values
// here cause the read to be rejected, even though clearly it was a good map.
//
// QUALITY_MARGIN is arbitrary, and we should verify whether or not it
// is any kind of an optimal choice.
//
#define QUALITY_MARGIN 60
#define CHECK_QUALITY_MARGIN if(mismatchCount > mismatchCutoff && (quality > bestMatchSumQ + QUALITY_MARGIN)) return -1
#else
#define CHECK_QUALITY_MARGIN
#endif

/*
 * return the sum of the quality at the mismatched sites for base space reads
 * @param matchPositionOffset   from where we start compare read and reference genome
 * @param mismatchCount         store results of how many bases are mismatched
 * @param bestMatchSumQ         so far the best sumQ value (TODO: probaby useless, may need to elimiate this param)
 * @param whichWord             indicate the position of the word index in the read
 * @return                      the quality which equals the sum of the base qualities at mismatched sites
 */
inline int ReadIndexer::getSumQ(
    genomeIndex_t matchPosition,
    int &mismatchCount,
    int bestMatchSumQ,
    int whichWord   // passed to getSumQSW if useGapped->true
)
{
    // fairly horrible hack, but cleanest way for now:
    if (useGapped) return getSumQSW(matchPosition, mismatchCount, bestMatchSumQ, whichWord);

    // if we have to check beyond the boundary of the maximum bases of reference genome,
    // or the lowest boundary 0, we will have to return UNSET_QUALITY, which is -1 ,
    // as the quality cannot be defined.
    if (matchPosition + this->packedRead.size() >= gs->getNumberBases() ||
            matchPosition < 0)
        return -1;//MatchedReadBase::UNSET_QUALITY;

    int quality = 0;
    mismatchCount = 0;

// XXX circular reference - fix it later...
//    if(bestMatchSumQ == MatchedReadBase::UNSET_QUALITY) bestMatchSumQ = 9999999;
// TODO by Xiaowei:
// bestMatchSumQ is never used in this function, so it maybe need to be deleted.
    if (bestMatchSumQ == -1) bestMatchSumQ = 9999999;

    vector<uint8_t>::iterator readIndexPosition = this->packedRead.packedBases.begin();

    int baseIndex = 0;

    // compare the firse base alone if the genome match position
    // is at an odd position.
    if (matchPosition&0x01)
    {
        if (read[baseIndex] != (*gs)[matchPosition++])
        {
            quality += binaryQuality[0];
            mismatchCount++;
            CHECK_QUALITY_MARGIN;
        }
        readIndexPosition = this->packedReadWithPad.packedBases.begin() + 1;
        baseIndex++;
    }

    uint8_t *referenceBytes = gs->getDataPtr(matchPosition);
    int byteCount = (this->packedRead.size() - baseIndex) / 2;


    // now compare two bases at a time, using 4 bit encoded
    // bases.
    while (byteCount--)
    {
        // two base compare:
        if (*readIndexPosition != *referenceBytes)
        {
            // at least one of the two bases is different,
            // so check them both:
            if (read[baseIndex] != (*gs)[matchPosition])
            {
                quality += binaryQuality[baseIndex];
                mismatchCount++;
            }
            if (read[baseIndex+1] != (*gs)[matchPosition+1])
            {
                quality += binaryQuality[baseIndex+1];
                mismatchCount++;
            }
            CHECK_QUALITY_MARGIN;
        }
        matchPosition+=2;
        baseIndex+=2;
        readIndexPosition++;
        referenceBytes++;
    }

    // if there is an odd left over base at the end, remember to
    // compre it as well:
    if (baseIndex < (int) read.size())
    {
        if (read[baseIndex] != (*gs)[matchPosition++])
        {
            quality += binaryQuality[baseIndex];
            mismatchCount++;
            CHECK_QUALITY_MARGIN;
        }
    }

    return quality;
}

#if 0
//
// this is a non-optimized version of getSumQ - I'd like to leave it here
// awhile longer for testing purposes in the future.
//
inline int ReadIndexer::getSumQOrig(
    genomeIndex_t matchPosition,
    int &mismatchCount,
    int bestMatchSumQ,
    int whichWord   // unused here
)
{
    // insert fast checking here.
    // XXX this is the slower old code... replace ASAP
    int quality = 0;
    mismatchCount = 0;

//    if(bestMatchSumQ == MatchedReadBase::UNSET_QUALITY) bestMatchSumQ = 9999999;
    if (bestMatchSumQ == -1) bestMatchSumQ = 9999999;

    for (int readIndexPosition=0; readIndexPosition<read.Length(); readIndexPosition++, matchPosition++)
    {
        if (read[readIndexPosition] != (*gs)[ matchPosition])
        {
            quality += binaryQuality[readIndexPosition];
#if 0
            if (quality > mapperOptions.qValueCutoff) return -1;
            if (readIndexPosition < (int) mapperOptions.readIndexCutoff)
            {
                mismatchCount++;
                if (mismatchCount > mismatchCutoff) return -1;
            }
#endif
            if (readIndexPosition < mapperOptions.readIndexCutoff)
            {
                mismatchCount++;
                if (mismatchCount > mismatchCutoff) return -1;

            }

            CHECK_QUALITY_MARGIN;
        }
    }

    return quality;
}
#endif

/*
 * return the sum of the quality at the mismatched sites for color space reads
 * NOTE: for base '5', the unknown read, its quality will not be counted
 * @param matchPositionOffset   from where we start compare read and reference genome
 * @param mismatchCount         store results of how many bases are mismatched
 * @param bestMatchSumQ         so far the best sumQ value (TODO: probaby useless, may need to elimiate this param)
 * @param whichWord             indicate the position of the word index in the read
 * @return                      the quality which equals the sum of the base qualities at mismatched sites
 */
inline int ReadIndexer::getColorSpaceSumQ(
    genomeIndex_t matchPosition,
    int &mismatchCount,
    int bestMatchSumQ,
    int whichWord
)
{
    // if we have to check beyond the boundary of the maximum bases of reference genome,
    // or the lowest boundary 0, we will have to return UNSET_QUALITY, which is -1 ,
    // as the quality cannot be defined.
    if (matchPosition + this->packedRead.size() >= gs->getNumberBases() ||
            matchPosition < 0)
        return -1;//MatchedReadBase::UNSET_QUALITY;

//#define DEBUG_GETCOLORSPACESUMQ
#ifdef DEBUG_GETCOLORSPACESUMQ
    std::cout<<"getColorSpaceSumQ() ref_pos="<<matchPosition<<std::endl;
    std::cout<<"read: "<<read<<std::endl;
    std::cout<<"ref : ";
    for (unsigned int i=0; i<read.size(); i++)
        std::cout<<(*gs)[matchPosition+i];
    std::cout<<std::endl;
#endif
    // insert fast checking here.
    // XXX this is the slower old code... replace ASAP
    int quality = 0;
    mismatchCount = 0;

//    if(bestMatchSumQ == MatchedReadBase::UNSET_QUALITY) bestMatchSumQ = 9999999;
    if (bestMatchSumQ == -1) bestMatchSumQ = 9999999;

    // here we will make it suitable for color space reads.
    unsigned int consecutiveMismatch = 0;
    // two consecutive mismatch could be a SNP ( freq: 0.1 % , so cap to 30)
    const unsigned int QUALITY_CAP = 30;

    int readLength= read.size();
    for (int readIndexPosition=0; readIndexPosition<readLength; readIndexPosition++, matchPosition++)
    {
        if (read[readIndexPosition] != (*gs)[ matchPosition] && read[readIndexPosition] != '5'
                && binaryQuality[readIndexPosition] > ('#' -'!') /* TODO Xiaowei: clear this up by using trimming!!*/
           )
        {
            consecutiveMismatch ++;
            quality += binaryQuality[readIndexPosition];
#if 0
            std::cout
                << index << ": want "
                << (*gs)[ matchPosition + index] << "  got " << read[index]
                << ", Q=" << readQuality[index]
                << ", sumQ=" << quality << std::endl;
#endif
#if 0
            if (quality > mapperOptions.qValueCutoff) return -1;
            if (readIndexPosition < (int) mapperOptions.readIndexCutoff)
            {
                mismatchCount++;
                if (mismatchCount > mismatchCutoff) return -1;
            }
#endif
            mismatchCount++;
            // mismatchCutoff is computed once per read, and used
            // here.
            CHECK_QUALITY_MARGIN;
        }
        else   // here the read matches the reference
        {
            if (consecutiveMismatch == 2 &&
                    checkColorSpaceSNP((*gs)[ matchPosition-1],
                                       (*gs)[ matchPosition-2],
                                       read[readIndexPosition-1],
                                       read[readIndexPosition-2]))
            {
                quality -= (binaryQuality[readIndexPosition-1] +
                            binaryQuality[readIndexPosition-2] -
                            QUALITY_CAP);
                mismatchCount --;
            }
            consecutiveMismatch = 0;
        }
    }

    // adjust quality if mismatch is in the last bp of the read.
    if (consecutiveMismatch == 2  &&
            checkColorSpaceSNP((*gs)[ matchPosition-1],
                               (*gs)[ matchPosition-2],
                               read[readLength-1],
                               read[readLength-2]))
        quality -= (binaryQuality[readLength-1] +
                    binaryQuality[readLength-2] -
                    QUALITY_CAP);
#ifdef DEBUG_GETCOLORSPACESUMQ
    std::cout<<"quality="<<quality<<std::endl;
#endif
    return quality;
}

// In color space read, two consecutive mismatches can mean a SNP position
// This function checks that:
// given 2 positions from color space reference and 2 positions from reads
// can these 2 position be a potential SNP ?
// e.g.
// If reference is 00, a potential reads representing a SNP can be : 00 11 22 33
// If reference is 01, a potential reads representing a SNP can be : 01 10 23 32
inline bool ReadIndexer::checkColorSpaceSNP(const char& reference1,
        const char& reference2,
        const char& read1,
        const char& read2)
{
    if (
        // condition 1: if reference is 00, 11, 22, 33, then
        // SNP mutation should also be 00, 11, 22, 33
        (
            (reference1 == reference2)
            && (read1 == read2)
        )
        ||
        // condition 2: if reference is 01, then SNP mutation must be
        // either 10, 23, 32
        (
            (reference1 != reference2) // e.g. reference is 01
            && (
                (// e.g. read is 01
                    (read1 == reference1)
                    && (read2 == reference2)
                )
                ||
                (// e.g. read is 10
                    (read1 == reference2)
                    && (read2 == reference1)
                )
                ||
                (// e.g. read is 23 or 32
                    (read1    != read2)
                    && (read1 != reference1)
                    && (read1 != reference2)
                    && (read2 != reference1)
                    && (read2 != reference2)
                )
            )
        )
    ) // end if
        return true;
    return false;
}

#endif
