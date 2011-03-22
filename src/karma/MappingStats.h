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

#ifndef _MAPPING_STATS
#define _MAPPING_STATS

#include "MapperPE.h"
#include "MapperSE.h"
#include "Performance.h"

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#include <iostream>

//
// keep and display genome matching statistics.
// this separates a bunch of nitpicking detail from the core loops,
// which helps clarify them.
//
// This is a simple class hierarchy: base, base->single_end, base->paired_end
//
class MappingStatsBase
{
protected:
    uint64_t    totalReads;
    uint64_t    totalMatches;
    uint64_t    badShortData;        // too short a read
    uint64_t    badUnequalData;      // read length!=quality length
    uint64_t    badIndexWords;       // too few index words (<2)

    uint64_t    repeatCountTooHigh;
    uint64_t    earlyDrops;
    uint64_t    noValidMatch;
    uint64_t    lowQualityDrops;
    uint64_t    invalidQualityDrops;
    uint64_t    totalBasesMapped;
    uint64_t    totalBasesMappedAndWritten;

public:
    MappingStatsBase();
    Timing  runTime;
    void    updateConsole(bool force = false);

    void    recordQualityInfo(MatchedReadBase &match);
    void    printStats(std::ostream &file, std::ostream &fileR);

    // simple functions for comparisons
    bool isTotalBasesMappedAndWrittenLessThan(const uint64_t& i)
    {
        return (totalBasesMappedAndWritten < i);
    };
    bool isTotalMatchesLessThan(const uint64_t& i)
    {
        return (totalMatches < i);
    };
    bool isTotalReadsLessThan(const uint64_t& i)
    {
        return (totalReads < i);
    }
    void updateBadInputStats(const int& rc)
    {
        badShortData += (rc==1);
        badUnequalData += (rc==2);
        badIndexWords += (rc==3);
    }
    void addTotalReadsByOne()
    {
        totalReads++;
    }
    void addTotalBasesMappedBy(const int& i)
    {
        totalBasesMapped += i;
    }
    void addTotalMatchesByOne()
    {
        totalMatches ++;
    }
    void addTotalBasesMappedAndWritten(const int& i)
    {
        totalBasesMappedAndWritten += i;
    }
};

class SingleEndStats: public MappingStatsBase
{
private:
    static const int qualityScoreBucketCount = 100;
    static const int qualityScoreBucketRange = 100 / qualityScoreBucketCount;
    uint64_t    qualityScoreHistogram[qualityScoreBucketCount];
public:
    SingleEndStats();
    void recordQualityInfo(MatchedReadSE &match);
    void printStats(std::ostream &file, std::ostream &fileR);
};

class PairedEndStats: public MappingStatsBase
{
private:
    uint64_t    noneMappable;
    uint64_t    oneMappable;

    static const int shortDistanceRange = 1000;    // choose a reasonable short distanace for histgrams
    static const int longDistanceBucketCount = 1000;    // choose a reasonable bucket count for whole range
    static const uint32_t bucketSize = (uint32_t)(((uint64_t) 1<<(sizeof(uint32_t)*8)) / longDistanceBucketCount);     // values per bucket

    struct
    {
        uint64_t    count;
        uint64_t    probesInOrder;
    } matchDirection[2][2];

    uint64_t    qValueBuckets[2][2][2][2][2][2][2][2];
    uint64_t    qValueBucketsDrop[2][2][2][2][2][2][2][2];
    uint64_t    qValueBucketsA[2][2][2][2];
    uint64_t    qValueBucketsB[2][2][2][2];

    uint64_t    shortDistanceHistogram[2*shortDistanceRange];
    uint64_t    longDistanceHistogram[2*longDistanceBucketCount];

public:
    PairedEndStats();
    void printHistograms(std::ostream &file, std::ostream &fileR);

    void addStats(MatchedReadPE &probeA, MatchedReadPE &probeB, int readLength);
    void addQValueBuckets(MatchedReadPE &probeA, MatchedReadPE &probeB);
    void addQValueBucketsDrop(MatchedReadPE &probeA, MatchedReadPE &probeB);
    void printQValueBuckets(std::ostream &file, std::ostream &fileR);
    void printStats(std::ostream &file, std::ostream &fileR);

};

#endif
