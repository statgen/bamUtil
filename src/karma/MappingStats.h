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

#include "MatchedReadSE.h"
#include "MatchedReadPE.h"
#include "Performance.h"

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <string>
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
    uint64_t    totalReads;          // number of Fastq reads that are from input
    uint64_t    totalMatches;        // number of Fastq reads that are mapped
    uint64_t    badShortData;        // too short a read
    uint64_t    badUnequalData;      // read length!=quality length
    uint64_t    badIndexWords;       // too few index words (<2)

    uint64_t    repeatCountTooHigh;
    uint64_t    earlyDrops;
    uint64_t    noValidMatch;
    uint64_t    lowQualityDrops;
    uint64_t    invalidQualityDrops;
    uint64_t    totalBasesMapped;

    void    recordQualityInfo(MatchedReadBase &match);
    void    printStats(std::ostream &file, std::ostream &fileR);

 public:
    Timing  runTime;
    MappingStatsBase();

    void    updateConsole(bool force = false);

    // update bad data stats, this is called after we prepared the mapper and before the alignment begin
    void updateBadInputStats(const int& rc)
    {
        badShortData += (rc==1);
        badUnequalData += (rc==2);
        badIndexWords += (rc==3);
    }

    // simple getter functions for comparisons
    uint64_t getTotalBasesMapped() 
    {
        return (totalBasesMapped);
    }
    uint64_t getTotalMatches() 
    {
        return (totalMatches );
    };
    uint64_t getTotalReads()
    {
        return (totalReads );
    }
};

class SingleEndStats: public MappingStatsBase
{
 private:
    static const int qualityScoreBucketCount = 100;
    static const int qualityScoreBucketRange = 100 / qualityScoreBucketCount;
    uint64_t    qualityScoreHistogram[qualityScoreBucketCount];

 protected:
    void recordQualityInfo(MatchedReadSE &match);
    void printStats(std::ostream &file, std::ostream &fileR);

 public:
    SingleEndStats();
    void recordMatchedRead(MatchedReadBase& matchedRead);

    // output to baseFileName.stat and baseFileName.R files including mapping related statistics
    void outputStatFile(std::string& baseFileName);
};

class PairedEndStats: public MappingStatsBase
{
 private:
    uint64_t    noneMappable;
    uint64_t    oneMappable;

    static const int shortDistanceRange = 100;    // choose a reasonable short distanace for histgrams
    static const int longDistanceBucketCount = 100;    // choose a reasonable bucket count for whole range
    static const uint32_t bucketSize = (uint32_t)(((uint64_t) 1<<(sizeof(uint32_t)*8)) / longDistanceBucketCount);     // values per bucket

    struct
    {
        uint64_t    count;
        uint64_t    probesInOrder;
    } matchDirection[2][2];

    // record probeA and probeB status in the 8 dimension below:
    // probeA:  [probeA.qualityIsValid()]
    //          [probeA.quality==MatchedReadBase::UNSET_QUALITY]
    //          [probeA.quality==MatchedReadBase::EARLYSTOP_QUALITY]
    //          [probeA.quality==MatchedReadBase::REPEAT_QUALITY]
    //          [probeB.qualityIsValid()]
    // probeB:  [probeB.quality==MatchedReadBase::UNSET_QUALITY]
    //          [probeB.quality==MatchedReadBase::EARLYSTOP_QUALITY]
    //          [probeB.quality==MatchedReadBase::REPEAT_QUALITY]
    uint64_t    qValueBuckets[2][2][2][2][2][2][2][2];
    uint64_t    qValueBucketsA[2][2][2][2];
    uint64_t    qValueBucketsB[2][2][2][2];

    uint64_t    shortDistanceHistogram[2*shortDistanceRange];
    uint64_t    longDistanceHistogram[2*longDistanceBucketCount];

 protected:
    void printHistograms(std::ostream &file, std::ostream &fileR);

    // record the match directions and outer distance distribution between two reads
    void addStats(MatchedReadBase &probeA, MatchedReadBase &probeB, int readLength);

    // record invalid quality (invalid, unset, early_stop, repeat) for each probe.
    void addQValueBuckets(MatchedReadBase &probeA, MatchedReadBase &probeB);

    void printQValueBuckets(std::ostream &file, std::ostream &fileR);

    void printStats(std::ostream &file, std::ostream &fileR);

 public:
    PairedEndStats();
    
    // record statistics from MatchedRead from the bestMatch of each aligner after alignment finished
    void recordMatchedRead(MatchedReadBase& matchedRead1, MatchedReadBase& matchedRead2);

    // output to baseFileName.stat and baseFileName.R files including mapping related statistics
    void outputStatFile(std::string& baseFileName);
    
#ifdef COMPILE_OBSOLETE_CODE
    uint64_t    qValueBucketsDrop[2][2][2][2][2][2][2][2];
    void addQValueBucketsDrop(MatchedReadPE &probeA, MatchedReadPE &probeB);
#endif
};

#endif
