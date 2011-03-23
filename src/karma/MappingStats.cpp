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

#include "MappingStats.h"

#include <iostream>
#include <iomanip>

using std::setw;
using std::endl;
using std::fixed;
using std::setprecision;

MappingStatsBase::MappingStatsBase()
{
    totalReads = 0;
    totalMatches = 0;
    badShortData = 0;
    badUnequalData = 0;
    badIndexWords = 0;

    repeatCountTooHigh = 0;
    earlyDrops = 0;
    noValidMatch = 0;
    lowQualityDrops = 0;
    invalidQualityDrops = 0;
    totalBasesMapped = 0;
    totalBasesMappedAndWritten = 0;
}

void MappingStatsBase::updateConsole(bool force)
{
    double pct = 100.0;

    if (totalReads) pct = 100.0 * totalMatches / totalReads;

    // 0x03ff = 1023,  (totalReads & 0x03ff) == 0 when totalReads is 0, 1024, 2048, .....
    if ((force || ((totalReads & 0x03ff) == 0)) && isatty(fileno(stdout)))
    {
        std::cout << fixed;
        std::cout << setprecision(2);
        std::cout
            << totalMatches
            << " matches of "
            << totalReads
            << " reads ("
            <<  pct
            << "%).\r"
            << std::flush;
    }
}



void MappingStatsBase::recordQualityInfo(MatchedReadBase &match)
{
    switch (match.quality)
    {
        case MatchedReadBase::UNSET_QUALITY:
            noValidMatch++;
            break;
        case MatchedReadBase::EARLYSTOP_QUALITY:
            earlyDrops++;
            break;
        case MatchedReadBase::REPEAT_QUALITY:
            repeatCountTooHigh++;
            break;
        default:
            // currently nothing...
            break;
    }

}

void MappingStatsBase::printStats(std::ostream &file, std::ostream &fileR)
{
    uint64_t    badData = badUnequalData + badShortData + badIndexWords;
    // execution time
    fileR << "totalUserSpaceExecutionTime=" << runTime.interval() << ", ";

    // reader stage inputs all records, drops a small number
    // due to basic format problems, such as mismatch read length/quality length
    fileR << "readerTotalInput=" << totalReads << ", ";
    fileR << "readerDropped=" << badData << ", ";
    fileR << "readerDroppedUnequal=" << badUnequalData << ", ";
    fileR << "readerDroppedBadIndex=" << badIndexWords << ", ";
    fileR << "readerDroppedShort=" << badShortData << ", ";

    // the mapper doesn't drop reads, but it inputs them and classifies their
    // quality scores.
    fileR << "mapperTotalInput=" << totalReads - badData << ", ";
    fileR << "mapperErrorNoValidMatch=" << noValidMatch << ", ";
    fileR << "mapperErrorHighRepeatRegion=" << repeatCountTooHigh << ", ";
    fileR << "mapperErrorEarlyDrops=" << earlyDrops << ", ";
    fileR << "mapperTotalBasesMapped=" << totalBasesMapped << ", ";

    //
    // the writer
    fileR << "writerTotalInput=" << totalReads - badData << ", ";
    fileR << "writerInvalidQualityDrops=" << invalidQualityDrops << ", ";
    fileR << "writerQualityDrops=" << lowQualityDrops << ", ";
    fileR << "writerTotalBasesMappedAndWritten=" << totalBasesMappedAndWritten << ", ";
    fileR << "writerTotalOutput=" << totalMatches << ", ";


    file << fixed;
    file << setprecision(2);

    file << std::endl;
    file << std::endl;

    file
    << totalReads
    << " total reads in "
    << runTime.interval()/3600.0
    << " hours at a rate of "
    << 3600.0 * totalReads / runTime.interval()
    << " probe reads per hour."
    << std::endl << std::endl;

    file
    << badUnequalData
    << " reads were ignored due to having unequal quality and data lengths."
    << std::endl << std::endl;

    file
    << badShortData
    << " reads were ignored due to being too short."
    << std::endl << std::endl;

    file
    << badIndexWords
    << " reads were ignored due to having too few valid index words."
    << std::endl << std::endl;


    file
    << " "
    << setw(8) << noValidMatch
    << " ("
    << setw(6) << 100.0 * noValidMatch / totalReads
    << "%) probe reads had no valid match."
    << std::endl;

    file
    << "+"
    << setw(8) << repeatCountTooHigh
    << " ("
    << setw(6) << 100.0 * repeatCountTooHigh / totalReads
    << "%) probe reads hit high repeat regions."
    << std::endl;

    file
    << "+"
    << setw(8) << earlyDrops
    << " ("
    << setw(6) <<  100.0 * earlyDrops / totalReads
    << "%) probe reads were early stops."
    << std::endl;

    file << "---------" << std::endl;

    //
    // sum of above three should == next one
    //
    file
    << "="
    << setw(8) << noValidMatch+repeatCountTooHigh + earlyDrops
    << " ("
    << setw(6) << 100.0 *(noValidMatch+repeatCountTooHigh+earlyDrops) / totalReads
    << "%) total probe reads that had an invalid quality."
    << std::endl;

    file << std::endl;
    file << std::endl;
    file << std::endl;

    file
    << " "
    << setw(8) << invalidQualityDrops
    << " ("
    << setw(6) << 100.0 * invalidQualityDrops / totalReads
    << "%) probe reads were dropped due to invalid quality (*)."
    << std::endl;

    file
    << "+"
    << setw(8) << lowQualityDrops
    << " ("
    << setw(6) << 100.0 * lowQualityDrops / totalReads
    << "%) probe reads were dropped for low quality."
    << std::endl;

    file
    << "+"
    << setw(8) << totalMatches
    << " ("
    << setw(6) << 100.0 * totalMatches / totalReads
    << "%) probe reads were successfully mapped."
    << std::endl;

    file << "---------" << std::endl;

    file
    << "="
    << setw(8) << totalReads
    << " ("
    << setw(6) << 100.0
    << "%) total probe reads were processed."
    << std::endl;

    file << "(*) this count may include good reads that were paired with a bad read," << std::endl;
    file << "where the pair were dropped together." << std::endl;
    file << std::endl;

    file << setprecision(3);

    file
    << totalBasesMapped
    << " ("
    << totalBasesMapped / 1000.0 / 1000.0 / 1000.0
    << "G) total bases mapped in "
    << runTime.interval()/3600.0
    << " hours."
    << std::endl
    << "This is a rate of "
    << 3600.0 * totalBasesMapped / 1000.0 / 1000.0 / 1000.0 / runTime.interval()
    << "G bases per hour."
    << std::endl
    << "This is a rate of "
    << 30*24*(3600.0 * totalBasesMapped / 1000.0 / 1000.0 / 1000.0 / runTime.interval())
    << "G bases per month per core."
    << std::endl
    << std::endl;

    file << std::endl;
    file << std::endl;
}

SingleEndStats::SingleEndStats()
{
    memset(qualityScoreHistogram, 0, sizeof(qualityScoreHistogram));
}

void SingleEndStats::recordQualityInfo(MatchedReadSE &match)
{
    MappingStatsBase::recordQualityInfo(match);
    double qs = match.getQualityScore();
    assert(qs>=0.0 && qs <= 100.0);

    if (qs >= 100.0)
    {
        qualityScoreHistogram[ qualityScoreBucketCount-1 ]++;
    }
    else
    {
        qualityScoreHistogram[(int)(qs/qualityScoreBucketRange)]++;
    }

}

void SingleEndStats::printStats(std::ostream &file, std::ostream &fileR)
{
    MappingStatsBase::printStats(file, fileR);
    file << std::endl;
    file << std::endl;
    file << "Quality score histogram for all reads:" << std::endl;
    fileR << "qsHistogram=c(";
    int totalInBuckets = 0;
    for (int i=0; i<qualityScoreBucketCount; i++)
    {
        file
        <<  i*qualityScoreBucketRange
        << " to "
        << (i+1)*qualityScoreBucketRange
        << " = "
        << qualityScoreHistogram[i]
        << std::endl;

        fileR << qualityScoreHistogram[i];
        if (i<qualityScoreBucketCount - 1) fileR << ", ";

        totalInBuckets += qualityScoreHistogram[i];
    }
    fileR << ")" << std::endl;
}

PairedEndStats::PairedEndStats()
{
    memset(matchDirection, 0, sizeof(matchDirection));
    memset(shortDistanceHistogram, 0, sizeof(shortDistanceHistogram));
    memset(longDistanceHistogram, 0, sizeof(longDistanceHistogram));
    memset(qValueBuckets, 0, sizeof(qValueBuckets));
    memset(qValueBucketsDrop, 0, sizeof(qValueBucketsDrop));
    memset(qValueBucketsA, 0, sizeof(qValueBucketsA));
    memset(qValueBucketsB, 0, sizeof(qValueBucketsB));
    oneMappable = 0;
    noneMappable = 0;
}

//
// F/R case:
//   (rightEndMatch2 - leftEndMatch1) = (geneMatchPos2 + READSIZE) - geneMatchPos1
//
//   R/F case:
//     (leftEndMatch2 - rightEndMatch1) = geneMatchPos2 - (geneMatchPos1 + READSIZE)
//
//
void PairedEndStats::addStats(MatchedReadPE &probeA, MatchedReadPE &probeB, int readLength)
{
    int64_t    widthRead, limitRead1, limitRead2 = 0;

    //
    // if either is unmapped, don't bother reporting insert
    // sizes
    //
    if (probeA.genomeMatchPosition == INVALID_GENOME_INDEX)
        return;

    if (probeB.genomeMatchPosition == INVALID_GENOME_INDEX)
        return;

    if (probeA.isForward())
    {
        limitRead1 = probeA.genomeMatchPosition;
    }
    else
    {
        limitRead1 = probeA.genomeMatchPosition + readLength;
    }

    if (probeB.isForward())
    {
        limitRead2 = probeB.genomeMatchPosition;
    }
    else
    {
        limitRead2 = probeB.genomeMatchPosition + readLength;
    }

    widthRead = (int64_t)  limitRead2 - limitRead1;

    widthRead *= (probeA.isForward()) ? 1 : -1;

    // count match directions
    matchDirection[probeA.isForward()][probeB.isForward()].count++;
    matchDirection[probeA.isForward()][probeB.isForward()].probesInOrder +=
        (!probeA.isForward() ? probeA.genomeMatchPosition < probeB.genomeMatchPosition :
         probeA.genomeMatchPosition > probeB.genomeMatchPosition);

    // histogram large range
    int longDistanceIndex = widthRead / bucketSize + longDistanceBucketCount;
    assert(longDistanceIndex < 2*longDistanceBucketCount);
    longDistanceHistogram[longDistanceIndex]++;

    // histogram short range
    if (widthRead >= shortDistanceRange)
    {
        widthRead = shortDistanceRange - 1;
    }
    else if (widthRead <= -shortDistanceRange)
    {
        widthRead = -shortDistanceRange;
    }
    int shortDistanceIndex = widthRead + shortDistanceRange;

    assert(shortDistanceIndex < 2*shortDistanceRange);
    shortDistanceHistogram[shortDistanceIndex]++;

    return;
}

void PairedEndStats::printStats(std::ostream &file, std::ostream &fileR)
{
    MappingStatsBase::printStats(file, fileR);

    fileR
    << "mapperPairedMatches = structure(.Dim = c(2,2), .Dimnames = list(c(\"Forward\", \"Reverse\"), c(\"Forward\", \"Reverse\")), c("
    << matchDirection[0][0].count << ", "
    << matchDirection[0][1].count << ", "
    << matchDirection[1][0].count << ", "
    << matchDirection[1][1].count
    << ")), ";

    fileR << "mapperOrderedPairedMatches = structure(.Dim = c(2,2), .Dimnames = list(c(\"Forward\", \"Reverse\"), c(\"Forward\", \"Reverse\")), c("
    << matchDirection[0][0].probesInOrder << ", "
    << matchDirection[0][1].probesInOrder << ", "
    << matchDirection[1][0].probesInOrder << ", "
    << matchDirection[1][1].probesInOrder
    << ")), ";

    file << fixed;
    file << setprecision(2);

    file
    << "Matched "
    << totalMatches
    << " paired ends out of "
    << totalReads
    << " paired ends read for a match rate of "
    << 100.0 * totalMatches/totalReads
    << "%."
    << std::endl << std::endl;

    file
    << " "
    << setw(8) << noneMappable
    << " ("
    << setw(6) << 100.0 * noneMappable/totalReads
    << "%) paired end reads were counted where neither end was mappable."
    << std::endl;

    file
    << " "
    << setw(8) << oneMappable
    << " ("
    << setw(6) << 100.0 * oneMappable/totalReads
    << "%) paired end reads were counted where only one end was mappable."
    << std::endl;

    file << std::endl << std::endl;

    file << "Forward/Forward matches: " <<  matchDirection[0][0].count << " (" << 100.0*matchDirection[0][0].count/totalMatches << "%) with " <<  matchDirection[0][0].probesInOrder << " in order" << std::endl;
    file << "Forward/Reverse matches: " <<  matchDirection[0][1].count << " (" << 100.0*matchDirection[0][1].count/totalMatches << "%) with " <<  matchDirection[0][1].probesInOrder << " in order" << std::endl;
    file << "Reverse/Forward matches: " <<  matchDirection[1][0].count << " (" << 100.0*matchDirection[1][0].count/totalMatches << "%) with " <<  matchDirection[1][0].probesInOrder << " in order" << std::endl;
    file << "Reverse/Reverse matches: " <<  matchDirection[1][1].count << " (" << 100.0*matchDirection[1][1].count/totalMatches << "%) with " <<  matchDirection[1][1].probesInOrder << " in order" << std::endl;

    file << std::endl;

    printQValueBuckets(file, fileR);

    printHistograms(file, fileR);
}

void PairedEndStats::printHistograms(std::ostream &file, std::ostream &fileR)
{
    // low value for R output at the moment, so no output

    // short range histogram:
    file << std::endl << std::endl;
    file << "Short Distance Historgram (" << 2*shortDistanceRange << " buckets)" << std::endl;
    file << "Outliers counted in Top/Bottom buckets." << std::endl;
    file << "Probe Separation\tFrequency Count" << std::endl;

    for (int i=0; i<2*shortDistanceRange; i++)
    {
        file << i - shortDistanceRange << "\t" <<  shortDistanceHistogram[i] << std::endl;
    }
    // long range histogram
    file << std::endl << std::endl;
    file << "Long Distance Historgram (" << longDistanceBucketCount << " buckets)" << std::endl;
    file << "Probe Separation by " << bucketSize << " low\tHigh\tFrequency Count" << std::endl;
    for (int i=0; i<2*longDistanceBucketCount; i++)
    {
        file << setprecision(0)
        << (1.0*i-longDistanceBucketCount-.5)*bucketSize
        << "\t"
        << (1.0*i-longDistanceBucketCount-.5)*bucketSize + bucketSize - 1
        << "\t"
        << longDistanceHistogram[i]
        << std::endl;
    }
}

// dropped cases:
void  PairedEndStats::addQValueBucketsDrop(MatchedReadPE &probeA, MatchedReadPE &probeB)
{
    qValueBucketsDrop
    [probeA.qualityIsValid()]
    [probeA.quality==MatchedReadBase::UNSET_QUALITY]
    [probeA.quality==MatchedReadBase::EARLYSTOP_QUALITY]
    [probeA.quality==MatchedReadBase::REPEAT_QUALITY]
    [probeB.qualityIsValid()]
    [probeB.quality==MatchedReadBase::UNSET_QUALITY]
    [probeB.quality==MatchedReadBase::EARLYSTOP_QUALITY]
    [probeB.quality==MatchedReadBase::REPEAT_QUALITY]
    ++;
}

// all cases before processing:
void  PairedEndStats::addQValueBuckets(MatchedReadPE &probeA, MatchedReadPE &probeB)
{
    qValueBuckets
    [probeA.qualityIsValid()]
    [probeA.quality==MatchedReadBase::UNSET_QUALITY]
    [probeA.quality==MatchedReadBase::EARLYSTOP_QUALITY]
    [probeA.quality==MatchedReadBase::REPEAT_QUALITY]
    [probeB.qualityIsValid()]
    [probeB.quality==MatchedReadBase::UNSET_QUALITY]
    [probeB.quality==MatchedReadBase::EARLYSTOP_QUALITY]
    [probeB.quality==MatchedReadBase::REPEAT_QUALITY]
    ++;
    qValueBucketsA
    [probeA.qualityIsValid()]
    [probeA.quality==MatchedReadBase::UNSET_QUALITY]
    [probeA.quality==MatchedReadBase::EARLYSTOP_QUALITY]
    [probeA.quality==MatchedReadBase::REPEAT_QUALITY]
    ++;
    qValueBucketsB
    [probeB.qualityIsValid()]
    [probeB.quality==MatchedReadBase::UNSET_QUALITY]
    [probeB.quality==MatchedReadBase::EARLYSTOP_QUALITY]
    [probeB.quality==MatchedReadBase::REPEAT_QUALITY]
    ++;
}


void  PairedEndStats::printQValueBuckets(std::ostream &file, std::ostream &fileR)
{
    // low value for R output at the moment, so no output

    const char *isValid[] = {"Invalid","Valid"};
    const char *isValidQuality[] = {"!UNSET_QUALITY","UNSET_QUALITY"};
    const char *isEarlyStopQuality[] = {"!EARLYSTOP_QUALITY","EARLYSTOP_QUALITY"};
    const char *isRepeatQuality[] = {"!REPEAT_QUALITY","REPEAT_QUALITY"};
    int i,j,k,l;
    int i2,j2,k2,l2;

    file << "PROBE STATES PRIOR TO PROCESSING:" << std::endl << std::endl;
    for (i=0; i<2; i++) for (j=0; j<2; j++) for (k=0; k<2; k++) for (l=0; l<2; l++)
                {
                    if (qValueBucketsA[i][j][k][l] == 0 && qValueBucketsB[i][j][k][l]==0) continue;
                    file
                    << isValid[i] << "\t"
                    << isValidQuality[j] << "\t"
                    << isEarlyStopQuality[k] << "\t"
                    << isRepeatQuality[l] << "\t"
                    << qValueBucketsA[i][j][k][l] << "\t"
                    << qValueBucketsB[i][j][k][l]
                    << std::endl;
                }
    file << std::endl << std::endl;

    file << "\t\t\tPROBE A\t\t\tPROBE B\tCOUNT" << std::endl << std::endl;

    for (i=0; i<2; i++) for (j=0; j<2; j++) for (k=0; k<2; k++) for (l=0; l<2; l++)
                {

                    for (i2=0; i2<2; i2++) for (j2=0; j2<2; j2++) for (k2=0; k2<2; k2++) for (l2=0; l2<2; l2++)
                                {
                                    if (qValueBuckets[i][j][k][l][i2][j2][k2][l2]==0) continue;

                                    file
                                    << isValid[i] << "\t"
                                    << isValidQuality[j] << "\t"
                                    << isEarlyStopQuality[k] << "\t"
                                    << isRepeatQuality[l] << "\t"
                                    << isValid[i2] << "\t"
                                    << isValidQuality[j2] << "\t"
                                    << isEarlyStopQuality[k2] << "\t"
                                    << isRepeatQuality[l2] << "\t"
                                    << qValueBuckets[i][j][k][l][i2][j2][k2][l2]
                                    << std::endl;
                                }
                }
    file << std::endl << std::endl;
    file << "PROBE STATES FOR DROPPED PAIRS:" << std::endl << std::endl;
    file << "\t\t\tPROBE A\t\t\tPROBE B\tCOUNT" << std::endl << std::endl;
    for (i=0; i<2; i++) for (j=0; j<2; j++) for (k=0; k<2; k++) for (l=0; l<2; l++)
                {

                    for (i2=0; i2<2; i2++) for (j2=0; j2<2; j2++) for (k2=0; k2<2; k2++) for (l2=0; l2<2; l2++)
                                {
                                    if (qValueBucketsDrop[i][j][k][l][i2][j2][k2][l2]==0) continue;
                                    file
                                    << isValid[i] << "\t"
                                    << isValidQuality[j] << "\t"
                                    << isEarlyStopQuality[k] << "\t"
                                    << isRepeatQuality[l] << "\t"
                                    << isValid[i2] << "\t"
                                    << isValidQuality[j2] << "\t"
                                    << isEarlyStopQuality[k2] << "\t"
                                    << isRepeatQuality[l2] << "\t"
                                    << qValueBucketsDrop[i][j][k][l][i2][j2][k2][l2]
                                    << std::endl;
                                }
                }
}
