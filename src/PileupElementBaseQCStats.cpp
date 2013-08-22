/*
 *  Copyright (C) 2011  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//////////////////////////////////////////////////////////////////////////
// This file contains the processing for the executable option "stats"
// which generates some statistics for SAM/BAM files.

#include <stdexcept>
#include "PileupElementBaseQCStats.h"
#include "SamFlag.h"

/////////////////////////////////////////////////////////////////////////////
//
// PileupElementBaseQCStats
//

bool PileupElementBaseQCStats::ourFilterDups = true;
bool PileupElementBaseQCStats::ourFilterQCFail = true;
int PileupElementBaseQCStats::ourMinMapQuality = 0;
IFILE PileupElementBaseQCStats::ourOutputFile = 0;
bool PileupElementBaseQCStats::ourPercentStats = false;
RunningStat PileupElementBaseQCStats::avgTotalReads; 
RunningStat PileupElementBaseQCStats::avgDups;
RunningStat PileupElementBaseQCStats::avgQCFail;
RunningStat PileupElementBaseQCStats::avgMapped;
RunningStat PileupElementBaseQCStats::avgPaired;
RunningStat PileupElementBaseQCStats::avgProperPaired;
RunningStat PileupElementBaseQCStats::avgZeroMapQ;
RunningStat PileupElementBaseQCStats::avgLT10MapQ;
RunningStat PileupElementBaseQCStats::avgMapQ255;
RunningStat PileupElementBaseQCStats::avgMapQPass;
RunningStat PileupElementBaseQCStats::avgAvgMapQ;
RunningStat PileupElementBaseQCStats::avgAvgMapQCount;
RunningStat PileupElementBaseQCStats::avgDepth;
RunningStat PileupElementBaseQCStats::avgQ20;
bool PileupElementBaseQCStats::ourBaseSum = false;

void PileupElementBaseQCStats::filterDups(bool filterDups)
{
    ourFilterDups = filterDups;
}

void PileupElementBaseQCStats::filterQCFail(bool filterQCFail)
{
    ourFilterQCFail = filterQCFail;
}

void PileupElementBaseQCStats::setMapQualFilter(int minMapQuality)
{
    ourMinMapQuality = minMapQuality;
}

void PileupElementBaseQCStats::setOutputFile(IFILE outputPtr)
{
    ourOutputFile = outputPtr;
}

void PileupElementBaseQCStats::printHeader()
{
    if(ourPercentStats)
    {
        ifprintf(ourOutputFile, "chrom\tchromStart\tchromEnd\tDepth\tQ20Bases\tQ20BasesPct(%)\tTotalReads\tMappedBases\tMappingRate(%)\tMapRate_MQPass(%)\tZeroMapQual(%)\tMapQual<10(%)\tPairedReads(%)\tProperPaired(%)\tDupRate(%)\tQCFailRate(%)\tAverageMapQuality\tAverageMapQualCount\n");
        //    ifprintf(ourOutputFile, "chrom\tchromStart\tchromEnd\tDepth\tQ20Bases(e9)\tQ20BasesPct(%)\tTotalReads(e6)\tMappedBases(e9)\tMappingRate(%)\tMapRate_MQPass(%)\tZeroMapQual(%)\tMapQual<10(%)\tPairedReads(%)\tProperPaired(%)\tDupRate(%)\tQCFailRate(%)\tAverageMapQuality\tAverageMapQualCount(e9)\n");
    }
    else
    {
        ifprintf(ourOutputFile, "chrom\tchromStart\tchromEnd\tTotalReads\tDups\tQCFail\tMapped\tPaired\tProperPaired\tZeroMapQual\tMapQual<10\tMapQual255\tPassMapQual\tAverageMapQuality\tAverageMapQualCount\tDepth\tQ20Bases\n");
    }
}


void PileupElementBaseQCStats::setPercentStats(bool percentStats)
{
    ourPercentStats = percentStats;
}


void PileupElementBaseQCStats::setBaseSum(bool baseSum)
{
    ourBaseSum = baseSum;
}


void PileupElementBaseQCStats::printSummary()
{
    if(ourBaseSum)
    {
        fprintf(stderr, "\nSummary of Pileup Stats (1st Row is Mean, 2nd Row is Standard Deviation)\nTotalReads\tDups\tQCFail\tMapped\tPaired\tProperPaired\tZeroMapQual\tMapQual<10\tMapQual255\tPassMapQual\tAverageMapQuality\tAverageMapQualCount\tDepth\tQ20Bases\n");
        
        fprintf(stderr, 
                 "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                 avgTotalReads.Mean(), avgDups.Mean(), avgQCFail.Mean(),
                 avgMapped.Mean(), avgPaired.Mean(), avgProperPaired.Mean(),
                 avgZeroMapQ.Mean(), avgLT10MapQ.Mean(), avgMapQ255.Mean(), 
                 avgMapQPass.Mean(), avgAvgMapQ.Mean(), avgAvgMapQCount.Mean(),
                 avgDepth.Mean(), avgQ20.Mean());
        fprintf(stderr, 
                 "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n\n",
                 avgTotalReads.StandardDeviation(), avgDups.StandardDeviation(), avgQCFail.StandardDeviation(),
                 avgMapped.StandardDeviation(), avgPaired.StandardDeviation(), avgProperPaired.StandardDeviation(),
                 avgZeroMapQ.StandardDeviation(), avgLT10MapQ.StandardDeviation(), avgMapQ255.StandardDeviation(), 
                 avgMapQPass.StandardDeviation(), avgAvgMapQ.StandardDeviation(), avgAvgMapQCount.StandardDeviation(),
                 avgDepth.StandardDeviation(), avgQ20.StandardDeviation());
    }
}


PileupElementBaseQCStats::PileupElementBaseQCStats()
    : PileupElement()
{
    initVars();
}


PileupElementBaseQCStats::~PileupElementBaseQCStats()
{
}


// Add an entry to this pileup element.  
void PileupElementBaseQCStats::addEntry(SamRecord& record)
{
    // Call the base class:
    PileupElement::addEntry(record);

    // Check to see if it is is a "match/mismatch"
    Cigar* cigar = record.getCigarInfo();
    
    if(cigar == NULL)
    {
        throw std::runtime_error("Failed to retrieve cigar info from the record.");
    }

    int32_t readIndex = 
        cigar->getQueryIndex(getRefPosition(), record.get0BasedPosition());

    // Increment the counts
    ++numEntries;
    uint16_t flag = record.getFlag();
    
    if(SamFlag::isDuplicate(flag))
    {
        ++numDups;
    }
    if(SamFlag::isQCFailure(flag))
    {
        ++numQCFail;
    }

    if((ourFilterDups && SamFlag::isDuplicate(flag)) ||
       (ourFilterQCFail && SamFlag::isQCFailure(flag)))
    {
        // Filtered for duplicate/QC
        return;
    }

    // Check if it mapped.
    if(!SamFlag::isMapped(flag))
    {
        // Not mapped, so filter.
        return;
    }

    ++numMapped;

    // It is mapped, so check pairing.
    if(SamFlag::isPaired(flag))
    {
        ++numPaired;
        if(SamFlag::isProperPair(flag))
        {
            ++numProperPaired;
        }
    }

    // It is mapped, so check the mapping quality.
    if(record.getMapQuality() < 10)
    {
        ++numLT10MapQ;
        if(record.getMapQuality() == 0)
        {
            ++numZeroMapQ;
        }
    }
    
    // Check for map filter.
    if(record.getMapQuality() >= ourMinMapQuality)
    {
        numMapQPass++;
    }

    if(record.getMapQuality() == 255)
    {
        // Do not include mapping quality greater than 255.
        ++numMapQ255;
        return;
    }
    else
    {
        // Store the mapping quality for calculating the average mapping quality
        // if it is not 255 (unknown).
        // Also store the number of entries used for calculating the mapping
        // quality.
        // Store this for overflow check.
        uint32_t prevMapQ = sumMapQ;
        sumMapQ += record.getMapQuality();
        // Increment the number of entries in the mapping quality sum.
        ++averageMapQCount;
        
        // Check for overflow.
        if(prevMapQ > sumMapQ)
        {
            std::cerr << "Mapping Quality Overflow for chromosome: "
                      << getChromosome() << ", Position: " << getRefPosition() 
                      << "\n";
            // So just calculate the previous average, then start adding to that.
            // This is not a good indicator, but it really shouldn't overflow.
            --averageMapQCount;
            sumMapQ = prevMapQ / averageMapQCount;
            sumMapQ += record.getMapQuality();
            averageMapQCount = 2;
        }
    }

    // Now that we have added the mapping quality for the average calculation,
    // filter out lower than minimum mapping quality,
    if(record.getMapQuality() < ourMinMapQuality)
    {
        return;
    }

    // Skip deletion for q20/depth since it does not have a base quality.
    if(readIndex == CigarRoller::INDEX_NA)
    {
        // filtered read, so do no more analysis on it.
        return;      
    }

    ++depth;
 
    // Check for Q20 base.
    if(record.getQuality(readIndex) >= Q20_CHAR_VAL)
    {
        // Greater than or equal to q20.
        ++numQ20;
    }    
}


// Perform the alalysis associated with this class.  May be a simple print, 
// a calculation, or something else.  Typically performed when this element
// has been fully populated by all records that cover the reference position.
void PileupElementBaseQCStats::analyze()
{
    // Only output if the position is covered.
    if(numEntries != 0)
    {
        if(ourPercentStats)
        {
            int32_t startPos = getRefPosition();
            myOutputString = getChromosome();
            myOutputString += "\t";
            myOutputString += startPos;
            myOutputString += "\t";
            myOutputString += startPos + 1;
            myOutputString += "\t";
            myOutputString += depth;
            myOutputString += "\t";
            myOutputString += numQ20;
            //        myOutputString += (double(numQ20))/E9_CALC;
            myOutputString += "\t";
            if(depth == 0)
            {
                myOutputString += (double)0;
            }
            else
            {
                myOutputString += 100 * (double(numQ20))/depth;
            }
            myOutputString += "\t";
            myOutputString += numEntries;
            //        myOutputString += numEntries/E6_CALC;
            myOutputString += "\t";
            myOutputString += numMapped;
            //        myOutputString += ((double)numMapped)/E9_CALC;
            myOutputString += "\t";
            if(numEntries == 0)
            {
                myOutputString += "0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000";
            }
            else
            {
                myOutputString += 100 * ((double)numMapped)/numEntries;
                myOutputString += "\t";
                myOutputString += 100 * ((double)numMapQPass)/numEntries;
                myOutputString += "\t";
                myOutputString += 100 * ((double)numZeroMapQ)/numEntries;
                myOutputString += "\t";
                myOutputString += 100 * ((double)numLT10MapQ)/numEntries;
                myOutputString += "\t";
                myOutputString += 100 * ((double)numPaired)/numEntries;
                myOutputString += "\t";
                myOutputString += 100 * ((double)numProperPaired)/numEntries;
                myOutputString += "\t";
                myOutputString += 100 * ((double)numDups)/numEntries;
                myOutputString += "\t";
                myOutputString += 100 * ((double)numQCFail)/numEntries;
            }

            myOutputString += "\t";
            if(averageMapQCount != 0)
            {
                myOutputString += ((double)sumMapQ)/averageMapQCount;
            }
            else
            {
                myOutputString += "0.000";
            }
            myOutputString += "\t";
            myOutputString += averageMapQCount;
            // myOutputString += ((double)averageMapQCount)/E9_CALC;
            myOutputString += "\n";
        
            ifprintf(ourOutputFile, myOutputString.c_str());
        }
        else
        {
            // Summary stats.
            int32_t startPos = getRefPosition();
            myOutputString = getChromosome();
            myOutputString += "\t";
            myOutputString += startPos;
            myOutputString += "\t";
            myOutputString += startPos + 1;
            myOutputString += "\t";
            myOutputString += numEntries;
            myOutputString += "\t";
            myOutputString += numDups;
            myOutputString += "\t";
            myOutputString += numQCFail;
            myOutputString += "\t";
            myOutputString += numMapped;
            myOutputString += "\t";
            myOutputString += numPaired;
            myOutputString += "\t";
            myOutputString += numProperPaired;
            myOutputString += "\t";
            myOutputString += numZeroMapQ;
            myOutputString += "\t";
            myOutputString += numLT10MapQ;
            myOutputString += "\t";
            myOutputString += numMapQ255;
            myOutputString += "\t";
            myOutputString += numMapQPass;
            myOutputString += "\t";
            if(averageMapQCount != 0)
            {
                myOutputString += ((double)sumMapQ)/averageMapQCount;
            }
            else
            {
                myOutputString += "0.000";
            }
            myOutputString += "\t";
            myOutputString += averageMapQCount;
            myOutputString += "\t";
            myOutputString += depth;
            myOutputString += "\t";
            myOutputString += numQ20;
            myOutputString += "\n";
        
            ifprintf(ourOutputFile, myOutputString.c_str());
        }

        if(ourBaseSum)
        {
            // Update the average values.
            avgTotalReads.Push(numEntries);
            avgDups.Push(numDups);
            avgQCFail.Push(numQCFail);
            avgMapped.Push(numMapped);
            avgPaired.Push(numPaired);
            avgProperPaired.Push(numProperPaired);
            avgZeroMapQ.Push(numZeroMapQ);
            avgLT10MapQ.Push(numLT10MapQ);
            avgMapQ255.Push(numMapQ255);
            avgMapQPass.Push(numMapQPass);
            avgAvgMapQ.Push( ((double)sumMapQ)/averageMapQCount);
            avgAvgMapQCount.Push(averageMapQCount);
            avgDepth.Push(depth);
            avgQ20.Push(numQ20);
        }
    }
}


// Resets the entry, setting the new position associated with this element.
void PileupElementBaseQCStats::reset(int32_t refPosition)
{
    // Call the base class.
    PileupElement::reset(refPosition);

    initVars();
}

void PileupElementBaseQCStats::initVars()
{
    numEntries = 0;
    numQ20 = 0;
    depth = 0;
    numDups = 0;
    numReads = 0;
    numMapped = 0;
    numMapQPass = 0;
    numZeroMapQ = 0;
    numLT10MapQ = 0;
    numPaired = 0;
    numProperPaired = 0;
    numQCFail = 0;
    numMapQ255 = 0;
    sumMapQ = 0;
    averageMapQCount = 0;
    myOutputString.Clear();
}

