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
BaseQCOutputFields PileupElementBaseQCStats::ourBaseQCOutputFields(true);


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
    std::string outputStr;
    getHeader(outputStr);
    ifprintf(ourOutputFile, outputStr.c_str());
}


void PileupElementBaseQCStats::getHeader(std::string& outputStr, bool summaryHdr)
{
    outputStr.clear();
    std::string prefix = "";

    if(!summaryHdr)
    {
        outputStr = "chrom\tchromStart";
        prefix = "\t";
        if(ourBaseQCOutputFields.chromEnd)
        {
            outputStr += "\tchromEnd";
        }
    }

    if(ourPercentStats && !summaryHdr)
    {
        if(ourBaseQCOutputFields.depth)        { outputStr += "\tDepth"; }
        if(ourBaseQCOutputFields.q20Bases)     { outputStr += "\tQ20Bases\tQ20BasesPct(%)"; }
        if(ourBaseQCOutputFields.totalReads)   { outputStr += "\tTotalReads"; }
        if(ourBaseQCOutputFields.mapped)       { outputStr += "\tMappedBases\tMappingRate(%)"; }
        if(ourBaseQCOutputFields.passMapQ)     { outputStr += "\tMapRate_MQPass(%)"; }
        if(ourBaseQCOutputFields.zeroMapQ)     { outputStr += "\tZeroMapQual(%)"; }
        if(ourBaseQCOutputFields.mapQlt10)     { outputStr += "\tMapQual<10(%)"; }
        if(ourBaseQCOutputFields.paired)       { outputStr += "\tPairedReads(%)"; }
        if(ourBaseQCOutputFields.properPaired) { outputStr += "\tProperPaired(%)"; }
        if(ourBaseQCOutputFields.dups)         { outputStr += "\tDupRate(%)"; }
        if(ourBaseQCOutputFields.qcFail)       { outputStr += "\tQCFailRate(%)"; }
        if(ourBaseQCOutputFields.avgMapQ)      { outputStr += "\tAverageMapQuality"; }
        if(ourBaseQCOutputFields.numAvgMapQ)   { outputStr += "\tAverageMapQualCount"; }
    }
    else
    {
        if(ourBaseQCOutputFields.totalReads)   { outputStr += prefix; outputStr += "TotalReads"; prefix = "\t"; }
        if(ourBaseQCOutputFields.dups)         { outputStr += prefix; outputStr += "Dups"; prefix = "\t"; }
        if(ourBaseQCOutputFields.qcFail)       { outputStr += prefix; outputStr += "QCFail"; prefix = "\t"; }
        if(ourBaseQCOutputFields.mapped)       { outputStr += prefix; outputStr += "Mapped"; prefix = "\t"; }
        if(ourBaseQCOutputFields.paired)       { outputStr += prefix; outputStr += "Paired"; prefix = "\t"; }
        if(ourBaseQCOutputFields.properPaired) { outputStr += prefix; outputStr += "ProperPaired"; prefix = "\t"; }
        if(ourBaseQCOutputFields.zeroMapQ)     { outputStr += prefix; outputStr += "ZeroMapQual"; prefix = "\t"; }
        if(ourBaseQCOutputFields.mapQlt10)     { outputStr += prefix; outputStr += "MapQual<10"; prefix = "\t"; }
        if(ourBaseQCOutputFields.mapQ255)      { outputStr += prefix; outputStr += "MapQual255"; prefix = "\t"; }
        if(ourBaseQCOutputFields.passMapQ)     { outputStr += prefix; outputStr += "PassMapQual"; prefix = "\t"; }
        if(ourBaseQCOutputFields.avgMapQ)      { outputStr += prefix; outputStr += "AverageMapQuality"; prefix = "\t"; }
        if(ourBaseQCOutputFields.numAvgMapQ)   { outputStr += prefix; outputStr += "AverageMapQualCount"; prefix = "\t"; }
        if(ourBaseQCOutputFields.depth)        { outputStr += prefix; outputStr += "Depth"; prefix = "\t"; }
        if(ourBaseQCOutputFields.q20Bases)     { outputStr += prefix; outputStr += "Q20Bases"; prefix = "\t"; };
    }
    outputStr += "\n";
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
        fprintf(stderr, "\nSummary of Pileup Stats (1st Row is Mean, 2nd Row is Standard Deviation)\n");
        
        std::string hdrStr;
        getHeader(hdrStr, true);
        fprintf(stderr, "%s", hdrStr.c_str());

        std::string prefix = "";

        if(ourBaseQCOutputFields.totalReads)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgTotalReads.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.dups)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgDups.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.qcFail)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgQCFail.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.mapped)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgMapped.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.paired)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgPaired.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.properPaired)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgProperPaired.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.zeroMapQ)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgZeroMapQ.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.mapQlt10)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgLT10MapQ.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.mapQ255)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgMapQ255.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.passMapQ)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgMapQPass.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.avgMapQ)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgAvgMapQ.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.numAvgMapQ)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgAvgMapQCount.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.depth)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgDepth.Mean());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.q20Bases)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgQ20.Mean());
            prefix = "\t";
        }

        fprintf(stderr, "\n");
        prefix = "";

        // Std Deviation
        if(ourBaseQCOutputFields.totalReads)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgTotalReads.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.dups)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgDups.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.qcFail)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgQCFail.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.mapped)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgMapped.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.paired)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgPaired.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.properPaired)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgProperPaired.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.zeroMapQ)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgZeroMapQ.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.mapQlt10)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgLT10MapQ.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.mapQ255)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgMapQ255.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.passMapQ)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgMapQPass.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.avgMapQ)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgAvgMapQ.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.numAvgMapQ)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgAvgMapQCount.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.depth)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgDepth.StandardDeviation());
            prefix = "\t";
        }
        if(ourBaseQCOutputFields.q20Bases)
        {
            fprintf(stderr, "%s%f", prefix.c_str(), avgQ20.StandardDeviation());
            prefix = "\t";
        }
        fprintf(stderr, "\n\n");
    }
}


void PileupElementBaseQCStats::setBaseQCOutputFields(const BaseQCOutputFields& baseQCFields)
{
    ourBaseQCOutputFields = baseQCFields;
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
            if(ourBaseQCOutputFields.chromEnd)
            { 
                myOutputString += "\t";
                myOutputString += startPos + 1;
            }
            if(ourBaseQCOutputFields.depth)
            {
                myOutputString += "\t";
                myOutputString += depth;
            }
            if(ourBaseQCOutputFields.q20Bases)
            {
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
            }
            if(ourBaseQCOutputFields.totalReads)
            {
                myOutputString += "\t";
                myOutputString += numEntries;
                //        myOutputString += numEntries/E6_CALC;
            }
            if(ourBaseQCOutputFields.mapped)
            {
                myOutputString += "\t";
                myOutputString += numMapped;
                //        myOutputString += ((double)numMapped)/E9_CALC;
                myOutputString += "\t";
                if(numEntries == 0)
                {
                    myOutputString += "0.000";
                }
                else
                {
                    myOutputString += 100 * ((double)numMapped)/numEntries;
                }
            }
            if(ourBaseQCOutputFields.passMapQ)
            {
                myOutputString += "\t";
                if(numEntries == 0)
                {
                    myOutputString += "0.000";
                }
                else
                {
                    myOutputString += 100 * ((double)numMapQPass)/numEntries;
                }
            }
            if(ourBaseQCOutputFields.zeroMapQ)
            {
                myOutputString += "\t";
                if(numEntries == 0)
                {
                    myOutputString += "0.000";
                }
                else
                {
                    myOutputString += 100 * ((double)numZeroMapQ)/numEntries;
                }
            }
            if(ourBaseQCOutputFields.mapQlt10)
            {
                myOutputString += "\t";
                if(numEntries == 0)
                {
                    myOutputString += "0.000";
                }
                else
                {
                    myOutputString += 100 * ((double)numLT10MapQ)/numEntries;
                }
            }
            if(ourBaseQCOutputFields.paired)
            {
                myOutputString += "\t";
                if(numEntries == 0)
                {
                    myOutputString += "0.000";
                }
                else
                {
                    myOutputString += 100 * ((double)numPaired)/numEntries;
                }
            }
            if(ourBaseQCOutputFields.properPaired)
            {
                myOutputString += "\t";
                if(numEntries == 0)
                {
                    myOutputString += "0.000";
                }
                else
                {
                    myOutputString += 100 * ((double)numProperPaired)/numEntries;
                }
            }
            if(ourBaseQCOutputFields.dups)
            {
                myOutputString += "\t";
                if(numEntries == 0)
                {
                    myOutputString += "0.000";
                }
                else
                {
                    myOutputString += 100 * ((double)numDups)/numEntries;
                }
            }
            if(ourBaseQCOutputFields.qcFail)
            {
                myOutputString += "\t";
                if(numEntries == 0)
                {
                    myOutputString += "0.000";
                }
                else
                {
                    myOutputString += 100 * ((double)numQCFail)/numEntries;
                }
            }

            if(ourBaseQCOutputFields.avgMapQ)
            {
                myOutputString += "\t";
                if(averageMapQCount != 0)
                {
                    myOutputString += ((double)sumMapQ)/averageMapQCount;
                }
                else
                {
                    myOutputString += "0.000";
                }
            }
            if(ourBaseQCOutputFields.numAvgMapQ)
            {
                myOutputString += "\t";
                myOutputString += averageMapQCount;
                // myOutputString += ((double)averageMapQCount)/E9_CALC;
            }
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
            
            if(ourBaseQCOutputFields.chromEnd)
            {
                myOutputString += "\t"; 
                myOutputString += startPos + 1;
            }
            if(ourBaseQCOutputFields.totalReads)
            {
                myOutputString += "\t";
                myOutputString += numEntries;
            }
            if(ourBaseQCOutputFields.dups)
            {
                myOutputString += "\t";
                myOutputString += numDups;
            }
            if(ourBaseQCOutputFields.qcFail)
            {
                myOutputString += "\t";
                myOutputString += numQCFail;
            }
            if(ourBaseQCOutputFields.mapped)
            {
                myOutputString += "\t";
                myOutputString += numMapped;
            }
            if(ourBaseQCOutputFields.paired)
            {
                myOutputString += "\t";
                myOutputString += numPaired;
            }
            if(ourBaseQCOutputFields.properPaired)
            {
                myOutputString += "\t";
                myOutputString += numProperPaired;
            }
            if(ourBaseQCOutputFields.zeroMapQ)
            {
                myOutputString += "\t";
                myOutputString += numZeroMapQ;
            }
            if(ourBaseQCOutputFields.mapQlt10)
            {
                myOutputString += "\t";
                myOutputString += numLT10MapQ;
            }
            if(ourBaseQCOutputFields.mapQ255)
            {
                myOutputString += "\t";
                myOutputString += numMapQ255;
            }
            if(ourBaseQCOutputFields.passMapQ)
            {
                myOutputString += "\t";
                myOutputString += numMapQPass;
            }
            if(ourBaseQCOutputFields.avgMapQ)
            {
                myOutputString += "\t";
                if(averageMapQCount != 0)
                {
                    myOutputString += ((double)sumMapQ)/averageMapQCount;
                }
                else
                {
                    myOutputString += "0.000";
                }
            }
            if(ourBaseQCOutputFields.numAvgMapQ)
            {
                myOutputString += "\t";
                myOutputString += averageMapQCount;
            }
            if(ourBaseQCOutputFields.depth)
            {
                myOutputString += "\t";
                myOutputString += depth;
            }
            if(ourBaseQCOutputFields.q20Bases)
            {
                myOutputString += "\t";
                myOutputString += numQ20;
            }
            myOutputString += "\n";
            ifprintf(ourOutputFile, myOutputString.c_str());
        }

        if(ourBaseSum)
        {
            // Update the average values.
            if(ourBaseQCOutputFields.totalReads)   { avgTotalReads.Push(numEntries); }
            if(ourBaseQCOutputFields.dups)         { avgDups.Push(numDups); }
            if(ourBaseQCOutputFields.qcFail)       { avgQCFail.Push(numQCFail); }
            if(ourBaseQCOutputFields.mapped)       { avgMapped.Push(numMapped); }
            if(ourBaseQCOutputFields.paired)       { avgPaired.Push(numPaired); }
            if(ourBaseQCOutputFields.properPaired) { avgProperPaired.Push(numProperPaired); }
            if(ourBaseQCOutputFields.zeroMapQ)     { avgZeroMapQ.Push(numZeroMapQ); }
            if(ourBaseQCOutputFields.mapQlt10)     { avgLT10MapQ.Push(numLT10MapQ); }
            if(ourBaseQCOutputFields.mapQ255)      { avgMapQ255.Push(numMapQ255); }
            if(ourBaseQCOutputFields.passMapQ)     { avgMapQPass.Push(numMapQPass); }
            if(ourBaseQCOutputFields.avgMapQ)      { avgAvgMapQ.Push( ((double)sumMapQ)/averageMapQCount); }
            if(ourBaseQCOutputFields.numAvgMapQ)   { avgAvgMapQCount.Push(averageMapQCount); }
            if(ourBaseQCOutputFields.depth)        { avgDepth.Push(depth); }
            if(ourBaseQCOutputFields.q20Bases)     { avgQ20.Push(numQ20); }
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

