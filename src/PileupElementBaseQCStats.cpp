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
    ifprintf(ourOutputFile, "chrom\tchromStart\tchromEnd\tDepth\tQ20Bases\tQ20BasesPct(%)\tTotalReads\tMappedBases\tMappingRate(%)\tMapRate_MQPass(%)\tZeroMapQual(%)\tMapQual<10(%)\tPairedReads(%)\tProperPaired(%)\tDupRate(%)\tQCFailRate(%)\n");
    //    ifprintf(ourOutputFile, "chrom\tchromStart\tchromEnd\tDepth\tQ20Bases(e9)\tQ20BasesPct(%)\tTotalReads(e6)\tMappedBases(e9)\tMappingRate(%)\tMapRate_MQPass(%)\tZeroMapQual(%)\tMapQual<10(%)\tPairedReads(%)\tProperPaired(%)\tDupRate(%)\tQCFailRate(%)\n");
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

    if(SamFlag::isPaired(flag))
    {
        ++numPaired;
        if(SamFlag::isProperPair(flag))
        {
            ++numProperPaired;
        }
    }
    if(SamFlag::isMapped(flag))
    {
        ++numMapped;
        // It is mapped, so check the mapping quality.
        if(record.getMapQuality() < 10)
        {
            ++numLT10MapQ;
            if(record.getMapQuality() == 0)
            {
                ++numZeroMapQ;
            }
        }

        // Check for map fileter.
        if(record.getMapQuality() < ourMinMapQuality)
        {
            numMapQFilter++;
        }
    }
    else
    {
        // Unmapped, so increment the unmapped/filtered reads
        ++numMapQFilter;
    }
    if(SamFlag::isDuplicate(flag))
    {
        ++numDups;
    }
    if(SamFlag::isQCFailure(flag))
    {
        ++numQCFail;
    }

    // Prior to doing any more analysis, check to see if it is filtered out.
    // Always filter out unmapped and mapped with lower than minimum mapping quality,
    // and based on the options, filter out duplicates and QC failures.
    if(!SamFlag::isMapped(flag) || (record.getMapQuality() < ourMinMapQuality) || (ourFilterDups && SamFlag::isDuplicate(flag)) ||
       (ourFilterQCFail && SamFlag::isQCFailure(flag)))
    {
        // filtered read, so do no more analysis on it.
        return;
    }

    // If the readPosition is N/A, this is a deletion.
    if(readIndex == CigarRoller::INDEX_NA)
    {
        // Do nothing unless this read has a base at this position.
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
            myOutputString += "0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\n";
        }
        else
        {
            myOutputString += 100 * ((double)numMapped)/numEntries;
            myOutputString += "\t";
            // calc % unmapped/filtered out due to mapquality
            // then subtract from 100 to get % that are mapped and pass
            // the quality filter.
            myOutputString += 100 - (100 * ((double)numMapQFilter)/numEntries);
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
            myOutputString += "\n";
        }
        ifprintf(ourOutputFile, myOutputString.c_str());
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
    numMapQFilter = 0;
    numZeroMapQ = 0;
    numLT10MapQ = 0;
    numPaired = 0;
    numProperPaired = 0;
    numQCFail = 0;
    myOutputString.Clear();
}

