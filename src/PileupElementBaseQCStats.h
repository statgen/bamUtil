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

#ifndef __PILEUP_ELEMENT_BASE_QC_STATS_H__
#define __PILEUP_ELEMENT_BASE_QC_STATS_H__

#include "PileupElement.h"
#include "SimpleStats.h"

class PileupElementBaseQCStats : public PileupElement
{
public:
    /// Set whether or not to filter duplicates (default is to filter them).
    static void filterDups(bool filterDups);
    /// Set whether or not to filter QC failures (default is to filter them).
    static void filterQCFail(bool filterQCFail);

    /// Set the minimum mapping quality required for a read to be 
    /// processed, and if it is below that, it will be filtered.
    static void setMapQualFilter(int minMapQuality);

    /// Set the output file to the already opened file.
    static void setOutputFile(IFILE outputPtr);

    // Print the output format, make sure you call after setSumStats
    // if you want summary statistics or your output file
    // will have the wrong header..
    static void printHeader();

    /// The default setting is to not do percentStats (percentStats = false)
    static void setPercentStats(bool percentStats);

    /// Set whether or not a summary of all bases should be collected.
    static void setBaseSum(bool baseSum);

    /// Prints a summary to stderr if setBaseSum was passed true.
    static void printSummary();

    PileupElementBaseQCStats();

    virtual ~PileupElementBaseQCStats();

    // Add an entry to this pileup element.  
    virtual void addEntry(SamRecord& record);

    // Perform the analysis associated with this class.
    virtual void analyze();

    // Resets the entry, setting the new position associated with this element.
    virtual void reset(int32_t refPosition);

private:
    PileupElementBaseQCStats(const PileupElement& q);

    void initVars();

    static bool ourFilterDups;
    static bool ourFilterQCFail;
    static int ourMinMapQuality;
    static IFILE ourOutputFile;
    static bool ourPercentStats;
    static const int Q20_CHAR_VAL = 53;
    static const int E9_CALC = 1000000000;
    static const int E6_CALC = 1000000;

    // These are for summary values.
    static RunningStat avgTotalReads; 
    static RunningStat avgDups;
    static RunningStat avgQCFail;
    static RunningStat avgMapped;
    static RunningStat avgPaired;
    static RunningStat avgProperPaired;
    static RunningStat avgZeroMapQ;
    static RunningStat avgLT10MapQ;
    static RunningStat avgMapQ255;
    static RunningStat avgMapQPass;
    static RunningStat avgAvgMapQ;
    static RunningStat avgAvgMapQCount;
    static RunningStat avgDepth;
    static RunningStat avgQ20;

    static bool ourBaseSum;

    int numEntries;
    int numQ20;
    int depth;
    int numDups;
    int numReads;
    int numMapped;
    int numMapQPass;
    int numZeroMapQ;
    int numLT10MapQ;
    int numPaired;
    int numProperPaired;
    int numQCFail;
    int numMapQ255;
    uint64_t sumMapQ;
    int averageMapQCount;

    String myOutputString;
};

#endif
