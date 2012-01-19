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

#ifndef __STATS_H__
#define __STATS_H__

#include "BamExecutable.h"
#include "SamFile.h"
#include "PosList.h"

class Stats : public BamExecutable
{
public:
    Stats();

    static void statsDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);

private:
    void reset();


    bool getNextSection(SamFile& samIn);

    // Read the specified dbsnp file using the specified sam/bam header for
    // determining the max reference length and for mapping the
    // reference name to the reference id.
    // Returns a pointer to a newly created dbsnp position list.  The caller of this
    // method should delete the list when done with it.
    PosList* readDBSnp(const String& dbsnpFileName, SamFileHeader& samHeader);


    // Calculate GC Content
    void calcGCContent(SamRecord& samRecord);

    // Check the sequence for the max number of repeats and return it.
    // Pass in repeatLen so it doesn't have to be recalculated for repeat everytime.
    int repeatCounter(const char *sequence, const char* repeat, int repeatLen);
    int repeatCounter(const char *sequence, const char* repeat, const char* reverse, int repeatLen);

    // Pointer to the region list file
    IFILE  myRegionList;

    int myStartPos;
    int myEndPos;

    String myRegBuffer;
    StringArray myRegColumn;


    struct gcContent
    {
        int allReads;
        int gcReads;
        int tmReads;

        void reset() { allReads = 0; gcReads = 0; tmReads = 0;}
    } myGcContent;
};

#endif
