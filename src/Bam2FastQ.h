/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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
// This file contains the processing for the executable option "bam2FastQ"
// which converts a bam file into fastq files.

#ifndef __BAM2FASTQ_H__
#define __BAM2FASTQ_H__

#include "BamExecutable.h"
#include "SamRecord.h"
#include "MateMapByCoord.h"
#include "SamCoordOutput.h"

class Bam2FastQ : public BamExecutable
{
public:
    Bam2FastQ();
    ~Bam2FastQ();

    static void bam2FastQDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);
    virtual const char* getProgramName() {return("bam:bam2FastQ");}

private:
    static const char* DEFAULT_FIRST_EXT;
    static const char* DEFAULT_SECOND_EXT;

    void handlePairedRN(SamRecord& samRec);
    void handlePairedCoord(SamRecord& samRec);
    // Handles a record, writing the fastq to the specified file.
    // Releases the record.
    void writeFastQ(SamRecord& samRec, IFILE filePtr,
                    const char* readNameExt = "");
    void cleanUpMateMap(uint64_t readPos, bool flushAll = false);

    void closeFiles();
    void getFileName(String& fn, const String& outBase, const char* ext);

    SamRecordPool myPool;
    MateMapByCoord myMateMap;

    IFILE myUnpairedFile;
    IFILE myFirstFile;
    IFILE mySecondFile;

    int myNumMateFailures;
    int myNumPairs;
    int myNumUnpaired;

    bool myReverseComp;
    bool myRNPlus;

    String myFirstRNExt;
    String mySecondRNExt;
};

#endif
