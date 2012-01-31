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
// This file contains the processing for the executable option "clipOverlap"
// which clips overlapping read pairs.

#ifndef __CLIP_OVERLAP_H__
#define __CLIP_OVERLAP_H__

#include "BamExecutable.h"
#include "SamFile.h"
#include "MateMapByCoord.h"
#include "SamCoordOutput.h"
#include "SimpleStats.h"

class ClipOverlap : public BamExecutable
{
public:
    ClipOverlap();

    static void clipOverlapDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);

private:
    static const int DEFAULT_POOL_SIZE = 1000000;

    int clipSortedByReadName(SamFile& samIn, SamFile& outFile);

    int clipSortedByCoord(SamFile& samIn, SamFile& outFile, int poolSize);
    
    void handleCoordRead(SamRecord& record, 
                         MateMapByCoord& mateMap,
                         SamCoordOutput& outputBuffer);

    // Flush up to the specified chromID/position, or if chromID is -1, flush
    // everything.
    bool coordFlush(int32_t chromID, int32_t position,
                    MateMapByCoord& mateMap, SamCoordOutput& outputBuffer);

    // Flush the first record from the mate map if there is one and flush the output buffer
    // up to and including that position (if there was nothing in the mateMap, flush everything).
    bool forceRecordFlush(MateMapByCoord& mateMap, SamCoordOutput& outputBuffer);

    // Clip overlapping reads and strands that extend past the other strand.
    void clip(SamRecord& firstRecord, SamRecord& secondRecord);

    // Clip an entire read.
    void clipEntire(SamRecord& record);

    // Calculate the average of the qualities at read positions starting at 
    // startPos and ending with endPos (included).  Stops reading if a 0 is
    // found, indicating the end of the string.
    double getAvgQual(SamRecord& record, int32_t startPos, int32_t endPos);

    String myStoreOrig;
    bool myStats;
    RunningStat myOverlaps;
    int myNumForwardClips;
    int myNumReverseClips;
    int myNumOrientationClips;

    int myNumMateFailures;
    int myNumPoolFail;
    bool myPoolSkipClip;
};

#endif
