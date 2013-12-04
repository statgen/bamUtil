/*
 *  Copyright (C) 2011-2012  Regents of the University of Michigan
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
#include "OverlapHandler.h"

class ClipOverlap : public BamExecutable
{
public:
    ClipOverlap();

    static void clipOverlapDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);
    virtual const char* getProgramName() {return("bam:clipOverlap");}

private:
    static const int DEFAULT_POOL_SIZE = 1000000;

    SamStatus::Status handleSortedByReadName(SamFile& samIn,
                                             SamFile* outFile);

    SamStatus::Status handleSortedByCoord(SamFile& samIn,
                                          SamCoordOutput* outputBufferPtr);
    
    ///////////////////////////////////////////////////////////////////
    // Methods to handle Coordinate Specific Processing.
    
    // Helper method to get a record ptr if needed and read the record.
    // returns success or the failure reason.
    SamStatus::Status readCoordRecord(SamFile& samIn,
                                      SamRecord** recordPtr,
                                      MateMapByCoord& mateMap, 
                                      SamCoordOutput* outputBufferPtr);

    // Flush the first record from the mate map if there is one and flush the output buffer
    // up to and including that position (if there was nothing in the mateMap, flush everything).
    bool forceRecordFlush(MateMapByCoord& mateMap, 
                          SamCoordOutput* outputBufferPtr);

    // Flush up to the first record in the mate map, or if it is empty,
    // flush up to and including the specified position.
    bool flushOutputBuffer(MateMapByCoord& mateMap,
                           SamCoordOutput& outputBuffer,
                           int32_t prevChrom,
                           int32_t prevPos);

    // Cleanup the mate map up to the specified record. 
    // If chrom is -1, empty the entire mateMap.
    void cleanupMateMap(MateMapByCoord& mateMap,
                        SamCoordOutput* outputBufferPtr,
                        int32_t chrom = -1, int32_t position = -1);

    ///////////////////////////////////////////////////////////////////
    // Private Member Data

    SamFileHeader mySamHeader;
    OverlapHandler* myOverlapHandler;
    SamRecordPool myPool; // used just for coord reads.
    bool myOverlapsOnly;
    uint16_t myIntExcludeFlags;

    uint32_t myNumMateFailures;
    uint32_t myNumPoolFail;
    uint32_t myNumPoolFailNoHandle;
    uint32_t myNumPoolFailHandled;
    uint32_t myNumOutOfOrder;
    bool myPoolSkipOverlap;
};

#endif
