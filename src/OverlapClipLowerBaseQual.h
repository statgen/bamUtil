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
// Paired Read Overlap Handler that clips the lower base quality read.

#ifndef __OVERLAP_CLIP_LOWER_BASE_QUAL_H__
#define __OVERLAP_CLIP_LOWER_BASE_QUAL_H__

#include "OverlapHandler.h"

class OverlapClipLowerBaseQual : public OverlapHandler
{
public:
    OverlapClipLowerBaseQual();

    virtual void printStats();

    virtual void handleOverlapPair(SamRecord& firstRecord,
                                   SamRecord& secondRecord);

    /// Handle the case where the mate was not found (maybe due to
    /// hitting a buffer limit).  Just clip this record at the overlap point
    /// or clip the entire read if it is the reverse strand and the mate is not.
    /// updateStats if set to update them (turn off for the 2nd read in the
    /// pair).  Return true, the record was updated.
    virtual bool handleOverlapWithoutMate(SamRecord& record, 
                                          bool updateStats = true);

private:
    // Calculate the average of the qualities at read positions starting at 
    // startPos and ending with endPos (included).  Stops reading if a 0 is
    // found, indicating the end of the string.
    double getAvgQual(SamRecord& record, int32_t startPos, int32_t endPos);


    int myNumForwardClips;
    int myNumReverseClips;
};

#endif
