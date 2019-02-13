/*
 *  Copyright (C) 2019  Regents of the University of Michigan
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

#ifndef __OVERLAP_SPLIT_CLIP_H__
#define __OVERLAP_SPLIT_CLIP_H__

#include "OverlapHandler.h"

class OverlapSplitClip : public OverlapHandler
{
public:
    OverlapSplitClip();

    virtual void printStats();

    virtual void handleOverlapPair(SamRecord& firstRecord,
                                   SamRecord& secondRecord);

    /// Handle the case where the pair does not overlap, but the
    /// reverse strand is completely prior to the forward strand.
    /// For this implementation, don't do anything, keep both reads.
    virtual void handleNoOverlapWrongOrientation(SamRecord& record,
                                                 bool updateStats = true,
                                                 bool mateUnmapped = true);

    /// Handle the case where the mate was not found (maybe due to
    /// hitting a buffer limit).  Just clip this record at the overlap point
    /// or clip the entire read if it is the reverse strand and the mate is not.
    /// updateStats if set to update them (turn off for the 2nd read in the
    /// pair).  Return true, the record was updated.
    virtual bool handleOverlapWithoutMate(SamRecord& record, 
                                          bool updateStats = true);

private:

    int myNumForwardClips;
    int myNumReverseClips;
};

#endif
