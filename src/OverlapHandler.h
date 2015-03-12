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
// Base Class Paired Read Overlap Handler.

#ifndef __OVERLAP_HANDLER_H__
#define __OVERLAP_HANDLER_H__

#include "SamRecord.h"
#include "SimpleStats.h"

class OverlapHandler
{
public:
    /// Enum for indicating whether or not there is an overlap based on
    /// a single read.  WRONG/unexpected orientation means that the reverse
    /// read starts before its mate or the forward read starts after its mate.
    enum OverlapInfo {
        NO_OVERLAP, ///< does not overlap its mate if it has one
        OVERLAP,    ///< read overlaps its mate
        NO_OVERLAP_WRONG_ORIENT,  ///< does not overlap its mate, but orientation is wrong
        SAME_START, ///< read & its mate start at the same position
        UNKNOWN_OVERLAP, ///< unkown if the read overlaps its mate
        UNKNOWN_OVERLAP_WRONG_ORIENT ///< unkown if the read overlaps its mate, but the orientation is unexpected
    };

    OverlapHandler() {init();}

    void keepStats(bool keepStats) {myStats = keepStats;}
    virtual void printStats();

    // Enable marking as unmapped if an entire read would be clipped
    // instead of the default option of marking the whole thing as clipped.
    void markAsUnmapped()
    { myUnmap = true; }

    // Specify storeOrigTag to be non-blank to store the original cigar
    // in the record under the specified tag.
    void storeOrigCigar(const String& storeOrigTag)
    { myStoreOrigCigar = storeOrigTag; }

    /// Get information on whether or not this record has an overlap
    /// with its mate (if it has one) (may be unknown if this is
    /// the 2nd read in the pair).
    OverlapInfo getOverlapInfo(SamRecord& record, uint16_t intExcludeFlags = 0);

    virtual void handleOverlapPair(SamRecord& firstRecord,
                                   SamRecord& secondRecord) = 0;

    /// Handle the case where the pair does not overlap, but the
    /// reverse strand is completely prior to the forward strand.
    /// updateStats specifies whether or not stats should be updated.
    /// They should only be updated once per pair.
    /// mateUnmapped indicates whether or not the mate unmapped flag
    /// should also be set.
    virtual void handleNoOverlapWrongOrientation(SamRecord& record,
                                                 bool updateStats = true,
                                                 bool mateUnmapped = true);

    /// Handle the case where the mate was not found (maybe due to
    /// hitting a buffer limit).  updateStats specifies whether or not
    /// stats should be updated. They should only be updated once per pair.
    /// Return whether or not the record was modified (not all handlers
    /// will change the record).  By default, just return false,
    /// no modifications.
    virtual bool handleOverlapWithoutMate(SamRecord& record,
                                          bool updateStats = true);

    virtual int numSteps() {return(1);}

    void incrementStep() { ++myStepNum; }

protected:
    // Will mark the mate as unmapped if myUnmap is set to true.
    void markMateUnmapped(SamRecord& record);

    int myStepNum;
    String myStoreOrigCigar;
    bool myStats;
    RunningStat myOverlaps;
    int myNumOrientationClips;
    // If set to true, mark entirely clipped records as unmapped.
    bool myUnmap;

private:
    void init()
    {
        myStepNum = 1;
        myStoreOrigCigar = "";
        myStats = false;
        myNumOrientationClips = 0;
        myUnmap = false;
    }
};

#endif
