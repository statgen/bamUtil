/*
 *  Copyright (C) 2010  Regents of the University of Michigan
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

#ifndef __SAM_FILTER_H__
#define __SAM_FILTER_H__

#include "SamRecord.h"
#include "GenomeSequence.h"

class SamFilter
{
public:
    enum FilterStatus {
        NONE, // The filter did not affect the read.
        CLIPPED, // Filtering clipped the read.
        FILTERED // Filtering caused the read to be modified to unmapped.
    };

    // Clip the read based on the specified mismatch threshold.
    // Returns how the read was affected, 
    //     NONE if the read was not modified,
    //     CLIPPED if the read was clipped,
    //     FILTERED if the whole read would have been clipped so instead the
    //              read was modified to unmapped.
    static FilterStatus clipOnMismatchThreshold(SamRecord& record, 
                                                GenomeSequence& refSequence,
                                                double mismatchThreshold);

    /// Soft clip the record from the front and/or the back.
    /// \param record record to be clipped (input/output parameter).
    /// \param numFrontClips number of bases that should be clipped from the
    /// front of the sequence read.  (total count, including any that are
    /// already clipped.)
    /// \param backClipPos number of bases that should be clipped from the
    /// back of the sequence read.  (total count, including any that are
    /// already clipped.)
    static FilterStatus softClip(SamRecord& record,
                                 int32_t numFrontClips, 
                                 int32_t numBackClips);

    /// Soft clip the cigar from the front and/or the back, writing the value
    /// into the new cigar, updatedCigar & startPos are only updated if
    /// the return FilterStatus is CLIPPED.
    /// \param oldCigar cigar prior to clipping
    /// \param numFrontClips number of bases that should be clipped from the
    /// front of the sequence read.  (total count, including any that are
    /// already clipped.)
    /// \param numBackClips number of bases that should be clipped from the
    /// back of the sequence read.  (total count, including any that are
    /// already clipped.)
    /// \param startPos 0-based start position associated with the
    /// cigar prior to updating (input) and set to the 0-based start position
    /// after updating (output) the cigar if it was CLIPPED.
    /// \param updatedCigar set to the clipped cigar if CLIPPED (output param).
    static FilterStatus softClip(Cigar& oldCigar, 
                                 int32_t numFrontClips,
                                 int32_t numBackClips,
                                 int32_t& startPos,
                                 CigarRoller& updatedCigar);

    // Filter the read based on the specified quality threshold.
    // Returns how the read was affected, 
    //     NONE if the read was not modified,
    //     FILTERED if the read was modified to unmapped because it was over
    //              the quality threshold.
    static FilterStatus filterOnMismatchQuality(SamRecord& record,
                                                GenomeSequence& refSequence,
                                                uint32_t qualityThreshold, 
                                                uint8_t defaultQualityInt);
    
    static uint32_t sumMismatchQuality(SamRecord& record, 
                                       GenomeSequence& refSequence,
                                       uint8_t defaultQualityInt);

    // Filter the read out (mark it as unmapped.
    static void filterRead(SamRecord& record);
};

#endif

