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

