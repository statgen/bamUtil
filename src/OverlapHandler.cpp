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
#include "OverlapHandler.h"
#include "SamFlag.h"
#include "CigarHelper.h"
#include "SamFilter.h"

void OverlapHandler::printStats()
{
    if(myStats)
    {
        std::cerr << "Overlap Statistics:" << std::endl;
        std::cerr << "Number of overlapping pairs: "
                  << myOverlaps.NumDataValues() << std::endl
                  << "Average # Reference Bases Overlapped: "
                  << myOverlaps.Mean() << std::endl
                  << "Variance of Reference Bases overlapped: "
                  << myOverlaps.Variance() << std::endl
                  << "Number of times orientation causes additional clipping: "
                  << myNumOrientationClips << std::endl;
    }
}


OverlapHandler::OverlapInfo OverlapHandler::getOverlapInfo(SamRecord& record, uint16_t intExcludeFlags)
{
    // Determine whether or not the reads overlap.
    int16_t flag = record.getFlag();
    // No overlap if:
    //  1) the read is not paired.
    //  2) read and its mate are on different chromosome ids
    //  3) read is unmapped
    //  4) mate is unmapped
    //  5) read contains specified exclusion flags.
    if(!SamFlag::isPaired(flag) || 
       (record.getMateReferenceID() != record.getReferenceID()) ||
       !SamFlag::isMapped(flag) || !SamFlag::isMateMapped(flag) ||
       ((flag & intExcludeFlags) != 0) )
    {
        return(NO_OVERLAP);
    }

    // Same chromosome and both reads are mapped,
    // so look at positions.
    int32_t readStart = record.get0BasedPosition();
    int32_t mateStart = record.get0BasedMatePosition();

    // If either position is unknown (-1), then just add
    // the record to the output buffer because no clipping
    // needs to be done.
    if((readStart == -1) || (mateStart == -1))
    {
        return(NO_OVERLAP);
    }

    if(readStart > mateStart)
    {
        // This read comes after the other, so can't tell if clipping
        // But check if this is forward and the mate is reverse.
        if(!SamFlag::isReverse(flag) && SamFlag::isMateReverse(flag))
        {
            return(UNKNOWN_OVERLAP_WRONG_ORIENT);
        }
        return(UNKNOWN_OVERLAP);
    }
    else if(readStart == mateStart)
    {
        // Start at the same spot.
        return(SAME_START);
    }
    
    // This read starts before the other.
    // No clipping if it ends before the other.
    int32_t readEnd = record.get0BasedAlignmentEnd();
    if(readEnd < mateStart)
    {
        // This read finishes before the mate starts so there is no overlap.
        // If this read is the reverse and the other read is the forward
        // strand, then they completely passed each other and both should
        // be clipped.
        if(SamFlag::isReverse(flag) && !SamFlag::isMateReverse(flag))
        {
            return(NO_OVERLAP_WRONG_ORIENT);
        }
        // First, but no clipping
        return(NO_OVERLAP);
    }

    // The reads overlap.
    return(OVERLAP);
}


void OverlapHandler::handleNoOverlapWrongOrientation(SamRecord& record,
                                                     bool updateStats,
                                                     bool mateUnmapped)
{
    static CigarRoller newCigar; // holds updated cigar.

    if(myStats && updateStats)
    {
        ++myNumOrientationClips;
    }

    if(!myUnmap)
    {
        // Clip the entire record.
        if(CigarHelper::softClipEndByRefPos(record, record.get0BasedPosition(),
                                            newCigar) != CigarHelper::NO_CLIP)
        { 
            // Write the original cigar into the specified tag.
            if(!myStoreOrigCigar.IsEmpty())
            {
                // Write original cigar.
                record.addTag(myStoreOrigCigar, 'Z', record.getCigar());
            }
            // Update the cigar.
            record.setCigar(newCigar);
        }
    }
    else
    {
        // Write the original cigar into the specified tag.
        if(!myStoreOrigCigar.IsEmpty())
        {
            // Write original cigar.
            record.addTag(myStoreOrigCigar, 'Z', record.getCigar());
        }
        // Filter read - mark it unmapped.
        SamFilter::filterRead(record);
        // If the mate will be unmapped, update the mate information.
        if(mateUnmapped)
        {
            markMateUnmapped(record);
        }
    }
}


bool OverlapHandler::handleOverlapWithoutMate(SamRecord& record,
                                              bool updateStats)
{
    return(false);
}


void OverlapHandler::markMateUnmapped(SamRecord& record)
{
    // Mark the mate as unmapped if we are doing unmapping
    if(myUnmap)
    {
        uint16_t flag = record.getFlag(); 
        flag |= SamFlag::MATE_UNMAPPED;
        // Also ensure proper pair is not set.
        flag &= ~SamFlag::PROPER_PAIR;
        record.setFlag(flag);
    }
}
