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

#include <assert.h>
#include "OverlapSplitClip.h"
#include "CigarHelper.h"
#include "SamFlag.h"
#include "SamFilter.h"

OverlapSplitClip::OverlapSplitClip()
    : OverlapHandler(),
      myNumForwardClips(0),
      myNumReverseClips(0)
{
}


void OverlapSplitClip:: printStats()
{
    if(myStats)
    {
        OverlapHandler::printStats();
        std::cerr << "Number of times the forward strand was clipped: "
                  << myNumForwardClips << std::endl
                  << "Number of times the reverse strand was clipped: "
                  << myNumReverseClips << std::endl;
    }
}


void OverlapSplitClip::handleOverlapPair(SamRecord& firstRecord,
                                         SamRecord& secondRecord)
{
    // Clip from 2nd record start to 1st record end
    int32_t overlapStart = secondRecord.get0BasedPosition();
    int32_t overlapEnd = firstRecord.get0BasedAlignmentEnd();
    int32_t secondEnd = secondRecord.get0BasedAlignmentEnd();
    uint16_t firstFlag = firstRecord.getFlag();
    uint16_t secondFlag = secondRecord.getFlag();

    // For now, just assert that first is forward & 2nd is reverse
    assert(!SamFlag::isReverse(firstFlag));
    assert(SamFlag::isReverse(secondFlag));
    
    // For now 2nd read must end last.
    assert(secondEnd>= overlapEnd);

    // Calculate the overlap length
    // (example: positions 3, 4, 5, 6, 7.  So 7-3+1=5 positions.)
    int32_t overlapLen = overlapEnd - overlapStart + 1;

    // If we are keeping stats, update the stats with the overlap length.
    if(myStats)
    {
        myOverlaps.Push(overlapLen);
    }

    // Clip half from 1st read & half from 2nd read.
    int32_t halfOverlap = overlapLen/2;
    // If the overlap is odd, clips will be different by 1.
    int32_t remainingOverlap = overlapLen - halfOverlap;
    
    static CigarRoller newFirstCigar; // holds updated cigar.
    static CigarRoller newSecondCigar; // holds updated cigar.

    // Clip the first record starting at the overlap start plus 1/2 the overlap
    // (example: overlap, positions 3,4,5,6,7.  halfOverlap is 2, so clip the
    // end starting at position 5 (3+2).
    int32_t firstClipPos = 
        CigarHelper::softClipEndByRefPos(firstRecord,
                                         overlapStart + halfOverlap, 
                                         newFirstCigar);

    // Clip the second record starting at the overlap end minus 1/2 the overlap
    // (example: overlap, positions 3,4,5,6,7.  halfOverlap is 2, so clip from
    // the beginning until overlapEnd - remainingOverlap (7-3)=4
    // Stop clipping at position 4.
    int32_t newSecondStartPos = 0;
    int32_t secondClipPos =
        CigarHelper::softClipBeginByRefPos(secondRecord,
                                           overlapEnd - remainingOverlap,
                                           newSecondCigar, newSecondStartPos);

    if(myStats)
    {
        // Since we are only handing where 1st is forward and 2nd is reverse
        // increment both types of clips
        ++myNumReverseClips;
        ++myNumForwardClips;
    }

    // Clip the overlaps.
    // Write the original cigars into the specified tag.
    if(!myStoreOrigCigar.IsEmpty())
    {
        // Write original cigar.
        firstRecord.addTag(myStoreOrigCigar, 'Z', 
                           firstRecord.getCigar());
        // Write original cigar.
        secondRecord.addTag(myStoreOrigCigar, 'Z', 
                            secondRecord.getCigar());
    }

    // Check if the entire first record will be clipped and should
    // instead be marked as unmapped.
    if(myUnmap && (firstClipPos == 0))
    {
        // Completely clipped, mark as unmapped
        SamFilter::filterRead(firstRecord);
        // Update the mate to indicate this record is unmapped
        markMateUnmapped(secondRecord);
    }
    else
    {
        // Write new Cigar
        firstRecord.setCigar(newFirstCigar);
    }

    // Check if the entire second record will be clipped and should
    // instead be marked as unmapped.
    if(myUnmap && (secondClipPos >= (secondRecord.getReadLength()-1)))
    {
        // Completely clipped, mark as unmapped
        SamFilter::filterRead(secondRecord);
        // Update the mate to indicate this record is unmapped
        markMateUnmapped(firstRecord);
    }
    else
    {
        // Update 2nd record's starting position
        secondRecord.set0BasedPosition(newSecondStartPos);
        firstRecord.set0BasedMatePosition(newSecondStartPos);
        secondRecord.setCigar(newSecondCigar);
    }
}


bool OverlapSplitClip::handleOverlapWithoutMate(SamRecord& record,
                                                bool updateStats)
{
    // Don't handle this right now
    assert(false);
    static CigarRoller newCigar;

    // Check if this is reverse and the mate is not.
    int16_t flag = record.getFlag();
    if(SamFlag::isReverse(flag) && !SamFlag::isMateReverse(flag))
    {
        if(myStats && updateStats)
        {
            ++myNumReverseClips;
            myOverlaps.Push(record.get0BasedAlignmentEnd() - 
                            record.get0BasedMatePosition() + 1);
        }
        handleNoOverlapWrongOrientation(record, updateStats, false);
    }
    else
    {
        // Just clip this one based on the overlap.
        int32_t clipPos = CigarHelper::softClipEndByRefPos(record, 
                                                           record.get0BasedMatePosition(), 
                                                           newCigar);
        if(clipPos != CigarHelper::NO_CLIP)
        {
            // Update the number of clips.
            if(myStats && updateStats)
            {
                // Clipping the first record's strand.
                if(SamFlag::isReverse(flag))
                {
                    ++myNumReverseClips;
                }
                else
                {
                    ++myNumForwardClips;
                }
                // Update the overlap length - the difference
                // between the mate start and this read's end.
                myOverlaps.Push(record.get0BasedAlignmentEnd() -
                                record.get0BasedMatePosition() + 1);
            }
            // Write the original cigar into the specified tag.
            if(!myStoreOrigCigar.IsEmpty())
            {
                // Write original cigar.
                record.addTag(myStoreOrigCigar, 'Z', record.getCigar());
            }

            if(myUnmap && (clipPos == 0))
            {
                // whole read got clipped
                SamFilter::filterRead(record);
            }
            else
            {
                record.setCigar(newCigar);
            }
        }
    }
    return(true);
}



