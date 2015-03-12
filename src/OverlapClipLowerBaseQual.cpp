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

#include "OverlapClipLowerBaseQual.h"
#include "CigarHelper.h"
#include "SamFlag.h"
#include "SamFilter.h"

OverlapClipLowerBaseQual::OverlapClipLowerBaseQual()
    : OverlapHandler(),
      myNumForwardClips(0),
      myNumReverseClips(0)
{
}


void OverlapClipLowerBaseQual:: printStats()
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


void OverlapClipLowerBaseQual::handleOverlapPair(SamRecord& firstRecord,
                                                 SamRecord& secondRecord)
{
    // Clip from 2nd record start to 1st record end
    int32_t overlapStart = secondRecord.get0BasedPosition();
    int32_t overlapEnd = firstRecord.get0BasedAlignmentEnd();

    uint16_t firstFlag = firstRecord.getFlag();
    uint16_t secondFlag = secondRecord.getFlag();

    // If we are keeping stats, update the stats with the overlap length.
    if(myStats)
    {
        // Both overlapStart & overlapEnd are 0-based inclusive.
        // Check to see which one ends last.
        int32_t secondEnd = secondRecord.get0BasedAlignmentEnd();
        int32_t end = overlapEnd;
        if(end > secondEnd)
        {
            // The first read ends after the 2nd one, so the
            // overlap ends at the end of the 2nd read.
            end = secondEnd;
        }
        // If overlapEnd = 5 & overlapStart = 3, they overlap at
        // positions 3, 4, & 5.  So 5-3+1=3 positions.
        myOverlaps.Push(end - overlapStart + 1);
    }

    static CigarRoller newFirstCigar; // holds updated cigar.
    static CigarRoller newSecondCigar; // holds updated cigar.

    // Determine which record will get clipped by determining
    // which record has a lower base quality in the overlapping region.
    // First check clipping the first region.
    // Determine the clipping on the 1st record at the start of the 2nd.
    int32_t firstClipPos = 
        CigarHelper::softClipEndByRefPos(firstRecord, overlapStart, 
                                         newFirstCigar);
    // Get the Quality of the clipped bases.
    // They run from the overlap start to the length of the read.
    double firstQualAvg = getAvgQual(firstRecord, firstClipPos, 
                                     firstRecord.getReadLength()-1);

    int32_t newPos = 0;
    int32_t secondClipPos =
        CigarHelper::softClipBeginByRefPos(secondRecord, overlapEnd,
                                           newSecondCigar, newPos);
    // Get the Quality of the clipped bases.
    // They run from the beginning until the overlap end(included).
    double secondQualAvg = getAvgQual(secondRecord, 0, secondClipPos);
        
    // Check to see whether the 1st record or the 2nd one should be clipped
    // based on which has the lower quality.
    if(firstQualAvg <= secondQualAvg)
    {
        // First read has lower or equal quality, so clip that.
        if(myStats)
        {
            if(SamFlag::isReverse(firstFlag))
            {
                ++myNumReverseClips;
            }
            else
            {
                ++myNumForwardClips;
            }
        }
        
        // Check to see if the entire read should be clipped by
        // checking forward/reverse.
        if(SamFlag::isReverse(firstFlag) && !SamFlag::isReverse(secondFlag))
        {
            // first record is reverse and 2nd is forward, so clip
            // those extending ends which means clipping all of the
            // reverse(first) strand since it's overlap has lower quality.
            if(myStats)
            {
                // Increment that we are doing an orientation clip.
                ++myNumOrientationClips;
            }

            // This will clip the entire read, do not update the stats.
            handleOverlapWithoutMate(firstRecord, false);
            markMateUnmapped(secondRecord);

            // Soft clip the end of the forward(second) strand that
            // extends past the reverse strand (overlapEnd + 1).
            if(CigarHelper::softClipEndByRefPos(secondRecord, overlapEnd+1,
                                                newSecondCigar)
               != CigarHelper::NO_CLIP)
            {
                // Write the original cigar into the specified tag.
                if(!myStoreOrigCigar.IsEmpty())
                {
                    // Write original cigar.
                    secondRecord.addTag(myStoreOrigCigar, 'Z', 
                                        secondRecord.getCigar());
                }
                secondRecord.setCigar(newSecondCigar);
            }
        }
        else
        {
            // No strand specific extended ends clipping is required,
            // so just clip the overlap.
            // Write the original cigar into the specified tag.
            if(!myStoreOrigCigar.IsEmpty())
            {
                // Write original cigar.
                firstRecord.addTag(myStoreOrigCigar, 'Z', 
                                   firstRecord.getCigar());
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
                firstRecord.setCigar(newFirstCigar);
            }
        }
    }
    else
    {
        // The 2nd clip has lower quality, so clip that.
        if(myStats)
        {
            if(SamFlag::isReverse(secondFlag))
            {
                ++myNumReverseClips;
            }
            else
            {
                ++myNumForwardClips;
            }
        }
        // Check to see if the entire read should be clipped by
        // checking forward/reverse.
        if(SamFlag::isReverse(firstFlag) && !SamFlag::isReverse(secondFlag))
        {
            // first record is reverse and 2nd is forward, so clip
            // those extending ends which means clipping all of the
            // forward(second) strand since it's overlap has lower quality.
            if(myStats)
            {
                // Increment that we are doing an orientation clip.
                ++myNumOrientationClips;
            }

            // This will clip the entire read
            handleOverlapWithoutMate(secondRecord, false);
            markMateUnmapped(firstRecord);

            // Soft clip the front of the reverse(first) strand that
            // extends past the forward strand (overlapStart - 1).
            if(CigarHelper::softClipBeginByRefPos(firstRecord, 
                                                  overlapStart-1,
                                                  newFirstCigar,
                                                  newPos)
               != CigarHelper::NO_CLIP)
            {
                // Write the original cigar into the specified tag.
                if(!myStoreOrigCigar.IsEmpty())
                {
                    // Write original cigar.
                    firstRecord.addTag(myStoreOrigCigar, 'Z', 
                                       firstRecord.getCigar());
                }
                firstRecord.setCigar(newFirstCigar);
                secondRecord.set0BasedMatePosition(newPos);
                firstRecord.set0BasedPosition(newPos);
                secondRecord.set0BasedMatePosition(newPos);
            }
        }
        else
        {
            // No strand specific extended ends clipping is required,
            // so just clip the overlap.
            // Write the original cigar into the specified tag.
            if(!myStoreOrigCigar.IsEmpty())
            {
                // Write original cigar.
                secondRecord.addTag(myStoreOrigCigar, 'Z', 
                                    secondRecord.getCigar());
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
                secondRecord.set0BasedPosition(newPos);
                firstRecord.set0BasedMatePosition(newPos);
                secondRecord.setCigar(newSecondCigar);
            }
        }
    }
}


bool OverlapClipLowerBaseQual::handleOverlapWithoutMate(SamRecord& record,
                                                        bool updateStats)
{
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


double OverlapClipLowerBaseQual::getAvgQual(SamRecord& record, 
                                            int32_t startPos, 
                                            int32_t endPos)
{
    int32_t qualSum = 0;
    const char* quality = record.getQuality();
    int numVals = 0;

    // Check for invalid start position.
    if((startPos < 0) || ((uint32_t)startPos > strlen(quality)))
    {
        // Invalid start position, just return 0.
        return(0);
    }

    for(int i = startPos; i <= endPos; i++)
    {
        // Check for the null terminator at the end of the quality string.
        if(quality[i] == 0)
        {
            // Qual is shorter than the read, so break.
            break;
        }
        qualSum += quality[i];
        // increment the number of values.
        ++numVals;
    }
    if(numVals != 0)
    {
        return(qualSum/(double)numVals);
    }
    else
    {
        return(0);
    }
}

