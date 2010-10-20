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

//////////////////////////////////////////////////////////////////////////


#include "SamFilter.h"

#include "SamQuerySeqWithRefHelper.h"
#include "BaseUtilities.h"
#include "SamFlag.h"

SamFilter::FilterStatus SamFilter::clipOnMismatchThreshold(SamRecord& record, 
                                                           GenomeSequence& refSequence,
                                                           double mismatchThreshold)
{
    // Read & clip from the left & right.    
    SamQuerySeqWithRefIter iterFromStart(record, refSequence, true);
    SamQuerySeqWithRefIter iterFromEnd(record, refSequence, false);

    SamSingleBaseMatchInfo baseMatchInfo;

    int32_t readLength = record.getReadLength();
    // Init start clip to be prior to the start index (0).
    const int32_t initialStartClipPos = -1;
    int32_t startClipPos = initialStartClipPos;
    // Init end clip to be past the last index (readLength).
    int32_t endClipPos = readLength;

    bool fromStartComplete = false;
    bool fromEndComplete = false;
    int32_t numBasesFromStart = 0;
    int32_t numBasesFromEnd = 0;
    int32_t numMismatchFromStart = 0;
    int32_t numMismatchFromEnd = 0;

    FilterStatus status = NONE;

    //////////////////////////////////////////////////////////
    // Determining the clip positions.
    while(!fromStartComplete || !fromEndComplete)
    {
        // Read from the start (left to right) on the read until
        // more have been read from that direction than the opposite direction.
        while(!fromStartComplete && 
              ((numBasesFromStart <= numBasesFromEnd) ||
               (fromEndComplete)))
        {
            if(iterFromStart.getNextMatchMismatch(baseMatchInfo) == false)
            {
                // Nothing more to read in this direction.
                fromStartComplete = true;
                break;
            }
            // Got a read.  Check to see if it is to or past the last clip.
            if(baseMatchInfo.getQueryIndex() >= endClipPos)
            {
                // This base is past where we are clipping, so we
                // are done reading in this direction.
                fromStartComplete = true;
                break;
            }
            // This is an actual base read from the left to the
            // right, so up the counter and determine if it was a mismatch.
            ++numBasesFromStart;

            if(baseMatchInfo.getType() == SamSingleBaseMatchInfo::MISMATCH)
            {
                // Mismatch
                ++numMismatchFromStart;
                // Check to see if it is over the threshold.
                double mismatchPercent = 
                    (double)numMismatchFromStart / numBasesFromStart;
                if(mismatchPercent > mismatchThreshold)
                {
                    // Need to clip.
                    startClipPos = baseMatchInfo.getQueryIndex();
                    // Reset the counters.
                    numBasesFromStart = 0;
                    numMismatchFromStart = 0;
                }
            }
        }

        // Now, read from right to left until more have been read
        // from the end than from the start.
        while(!fromEndComplete && 
              ((numBasesFromEnd <= numBasesFromStart) ||
               (fromStartComplete)))
        {
            if(iterFromEnd.getNextMatchMismatch(baseMatchInfo) == false)
            {
                // Nothing more to read in this direction.
                fromEndComplete = true;
                break;
            }
            // Got a read.  Check to see if it is to or past the first clip.
            if(baseMatchInfo.getQueryIndex() <= startClipPos)
            {
                // This base is past where we are clipping, so we
                // are done reading in this direction.
                fromEndComplete = true;
                break;
            }
            // This is an actual base read from the right to the
            // left, so up the counter and determine if it was a mismatch.
            ++numBasesFromEnd;

            if(baseMatchInfo.getType() == SamSingleBaseMatchInfo::MISMATCH)
            {
                // Mismatch
                ++numMismatchFromEnd;
                // Check to see if it is over the threshold.
                double mismatchPercent = 
                    (double)numMismatchFromEnd / numBasesFromEnd;
                if(mismatchPercent > mismatchThreshold)
                {
                    // Need to clip.
                    endClipPos = baseMatchInfo.getQueryIndex();
                    // Reset the counters.
                    numBasesFromEnd = 0;
                    numMismatchFromEnd = 0;
                }
            }
        }
    }

    //////////////////////////////////////////////////////////
    // Done determining the clip positions.
    // Check to see if clipping needs to be done.
    Cigar* cigar = record.getCigarInfo();
    if((startClipPos != initialStartClipPos) || (endClipPos != readLength))
    {
        // Clipping from front and/or from the end.
        // Check to see if the entire read was clipped.
        if((startClipPos + 1) >= endClipPos)
        {
            /////////////////////////////
            // The entire read was clipped.
            filterRead(record);
            status = FILTERED;
        }
        else
        {
            // Part of the read was clipped.
            status = CLIPPED;
            
            // Need to create a new cigar.
            CigarRoller newCigar;
            
            // Loop through, creating an updated cigar.

            int origCigarOpIndex = 0;
            // Track how many read positions are covered up to this
            // point by the cigar to determine up to up to what
            // point in the cigar is affected by this clipping.
            int32_t numPositions = 0;

            // Add 1 to the startClipPos to determine the number of
            // softclips since startClipPos starts at 0.
            int32_t numSoftClips = startClipPos + 1;

            // Track if any non-clips are in the new cigar.
            bool onlyClips = true;

            const Cigar::CigarOperator* op = NULL;

            //////////////////
            // Clip from front
            // keep any hard clips that are prior to the soft clip that
            // is being performed..
            while((origCigarOpIndex < cigar->size()) &&
                  (numPositions < numSoftClips))
            {
                op = &(cigar->getOperator(origCigarOpIndex));
                switch(op->operation)
                {
                    case Cigar::hardClip:
                        // Keep this operation as the new clips do not
                        // affect other clips.
                        newCigar += *op;
                        break;
                    case Cigar::del:
                    case Cigar::skip:
                        // Skip and delete are going to be dropped, and
                        // are not in the read, so the read index doesn't
                        // need to be updated
                        break;
                    case Cigar::insert:
                    case Cigar::match:
                    case Cigar::mismatch:
                    case Cigar::softClip:
                        // Update the read index as these types
                        // are found in the read.
                        numPositions += op->count;
                        break;
                    case Cigar::none:
                    default:
                        // Nothing to do for none.
                        break;
                };
                ++origCigarOpIndex;
            }
             
            if(numSoftClips != 0)
            {
                // Add the softclip from the front of the read.
                newCigar.Add(Cigar::softClip, numSoftClips);

                // Add the rest of the last Cigar operation if
                // it is not entirely clipped.
                if(numPositions > numSoftClips)
                {
                    // Before adding it, check to see if the same
                    // operation is clipped from the end.
                    // numPositions greater than the endClipPos
                    // means that it is equal or past that position,
                    // so shorten the number of positions.
                    int32_t newCount = numPositions - numSoftClips;
                    if(numPositions > endClipPos)
                    {
                        newCount -= (numPositions - endClipPos);
                    }
                    if(newCount > 0)
                    {
                        newCigar.Add(op->operation, newCount);
                        if(!Cigar::isClip(op->operation))
                        {
                            onlyClips = false;
                        }
                    }
                }
            }

            // Add operations until the point of the end clip is reached.
            // For example...
            //   2M1D3M = MMDMMM  readLength = 5
            // readIndex: 01 234
            //   at cigarOpIndex 0 (2M), numPositions = 2.
            //   at cigarOpIndex 1 (1D), numPositions = 2.
            //   at cigarOpIndex 2 (3M), numPositions = 5.
            // if endClipPos = 2, we still want to consume the 1D, so
            // need to keep looping until numPositions > endClipPos
            while((origCigarOpIndex < cigar->size()) &&
                  (numPositions <= endClipPos))
            {
                op = &(cigar->getOperator(origCigarOpIndex));
                
                // Update the numPositions count if the operations indicates
                // bases within the read.
                if(!Cigar::foundInQuery(op->operation))
                {
                    // This operation is not in the query read sequence,
                    // so it is not yet to the endClipPos, just add the
                    // operation do not increment the number of positions.
                    newCigar += *op;
                    if(!Cigar::isClip(op->operation))
                    {
                        onlyClips = false;
                    }
                }
                else
                {
                    // This operation appears in the query sequence, so
                    // check to see if the clip occurs in this operation.

                    // endClipPos is 0 based & numPositions is a count.
                    // If endClipPos is 4, then it is the 5th position.
                    // If 4 positions are covered so far (numPositions = 4), 
                    // then we are right at endCLipPos: 4-4 = 0, none of 
                    // this operation should be kept. 
                    // If only 3 positions were covered, then we are at offset
                    // 3, so offset 3 should be added: 4-3 = 1.
                    uint32_t numPosTilClip = endClipPos - numPositions;

                    if(numPosTilClip < op->count)
                    {
                        // this operation is partially clipped, write the part
                        // that was not clipped if it is not all clipped.
                        if(numPosTilClip != 0)
                        {
                            newCigar.Add(op->operation,
                                         numPosTilClip);
                            if(!Cigar::isClip(op->operation))
                            {
                                onlyClips = false;
                            }
                        }
                    }
                    else
                    {
                        // This operation is not clipped, so add it
                        newCigar += *op;
                        if(!Cigar::isClip(op->operation))
                        {
                            onlyClips = false;
                        }
                    }
                    // This operation occurs in the query sequence, so 
                    // increment the number of positions covered.
                    numPositions += op->count;
                }

                // Move to the next cigar position.
                ++origCigarOpIndex;
            }
            
            //////////////////
            // Add the softclip to the back.
            numSoftClips = readLength - endClipPos;
            if(numSoftClips != 0)
            {
                // Add the softclip to the end
                newCigar.Add(Cigar::softClip, numSoftClips);
            }

            //////////////////
            // Add any hardclips remaining in the original cigar to the back.
            while(origCigarOpIndex < cigar->size())
            {
                op = &(cigar->getOperator(origCigarOpIndex));
                if(op->operation == Cigar::hardClip)
                {
                    // Keep this operation as the new clips do not
                    // affect other clips.
                    newCigar += *op;
                }
                ++origCigarOpIndex;
            }

            // Check to see if the new cigar is only clips.
            if(onlyClips)
            {
                // Only clips in the new cigar, so mark the read as filtered
                // instead of updating the cigar.
                /////////////////////////////
                // The entire read was clipped.
                filterRead(record);
                status = FILTERED;
            }
            else
            {
                // Part of the read was filtered, and now that we have
                // an updated cigar, update the read.
                record.setCigar(newCigar);

                // Update the starting position if a clip was added to
                // the front.
                if(startClipPos != initialStartClipPos)
                {
                    // Convert from query index to reference position (from the
                    // old cigar)
                    // Get the position for the start clipped position by
                    // getting the position associated with the clipped base on
                    // the reference.  Then add one to get to the first
                    // non-clipped position.
                    uint32_t newStartPos = 
                        cigar->getRefPosition(startClipPos, 
                                              record.get0BasedPosition()) + 1;
                    record.set0BasedPosition(newStartPos);
                }
            }
        }
    }
    return(status);
}


SamFilter::FilterStatus SamFilter::filterOnMismatchQuality(SamRecord& record, 
                                                           GenomeSequence& refSequence,
                                                           uint32_t qualityThreshold, 
                                                           uint8_t defaultQualityInt)
{
    uint32_t totalMismatchQuality = 
        sumMismatchQuality(record, refSequence, defaultQualityInt); 
    
    // If the total mismatch quality is over the threshold, 
    // filter the read.
    if(totalMismatchQuality > qualityThreshold)
    {
        filterRead(record);
        return(FILTERED);
    }
    return(NONE);
}


// NOTE: Only positions where the reference and read both have bases that
//       are different and not 'N' are considered mismatches.
uint32_t SamFilter::sumMismatchQuality(SamRecord& record, 
                                       GenomeSequence& refSequence,
                                       uint8_t defaultQualityInt)
{
    // Track the mismatch info.
    int mismatchQual = 0;
    int numMismatch = 0;

    SamQuerySeqWithRefIter sequenceIter(record, refSequence);

    SamSingleBaseMatchInfo baseMatchInfo;
    while(sequenceIter.getNextMatchMismatch(baseMatchInfo))
    {
        if(baseMatchInfo.getType() == SamSingleBaseMatchInfo::MISMATCH)
        {
            // Got a mismatch, get the associated quality.
            char readQualityChar = 
                record.getQuality(baseMatchInfo.getQueryIndex());
            uint8_t readQualityInt = 
                BaseUtilities::getPhredBaseQuality(readQualityChar);
            
            if(readQualityInt == BaseUtilities::UNKNOWN_QUALITY_INT)
            {
                // Quality was not specified, so use the configured setting.
                readQualityInt = defaultQualityInt;
            }
            mismatchQual += readQualityInt;
            ++numMismatch;
        }
    }

    return(mismatchQual);
}


void SamFilter::filterRead(SamRecord& record)
{
    // Filter the read by marking it as unmapped.
    uint16_t flag = record.getFlag(); 
    SamFlag::setUnmapped(flag);
    record.setFlag(flag);
}
