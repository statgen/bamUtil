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

#include <stdexcept>

#include "SamQuerySeqWithRefHelper.h"
#include "BaseUtilities.h"
#include "SamFlag.h"

SamQuerySeqWithRefIter::SamQuerySeqWithRefIter(SamRecord& record,
                                               GenomeSequence& refSequence,
                                               bool forward)
    : myRecord(record),
      myRefSequence(refSequence),
      myCigar(NULL),
      myStartOfReadOnRefIndex(0),
      myQueryIndex(0),
      myForward(forward)
{
    myCigar = myRecord.getCigarInfo();   
    myStartOfReadOnRefIndex = 
        myRefSequence.getGenomePosition(myRecord.getReferenceName()) +
        myRecord.get0BasedPosition();

    if(!forward)
    {
        myQueryIndex = myRecord.getReadLength() - 1;
    }
}


SamQuerySeqWithRefIter::~SamQuerySeqWithRefIter()
{
}



bool SamQuerySeqWithRefIter::reset(bool forward)
{
    myCigar = myRecord.getCigarInfo();
    if(myCigar == NULL)
    {
        // Failed to get Cigar.
        return(false);
    }

    // Get where the position of where this read starts as mapped to the 
    // reference.
    myStartOfReadOnRefIndex = 
        myRefSequence.getGenomePosition(myRecord.getReferenceName()) +
        myRecord.get0BasedPosition();

    forward = myForward;

    if(myForward)
    {
        myQueryIndex = 0;
    }
    else
    {
        // reverse, so start at the last entry.
        myQueryIndex = myRecord.getReadLength() - 1;
    }
    return(true);
}


// Returns information for the next position where the query and the 
// reference match or mismatch.  To be a match or mismatch, both the query
// and reference must have a base that is not 'N'.
// This means:
//    insertions and deletions are not mismatches or matches.
//    'N' bases are not matches or mismatches
// NOTE: Only returns positions where both the query and the reference
// Returns true if an entry was found, false if there are no more matches or
// mismatches.
bool SamQuerySeqWithRefIter::getNextMatchMismatch(SamSingleBaseMatchInfo& matchMismatchInfo)
{
    // Check whether or not this read is mapped. 
    // If the read is not mapped, return no matches.
    if(!SamFlag::isMapped(myRecord.getFlag()))
    {
        // Not mapped.
        return(false);
    }

    // Check that the Cigar is set.
    if(myCigar == NULL)
    {
        // Error.
        throw(std::runtime_error("Cannot determine matches/mismatches since failed to retrieve the cigar"));
        return(false);
    }

    // Repull the read length from the record to check just in case the
    // record has changed length.
    // Loop until a match or mismatch is found as long as query index
    // is still within the read  (Loop is broken by a return.
    while((myQueryIndex < myRecord.getReadLength()) &&
          (myQueryIndex >= 0))
    {
        // Still more bases, look for a match/mismatch.

        // Get the reference offset for this read position.
        int32_t refOffset = myCigar->getRefOffset(myQueryIndex);
        if(refOffset == Cigar::INDEX_NA)
        {
            // This is either a softclip or an insertion
            // which do not count as a match or a mismatch, so
            // go to the next index.
            nextIndex();
            continue;
        }
        
        // Both the reference and the read have a base, so get the bases.
        char readBase = myRecord.getSequence(myQueryIndex);
        char refBase = myRefSequence[myStartOfReadOnRefIndex + refOffset];
       
        // If either the read or the reference bases are unknown, then
        // it does not count as a match or a mismatch.
        if(BaseUtilities::isAmbiguous(readBase) || 
           BaseUtilities::isAmbiguous(refBase))
        {
            // Either the reference base or the read base are unknown,
            // so skip this position.
            nextIndex();
            continue;
        }

        // Both the read & the reference have a known base, so it is either
        // a match or a mismatch.
        matchMismatchInfo.setQueryIndex(myQueryIndex);

        // Check if they are equal.
        if(BaseUtilities::areEqual(readBase, refBase))
        {
            // Match.
            matchMismatchInfo.setType(SamSingleBaseMatchInfo::MATCH);
            // Increment the query index to the next position.
            nextIndex();
            return(true);
        }
        else
        {
            // Mismatch
            matchMismatchInfo.setType(SamSingleBaseMatchInfo::MISMATCH);
            // Increment the query index to the next position.
            nextIndex();
            return(true);
        }
    }

    // No matches or mismatches were found, so return false.
    return(false);
}


void SamQuerySeqWithRefIter::nextIndex()
{
    if(myForward)
    {
        ++myQueryIndex;
    }
    else
    {
        --myQueryIndex;
    }
}


SamSingleBaseMatchInfo::SamSingleBaseMatchInfo()
    : myType(UNKNOWN),
      myQueryIndex(0)
{
}


SamSingleBaseMatchInfo::~SamSingleBaseMatchInfo()
{
}


SamSingleBaseMatchInfo::Type SamSingleBaseMatchInfo::getType()
{
    return(myType);
}


int32_t SamSingleBaseMatchInfo::getQueryIndex()
{
    return(myQueryIndex);
}


void SamSingleBaseMatchInfo::setType(Type newType)
{
    myType = newType;
}


void SamSingleBaseMatchInfo::setQueryIndex(int32_t queryIndex)
{
    myQueryIndex = queryIndex;
}
