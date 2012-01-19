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
#include "PileupElementBaseInfo.h"
#include "SamFlag.h"
#include <stdexcept>

/////////////////////////////////////////////////////////////////////////////
//
// PileupElementBaseInfo
//

IFILE PileupElementBaseInfo::ourOutputFile = NULL;
int PileupElementBaseInfo::ourGapSize = 0;
bool PileupElementBaseInfo::ourIgnoreDeletion = false;

int PileupElementBaseInfo::ourPrevPos = PileupElement::UNSET_POSITION;
int PileupElementBaseInfo::ourPrevChromID = INVALID_CHROMOSOME_INDEX;

void PileupElementBaseInfo::setOutputFile(const char* outputFile)
{
    ourOutputFile = ifopen(outputFile, "w");
}

void PileupElementBaseInfo::closeOutputFile()
{
    ifclose(ourOutputFile);
}


void PileupElementBaseInfo::setGapSize(int gapSize)
{
    ourGapSize = gapSize;
}


void PileupElementBaseInfo::setIgnoreDeletion(bool ignoreDeletion)
{
    ourIgnoreDeletion = ignoreDeletion;
}


void PileupElementBaseInfo::printHeader()
{
}


PileupElementBaseInfo::PileupElementBaseInfo()
    : PileupElement()
{
    initVars();
}


PileupElementBaseInfo::~PileupElementBaseInfo()
{
}


// Add an entry to this pileup element.  
void PileupElementBaseInfo::addEntry(SamRecord& record)
{
    if(record.getReadLength() == 0)
    {
        // There is no specific position information for this record,
        // so just return without processing.
        return;
    }

    // Call the base class:
    PileupElement::addEntry(record);

    if(myChromID == INVALID_CHROMOSOME_INDEX)
    {
        myChromID = record.getReferenceID();
    }

    // Get the position within the read of this entry.
    Cigar* cigar = record.getCigarInfo();
    if(cigar == NULL)
    {
        throw std::runtime_error("Failed to retrieve cigar info from the record.");
    }

    // The cycle is the query index.
    int32_t cycle = 
        cigar->getQueryIndex(getRefPosition(), record.get0BasedPosition());
    
    char base;
    char qual;
    if(cycle == Cigar::INDEX_NA)
    {
        // This position was not found in the read, so is a deletion.
        if(ourIgnoreDeletion)
        {
            // We are set to ignore deletions, so just return without processing
            return;
        }
        // set the base for the deletion.
        base = DELETION_BASE;
        // No quality, so set to 0.
        qual = 0;
    }
    else
    {
        // Get the base & quality.
        base = record.getSequence(cycle);
        qual = record.getQuality(cycle);
    }
    bool strand = !SamFlag::isFirstFragment(record.getFlag());
    uint8_t mq = record.getMapQuality();

    // Check to see if this entry matches the reference base.
    if(base != myRefBase)
    {
        myAllRef = false;
    }

    myBaseInfoRecord.add(base, qual, cycle, strand, mq);
}


// Perform the analysis associated with this class.  May be a simple print, 
// a calculation, or something else.  Typically performed when this element
// has been fully populated by all records that cover the reference position.
void PileupElementBaseInfo::analyze()
{
    // Check the previous position.
    int posDiff = getRefPosition() - ourPrevPos;
    if(ourPrevChromID != myChromID || posDiff > ourGapSize)
    {
        // Write a new position record.
        // TODO        BaseInfoRecord::writePos(getChromosome(), getRefPosition(), ourOutputFile);
        BaseInfoRecord::writePos(myChromID, getRefPosition(), ourOutputFile);
    }
    else if(posDiff < 0)
    {
        // Position is going back.  This is an error.
        std::cerr << "PileupElementBaseInfo is being analyzed out of order!\n";
        return;
    }
    else
    {
        // Write empty records to fill in the gap to this position.
        while(posDiff > 0)
        {
            BaseInfoRecord::writeEmpty(ourOutputFile);
            --posDiff;
        }
    }

    // If all positions match the reference, write the short record.
    if(myAllRef)
    {
        myBaseInfoRecord.writeRefOnly(ourOutputFile);
    }
    else
    {
        // Write the full record.
        myBaseInfoRecord.writeDetailed(ourOutputFile);
    }

    // Update the last position written.
    ourPrevChromID = myChromID;
    ourPrevPos = getRefPosition();
}


// Resets the entry, setting the new position associated with this element.
void PileupElementBaseInfo::reset(int32_t refPosition)
{
    // Call the base class.
    PileupElement::reset(refPosition);

    initVars();

    // get the reference base.
    myRefBase = getRefBase();

    myBaseInfoRecord.reset();

    myChromID = INVALID_CHROMOSOME_INDEX;
}

void PileupElementBaseInfo::initVars()
{
    myOutputString.Clear();
    // Default that all bases match the reference.
    // This will be changed when one does not.
    myAllRef = true;
}

