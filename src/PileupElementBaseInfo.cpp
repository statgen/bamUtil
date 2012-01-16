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

/////////////////////////////////////////////////////////////////////////////
//
// PileupElementBaseInfo
//

IFILE PileupElementBaseInfo::ourOutputFile = NULL;
int PileupElementBaseInfo::ourGapSize = 0;
int PileupElementBaseInfo::ourPrevPos = -1;

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
    // Call the base class:
    PileupElement::addEntry(record);

    // Get the position within the read of this entry.
    Cigar* cigar = record.getCigarInfo();
    if(cigar == NULL)
    {
        throw std::runtime_error("Failed to retrieve cigar info from the record.");
    }

    // The cycle is the query index.
    int32_t cycle = 
        cigar->getQueryIndex(getRefPosition(), record.get0BasedPosition());
    
    // Get the base.
    char base = record.getSequence(cycle);
    char qual = record.getQuality(cycle);
    bool strand = !SamFlag::isFirstFragment(record.getFlag());
    uint8_t mq = record.getMapQuality();

    // Check to see if this entry matches the reference base.
    myBaseInfoRecord.add(base, qual, cycle, strand, mq);
}


// Perform the alalysis associated with this class.  May be a simple print, 
// a calculation, or something else.  Typically performed when this element
// has been fully populated by all records that cover the reference position.
void PileupElementBaseInfo::analyze()
{

}


// Resets the entry, setting the new position associated with this element.
void PileupElementBaseInfo::reset(int32_t refPosition)
{
    // Call the base class.
    PileupElement::reset(refPosition);

    initVars();

    // get the reference base.
    myRefBase = getRefBase();
}

void PileupElementBaseInfo::initVars()
{
    myOutputString.Clear();
}

