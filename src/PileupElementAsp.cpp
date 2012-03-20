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
#include "PileupElementAsp.h"
#include <stdexcept>

#include "SamFlag.h"
#include "BaseUtilities.h"

/////////////////////////////////////////////////////////////////////////////
//
// PileupElementAsp
//

IFILE PileupElementAsp::ourOutputFile = NULL;
int PileupElementAsp::ourGapSize = 0;
bool PileupElementAsp::ourIgnoreDeletion = false;
bool PileupElementAsp::ourReportOverMax = true;

int PileupElementAsp::ourPrevPos = PileupElement::UNSET_POSITION;
int PileupElementAsp::ourPrevChromID = INVALID_CHROMOSOME_INDEX;

int PileupElementAsp::ourNumPosRecs = 0;
int PileupElementAsp::ourNumEmptyRecs = 0;
int PileupElementAsp::ourNumRefOnlyRecs = 0;
int PileupElementAsp::ourNumDetailedRecs = 0;

void PileupElementAsp::setOutputFile(const char* outputFile)
{
    ourOutputFile = ifopen(outputFile, "w", InputFile::BGZF);
}

void PileupElementAsp::closeOutputFile()
{
    ifclose(ourOutputFile);
}


void PileupElementAsp::setGapSize(int gapSize)
{
    ourGapSize = gapSize;
}


void PileupElementAsp::setIgnoreDeletion(bool ignoreDeletion)
{
    ourIgnoreDeletion = ignoreDeletion;
}


void PileupElementAsp::printHeader()
{
}


PileupElementAsp::PileupElementAsp()
    : PileupElement()
{
    initVars();
}


PileupElementAsp::~PileupElementAsp()
{
}


// Add an entry to this pileup element.  
void PileupElementAsp::addEntry(SamRecord& record)
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
        // This means the reference base wasn't set either, so set it now.
        myRefBase = getRefBase();
    }

    // If the reference is 'N', skip this position.
    if(myRefBase == 'N')
    {
        return;
    }

    // Get the position within the read of this entry.
    Cigar* cigar = record.getCigarInfo();
    if(cigar == NULL)
    {
        throw std::runtime_error("Failed to retrieve cigar info from the record.");
    }
    if(cigar->size() == 0)
    {
        // There is no cigar, so return.
        return;
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
        base = AspRecord::DELETION_BASE;
        // No quality, so set to unknown.
        qual = BaseUtilities::UNKNOWN_QUALITY_CHAR;
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
    if(base == '=')
    {
        base = myRefBase;
    }
    else if(base != myRefBase)
    {
        myAllRef = false;
    }
    if(!myAspRecord.add(base, qual, cycle, strand, mq))
    {
        // Count the number of entries over the max.
        ++myOverMax;
    }
}


// Perform the analysis associated with this class.  May be a simple print, 
// a calculation, or something else.  Typically performed when this element
// has been fully populated by all records that cover the reference position.
void PileupElementAsp::analyze()
{
    // If the reference is 'N', skip this position.
    if(myRefBase == 'N')
    {
        return;
    }
    if((myOverMax != 0) && (ourReportOverMax))
    {
        // This position had more than the max allowed, so not
        // all got added.
        std::cerr << "Went over the maximum number of records at a position, "
                  << getChromosome() << ":"
                  << getRefPosition() << "by " << myOverMax << " records\n";
    }

    // Check the previous position.
    int posDiff = getRefPosition() - ourPrevPos;
    if(ourPrevChromID != myChromID || posDiff > ourGapSize)
    {
        // Write a new position record.
        AspRecord::writePos(myChromID, getRefPosition(), ourOutputFile);
        ++ourNumPosRecs;
    }
    else if(posDiff < 0)
    {
        // Position is going back.  This is an error.
        std::cerr << "PileupElementAsp is being analyzed out of order!\n";
        return;
    }
    else
    {
        // Write empty records to fill in the gap to this position.
        while(posDiff > 1)
        {
            AspRecord::writeEmpty(ourOutputFile);
            --posDiff;
            ++ourNumEmptyRecs;
        }
    }

    // If all positions match the reference, write the short record.
    if(myAllRef)
    {
        myAspRecord.setRefOnlyType();
        myAspRecord.write(ourOutputFile);
        ++ourNumRefOnlyRecs;
    }
    else
    {
        // Write the full record.
        myAspRecord.write(ourOutputFile);
        ++ourNumDetailedRecs;
    }

    // Update the last position written.
    ourPrevChromID = myChromID;
    ourPrevPos = getRefPosition();
}


// Resets the entry, setting the new position associated with this element.
void PileupElementAsp::reset(int32_t refPosition)
{
    // Call the base class.
    PileupElement::reset(refPosition);

    initVars();

    myAspRecord.reset();

    myChromID = INVALID_CHROMOSOME_INDEX;
}

void PileupElementAsp::initVars()
{
    myOutputString.Clear();
    // Default that all bases match the reference.
    // This will be changed when one does not.
    myAllRef = true;

    // Did not go over the max allowed records.
    myOverMax = 0;
}

