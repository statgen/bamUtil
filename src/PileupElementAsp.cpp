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

AspFileWriter* PileupElementAsp::ourOutputFile = NULL;
bool PileupElementAsp::ourIgnoreDeletion = false;
bool PileupElementAsp::ourReportOverMax = true;

void PileupElementAsp::setOutputFile(AspFileWriter* outputFile)
{
    ourOutputFile = outputFile;
}


void PileupElementAsp::setIgnoreDeletion(bool ignoreDeletion)
{
    ourIgnoreDeletion = ignoreDeletion;
}


PileupElementAsp::PileupElementAsp()
    : PileupElement(),
      myAspRecord()
{
    initVars();
}


// NOTE that this method does not actually copy, it just resets.
PileupElementAsp::PileupElementAsp(const PileupElementAsp& q)
    : PileupElement(),
      myAspRecord()
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

    // Call the base class since there is data to process.
    PileupElement::addEntry(record);

    if(myRefBase == UNKNOWN_REF_BASE)
    {
        // This means the reference base wasn't set, so set it now.
        myRefBase = getRefBase();
    }

    // If the reference is 'N', skip this position.
    if(myRefBase == 'N')
    {
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
    bool strand = SamFlag::isReverse(record.getFlag());
    uint8_t mq = record.getMapQuality();

    // Check to see if this entry matches the reference base.
    if(base == '=')
    {
        base = myRefBase;
    }
    // Check if it does not match the reference base and is not ambiguous.
    else if(!BaseUtilities::isAmbiguous(base) && (base != myRefBase))
    {
        myAllRef = false;
    }
    myAspRecord.add(myRefBase, base, qual, cycle, strand, mq);
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
    if((myAspRecord.getNumOverMaxBases() != 0) && (ourReportOverMax))
    {
        // This position had more than the max allowed, so not
        // all got added.
        std::cerr << "Went over the maximum number of records at a position, "
                  << getChromosome() << ":"
                  << getRefPosition() << "by "
                  << myAspRecord.getNumOverMaxBases() << " records\n";
    }

    // If all positions match the reference, write the short record.
    if(myAllRef)
    {
        // Set to ref only (default is full record).
        myAspRecord.setRefOnlyType();
    }
    else
    {
        // Write the full record.
        myAspRecord.setDetailedType();
    }
    // Write the record.
    if(ourOutputFile != NULL)
    {
        ourOutputFile->write(myAspRecord, getChromosome(), getRefPosition());
    }
}


// Resets the entry, setting the new position associated with this element.
void PileupElementAsp::reset(int32_t refPosition)
{
    // Call the base class.
    PileupElement::reset(refPosition);

    initVars();

    myAspRecord.reset();
}

void PileupElementAsp::initVars()
{
    myOutputString.Clear();
    // Default that all bases match the reference.
    // This will be changed when one does not.
    myAllRef = true;
    myRefBase = UNKNOWN_REF_BASE;
}

