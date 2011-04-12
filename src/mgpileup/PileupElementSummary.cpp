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

#include "PileupElementSummary.h"
#include "BaseAsciiMap.h"

PileupElementSummary::PileupElementSummary()
    : PileupElement(),
      myRefAllele(' ')
{
    reset();
}


// NOTE that this method does not actually copy, it just resets.
PileupElementSummary::PileupElementSummary(const PileupElementSummary& q)
    : PileupElement(),
      myRefAllele(' ')
{
    reset();
}


PileupElementSummary::~PileupElementSummary()
{
}


// Add an entry to this pileup element.  
void PileupElementSummary::addEntry(SamRecord& record)
{
    // Call the base class:
    PileupElement::addEntry(record);

    // Must wait til here to get positions since PileupElement::addEntry is what sets the chromosome.
    GenomeSequence* refPtr = getReference();
    if((myRefAllele == ' ') && (refPtr != NULL))
    {
        genomeIndex_t markerIndex = refPtr->getGenomePosition(getChromosome(), static_cast<uint32_t>(getRefPosition()+1));
        myRefAllele = (*refPtr)[markerIndex];
    }

    if(strcmp(record.getSequence(), "*") == 0)
    {
        // no sequence, so just return.
        return;
    }

    Cigar* cigar = record.getCigarInfo();
    if(cigar == NULL)
    {
        throw std::runtime_error("Failed to retrieve cigar info from the record.");
    }
    
    ++myDepth;

    int32_t readIndex = 
        cigar->getQueryIndex(getRefPosition(), record.get0BasedPosition());

    // If the readPosition is N/A, this is a deletion.
    // Nothing is done to identify an insertion, just pull the "matching" base.
    if(readIndex != CigarRoller::INDEX_NA)
    {
        char base = record.getSequence(readIndex);

        myNumAlleles[BaseAsciiMap::baseColor2int[base]]++;
    }
    else
    {
        ++myNumDeletes;
    }
}

void PileupElementSummary::analyze()
{
    if(getRefPosition() != UNSET_POSITION)
    {
        std::cout << getChromosome() << "\t" << getRefPosition() + 1
                  << "\tRefBase:" << myRefAllele
                  << "\tDepth:" << myDepth;
        for(int i = 0; i <= BaseAsciiMap::baseXIndex; i++)
        {
            std::cout << "\t" << BaseAsciiMap::int2base[i] 
                      << ":" << myNumAlleles[i];
        }
        std::cout << "\tDel:" << myNumDeletes << "\n";
    }
}

void PileupElementSummary::reset(int refPosition)
{
    // Call the base class.
    PileupElement::reset(refPosition);
    reset();
 }

void PileupElementSummary::reset()
{
    for(int i = 0; i <= BaseAsciiMap::baseXIndex; i++)
    {
        myNumAlleles[i] = 0;
    }
    myDepth = 0;
    myNumDeletes = 0;
    myRefAllele = ' ';
}
