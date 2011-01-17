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

#ifndef __PILEUP_ELEMENT_BASE_QUAL_H__
#define __PILEUP_ELEMENT_BASE_QUAL_H__

#include <stdint.h>
#include "PileupElement.h"
#include "GenomeSequence.h"

/// This class inherits from the base class and stores base and qualities.
class PileupElementBaseQual : public PileupElement
{
public:
    PileupElementBaseQual();
    PileupElementBaseQual(bool addDelAsBase);
    PileupElementBaseQual(const PileupElementBaseQual& q);
    virtual ~PileupElementBaseQual();
 
    // Add an entry to this pileup element.  
    virtual void addEntry(SamRecord& record);

    // Perform the alalysis associated with this class.  In this case, it is
    // a print of the base & quality information associated with this position.
    virtual void analyze();

	//computes GL Scores
	void computeGLScores(int index, int16_t* GLScores, char* bases, int8_t* baseQualities);

    // Resets the entry, setting the new position associated with this element.
    virtual void reset(int refPosition, GenomeSequence* refSeq, InputFile* vcfOutFile, bool addDelAsBase, double*** logGLMatrix);
    virtual void reset(int refPosition);
    		
    // Allows for repeat polymorphisms
    virtual const char* getRefAllele();
    
private:
    static const char UNSET_QUAL = 0xFF;

    char* myBases;
    int8_t* myMapQualities;
    int8_t* myQualities;
    char* myStrands;
    int8_t* myCycles;
    int16_t* myGLScores;
    int myAllocatedSize;
    int myIndex;
    bool myAddDelAsBase;
    std::string myRefAllele;
    GenomeSequence* myRefSeq;
	InputFile* myVcfOutFile;
	double*** myLogGLMatrix;
};

#endif
