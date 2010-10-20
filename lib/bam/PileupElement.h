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

#ifndef __PILEUP_ELEMENT_H__
#define __PILEUP_ELEMENT_H__

#include "SamRecord.h"

/// This is a base class pileup component, representing the information
/// for one reference position.  Child classes will be defined to detail more
/// information that needs to be saved and how it should be analyzed.
class PileupElement
{
public:
    PileupElement();

    PileupElement(const PileupElement& q);

    virtual ~PileupElement();


    // Add an entry to this pileup element.  
    virtual void addEntry(SamRecord& record);

    // Perform the analysis associated with this class.
    virtual void analyze();

    // Resets the entry, setting the new position associated with this element.
    virtual void reset(int32_t refPosition);
    
    const char* getChromosome() const { return(myChromosome.c_str()); }

    int32_t getRefPosition()  const { return(myRefPosition); }

protected:
    static const int32_t UNSET_POSITION = -1;

private:
    int32_t myRefPosition;
    std::string myChromosome;
};


#endif
