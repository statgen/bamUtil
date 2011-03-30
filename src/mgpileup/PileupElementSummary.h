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

#ifndef __PILEUP_ELEMENT_SUMMARY_H__
#define __PILEUP_ELEMENT_SUMMARY_H__

#include <stdint.h>
#include "PileupElement.h"

/// This class inherits from the base class and stores base count summaries.
class PileupElementSummary : public PileupElement
{
public:
    PileupElementSummary();
    // NOTE that this method does not actually copy, it just resets.
    PileupElementSummary(const PileupElementSummary& q);
    virtual ~PileupElementSummary();

    // Add an entry to this pileup element.
    virtual void addEntry(SamRecord& record);

    // Perform the alalysis associated with this class.  In this case, it is
    // a print of the base & quality information associated with this position.
    virtual void analyze();

    // Resets the entry, setting the new position associated with this element.
    virtual void reset(int32_t refPosition);

private:
    // Reset the allele counts.
    void reset();

    int myNumAlleles[6];
    int myDepth;
    int myNumDeletes;
    char myRefAllele;
};

#endif
