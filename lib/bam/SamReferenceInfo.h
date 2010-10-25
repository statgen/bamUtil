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

#ifndef __SAM_REFERENCE_INFO_H__
#define __SAM_REFERENCE_INFO_H__

#include "StringArray.h"
#include "StringHash.h"
#include "IntArray.h"

class SamReferenceInfo
{
public:
    SamReferenceInfo();
    ~SamReferenceInfo();
    // Add reference sequence name and reference sequence length.
    void add(const char* referenceSequenceName, 
             int32_t referenceSequenceLength);

    int getReferenceID(const String & referenceName);
    int getReferenceID(const char* referenceName);
    const String & getReferenceLabel(int id) const;

    // Get the number of entries contained here.
    int32_t getNumEntries() const;

    // Return the reference name at the specified index.
    // Returns "" if index is out of bounds
    const char* getReferenceName(int index) const;
    
    // Return the reference length at the specified index.
    // Returns 0 if index is out of bounds
    int32_t getReferenceLength(int index) const;

    // Reset this reference info.
    void clear();

    SamReferenceInfo & operator = (const SamReferenceInfo & rhs);

private:
    // Reference Name information
    StringArray    myReferenceContigs;
    StringIntHash  myReferenceHash;
    IntArray       myReferenceLengths;
};

#endif

