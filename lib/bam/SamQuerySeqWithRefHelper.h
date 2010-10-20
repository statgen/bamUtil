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

#ifndef __SAM_QUERY_SEQ_WITH_REF_HELPER_H__
#define __SAM_QUERY_SEQ_WITH_REF_HELPER_H__

#include <stdint.h>

#include "SamRecord.h"
#include "GenomeSequence.h"

// This class contains the match/mismatch information
// between the reference and a read for a single base.
class SamSingleBaseMatchInfo
{
public:
    // More types can be added later as needed.
    enum Type {UNKNOWN, MATCH, MISMATCH};

    SamSingleBaseMatchInfo();
    ~SamSingleBaseMatchInfo();


    // Get info from this class.
    Type getType();
    int32_t getQueryIndex();

    // Set info in this class.
    void setType(Type newType);
    void setQueryIndex(int32_t queryIndex);

private:
    Type myType;
    int32_t myQueryIndex;
};

// Iterates through the query and compare with reference.
// NOTE: References to the GenomeSequence and SamRecord are stored, the objects
// are not copied, so they must remain valid as long as this class is used.
class SamQuerySeqWithRefIter
{
public:
    SamQuerySeqWithRefIter(SamRecord& record, GenomeSequence& refSequence,
                           bool forward = true);
    virtual ~SamQuerySeqWithRefIter();
    
    // Reset to start at the beginning of the record.
    // This will re-read values from SamRecord, so can be used if it has
    // changed to contain information for a new record.
    // forward = true means to start from the beginning and go to the end.
    // forward = false means to start from the end and go to the beginning.
    bool reset(bool forward = true);
    
    bool getNextMatchMismatch(SamSingleBaseMatchInfo& matchMismatchInfo);
    
private:

    SamQuerySeqWithRefIter();
    
    void nextIndex();

    // Tells whether or not the two bases are equal
    bool areEqual(char base1, char base2);

    SamRecord& myRecord;
    GenomeSequence& myRefSequence;
    Cigar* myCigar;
    uint32_t myStartOfReadOnRefIndex;
    int32_t myQueryIndex;
    bool myForward;
};

#endif
