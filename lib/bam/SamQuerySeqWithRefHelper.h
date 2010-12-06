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

/// This class contains the match/mismatch information
/// between the reference and a read for a single base.
class SamSingleBaseMatchInfo
{
public:
    /// More types can be added later as needed.
    enum Type {UNKNOWN, MATCH, MISMATCH};

    SamSingleBaseMatchInfo();
    ~SamSingleBaseMatchInfo();


    /// Get the type (match/mismatch/unknown) for this object.
    Type getType();

    /// Get the query index for this object.
    int32_t getQueryIndex();

    /// Set the type (match/mismatch/unkown) for this object.
    void setType(Type newType);

    /// Set the query index for this object.
    void setQueryIndex(int32_t queryIndex);

private:
    Type myType;
    int32_t myQueryIndex;
};

/// Iterates through the query and compare with reference.
/// NOTE: References to the GenomeSequence and SamRecord are stored, the objects
/// are not copied, so they must remain valid as long as this class is used.
class SamQuerySeqWithRefIter
{
public:
    SamQuerySeqWithRefIter(SamRecord& record, GenomeSequence& refSequence,
                           bool forward = true);
    virtual ~SamQuerySeqWithRefIter();
    
    /// Reset to start at the beginning of the record.
    /// This will re-read values from SamRecord, so can be used if it has
    /// changed to contain information for a new record.
    /// \param forward true means to start from the beginning and go to the end;
    /// false means to start from the end and go to the beginning.
    /// \return true if successfully reset; false if failed to read the Cigar.
    bool reset(bool forward = true);
    
    /// Returns information for the next position where the query and the 
    /// reference match or mismatch.  To be a match or mismatch, both the query
    /// and reference must have a base that is not 'N'.
    /// This means:
    ///    insertions and deletions are not mismatches or matches.
    ///    'N' bases are not matches or mismatches
    /// \param matchMismatchInfo return parameter with the information about
    /// the matching/mismatching base.
    /// \return true if there was another match/mismatch
    /// (matchMismatchInfo was set); false if not.
    bool getNextMatchMismatch(SamSingleBaseMatchInfo& matchMismatchInfo);
    
private:

    SamQuerySeqWithRefIter();
    
    void nextIndex();

    SamRecord& myRecord;
    GenomeSequence& myRefSequence;
    Cigar* myCigar;
    uint32_t myStartOfReadOnRefIndex;
    int32_t myQueryIndex;
    bool myForward;
};


/// Contains methods for converting between query sequence and reference.
/// Sequence.
/// NOTE: References to the GenomeSequence and SamRecord are stored, the objects
/// are not copied, so they must remain valid as long as this class is used.
class SamQuerySeqWithRef
{
public:
    static void seqWithEquals(const char* currentSeq,
                              int32_t seq0BasedPos,
                              Cigar& cigar, 
                              const char* referenceName,
                              GenomeSequence& refSequence,
                              std::string& updatedSeq);

    static void seqWithoutEquals(const char* currentSeq,
                                 int32_t seq0BasedPos,
                                 Cigar& cigar, 
                                 const char* referenceName,
                                 GenomeSequence& refSequence,
                                 std::string& updatedSeq);

private:
    SamQuerySeqWithRef();
};
#endif
