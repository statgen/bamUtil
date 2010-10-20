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

#ifndef _PILEUP_READER_H
#define _PILEUP_READER_H

#include <list>
#include <stdint.h>
#include <string>
#include <queue>

#include "SamFile.h"

///
/// \file this file implements a means of reading the
/// bases and structural variations of interest that are
/// present in any BAM file.
/// 
/// Two data structures contain the pertinent information:
///
/// - ReadBase contains information about a single base in
///   a single sequence, and is a copy of the useful fields in
///   the BAM record that it originated from.
///
/// - ReadInsertion contains the information that can't
///   be represented as a base itself, which is material
///   inserted at the genome at a given location.
///
/// The class ReadBase represents either bases that are directly
/// present in the class or it represents deletions of those
/// bases.  So for example, if you ask for chromosome 1, position
/// 25 and the corresponding aligned read indicates in the CIGAR
/// string that there is a deletion there, the ReadBase object
/// will not show a valid base, but rather indicate that this
/// base and some number of subsequent bases were deleted.
///
/// Insertions cannot be represented directly as a single added
/// base nor as a set of bases that were deleted (which ReadBase
/// does) but rather a set of bases that appears in between two
/// positions on the reference.
///
/// The PileupReader class provides the mechanism for attaching
/// to a BAM file and retrieving the above data.
///
/// Use cases:
///
/// - read the BAM data associated with specific known marker
///   positions, usually those read locations that are provided
///   by a sequencing chip.  In this case, PileupReader will by
///   default use an index if present to rapidly skip from marker
///   to marker.  This behavior may be controlled by using the
///   useIndex() method.
///
/// - scan the entire genome, base by base, reading bases and
///   insertions and evaluating them, for example to do novel
///   SNP detection.
/// 

///
/// When we call PileupReader::getPileup, we return a set
/// of these objects representing the bases found in all
/// overlapping reads in the given BAM file.
///
/// The information in this class is the pertinent information
/// from the SAM format record.
///
class ReadBase {
    SamRecord   *_samRecord;
    int         _deletedBaseCount;   // how many bases, including this one, are deleted
    int         _insertedBasesOverlapping;
    bool        _InsertionFollowing;        // THIS read contains an insertion immediately after this base

    // More private state is needed - e.g. offset and length into read.
    // But also cigar string parsing state information - probably need
    // a small common state class used here and in ReadInsertion
public:
    ReadBase() : _deletedBaseCount(0), _insertedBasesOverlapping(0), _InsertionFollowing(false) {;}

    bool set(SamRecord &, int basePosition);

    char getBase();
    char getPhredQualityChar();
    int  getPhredQualityScore();

    int  getDeletedBaseCount();
    int  getInsertedBasesOverlapping();

    bool hasInsertionFollowing();
};

///
/// Contain all pertinent information from overlapping
/// SAM records that START an insertion at this point.
///
class ReadInsertion {
    SamRecord *_samRecord;
    // More private state is needed - e.g. offset and length into read.
    // But also cigar string parsing state information - probably need
    // a small common state class used here and in ReadBase
public:
    bool set(SamRecord &r);

    /// STL single copy access:
    bool        getQualities(std::string &phredQualities);
    bool        getBases(std::string &bases);

    /// zero copy access:
    int         getQualities(char *phredQualities, int maxLength);
    int         getBases(char *bases, int maxLength);

    uint8_t     getMapQuality();
    uint16_t    getSAMFlag();
    SamRecord & getSamRecord();
};

///
/// class for tracking BAM file sequence data that covers markers of interest.
///
/// Given a sorted BAM file, allow us to scan over the file, either base by
/// base or marker by marker in ascending order, and obtain a list of
/// bases from the BAM file at the corresponding position.
///
/// In the event that an insertion is noticed after the marker of interest,
/// signal that in getPileup, and allow the data to be returned by calling
/// getInsertions().
///
///
class SequenceCoverageReader {

private:
    SamFile&    _samFile;
    int         _referenceID;
    uint32_t    _position;
    bool        _useIndex;
    std::list<SamRecord *>    _samRecords;    // list of ptrs to records that overlap the last read pileup position
    std::queue<SamRecord *> _unused;         // fifo to recliam pre-allocated records from

public:

    SequenceCoverageReader(SamFile &samFile) :
        _samFile(samFile),
        _referenceID(0),
        _position(0),
        _useIndex(true)
        { ; }

    void useIndex(bool u=true) { _useIndex = u; }

    /// \brief getBases returns a pileup vector of the same
    /// length as the number of reads that overlaps the given
    /// chromosome and position
    /// 
    /// \param referenceID  identifies the chromosome in which
    /// to retrieve bases from.
    /// 
    /// \param position  indicates the 1-based offset of the 
    /// base of interest on that chromosome.  This number is
    /// relative to the genome reference that the BAM file
    /// was aligned against.
    /// 
    /// \param readBases a caller allocated vector of of readBase
    /// objects
    /// 
    /// \param basesInserted bool that indicates if getInsertions() will
    /// return a non-zero vector of bases that were inserted AFTER
    /// this position.  Since reads may truncate this insertion,
    /// they may be different lengths.
    ///
    /// \return true if more data exists in the BAM file, false otherwise
    ///
    bool getBases(int32_t referenceID, uint32_t position, std::vector<ReadBase> &readBases, bool &basesInserted);

    bool getInsertedBases(int referenceID, uint32_t position, std::vector<ReadInsertion> &readInsertions);

};


#endif
