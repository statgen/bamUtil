/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan
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

#ifndef __DE_DUP_H
#define __DE_DUP_H

#include "BamExecutable.h"
#include <vector>
#include <set>
#include <map>

/*---------------------------------------------------------------/
  /
  / This class will remove or mark duplicate reads from a BAM file.
  /
  /---------------------------------------------------------------*/
class Dedup : public BamExecutable
{
public:
    static void dedupDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);

    // Each read is assigned a key based on its referenceID, coordinate, orientation, and libraryID
    // This structure stores the two keys in a paired end read.
    struct PairedKey {
        uint64_t key1;
        uint64_t key2;
        PairedKey(uint64_t k1, uint64_t k2) : key1(k1), key2(k2) {}
    };

    // Paired key comparison operator used for sorting paired end reads.
    struct PairedKeyComparator {
        inline bool operator() (const PairedKey& lhs, const PairedKey& rhs) const {
            if (lhs.key2 < rhs.key2) return true;
            if (lhs.key2 > rhs.key2) return false;
            return lhs.key1 < rhs.key1;
        }
    };

    // When we have an overlapping paired end read, we store necessary data for updating the
    // records in this structure.
    struct UpdateData {
        Cigar cigar;
        std::string sequence;
        std::string baseQualities;
    };

    // This structure stores information needed to determine when two ends in a pair overlap
    struct OverlapData {
        Cigar cigar;
        int32_t clippedEnd;
        int32_t start;
        std::string sequence;
        std::string baseQualities;
    };

    // This structure stores some basic information from either a single read or paired read
    struct ReadData {
        int baseQuality;
        uint64_t key1, key2;
        uint32_t recordCount1, recordCount2;
        bool paired;
        OverlapData * overlapData;
        std::string readName;
        inline int getPairedBaseQuality() {
            return baseQuality + ( paired ? PAIRED_QUALITY_OFFSET : 0 );
        }
    };

    // A map from read group IDs to its libraryID
    typedef std::map< std::string, uint32_t, std::less<std::string> > StringToInt32Map;
    StringToInt32Map rgidLibMap;

    // A map from the read name to its read data
    typedef std::map< std::string, ReadData*, std::less<std::string> > StringToReadDataPointerMap;
    typedef std::map< std::string, ReadData*, std::less<std::string> >::iterator StringToReadDataPointerMapIterator;
    StringToReadDataPointerMap readDataMap; // key -> ReadData*

    // A vector of ReadData's
    typedef std::vector<ReadData*> ReadDataPointerVector;
    typedef std::vector<ReadData*>::iterator ReadDataPointerVectorIterator;

    // A map from the key of a single read to its read data
    typedef std::map< uint64_t, ReadDataPointerVector, std::less<uint64_t> > SingleKeyToReadDataPointerMap;
    typedef std::map< uint64_t, ReadDataPointerVector, std::less<uint64_t> >::iterator SingleKeyToReadDataPointerMapIterator;
    SingleKeyToReadDataPointerMap fragmentMap;

    // A map from the key of a paired read to its read data
    typedef std::map< PairedKey, ReadDataPointerVector, PairedKeyComparator > PairedKeyToReadDataPointerMap;
    typedef std::map< PairedKey, ReadDataPointerVector, PairedKeyComparator >::iterator PairedKeyToReadDataPointerMapIterator;
    PairedKeyToReadDataPointerMap pairedMap;

    // Stores the record counts of duplicates reads
    typedef std::vector< uint32_t > Int32Vector;
    typedef std::vector< uint32_t >::iterator Int32VectorIterator;
    Int32Vector duplicateIndices;

    // A map from record number to updated data for that record
    typedef std::map< uint32_t, UpdateData*, std::less<uint32_t> > UInt32ToUpdateDataPointerMap;
    typedef std::map< uint32_t, UpdateData*, std::less<uint32_t> >::iterator UInt32ToUpdateDataPointerIterator;
    UInt32ToUpdateDataPointerMap recordUpdates;

    int lastCoordinate;
    int lastReference;
    uint32_t numLibraries;
    uint32_t singleDuplicates, pairedDuplicates, overlappingPairs;
    bool removeFlag, forceFlag, verboseFlag, overlapFlag, overlapReplacementFlag;

    static const uint32_t CLIP_OFFSET;
    static const uint64_t UNMAPPED_SINGLE_KEY;
    static const uint64_t UNMAPPED_PAIRED_KEY;
    static const uint32_t EMPTY_RECORD_COUNT;
    static const int PAIRED_QUALITY_OFFSET;
    static const int MAX_REF_ID;
    static const int LOOK_BACK;

    // phased out
    int unmapped;
    int singleRead;
    int firstPair;
    int foundPair;
    int properPair;

    // builds the read group library map
    void buildReadGroupLibraryMap(SamFileHeader& header);

    // When a record is read, we put it into the appropriate maps
    void placeRecordInMaps(SamRecord & record, uint32_t recordCount);

    // Returns the libraryID of a record
    uint32_t getLibraryID(SamRecord& record, bool checkTags = false);

    // Add the base qualities in a read
    int getBaseQuality(SamRecord& record);

    // Returns the key constructed from those four pieces of information
    uint64_t makeKey(uint32_t reference, uint32_t coordinate, bool orientation, uint32_t libraryID);

    // Returns the key constructed from a record
    uint64_t makeKeyFromRecord(SamRecord & record);

    // Forms a cigar from an expanded string
    Cigar * rollupCigar(std::string cigarString);

    // Inserts soft clips into a cigar
    Cigar * insertClipsIntoCigar(Cigar * cigar, int32_t begin, int32_t end);
    
    // Determines if the current position has changed when we read record
    bool hasPositionChanged(SamRecord & record);

    // Once record is read, look back at previous reads and determine duplicates
    void markDuplicatesBefore(SamRecord & record);

    // Same as above, but it uses record's referenceID and coordinate
    void markDuplicatesBefore(uint32_t referenceID, uint32_t coordinate);

    // Handle the overlap between the reads in readData and record
    void processOverlap(ReadData * readData, SamRecord & record, uint32_t recordCount);

public:
    Dedup(): lastCoordinate(-1), lastReference(-1), numLibraries(0), 
             singleDuplicates(0),
             pairedDuplicates(0),
             overlappingPairs(0),
             removeFlag(false),
             forceFlag(false),
             verboseFlag(false),
             overlapFlag(false),
             overlapReplacementFlag(false),
             unmapped(0),
             singleRead(0),
             firstPair(0),
             foundPair(0),
             properPair(0) {}
};

#endif // __DE_DUP_H
