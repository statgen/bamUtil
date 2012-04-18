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
#include "SamRecordPool.h"

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

    Dedup():
        mySamPool(),
        lastCoordinate(-1), lastReference(-1), numLibraries(0), 
        singleDuplicates(0),
        pairedDuplicates(0),
        removeFlag(false),
        forceFlag(false),
        verboseFlag(false),
        singleRead(0),
        firstPair(0),
        foundPair(0)
    {}

    ~Dedup();

private:
    struct ReadData
    {
        int sumBaseQual;
        SamRecord* recordPtr;
        int recordIndex;
        bool paired;
        ReadData()
            : sumBaseQual(0), recordPtr(NULL), recordIndex(0), paired(false) {}
    };
    struct PairedData
    {
        int sumBaseQual;
        SamRecord* record1Ptr;
        SamRecord* record2Ptr;
        int record1Index;
        int record2Index;
        PairedData()
            : sumBaseQual(0), record1Ptr(NULL), record2Ptr(NULL),
              record1Index(0), record2Index(0) {}
    };

    // Each read is assigned a key based on its referenceID, coordinate, orientation, and libraryID
    // This structure stores the two keys in a paired end read.
    struct PairedKey {
        uint64_t key1;
        uint64_t key2;
        PairedKey(uint64_t k1, uint64_t k2)
        {
            if(k1 <= k2)
            {
                key1 = k1;
                key2 = k2;
            }
            else
            {
                key1 = k2;
                key2 = k1;
            }
        }
    };

    // Paired key comparison operator used for sorting paired end reads.
    struct PairedKeyComparator {
        inline bool operator() (const PairedKey& lhs, const PairedKey& rhs) const {
            if (lhs.key2 < rhs.key2) return true;
            if (lhs.key2 > rhs.key2) return false;
            return lhs.key1 < rhs.key1;
        }
    };

    // A map from read group IDs to its libraryID
    typedef std::map< std::string, uint32_t, std::less<std::string> > StringToInt32Map;
    StringToInt32Map rgidLibMap;

    // A map from the key of a single read to its read data
    typedef std::map< uint64_t, ReadData, std::less<uint64_t> > FragmentMap;
    typedef std::pair<FragmentMap::iterator,bool> FragmentMapInsertReturn;
    FragmentMap myFragmentMap;

    // Map for storing reads until the mate is found.
    typedef std::multimap<uint64_t, ReadData> MateMap;
    typedef std::pair<MateMap::iterator,bool> MateMapInsertReturn;
    MateMap myMateMap;

    // A map from the key of a paired read to its read data
    typedef std::map<PairedKey, PairedData, PairedKeyComparator> PairedMap;
    typedef std::pair<PairedMap::iterator,bool> PairedMapInsertReturn;
    PairedMap myPairedMap;

    // Stores the record counts of duplicates reads
    typedef std::vector< uint32_t > Int32Vector;
    typedef std::vector< uint32_t >::iterator Int32VectorIterator;
    Int32Vector myDupList;

    // Pool of sam records.
    SamRecordPool mySamPool;

    int lastCoordinate;
    int lastReference;
    uint32_t numLibraries;
    uint32_t singleDuplicates, pairedDuplicates;
    bool removeFlag, forceFlag, verboseFlag;

    static const uint32_t CLIP_OFFSET;
    static const uint64_t UNMAPPED_SINGLE_KEY;
    static const uint64_t UNMAPPED_PAIRED_KEY;
    static const uint32_t EMPTY_RECORD_COUNT;
    static const int PAIRED_QUALITY_OFFSET;
    static const int MAX_REF_ID;
    static const int LOOK_BACK;

    // phased out
    int singleRead;
    int firstPair;
    int foundPair;

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

    // Determines if the current position has changed when we read record
    bool hasPositionChanged(SamRecord & record);

    // Once record is read, look back at previous reads and determine duplicates
    void markDuplicatesBefore(SamRecord & record);

    // Same as above, but it uses record's referenceID and coordinate
    void markDuplicatesBefore(uint32_t referenceID, uint32_t coordinate);

};

#endif // __DE_DUP_H
