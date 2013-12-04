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
#include "Recab.h"
#include "SamFlag.h"

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
    virtual const char* getProgramName() {return("bam:dedup");}

    Dedup():
        myRecab(),
        myDoRecab(false),
        myOneChrom(false),
        mySamPool(),
        lastCoordinate(-1), lastReference(-1), numLibraries(0), 
        myNumMissingMate(0),
        myForceFlag(false),
        myMinQual(15)
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


    struct DupKey
    {
        int32_t reference;
        int32_t coordinate;
        bool orientation;
        uint32_t libraryID;
        DupKey()
            : reference(0), coordinate(0), orientation(false), libraryID(0) {}
        DupKey(int32_t ref, int32_t coord, bool orient, uint32_t lib)
            : reference(ref), coordinate(coord), orientation(orient), libraryID(lib) {}
        inline void updateKey(SamRecord& record, uint32_t libID)
        {
            reference = record.getReferenceID();
            coordinate = record.get0BasedUnclippedStart();
            orientation = SamFlag::isReverse(record.getFlag());
            if(orientation)
            {
                // Reverse, so get the unclipped end.
                coordinate = record.get0BasedUnclippedEnd();
            }
            libraryID = libID;
        }
        inline void cleanupKey(int32_t referenceID, int32_t coord)
        {
            reference = referenceID;
            coordinate = coord - CLIP_OFFSET;
            orientation = false;
            libraryID = 0;
        }
        inline bool operator <(const DupKey& key) const
        {
            if(reference == key.reference)
            {
                // same reference, so check the coordinate.
                if(coordinate == key.coordinate)
                {
                    // Same coordinate, so check orientation.
                    if(orientation == key.orientation)
                    {
                        // Same orientation, so check library id.
                        return(libraryID < key.libraryID);
                    }
                    // Different orientations, so this is less than the
                    // other if this orientation is false.
                    return(!orientation);
                }
                // Same Ref, & different coordinates, so just check that.
                return(coordinate < key.coordinate);
            }
            // Different references, so less than is determined by 
            // checking the reference only.
            return(reference < key.reference);
        }
        inline DupKey&  operator = (const DupKey& key)
        {
            reference = key.reference;
            coordinate = key.coordinate;
            orientation = key.orientation;
            libraryID = key.libraryID;
            return(*this);
        }
    };
    
    // Each read is assigned a key based on its referenceID, coordinate, orientation, and libraryID
    // This structure stores the two keys in a paired end read.
    struct PairedKey {
        DupKey key1;
        DupKey key2;
        PairedKey(DupKey k1, DupKey k2)
        {
            if(k2 < k1)
            {
                key1 = k2;
                key2 = k1;
            }
            else
            {
                key1 = k1;
                key2 = k2;
            }
        }
    };


    // Paired key comparison operator used for sorting paired end reads.
    struct PairedKeyComparator {
        inline bool operator() (const PairedKey& lhs, const PairedKey& rhs) const {
            if (lhs.key2 < rhs.key2) return true;
            if (rhs.key2 < lhs.key2) return false;
            return lhs.key1 < rhs.key1;
        }
    };

    // A map from read group IDs to its libraryID
    typedef std::map< std::string, uint32_t, std::less<std::string> > StringToInt32Map;
    StringToInt32Map rgidLibMap;

    // A map from the key of a single read to its read data
    typedef std::map< DupKey, ReadData > FragmentMap;
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

    // Recalibrator logic.
    Recab myRecab;
    bool myDoRecab;
    bool myOneChrom;
    
    // Pool of sam records.
    SamRecordPool mySamPool;
    
    int lastCoordinate;
    int lastReference;
    uint32_t numLibraries;
    int myNumMissingMate;
    bool myForceFlag;
    int myMinQual;

    static const int DEFAULT_MIN_QUAL;
    static const uint32_t CLIP_OFFSET;

    // Once record is read, look back at previous reads and determine 
    // if any no longer need to be kept for duplicate checking.
    // Call with NULL to cleanup all records.
    void cleanupPriorReads(SamRecord* record);

    // Determines if the current position has changed when we read record
    bool hasPositionChanged(SamRecord & record);
    
    // When a record is read, check if it is a duplicate or
    // store for future checking.
    void checkDups(SamRecord & record, uint32_t recordCount);

    // Add the base qualities in a read
    int getBaseQuality(SamRecord& record);

    // builds the read group library map
    void buildReadGroupLibraryMap(SamFileHeader& header);

    // Returns the libraryID of a record
    uint32_t getLibraryID(SamRecord& record, bool checkTags = false);

    // Handle records that are not to be marked duplicates.
    // Performs any additional processing the first time through the file.
    void handleNonDuplicate(SamRecord* recordPtr);

    // Handle records whose mate was not found.  This will handle all processing
    // including calling handleNonDuplicate. 
    void handleMissingMate(SamRecord* recordPtr);

    void handleDuplicate(uint32_t index, SamRecord* recordPtr);

    inline int getFirstIndex(const DupKey& key1, 
                             int key1Index,
                             const DupKey& key2,
                             int key2Index)
    {
        if(key1.reference < key2.reference)
        {
            // key1 has a smaller chromosome
            return(key1Index);
        }
        if(key1.reference > key2.reference)
        {
            // key2 has a smaller chromosome
            return(key2Index);
        }
        // Same chromosome, so check the coordinate.
        if(key1.coordinate < key2.coordinate)
        {
            // key1 has a smaller coordinate
            return(key1Index);
        }
        if(key1.coordinate < key2.coordinate)
        {
            // key2 has a smaller coordinate
            return(key2Index);
        }

        // Same chromosome & Coordinate, so return the
        // smaller index.
        if(key1Index < key2Index)
        {
            return(key1Index);
        }
        return(key2Index);
    }
};

#endif // __DE_DUP_H
