/*
 *  Copyright (C) 2016  Regents of the University of Michigan
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

#ifndef __DE_DUP_BASE_H
#define __DE_DUP_BASE_H

#include "BamExecutable.h"
#include <vector>
#include <map>
#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <unordered_map>
#endif
#include "SamRecordPool.h"
#include "Recab.h"

/*---------------------------------------------------------------/
  /
  / This class is the base class for duplicate marking.
  /
  /---------------------------------------------------------------*/
class DedupBase : public BamExecutable
{
protected:
    DedupBase():
        myForceFlag(false),
        myMarkNonPrimary(true),
        mySamPool(),
        mySecondarySupplementaryMap()
    {}
 
    ~DedupBase();


    // Second pass through the file, marking duplicates
    // and calling recalibration if necessary
    void markDuplicateLoop(bool verboseFlag, bool removeFlag, const String& inFile, const String& outFile, Recab* recabPtr);

    void setDuplicate(SamRecord* recordPtr, bool duplicate);

    void markDuplicateInNonPrimaryMaps(SamRecord* recordPtr);

    bool myForceFlag;
    bool myMarkNonPrimary;

    // Pool of sam records.
    SamRecordPool mySamPool;
    
    // Stores the record counts of duplicates reads
    typedef std::vector< uint32_t > Int32Vector;
    typedef std::vector< uint32_t >::iterator Int32VectorIterator;
    Int32Vector myDupList;
    
    // Secondary and Supplementary Read Map
    // True means it should be marked duplicate, 
    // false means it should not be marked duplicate.
    typedef std::pair<std::string,bool> NonPrimaryPair;
    #ifdef __GXX_EXPERIMENTAL_CXX0X__
    typedef std::unordered_map<std::string,bool> NonPrimaryMap;
    #else
    typedef std::map<std::string,bool> NonPrimaryMap;
    #endif
    NonPrimaryMap mySecondarySupplementaryMap;

};

#endif // __DE_DUP_H
