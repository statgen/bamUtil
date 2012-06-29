/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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

#ifndef __MATE_MAP_BY_COORD_H__
#define __MATE_MAP_BY_COORD_H__

#include "SamFile.h"


/// Class for buffering up reads waiting to wait for the mate's match.
/// Assumes the records are added in a coordinate sorted order.
/// Assumes the mate chromosome/position information on reads are accurate
/// otherwise the mates will not be found.
class MateMapByCoord
{
public:

    /// mateCoord set to false (default) means to store the records in the
    /// map based on the record's coordinate.  That way the first record in
    /// the map has the earliest position.
    /// mateCoord set to true means to store the records in the map based on
    /// the mate's coordinate.  That way the first record in the map has 
    /// the earliest mate position.
    MateMapByCoord(bool mateCoord = false);
    ~MateMapByCoord();

    /// Return the mate of the specified record, or NULL if mate is not found.
    /// The mate is removed from this map if it is found.
    SamRecord* getMate(SamRecord& record);

    /// Add the specified record to this mate map.
    void add(SamRecord& record);

    /// Return the first record (position-wise).
    /// If there are no records, NULL is returned.
    /// The record is NOT removed from the map, call popFirst to remove the 
    /// record.
    SamRecord* first();

    /// Remove the first record from the map.
    void popFirst();

protected:

private:
    typedef std::pair<uint64_t, SamRecord*> MATE_MAP_PAIR;
    typedef std::multimap<uint64_t, SamRecord*> MATE_MAP;
    
    MATE_MAP myMateBuffer;
    bool myMateCoord;
};


#endif
