/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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

//////////////////////////////////////////////////////////////////////////
#include "AspReadNameID.h"

AspReadNameID::AspReadNameID()
    :  myNextID(0),
       myReadNameVector(),
       myReadNameMap()
{
    myReadNameVector.resize(MAX_ID + 1);
}

AspReadNameID::~AspReadNameID()
{
    myReadNameVector.clear();
    myReadNameMap.clear();
}


uint16_t AspReadNameID::getReadNameID(const char* readName)
{
    // Check to see if the readName is already in the map.
    AspReadNameMapPair mapPair(readName, myNextID);
    std::pair<AspReadNameMapType::iterator, bool> insertReturn = 
        myReadNameMap.insert(mapPair);
    // If the read name is new, the id is the next id.
    uint16_t returnID = myNextID;
    if(insertReturn.second)
    {
        // Inserted a new read name, so update the nextID info.
        if(!myReadNameVector[myNextID].empty())
        {
            // Already an element with this id, so remove it from
            // the read name map.
            myReadNameMap.erase(myReadNameVector[myNextID]);
        }
        // Update this ID to be associated with this read name.
        myReadNameVector[myNextID] = readName;
        // Update the next id.
        if(myNextID < MAX_ID)
        {
            ++myNextID;
        }
        else
        {
            // Hit the maximum id, so wrap around to 0.
            myNextID = 0;
        }
    }
    else
    {
        // This read name has already been assigned an id, so return that.
        returnID = (insertReturn.first)->second;
    }
    return(returnID);
}
