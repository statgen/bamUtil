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

#include "MateMapByCoord.h"
#include "SamHelper.h"

MateMapByCoord::MateMapByCoord(bool mateCoord)
    : myMateBuffer(),
      myMateCoord(mateCoord)
{
}


MateMapByCoord::~MateMapByCoord()
{
    myMateBuffer.clear();
}


SamRecord* MateMapByCoord::getMate(SamRecord& record)
{
    SamRecord* mate = NULL;
    const char* readName = record.getReadName();

    // Get the key for finding this mate (its chrom/pos).
    uint64_t mateKey;
    if(myMateCoord)
    {
        // Since the record's are stored by the mate coordinate, 
        // search the map using the passed in record's coordinate.
        mateKey = SamHelper::combineChromPos(record.getReferenceID(), 
                                             record.get0BasedPosition());
    }
    else
    {
        // Since the record's are stored by the record's coordinate, 
        // search the map using the passed in record's mate coordinate.
        mateKey = SamHelper::combineChromPos(record.getMateReferenceID(), 
                                         record.get0BasedMatePosition());
    }

    std::pair<MATE_MAP::iterator,MATE_MAP::iterator> matches;

    // Find for the elements with this key.
    matches = myMateBuffer.equal_range(mateKey);

    // Loop through the elements that matched the key looking for the mate.
    for(MATE_MAP::iterator iter = matches.first; 
        iter != matches.second; iter++)
    {
        if(strcmp(((*iter).second)->getReadName(), readName) == 0)
        {
            // Found the match.
            mate = (*iter).second;
            // Remove the entry from the map.
            myMateBuffer.erase(iter);
            break;
        }
    }
    return(mate);
}


void MateMapByCoord::add(SamRecord& record)
{
    uint64_t chromPos;

    if(myMateCoord)
    {
        chromPos = SamHelper::combineChromPos(record.getMateReferenceID(),
                                              record.get0BasedMatePosition());
    }
    else
    {
        chromPos = SamHelper::combineChromPos(record.getReferenceID(),
                                              record.get0BasedPosition());
    }

    myMateBuffer.insert(MATE_MAP_PAIR(chromPos, &record));
}


SamRecord* MateMapByCoord::first()
{
    MATE_MAP::iterator first = myMateBuffer.begin();
    if(first == myMateBuffer.end())
    {
        return(NULL);
    }

    // Return the record from the first element.
    return(first->second);
}


void MateMapByCoord::popFirst()
{
    MATE_MAP::iterator first = myMateBuffer.begin();
    if(first != myMateBuffer.end())
    {
        // There is a first element, so remove it.
        myMateBuffer.erase(first);
    }
    return;
}
