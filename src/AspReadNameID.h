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


#ifndef __ASP_READNAMEID_H__
#define __ASP_READNAMEID_H__

#include "SamFile.h"
#include "BaseAsciiMap.h"
#include <unordered_map>

class AspReadNameID
{
public:
    AspReadNameID();
    ~AspReadNameID();

    /// Get the ID associated with the specified read name.
    uint16_t getReadNameID(const char* readName);

private:
    AspReadNameID(const AspReadNameID & readNameID);

    static const unsigned int MAX_ID = 0xFFFF;

    uint16_t myNextID;

    std::vector<std::string> myReadNameVector;
    typedef std::unordered_map<std::string, uint16_t> AspReadNameMapType;
    //    typedef std::map<std::string, uint16_t> AspReadNameMapType;
    typedef std::pair<std::string, uint16_t> AspReadNameMapPair;
    AspReadNameMapType myReadNameMap;
};

#endif
