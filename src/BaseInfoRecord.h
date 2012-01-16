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


#ifndef __BASE_INFO_RECORD_H__
#define __BASE_INFO_RECORD_H__

#include "StringArray.h"
#include "BamExecutable.h"
#include "SamFile.h"

class BaseInfoRecord
{
public:
    BaseInfoRecord();
    
    ~BaseInfoRecord();
    
    void add(char base, int qual, int cycle, bool strand, int mq);

    void reset();

private:
    bool myAllRefBase;


    TODO - these are max size of 255, so make them arrays of 255.
        Maybe even as possible such that they can be directly written out to avoid a 2nd loop?

    std::vector<int> myBases;
    std::vector<int> myQuals;
    std::vector<int> myCycles;
    std::vector<bool> myStrands;
    std::vector<uint8_t> myMQs;
};


#endif
