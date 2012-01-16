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
//////////////////////////////////////////////////////////////////////////
#include "BaseInfoRecord.h"

/////////////////////////////////////////////////////////////////////////////
//
// BaseInfoRecord
//
BaseInfoRecord::BaseInfoRecord()
{
    reset();
}


BaseInfoRecord::~BaseInfoRecord()
{
}


// Add an entry
void BaseInfoRecord::add(char base, int qual, int cycle, bool strand, int mq)
{
    myBases.push_back(base);
    myQuals.push_back(qual);
    myCycles.push_back(cycle);
    myStrands.push_back(strand);
    myMQs.push_back(mq);
}


void BaseInfoRecord::reset()
{
    myBases.clear();
    myQuals.clear();
    myCycles.clear();
    myStrands.clear();
    myMQs.clear();
    myAllRefBase = true;
}

