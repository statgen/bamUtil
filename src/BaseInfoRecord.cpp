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

IFILE BaseInfoRecord::ourOutput = NULL;
const unsigned int BaseInfoRecord::REC_TYPE_LEN = 1;
const uint8_t BaseInfoRecord::EMPTY_REC = 0x0;
const uint8_t BaseInfoRecord::POS_REC = 0x1;
const uint8_t BaseInfoRecord::REF_ONLY_REC = 0x2;
const uint8_t BaseInfoRecord::DETAILED_REC = 0x3;

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
    int basesIndex = (myNumBases + 1)/2;
    if((myNumBases % 2) == 0)
    {
        // When myNumBases is a multiple of 2, 
        // we want the upper bits of the index.
        myBases[basesIndex] = (base && 0xF) << 4;
    }
    else
    {
        // Not a multiple of 2, so add it to the lower bits.
        myBases[basesIndex] |= (base && 0xF);
    }
    myQuals[myNumBases] = qual;
    myCycles[myNumBases] = cycle;

    if((myNumBases % 8) == 0)
    {
        // Because they will be written on a byte boundary, 
        // clear out the other strands in this boundary.
        myStrands[myNumBases+1] = 0;
        myStrands[myNumBases+2] = 0;
        myStrands[myNumBases+3] = 0;
        myStrands[myNumBases+4] = 0;
        myStrands[myNumBases+5] = 0;
        myStrands[myNumBases+6] = 0;
        myStrands[myNumBases+7] = 0;
    }
    myStrands[myNumBases] = strand;
    myMQs[myNumBases] = mq;
    ++myNumBases;
}


void BaseInfoRecord::reset()
{
    myNumBases = 0;
    myChromID = INVALID_CHROMOSOME_INDEX;
    my0BasedPos = -1;
}


bool read(filePtr)
{
    reset();
    // Read the first byte that contains the type.
    uint8_t type;
    if(ifread(filePtr, type, REC_TYPE_LEN) != REC_TYPE_LEN)
    {
        TODO error;
        return(false);
    }

    switch(type)
    {
        case EMPTY_REC:
            // This is an empty rec which has been entirely read, so just retun.
            return(true);
            break;
        case POS_REC:
            // Position record type.
            return(readPosRecord());
            break;
        case REF_ONLY_REC:
            return(readRefOnlyRecord());
            break;
        case DETAILED_REC:
            return(readDetailedRecord());
            break;
        default:
            TODO error;
            return(false);
    }
}


bool BaseInfoRecord::isEmpty()
{
    if((myNumBases == 0) && (myChromID == INVALID_CHROMOSOME_INDEX))
    {
        return(true);
    }
    return(false);
}


bool BaseInfoRecord::isPos()
{
    if((myNumBases == 0) && (myChromID != INVALID_CHROMOSOME_INDEX))
    {
        return(true);
    }
    return(false);
}
bool BaseInfoRecord::isRefOnly()
{
    if((myNumBases != 0) && ())
    {
        return(true);
    }
    return(false);
}
bool BaseInfoRecord::isDetailed()
{
    if((myNumBases != 0) && ())
    {
        return(true);
    }
    return(false);
}


void BaseInfoRecord::writeEmpty(IFILE outputFile)
{
    ifwrite(outputFile, &EMPTY_REC, REC_TYPE_LEN);
}


void BaseInfoRecord::writePos(int32_t chrom, int32_t pos, IFILE outputFile)
{
    ifwrite(outputFile, &POS_REC, REC_TYPE_LEN);
    ifwrite(outputFile, &chrom, sizeof(int32_t));
    ifwrite(outputFile, &pos, sizeof(int32_t));
}


void BaseInfoRecord::writeRefOnly(IFILE outputFile)
{
    ifwrite(outputFile, &REF_ONLY_REC, REC_TYPE_LEN);
    int32_t glh = 0;
    int32_t gla = 0;
    
    ifwrite(outputFile, &myNumBases, sizeof(myNumBases));
    ifwrite(outputFile, &glh, sizeof(int32_t));
    ifwrite(outputFile, &gla, sizeof(int32_t));
}


void BaseInfoRecord::writeDetailed(IFILE outputFile)
{
    ifwrite(outputFile, &DETAILED_REC, REC_TYPE_LEN);

    // Write the number of bases.
    ifwrite(outputFile, &myNumBases, sizeof(myNumBases));
    int basesSize = (myNumBases+1)/2;
    ifwrite(outputFile, myBases, basesSize);

    ifwrite(outputFile, myQuals, myNumBases);
    ifwrite(outputFile, myCycles, myNumBases);
    ifwrite(outputFile, myStrands, (myNumBases+7)/8);
    ifwrite(outputFile, myMQs, myNumBases);
}
