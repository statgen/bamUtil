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

#include "SamFile.h"

class BaseInfoRecord
{
public:
    BaseInfoRecord();
    
    ~BaseInfoRecord();
    
    static bool setOutputFile(const char* fileName);

    void add(char base, int qual, int cycle, bool strand, int mq);

    void reset();

    //////////////////////////////////
    // Reading methods.

    // Read a record from the file.  It is assumed the file is in the 
    // correct position.
    bool read(IFILE filePtr);

    // 
    bool isEmpty();
    bool isPos();
    bool isRefOnly();
    bool isDetailed();

    /////////////////////////////////
    // Accessor methods
    

    ////////////////////////////////
    // Writing methods.

    static void writeEmpty(IFILE outputFile);
    static void writePos(int32_t chrom, int32_t pos, IFILE outputFile);

    void writeRefOnly(IFILE outputFile);
    void writeDetailed(IFILE outputFile);

private:
    static IFILE ourOutput;
    static const unsigned int REC_TYPE_LEN;
    static const uint8_t EMPTY_REC;
    static const uint8_t POS_REC;
    static const uint8_t REF_ONLY_REC;
    static const uint8_t DETAILED_REC;
    static const int MAX_NUM_BASES = 255;

    int myNumBases;

    int myGLH;
    int myGLA;

    // Since each base is only 4 bits, each index holds two bases.
    // The earlier base is in the upper bits.
    int8_t myBases[(MAX_NUM_BASES+1)/2];
    int8_t myQuals[MAX_NUM_BASES];
    int8_t myCycles[MAX_NUM_BASES];

    // Since strands are only 1 bit, each index holds 4 strands
    // with the first strand in the uppermost bit.
    bool myStrands[MAX_NUM_BASES];
    uint8_t myMQs[MAX_NUM_BASES];

    // The choromsome & position for this record.
    int32_t myChormID;
    int32_t my0BasedPos;
};


#endif
