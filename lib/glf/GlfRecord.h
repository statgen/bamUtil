/*
 *  Copyright (C) 2010  Regents of the University of Michigan
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

#ifndef __GLF_RECORD_H__
#define __GLF_RECORD_H__

#include <map>
#include <stdint.h>

#include "InputFile.h" 
#include "CharBuffer.h"

class GlfRecord
{
public:
    GlfRecord();
    ~GlfRecord();

//     // Copy Constructor   
//     GlfRecord(const GlfRecord& record);

//     // Overload operator = to copy the passed in record into this record.
//     GlfRecord & operator = (const GlfRecord& record);

//     // Overload operator = to copy the passed in record into this record.
//     bool copy(const GlfRecord& record);

    void reset();
   
    // Read the record from the specified file.  Assumes the file is in
    // the correct position for reading the record.
    bool read(IFILE filePtr);

    // Write the record to the specified file.
    bool write(IFILE filePtr) const;

    void print() const;

    // Accessors to set the generic values.
    bool setRtypeRef(uint8_t rtypeRef);
    bool setRecordType(uint8_t recType);
    bool setRefBaseInt(uint8_t refBase);
    // TODO   bool setRefBaseChar(char refBase);
    bool setOffset(uint32_t offset);
    bool setMinDepth(uint32_t minDepth);
    bool setMinLk(uint8_t minLk);
    bool setReadDepth(uint32_t readDepth);
    bool setRmsMapQ(uint8_t rmsMapQ);
    
    // Accessors for setting record type 1;
    //bool setType1(all fields for type 1);
    bool setLk(int index, uint8_t value);

    // Accessors for setting record type 2
//     bool setType2(all fields for type 2);
    bool setLkHom1(uint8_t lk);
    bool setLkHom2(uint8_t lk);
    bool setLkHet(uint8_t lk);
    bool setInsertionIndel1(const std::string& indelSeq);
    bool setDeletionIndel1(const std::string& indelSeq);
    bool setInsertionIndel2(const std::string& indelSeq);
    bool setDeletionIndel2(const std::string& indelSeq);

    // Accessors to get the generic values.
    inline int getRecordType() const
    {
        return(myRecTypeRefBase >> REC_TYPE_SHIFT);
    }
    inline int getRefBase() const
    {
        return(myRecTypeRefBase & REF_BASE_MASK);
    }

    char getRefBaseChar() const;
    uint32_t getOffset();
    uint8_t getMinDepth();
    uint8_t getMinLk();
    uint32_t getReadDepth();
    uint8_t getRmsMapQ();
    
    // Accessors for getting record type 1;
    //bool getType1(all fields for type 1);
    uint8_t getLk(int index);
    
    //     // Accessors for setting record type 1;
    //     bool setType2(all fields for type 2);
    uint8_t getLkHom1();
    uint8_t getLkHom2();
    uint8_t getLkHet();
    int16_t getIndel1(std::string& indelSeq);
    int16_t getIndel2(std::string& indelSeq);
    

private:
    // Read a record of record type 1.
    void readType1(IFILE filePtr);

    // Read a record of record type 2.
    void readType2(IFILE filePtr);


    // Write the rtyperef field.
    void writeRtypeRef(IFILE filePtr) const;


    // Write a record of record type 1.
    void writeType1(IFILE filePtr) const;

    // Write a record of record type 2.
    void writeType2(IFILE filePtr) const;

    // Contains record_type and ref_base.
    uint8_t myRecTypeRefBase;

    static const uint8_t REC_TYPE_SHIFT = 4;
    static const uint8_t REF_BASE_MASK = 0xF;
    static const uint8_t REC_TYPE_MASK = 0xF0;

    static const uint32_t MIN_LK_SHIFT = 24;
    static const uint32_t READ_DEPTH_MASK = 0xFFFFFF;
    static const uint32_t MIN_LK_MASK = 0xFF000000;

    static const char REF_BASE_MAX = 15;
    static std::string REF_BASE_CHAR;

    static const int NUM_REC1_LIKELIHOOD = 10;

    struct
    {
        uint32_t offset;
        uint32_t min_depth;
        uint8_t rmsMapQ;
        uint8_t lk[NUM_REC1_LIKELIHOOD];
    } myRec1Base;

    static const int REC1_BASE_SIZE = 19;

    struct
    {
        uint32_t offset;
        uint32_t min_depth;
        uint8_t rmsMapQ;
        uint8_t lkHom1;
        uint8_t lkHom2;
        uint8_t lkHet;
        int16_t indelLen1;
        int16_t indelLen2;
    } myRec2Base;

    // TODO rest of rec 2.
    CharBuffer myIndelSeq1;
    CharBuffer myIndelSeq2;

    static const int REC2_BASE_SIZE = 16;

};

#endif
