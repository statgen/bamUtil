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

#ifndef __SAM_RECORD_H__
#define __SAM_RECORD_H__

#include <stdint.h>

#include "SamStatus.h"
#include "LongHash.h"
#include "MathVector.h"
#include "StringArray.h"
#include "IntArray.h"
#include "SamFileHeader.h"
#include "CigarRoller.h"

struct bamRecordStruct
{
public:
    int32_t      myBlockSize;
    int32_t      myReferenceID;
    int32_t      myPosition;
    uint32_t     myReadNameLength : 8, myMapQuality : 8, myBin : 16;
    uint32_t     myCigarLength : 16, myFlag : 16;
    int32_t      myReadLength;
    int32_t      myMateReferenceID;
    int32_t      myMatePosition;
    int32_t      myInsertSize;             // Outer fragment length
    char  myData[1];
};

class SamRecord
{
public:
   
    /// Default Constructor.
    SamRecord();

    /// Constructor that sets the error handling type.
    /// \param errorHandlingType how to handle errors.
    SamRecord(ErrorHandler::HandlingType errorHandlingType);

    ~SamRecord();

    // Reset the fields of the record to a default value.
    void resetRecord();
    // Reset the tag iterator to the beginning of the tags.
    void resetTagIter();
 
    // Returns whether or not the record is valid.
    // Header is needed to perform some validation against it.
    // Sets the status to indicate success or failure.
    bool isValid(SamFileHeader& header);

    ///////////////////////
    // Set alignment data
    ///////////////////////
    // Set methods for record fields.  All of the "set" methods set the
    // status to indicate success or the failure reason.
    bool setReadName(const char* readName);
    bool setFlag(uint16_t flag);
    bool setReferenceName(SamFileHeader& header, 
                          const char* referenceName);
    bool set1BasedPosition(int32_t position);
    bool set0BasedPosition(int32_t position);
    bool setMapQuality(uint8_t mapQuality);
    bool setCigar(const char* cigar);
    bool setCigar(const Cigar& cigar);
    bool setMateReferenceName(SamFileHeader& header,
                              const char* mateReferenceName);
    bool set1BasedMatePosition(int32_t matePosition);
    bool set0BasedMatePosition(int32_t matePosition);
    bool setInsertSize(int32_t insertSize);
    bool setSequence(const char* seq);
    bool setQuality(const char* quality);

    // Read the BAM record from a file.
    SamStatus::Status setBufferFromFile(IFILE filePtr, SamFileHeader& header);

    // Read the BAM record from a file.
    SamStatus::Status setBuffer(const char* fromBuffer, uint32_t fromBufferSize,
                                SamFileHeader& header);

    // Add the specified tag to the record.
    // Returns true if the tag was successfully added, false otherwise.
    // Sets the status.
    bool addTag(const char* tag, char vtype, const char* value);

    // Get methods for record fields.  All of the "get" methods set the
    // status to indicate success or the failure reason.
    const void* getRecordBuffer();
    SamStatus::Status writeRecordBuffer(IFILE filePtr);
    int32_t getBlockSize();
    const char* getReferenceName();
    int32_t getReferenceID();
    int32_t get1BasedPosition();
    int32_t get0BasedPosition();
    uint8_t getReadNameLength();
    uint8_t getMapQuality();
    uint16_t getBin();
    uint16_t getCigarLength();
    uint16_t getFlag();
    int32_t getReadLength();

    // This method returns the mate reference name.  If it is equal to the
    // reference name, it still returns the reference name.
    const char* getMateReferenceName();

    // This method returns the mate reference name.  If it is equal to the
    // reference name, it returns "=", unless they are both "*" in which case
    // "*" is returned.
    const char* getMateReferenceNameOrEqual();
    int32_t getMateReferenceID();
    int32_t get1BasedMatePosition();
    int32_t get0BasedMatePosition();
    int32_t getInsertSize();

    // Returns the inclusive rightmost position of the clipped sequence.
    int32_t get0BasedAlignmentEnd();
    int32_t get1BasedAlignmentEnd();
   
    // Return the length of the alignment.
    int32_t getAlignmentLength();

    // Returns the inclusive left-most position adjust for clipped bases.
    int32_t get0BasedUnclippedStart();
    int32_t get1BasedUnclippedStart();
    // Returns the inclusive right-most position adjust for clipped bases.
    int32_t get0BasedUnclippedEnd();
    int32_t get1BasedUnclippedEnd();

    const char* getReadName();
    const char* getCigar();
    const char* getSequence();
    const char* getQuality();

    // Get the sequence base at the specified index into this sequence 0 to
    // readLength - 1.
    char getSequence(int index);

    // Get the quality char at the specified index into this quality 0 to
    // readLength - 1.
    char getQuality(int index);
   
    // TODO - want this to be getCigar
    Cigar* getCigarInfo();

    uint32_t getTagLength();

    // Sets the Status to SUCCESS when a tag is successfully returned or
    // when there are no more tags.  Otherwise the status is set to describe
    // why it failed (parsing, etc).
    bool getNextSamTag(char* tag, char& vtype, void** value);

    // Returns the values of all fields except the tags.
    bool getFields(bamRecordStruct& recStruct, String& readName, 
                   String& cigar, String& sequence, String& quality);

    // The following set of methods do not set the status.
    bool isIntegerType(char vtype) const;
    bool isDoubleType(char vtype) const;
    bool isCharType(char vtype) const;
    bool isStringType(char vtype) const;

    // The following set of methods do not set the status.
    void clearTags();
   
    // Returns the status associated with the last method
    // that sets the status.
    const SamStatus& getStatus();
    
    // The following set of methods do not set the status.
    String & getString(const char * tag);
    int &    getInteger(const char * tag);
    double & getDouble(const char * tag);


//     void getSamExtraFieldFromKey(int key, String& extraField);
    
    // The following set of methods do not set the status.
    bool checkString(const char * tag)    { return checkTag(tag, 'Z'); }
    bool checkInteger(const char * tag)   { return checkTag(tag, 'i'); }
    bool checkDouble(const char * tag)    { return checkTag(tag, 'f'); }
    bool checkTag(const char * tag, char type);


    
    // Return the number of bases in this read that overlap the passed in
    // region.
    // start : inclusive 0-based start position (reference position) of the
    //         region to check for overlaps in.
    //         (-1 indicates to start at the beginning of the reference.)
    // end   : exclusive 0-based end position (reference position) of the
    //          region to check for overlaps in.
    //         (-1 indicates to go to the end of the reference.)
    // Returns the number of overlapping bases
    // (matches in the cigar - not skips/deletions)
    uint32_t getNumOverlaps(int32_t start, int32_t end);


private:
    static int MAKEKEY(char ch1, char ch2, char type)
    { return (type << 16) + (ch2 << 8) + ch1; }

    // Allocate space for the record - does a realloc.  
    // The passed in size is the size of the entire record including the
    // block size field.
    // Adds any errors to myStatus.
    bool allocateRecordStructure(int size);


    void* getStringPtr(int offset);
    void* getIntegerPtr(int offset);
    void* getDoublePtr(int offset);

    // Fixes the buffer to match the variable length fields.
    // Adds any errors to myStatus.
    bool fixBuffer();

    // Sets the Sequence and Quality strings from the buffer.
    // They are done together in one method because they require the same
    // loop, so might as well be done at the same time.
    // Adds any errors to myStatus.
    void setSequenceAndQualityFromBuffer();

    // Parse the cigar to calculate the alignment/unclipped ends and convert
    // to SAM/BAM format.
    // Adds any errors to myStatus.
    bool parseCigar();
    // Parse the cigar string to calculate the cigar length and alignment end
    // and convert to SAM format.
    // Adds any errors to myStatus.
    bool parseCigarBinary();
    // Parse the cigar string to calculate the cigar length and alignment end
    // and convert to BAM format.
    // Adds any errors to myStatus.
    bool parseCigarString();

    // Set the tags from the buffer.
    // Adds any errors to myStatus.
    bool setTagsFromBuffer();

    // Set the tags in the buffer.
    // Adds any errors to myStatus.
    bool setTagsInBuffer();

    void setVariablesForNewBuffer(SamFileHeader& header);

    void getVtype(int key, char& vtype) const;
    void getTag(int key, char* tag) const;

    String & getString(int offset);
    int &    getInteger(int offset);
    double & getDouble(int offset);

    static const int DEFAULT_BLOCK_SIZE = 40;
    static const int DEFAULT_BIN = 4680;
    static const int DEFAULT_READ_NAME_LENGTH = 8;
    static const char* DEFAULT_READ_NAME;
    static const char* FIELD_ABSENT_STRING;

    bamRecordStruct * myRecordPtr;
    int allocatedSize;

    // Pointer to a temporary cigar buffer that can be used during string
    // parsing before it is ready to be copied into the actual record.
    uint32_t* myCigarTempBuffer;

    // Size of the currently allocated temporary cigar buffer.
    int myCigarTempBufferAllocatedSize;

    // Length of the cigar currently contained in the temporary buffer.
    int myCigarTempBufferLength;

    // Track if the buffer is in sync with the Strings/Tags.
    // Set to false if any of the variable length fields are modified.
    // Set to true when the buffer is updated to match the variable length
    // fields.
    bool myIsBufferSynced;

    // Track if the tags need to be set from the buffer.
    bool myNeedToSetTagsFromBuffer;

    // Trag if the tags need to be set in the buffer.
    // Allows you to set just the tags if they are the only thing that changed
    // in the buffer.
    bool myNeedToSetTagsInBuffer;

    int myTagBufferSize;
    int myLastTagIndex;

    String myReadName;
    String myReferenceName;
    String myMateReferenceName;
    String myCigar;
    String mySequence;
    String myQuality;

    // The length of the alignment.
    int32_t myAlignmentLength;
    // Unclipped alignment positions.
    int32_t myUnclippedStartOffset;
    int32_t myUnclippedEndOffset;
    
    CigarRoller myCigarRoller;

    LongHash<int>  extras;
    StringArray    strings;
    IntArray       integers;
    Vector         doubles;


    // Track whether or not the buffer values are correct for
    // each setting.
    bool myIsReadNameBufferValid;
    bool myIsCigarBufferValid;
    bool myIsSequenceBufferValid;
    bool myIsQualityBufferValid;
    bool myIsTagsBufferValid;
    bool myIsBinValid;

    SamStatus myStatus;

    String NOT_FOUND_TAG_STRING;
    int NOT_FOUND_TAG_INT;
    double NOT_FOUND_TAG_DOUBLE;
};

#endif
