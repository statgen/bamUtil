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

#include <stdlib.h>
#include <limits>
#include <stdexcept>

#include "bam.h"

#include "SamRecord.h"
#include "SamValidation.h"

#include "BaseUtilities.h"
#include "SamQuerySeqWithRefHelper.h"

const char* SamRecord::DEFAULT_READ_NAME = "UNKNOWN";
const char* SamRecord::FIELD_ABSENT_STRING = "=";


SamRecord::SamRecord()
    : myStatus(),
      myRefPtr(NULL),
      mySequenceTranslation(NONE)
{
    int32_t defaultAllocSize = DEFAULT_BLOCK_SIZE + sizeof(int32_t);

    myRecordPtr = 
        (bamRecordStruct *) malloc(defaultAllocSize);

    myCigarTempBuffer = NULL;
    myCigarTempBufferAllocatedSize = 0;

    allocatedSize = defaultAllocSize;

    resetRecord();
}


SamRecord::SamRecord(ErrorHandler::HandlingType errorHandlingType)
    : myStatus(errorHandlingType),
      myRefPtr(NULL),
      mySequenceTranslation(NONE)
{
    int32_t defaultAllocSize = DEFAULT_BLOCK_SIZE + sizeof(int32_t);

    myRecordPtr = 
        (bamRecordStruct *) malloc(defaultAllocSize);

    myCigarTempBuffer = NULL;
    myCigarTempBufferAllocatedSize = 0;

    allocatedSize = defaultAllocSize;

    resetRecord();
}


SamRecord::~SamRecord()
{
    resetRecord();

    if(myRecordPtr != NULL)
    {
        free(myRecordPtr);
        myRecordPtr = NULL;
    }
    if(myCigarTempBuffer != NULL)
    {
        free(myCigarTempBuffer);
        myCigarTempBuffer = NULL;
        myCigarTempBufferAllocatedSize = 0;
    }
}


// Resets the fields of the record to a default value.
void SamRecord::resetRecord()
{
    myIsBufferSynced = true;

    myRecordPtr->myBlockSize = DEFAULT_BLOCK_SIZE;
    myRecordPtr->myReferenceID = -1;
    myRecordPtr->myPosition = -1;
    myRecordPtr->myReadNameLength = DEFAULT_READ_NAME_LENGTH;
    myRecordPtr->myMapQuality = 0;
    myRecordPtr->myBin = DEFAULT_BIN;
    myRecordPtr->myCigarLength = 0;
    myRecordPtr->myFlag = 0;
    myRecordPtr->myReadLength = 0;
    myRecordPtr->myMateReferenceID = -1;
    myRecordPtr->myMatePosition = -1;
    myRecordPtr->myInsertSize = 0;
   
    // Set the sam values for the variable length fields.
    // TODO - one way to speed this up might be to not set to "*" and just
    // clear them, and write out a '*' for SAM if it is empty.
    myReadName = DEFAULT_READ_NAME;
    myReferenceName = "*";
    myMateReferenceName = "*";
    myCigar = "*";
    mySequence = "*";
    mySeqWithEq.clear();
    mySeqWithoutEq.clear();
    myQuality = "*";
    myNeedToSetTagsFromBuffer = false;
    myNeedToSetTagsInBuffer = false;

    // Initialize the calculated alignment info to the uncalculated value.
    myAlignmentLength = -1;
    myUnclippedStartOffset = -1;
    myUnclippedEndOffset = -1;

    clearTags();

    // Set the bam values for the variable length fields.
    // Only the read name needs to be set, the others are a length of 0.
    // Set the read name.  The min size of myRecordPtr includes the size for
    // the default read name.
    memcpy(&(myRecordPtr->myData), myReadName.c_str(), 
           myRecordPtr->myReadNameLength);

    // Set that the variable length buffer fields are valid.
    myIsReadNameBufferValid = true;
    myIsCigarBufferValid = true;
    myIsSequenceBufferValid = true;
    myBufferSequenceTranslation = NONE;
    myIsQualityBufferValid = true;
    myIsTagsBufferValid = true;
    myIsBinValid = true;

    myCigarTempBufferLength = -1;

    myStatus = SamStatus::SUCCESS;

    NOT_FOUND_TAG_STRING = "";
    NOT_FOUND_TAG_INT = -1;
    NOT_FOUND_TAG_DOUBLE = -1;
}


// Reset the tag iterator to the beginning of the tags.
void SamRecord::resetTagIter()
{
    myLastTagIndex = -1;
}


// Returns whether or not the record is valid.
// Header is needed to perform some validation against it.
bool SamRecord::isValid(SamFileHeader& header)
{
    myStatus = SamStatus::SUCCESS;
    SamValidationErrors invalidSamErrors;
    if(!SamValidator::isValid(header, *this, invalidSamErrors))
    {
        // The record is not valid.
        std::string errorMessage = "";
        invalidSamErrors.getErrorString(errorMessage);
        myStatus.setStatus(SamStatus::INVALID, errorMessage.c_str());
        return(false);
    }
    // The record is valid.
    return(true);
}


// Read the BAM record from a file.
SamStatus::Status SamRecord::setBufferFromFile(IFILE filePtr, 
                                               SamFileHeader& header)
{
    myStatus = SamStatus::SUCCESS;
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return failure.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Can't read from an unopened file.");
        return(SamStatus::FAIL_ORDER);
    }

    // Clear the record.
    resetRecord();

    // read the record size.
    int numBytes = 
        ifread(filePtr, &(myRecordPtr->myBlockSize), sizeof(int32_t));

    if(ifeof(filePtr))
    {
        if(numBytes == 0)
        {
            // End of file, nothing was read, no more records.
            myStatus.setStatus(SamStatus::NO_MORE_RECS,
                               "No more records left to read.");
            return(SamStatus::NO_MORE_RECS);
        }
        else
        {
            // Error: end of the file reached prior to reading the rest of the
            // record.
            myStatus.setStatus(SamStatus::FAIL_PARSE, 
                               "EOF reached in the middle of a record.");
            return(SamStatus::FAIL_PARSE);
        }
    }

    // allocate space for the record size.
    if(!allocateRecordStructure(myRecordPtr->myBlockSize + sizeof(int32_t)))
    {
        // Failed to allocate space.
        // Status is set by allocateRecordStructure.
        return(SamStatus::FAIL_MEM);
    }

    // Read the rest of the alignment block, starting at the reference id.
    if(ifread(filePtr, &(myRecordPtr->myReferenceID), myRecordPtr->myBlockSize)
       != (unsigned int)myRecordPtr->myBlockSize)
    {
        // Error reading the record.  Reset it and return failure.
        resetRecord();
        myStatus.setStatus(SamStatus::FAIL_IO,
                           "Failed to read the record");
        return(SamStatus::FAIL_IO);
    }

    setVariablesForNewBuffer(header);

    // Return the status of the record.
    return(SamStatus::SUCCESS);
}


void SamRecord::setReference(GenomeSequence* reference)
{
    myRefPtr = reference;
}


// Set the type of sequence translation to use when getting
// the sequence.  The default type (if this method is never called) is
// NONE (the sequence is left as-is).  This is used 
void SamRecord::setSequenceTranslation(SequenceTranslation translation)
{
    mySequenceTranslation = translation;
}


bool SamRecord::setReadName(const char* readName) 
{
    myReadName = readName;
    myIsBufferSynced = false;
    myIsReadNameBufferValid = false;
    myStatus = SamStatus::SUCCESS;

    // The read name must at least have some length, otherwise this is a parsing
    // error.
    if(myReadName.Length() == 0)
    {
        // Invalid - reset ReadName return false.
        myReadName = DEFAULT_READ_NAME;
        myRecordPtr->myReadNameLength = DEFAULT_READ_NAME_LENGTH;
        myStatus.setStatus(SamStatus::INVALID, "0 length Query Name.");
        return(false);
    }

    return true;
}


bool SamRecord::setFlag(uint16_t flag)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myFlag = flag;
    return true;
}


bool SamRecord::setReferenceName(SamFileHeader& header,
                                 const char* referenceName) 
{
    myStatus = SamStatus::SUCCESS;

    myReferenceName = referenceName;
    myRecordPtr->myReferenceID = header.getReferenceID(referenceName);

    return true;
}


bool SamRecord::set1BasedPosition(int32_t position) 
{
    return(set0BasedPosition(position - 1));
}


bool SamRecord::set0BasedPosition(int32_t position)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myPosition = position;
    myIsBinValid = false;
    return true;
}


bool SamRecord::setMapQuality(uint8_t mapQuality)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myMapQuality = mapQuality;
    return true;
}


bool SamRecord::setCigar(const char* cigar)
{
    myStatus = SamStatus::SUCCESS;
    myCigar = cigar;
 
    myIsBufferSynced = false;
    myIsCigarBufferValid = false;
    myCigarTempBufferLength = -1;
    myIsBinValid = false;

    // Initialize the calculated alignment info to the uncalculated value.
    myAlignmentLength = -1;
    myUnclippedStartOffset = -1;
    myUnclippedEndOffset = -1;

    return true;
}


bool SamRecord::setCigar(const Cigar& cigar)
{
    myStatus = SamStatus::SUCCESS;
    cigar.getCigarString(myCigar);
 
    myIsBufferSynced = false;
    myIsCigarBufferValid = false;
    myCigarTempBufferLength = -1;
    myIsBinValid = false;

    // Initialize the calculated alignment info to the uncalculated value.
    myAlignmentLength = -1;
    myUnclippedStartOffset = -1;
    myUnclippedEndOffset = -1;

    return true;
}


bool SamRecord::setMateReferenceName(SamFileHeader& header,
                                     const char* mateReferenceName) 
{
    myStatus = SamStatus::SUCCESS;
    // Set the mate reference, if it is "=", set it to be equal
    // to myReferenceName.  This assumes that myReferenceName has already
    // been called.
    if(strcmp(mateReferenceName, FIELD_ABSENT_STRING) == 0)
    {
        myMateReferenceName = myReferenceName;
    }
    else
    {
        myMateReferenceName = mateReferenceName;
    }

    // Set the Mate Reference ID.
    myRecordPtr->myMateReferenceID = 
        header.getReferenceID(myMateReferenceName);

    return true;
}


bool SamRecord::set1BasedMatePosition(int32_t matePosition)
{
    return(set0BasedMatePosition(matePosition - 1));
}


bool SamRecord::set0BasedMatePosition(int32_t matePosition)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myMatePosition = matePosition;
    return true;
}


bool SamRecord::setInsertSize(int32_t insertSize)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myInsertSize = insertSize;
    return true;
}


bool SamRecord::setSequence(const char* seq) 
{
    myStatus = SamStatus::SUCCESS;
    mySequence = seq;
    mySeqWithEq.clear();
    mySeqWithoutEq.clear();
   
    myIsBufferSynced = false;
    myIsSequenceBufferValid = false;
    return true;
}


bool SamRecord::setQuality(const char* quality) 
{
    myStatus = SamStatus::SUCCESS;
    myQuality = quality;
    myIsBufferSynced = false;
    myIsQualityBufferValid = false;
    return true;
}


// Set the BAM record from the passeed in buffer of the specified size.
// Note: The size includes the block size.
SamStatus::Status SamRecord::setBuffer(const char* fromBuffer,
                                       uint32_t fromBufferSize,
                                       SamFileHeader& header)
{
    myStatus = SamStatus::SUCCESS;
    if((fromBuffer == NULL) || (fromBufferSize == 0))
    {
        // Buffer is empty.
        myStatus.setStatus(SamStatus::FAIL_PARSE,
                           "Cannot parse an empty file.");
        return(SamStatus::FAIL_PARSE);
    }

    // Clear the record.   
    resetRecord();

    // allocate space for the record size.
    if(!allocateRecordStructure(fromBufferSize))
    {
        // Failed to allocate space.
        return(SamStatus::FAIL_MEM);
    }
   
    memcpy(myRecordPtr, fromBuffer, fromBufferSize);

    setVariablesForNewBuffer(header);

    // Return the status of the record.
    return(SamStatus::SUCCESS);
}


// Add the specified tag to the record.
// Returns true if the tag was successfully added, false otherwise.
bool SamRecord::addTag(const char* tag, char vtype, const char* valuePtr)
{
    myStatus = SamStatus::SUCCESS;
    bool status = true; // default to successful.
    int key = 0;
    int intVal = 0;
    int index = 0;
    char bamvtype;

    int tagBufferSize = 0;

    // First check to see if the tags need to be synced to the buffer.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read tags from the buffer, so cannot add new ones.
            return(false);
        }
    }

    switch (vtype)
    {
        case 'A' :
            index = integers.Length();
            bamvtype = vtype;
            integers.Push((const int)*(valuePtr));
            tagBufferSize += 4;
            break;
        case 'i' :
            index = integers.Length();
            intVal = atoi((const char *)valuePtr);
            // Ints come in as int.  But it can be represented in fewer bits.
            // So determine a more specific type that is in line with the
            // types for BAM files.
            // First check to see if it is a negative.
            if(intVal < 0)
            {
                // The int is negative, so it will need to use a signed type.
                // See if it is greater than the min value for a char.
                if(intVal > std::numeric_limits<char>::min())
                {
                    // It can be stored in a signed char.
                    bamvtype = 'c';
                    tagBufferSize += 4;
                }
                else if(intVal > std::numeric_limits<short>::min())
                {
                    // It fits in a signed short.
                    bamvtype = 's';
                    tagBufferSize += 5;
                }
                else
                {
                    // Just store it as a signed int.
                    bamvtype = 'i';
                    tagBufferSize += 7;
                }
            }
            else
            {
                // It is positive, so an unsigned type can be used.
                if(intVal < std::numeric_limits<unsigned char>::max())
                {
                    // It is under the max of an unsigned char.
                    bamvtype = 'C';
                    tagBufferSize += 4;
                }
                else if(intVal < std::numeric_limits<unsigned short>::max())
                {
                    // It is under the max of an unsigned short.
                    bamvtype = 'S';
                    tagBufferSize += 5;
                }
                else
                {
                    // Just store it as an unsigned int.
                    bamvtype = 'I';
                    tagBufferSize += 7;
                }
            }
            integers.Push(intVal);
            break;
        case 'Z' :
            index = strings.Length();
            bamvtype = vtype;
            strings.Push((const char *)valuePtr);
            tagBufferSize += 4 + strings.Last().Length();
            break;
        case 'f' :
            index = doubles.Length();
            bamvtype = vtype;
            doubles.Push(atof((const char *)valuePtr));
            tagBufferSize += 7;
            break;
        default :
            fprintf(stderr,
                    "samFile::ReadSAM() - Unknown custom field of type %c\n",
                    vtype);
            myStatus.setStatus(SamStatus::FAIL_PARSE, 
                               "Unknown custom field in a tag");
            status = false;
    }

    // Only add the tag if it has so far been successfully processed.
    if(status)
    {
        // The buffer tags are now out of sync.
        myNeedToSetTagsInBuffer = true;
        myIsTagsBufferValid = false;
        myIsBufferSynced = false;

        key = MAKEKEY(tag[0], tag[1], bamvtype);
        extras.Add(key, index);
        myTagBufferSize += tagBufferSize;
    }
    return(status);
}


// Get methods for record fields.
const void* SamRecord::getRecordBuffer()
{
    return(getRecordBuffer(mySequenceTranslation));
}


// Get methods for record fields.
const void* SamRecord::getRecordBuffer(SequenceTranslation translation) 
{
    myStatus = SamStatus::SUCCESS;
    bool status = true;
    // If the buffer is not synced or the sequence in the buffer is not
    // properly translated, fix the buffer.
    if((myIsBufferSynced == false) ||
       (myBufferSequenceTranslation != translation))
    {
        status &= fixBuffer(translation);
    }
    // If the buffer is synced, check to see if the tags need to be synced.
    if(myNeedToSetTagsInBuffer)
    {
        status &= setTagsInBuffer();
    }
    if(!status)
    {
        return(NULL);
    }
    return (const void *)myRecordPtr;
}


// Write the record as a buffer into the file using the class's 
// sequence translation setting.
SamStatus::Status SamRecord::writeRecordBuffer(IFILE filePtr)
{
    return(writeRecordBuffer(filePtr, mySequenceTranslation));
}


// Write the record as a buffer into the file using the specified translation.
SamStatus::Status SamRecord::writeRecordBuffer(IFILE filePtr, 
                                               SequenceTranslation translation)
{
    myStatus = SamStatus::SUCCESS;
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return failure.
        myStatus.setStatus(SamStatus::FAIL_ORDER,
                           "Can't write to an unopened file.");
        return(SamStatus::FAIL_ORDER);
    }

    if((myIsBufferSynced == false) ||
       (myBufferSequenceTranslation != translation))
    {
        if(!fixBuffer(translation))
        {
            return(myStatus.getStatus());
        }
    }

    // Write the record.
    unsigned int numBytesToWrite = myRecordPtr->myBlockSize + sizeof(int32_t);
    unsigned int numBytesWritten = 
        ifwrite(filePtr, myRecordPtr, numBytesToWrite);

    // Return status based on if the correct number of bytes were written.
    if(numBytesToWrite == numBytesWritten)
    {
        return(SamStatus::SUCCESS);
    }
    // The correct number of bytes were not written.
    myStatus.setStatus(SamStatus::FAIL_IO, "Failed to write the entire record.");
    return(SamStatus::FAIL_IO);
}


int32_t SamRecord::getBlockSize() 
{
    myStatus = SamStatus::SUCCESS;
    // If the buffer isn't synced, sync the buffer to determine the
    // block size.
    if(myIsBufferSynced == false)
    {
        // Since this just returns the block size, the translation of
        // the sequence does not matter, so just use the currently set
        // value.
        fixBuffer(myBufferSequenceTranslation);
    }
    return myRecordPtr->myBlockSize;
}


// This method returns the reference name.
const char* SamRecord::getReferenceName()
{
    myStatus = SamStatus::SUCCESS;
    return myReferenceName.c_str();
}


int32_t SamRecord::getReferenceID()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myReferenceID;
}


int32_t SamRecord::get1BasedPosition()
{
    myStatus = SamStatus::SUCCESS;
    return (myRecordPtr->myPosition + 1);
}


int32_t SamRecord::get0BasedPosition()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myPosition;
}


uint8_t SamRecord::getReadNameLength() 
{
    myStatus = SamStatus::SUCCESS;
    // If the buffer is valid, return the size from there, otherwise get the 
    // size from the string length + 1 (ending null).
    if(myIsReadNameBufferValid)
    {
        return(myRecordPtr->myReadNameLength);
    }
   
    return(myReadName.Length() + 1);
}


uint8_t SamRecord::getMapQuality()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myMapQuality;
}


uint16_t SamRecord::getBin()
{
    myStatus = SamStatus::SUCCESS;
    if(!myIsBinValid)
    {
        // The bin that is set in the record is not valid, so
        // reset it.
        myRecordPtr->myBin = 
            bam_reg2bin(myRecordPtr->myPosition, get1BasedAlignmentEnd());      
        myIsBinValid = true;
    }
    return(myRecordPtr->myBin);
}


uint16_t SamRecord::getCigarLength()
{
    myStatus = SamStatus::SUCCESS;
    // If the cigar buffer is valid
    // then get the length from there.
    if(myIsCigarBufferValid)
    {
        return myRecordPtr->myCigarLength;      
    }

    if(myCigarTempBufferLength == -1)
    {
        // The cigar buffer is not valid and the cigar temp buffer is not set,
        // so parse the string.
        parseCigarString();
    }
   
    // The temp buffer is now set, so return the size.
    return(myCigarTempBufferLength);
}


uint16_t SamRecord::getFlag()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myFlag;
}


int32_t SamRecord::getReadLength() 
{
    myStatus = SamStatus::SUCCESS;
    if(myIsSequenceBufferValid == false)
    {
        // If the sequence is "*", then return 0.
        if((mySequence.Length() == 1) && (mySequence[0] == '*'))
        {
            return(0);
        }
        // Do not add 1 since it is not null terminated.
        return(mySequence.Length());
    }
    return(myRecordPtr->myReadLength);
}


// This method returns the mate reference name.  If it is equal to the
// reference name, it still returns the reference name.
const char* SamRecord::getMateReferenceName()
{
    myStatus = SamStatus::SUCCESS;
    return myMateReferenceName.c_str();
}


// This method returns the mate reference name.  If it is equal to the
// reference name, it returns "=", unless they are both "*" in which case
// "*" is returned.
const char* SamRecord::getMateReferenceNameOrEqual()
{
    myStatus = SamStatus::SUCCESS;
    if(myMateReferenceName == "*")
    {
        return(myMateReferenceName);
    }
    if(myMateReferenceName == getReferenceName())
    {
        return(FIELD_ABSENT_STRING);
    }
    else
    {
        return(myMateReferenceName);
    }
}


int32_t SamRecord::getMateReferenceID()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myMateReferenceID;
}


int32_t SamRecord::get1BasedMatePosition()
{
    myStatus = SamStatus::SUCCESS;
    return (myRecordPtr->myMatePosition + 1);
}


int32_t SamRecord::get0BasedMatePosition()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myMatePosition;
}


int32_t SamRecord::getInsertSize()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myInsertSize;
}


// Returns the inclusive rightmost position of the clipped sequence.
int32_t SamRecord::get0BasedAlignmentEnd()
{
    myStatus = SamStatus::SUCCESS;
    if(myAlignmentLength == -1)
    {
        // Alignment end has not been set, so calculate it.
        parseCigar();
    }
    // If alignment length > 0, subtract 1 from it to get the end.
    if(myAlignmentLength == 0)
    {
        // Length is 0, just return the start position.
        return(myRecordPtr->myPosition);
    }
    return(myRecordPtr->myPosition + myAlignmentLength - 1);
}


// Returns the inclusive rightmost position of the clipped sequence.
int32_t SamRecord::get1BasedAlignmentEnd()
{
    return(get0BasedAlignmentEnd() + 1);
}

   
// Return the length of the alignment.
int32_t SamRecord::getAlignmentLength()
{
    myStatus = SamStatus::SUCCESS;
    if(myAlignmentLength == -1)
    {
        // Alignment end has not been set, so calculate it.
        parseCigar();
    }
    // Return the alignment length.
    return(myAlignmentLength);
}

// Returns the inclusive left-most position adjust for clipped bases.
int32_t SamRecord::get0BasedUnclippedStart()
{
    myStatus = SamStatus::SUCCESS;
    if(myUnclippedStartOffset == -1)
    {
        // Unclipped has not yet been calculated, so parse the cigar to get it
        parseCigar();
    }
    return(myRecordPtr->myPosition - myUnclippedStartOffset);
}


// Returns the inclusive left-most position adjust for clipped bases.
int32_t SamRecord::get1BasedUnclippedStart()
{
    return(get0BasedUnclippedStart() + 1);
}


// Returns the inclusive right-most position adjust for clipped bases.
int32_t SamRecord::get0BasedUnclippedEnd()
{
    // myUnclippedEndOffset will be set by get0BasedAlignmentEnd if the 
    // cigar has not yet been parsed, so no need to check it here.
    return(get0BasedAlignmentEnd() + myUnclippedEndOffset);
}


// Returns the inclusive right-most position adjust for clipped bases.
int32_t SamRecord::get1BasedUnclippedEnd()
{
    return(get0BasedUnclippedEnd() + 1);
}


// Get the read name.
const char* SamRecord::getReadName()
{
    myStatus = SamStatus::SUCCESS;
    if(myReadName.Length() == 0)
    {
        // 0 Length, means that it is in the buffer, but has not yet
        // been synced to the string, so do the sync.
        myReadName = (char*)&(myRecordPtr->myData);
    }
    return myReadName.c_str();
}


const char* SamRecord::getCigar()
{
    myStatus = SamStatus::SUCCESS;
    if(myCigar.Length() == 0)
    {
        // 0 Length, means that it is in the buffer, but has not yet
        // been synced to the string, so do the sync.
        parseCigarBinary();
    }
    return myCigar.c_str();
}


const char* SamRecord::getSequence()
{
    return(getSequence(mySequenceTranslation));
}


const char* SamRecord::getSequence(SequenceTranslation translation)
{
    myStatus = SamStatus::SUCCESS;
    if(mySequence.Length() == 0)
    {
        // 0 Length, means that it is in the buffer, but has not yet
        // been synced to the string, so do the sync.
        setSequenceAndQualityFromBuffer();
    }

    // Determine if translation needs to be done.
    if((translation == NONE) || (myRefPtr == NULL))
    {
        return mySequence.c_str();
    }
    else if(translation == EQUAL)
    {
        if(mySeqWithEq.length() == 0)
        {
            // Check to see if the sequence is defined.
            if(mySequence == "*")
            {
                // Sequence is undefined, so no translation necessary.
                mySeqWithEq = '*';
            }
            else
            {
                // Sequence defined, so translate it.
                SamQuerySeqWithRef::seqWithEquals(mySequence.c_str(), 
                                                  myRecordPtr->myPosition,
                                                  *(getCigarInfo()),
                                                  getReferenceName(),
                                                  *myRefPtr,
                                                  mySeqWithEq);
            }
        }
        return(mySeqWithEq.c_str());
    }
    else
    {
        // translation == BASES
        if(mySeqWithoutEq.length() == 0)
        {
            if(mySequence == "*")
            {
                // Sequence is undefined, so no translation necessary.
                mySeqWithoutEq = '*';
            }
            else
            {
                // Sequence defined, so translate it.
                SamQuerySeqWithRef::seqWithoutEquals(mySequence.c_str(), 
                                                     myRecordPtr->myPosition,
                                                     *(getCigarInfo()),
                                                     getReferenceName(),
                                                     *myRefPtr,
                                                     mySeqWithoutEq);
            }
        }
        return(mySeqWithoutEq.c_str());
    }
}


const char* SamRecord::getQuality() 
{
    myStatus = SamStatus::SUCCESS;
    if(myQuality.Length() == 0)
    {
        // 0 Length, means that it is in the buffer, but has not yet
        // been synced to the string, so do the sync.
        setSequenceAndQualityFromBuffer();      
    }
    return myQuality.c_str();
}


char SamRecord::getSequence(int index)
{
    return(getSequence(index, mySequenceTranslation));
}


char SamRecord::getSequence(int index, SequenceTranslation translation)
{
    static const char * asciiBases = "=AC.G...T......N";

    // Determine the read length.
    int32_t readLen = getReadLength();

    // If the read length is 0, this method should not be called.
    if(readLen == 0)
    {
        String exceptionString = "SamRecord::getSequence(";
        exceptionString += index;
        exceptionString += ") is not allowed since sequence = '*'";
        throw std::runtime_error(exceptionString.c_str());
    }
    else if((index < 0) || (index >= readLen))
    {
        // Only get here if the index was out of range, so thow an exception.
        String exceptionString = "SamRecord::getSequence(";
        exceptionString += index;
        exceptionString += ") is out of range. Index must be between 0 and ";
        exceptionString += (readLen - 1);
        throw std::runtime_error(exceptionString.c_str());
    }

    // Determine if translation needs to be done.
    if((translation == NONE) || (myRefPtr == NULL))
    {
        // No translation needs to be done.
        if(mySequence.Length() == 0)
        {
            // Parse BAM sequence.
            // TODO - maybe store this pointer - and use that to track when
            // valid?
            unsigned char * packedSequence = 
                (unsigned char *)myRecordPtr->myData + 
                myRecordPtr->myReadNameLength +
                myRecordPtr->myCigarLength * sizeof(int);
            
            return(index & 1 ?
                   asciiBases[packedSequence[index / 2] & 0xF] :
                   asciiBases[packedSequence[index / 2] >> 4]);
        }
        // Already have string.
        return(mySequence[index]);
    }
    else
    {
        // Need to translate the sequence either to have '=' or to not
        // have it.
        // First check to see if the sequence has been set.
        if(mySequence.Length() == 0)
        {
            // 0 Length, means that it is in the buffer, but has not yet
            // been synced to the string, so do the sync.
            setSequenceAndQualityFromBuffer();
        }

        // Check the type of translation.
        if(translation == EQUAL)
        {
            // Check whether or not the string has already been 
            // retrieved that has the '=' in it.
            if(mySeqWithEq.length() == 0)
            {
                // The string with '=' has not yet been determined,
                // so get the string.
                // Check to see if the sequence is defined.
                if(mySequence == "*")
                {
                    // Sequence is undefined, so no translation necessary.
                    mySeqWithEq = '*';
                }
                else
                {
                    // Sequence defined, so translate it.
                    SamQuerySeqWithRef::seqWithEquals(mySequence.c_str(), 
                                                      myRecordPtr->myPosition, 
                                                      *(getCigarInfo()),
                                                      getReferenceName(),
                                                      *myRefPtr,
                                                      mySeqWithEq);
                }
            }
            // Sequence is set, so return it.
            return(mySeqWithEq[index]);
        }
        else
        {
            // translation == BASES
            // Check whether or not the string has already been 
            // retrieved that does not have the '=' in it.
            if(mySeqWithoutEq.length() == 0)
            {
                // The string with '=' has not yet been determined,
                // so get the string.
                // Check to see if the sequence is defined.
                if(mySequence == "*")
                {
                    // Sequence is undefined, so no translation necessary.
                    mySeqWithoutEq = '*';
                }
                else
                {
                    // Sequence defined, so translate it.
                    // The string without '=' has not yet been determined,
                    // so get the string.
                    SamQuerySeqWithRef::seqWithoutEquals(mySequence.c_str(), 
                                                         myRecordPtr->myPosition, 
                                                         *(getCigarInfo()),
                                                         getReferenceName(),
                                                         *myRefPtr,
                                                         mySeqWithoutEq);
                }
            }
            // Sequence is set, so return it.
            return(mySeqWithoutEq[index]);
        }
    }
}


char SamRecord::getQuality(int index)
{
    // Determine the read length.
    int32_t readLen = getReadLength();

    // If the read length is 0, return ' ' whose ascii code is below
    // the minimum ascii code for qualities.
    if(readLen == 0)
    {
        return(BaseUtilities::UNKNOWN_QUALITY_CHAR);
    }
    else if((index < 0) || (index >= readLen))
    {
        // Only get here if the index was out of range, so thow an exception.
        String exceptionString = "SamRecord::getQuality(";
        exceptionString += index;
        exceptionString += ") is out of range. Index must be between 0 and ";
        exceptionString += (readLen - 1);
        throw std::runtime_error(exceptionString.c_str());
    }

    if(myQuality.Length() == 0)
    {
        // Parse BAM Quality.
        unsigned char * packedQuality = 
            (unsigned char *)myRecordPtr->myData +
            myRecordPtr->myReadNameLength + 
            myRecordPtr->myCigarLength * sizeof(int) + 
            (myRecordPtr->myReadLength + 1) / 2;
        return(packedQuality[index] + 33);
    }
    else
    {
        // Already have string.
        if((myQuality.Length() == 1) && (myQuality[0] == '*'))
        {
            // Return 0xFF like it does for BAM.
            return(0xFF);
        }
        else
        {
            return(myQuality[index]);
        }
    }
}

   
Cigar* SamRecord::getCigarInfo()
{
    // Check to see whether or not the Cigar has already been
    // set - this is determined by checking if alignment length
    // is set since alignment length and the cigar are set
    // at the same time.
    if(myAlignmentLength == -1)
    {
        // Not been set, so calculate it.
        parseCigar();
    }
    return(&myCigarRoller);
}


uint32_t SamRecord::getTagLength()
{
    myStatus = SamStatus::SUCCESS;
    if(myNeedToSetTagsFromBuffer)
    {
        // Tags are only set in the buffer, so the size of the tags is 
        // the length of the record minus the starting location of the tags.
        unsigned char * tagStart = 
            (unsigned char *)myRecordPtr->myData 
            + myRecordPtr->myReadNameLength 
            + myRecordPtr->myCigarLength * sizeof(int)
            + (myRecordPtr->myReadLength + 1) / 2 + myRecordPtr->myReadLength;
      
        // The non-tags take up from the start of the record to the tag start.
        // Do not include the block size part of the record since it is not
        // included in the size.
        uint32_t nonTagSize = 
            tagStart - (unsigned char*)&(myRecordPtr->myReferenceID);
        // Tags take up the size of the block minus the non-tag section.
        uint32_t tagSize = myRecordPtr->myBlockSize - nonTagSize;
        return(tagSize);
    }

    // Tags are stored outside the buffer, so myTagBufferSize is set.
    return(myTagBufferSize);
}


// Returns true if there is another tag and sets tag and vtype to the
// appropriate values, and returns a pointer to the value.
// Sets the Status to SUCCESS when a tag is successfully returned or
// when there are no more tags.  Otherwise the status is set to describe
// why it failed (parsing, etc).
bool SamRecord::getNextSamTag(char* tag, char& vtype, void** value)
{
    myStatus = SamStatus::SUCCESS;
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            return(false);
        }
    }

    // Increment the tag index to start looking at the next tag.
    // At the beginning, it is set to -1.
    myLastTagIndex++;
    int maxTagIndex = extras.Capacity();
    if(myLastTagIndex >= maxTagIndex)
    {
        // Hit the end of the tags, return false, no more tags.
        // Status is still success since this is not an error, 
        // it is just the end of the list.
        return(false);
    }

    bool tagFound = false;
    // Loop until a tag is found or the end of extras is hit.
    while((tagFound == false) && (myLastTagIndex < maxTagIndex))
    {
        if(extras.SlotInUse(myLastTagIndex))
        {
            // Found a slot to use.
            int key = extras.GetKey(myLastTagIndex);
            getTag(key, tag);
            getVtype(key, vtype);
            tagFound = true;
            // Get the value associated with the key based on the vtype.
            switch (vtype)
            {
                case 'A' :
                    *value = getIntegerPtr(myLastTagIndex);
                    break;
                case 'f' :
                    *value = getDoublePtr(myLastTagIndex);
                    break;
                case 'c' :
                case 'C' :
                case 's' :
                case 'S' :
                case 'i' :
                case 'I' :
                    vtype = 'i';
                    *value = getIntegerPtr(myLastTagIndex);
                    break;
                case 'Z' :
                    *value = getStringPtr(myLastTagIndex);
                    break;
                default:
                    myStatus.setStatus(SamStatus::FAIL_PARSE,
                                       "Unknown tag type");
                    tagFound = false;
                    break;
            }
        }
        if(!tagFound)
        {
            // Increment the index since a tag was not found.
            myLastTagIndex++;
        }
    }
    return(tagFound);
}


// Returns the values of all fields except the tags.
bool SamRecord::getFields(bamRecordStruct& recStruct, String& readName, 
                          String& cigar, String& sequence, String& quality)
{
    return(getFields(recStruct, readName, cigar, sequence, quality,
                     mySequenceTranslation));
}


// Returns the values of all fields except the tags.
bool SamRecord::getFields(bamRecordStruct& recStruct, String& readName, 
                          String& cigar, String& sequence, String& quality,
                          SequenceTranslation translation)
{
    myStatus = SamStatus::SUCCESS;
    if(myIsBufferSynced == false)
    {
        if(!fixBuffer(translation))
        {
            // failed to set the buffer, return false.
            return(false);
        }
    }
    memcpy(&recStruct, myRecordPtr, sizeof(bamRecordStruct));

    readName = getReadName();
    // Check the status.
    if(myStatus != SamStatus::SUCCESS)
    {
        // Failed to set the fields, return false.
        return(false);
    }
    cigar = getCigar();
    // Check the status.
    if(myStatus != SamStatus::SUCCESS)
    {
        // Failed to set the fields, return false.
        return(false);
    }
    sequence = getSequence(translation);
    // Check the status.
    if(myStatus != SamStatus::SUCCESS)
    {
        // Failed to set the fields, return false.
        return(false);
    }
    quality = getQuality();
    // Check the status.
    if(myStatus != SamStatus::SUCCESS)
    {
        // Failed to set the fields, return false.
        return(false);
    }
    return(true);
}


bool SamRecord::isIntegerType(char vtype) const
{
    if((vtype == 'c') || (vtype == 'C') ||
       (vtype == 's') || (vtype == 'S') ||
       (vtype == 'i') || (vtype == 'I'))
    {
        return(true);
    }
    return(false);
}


bool SamRecord::isDoubleType(char vtype) const
{
    if(vtype == 'f')
    {
        return(true);
    }
    return(false);
}


bool SamRecord::isCharType(char vtype) const
{
    if(vtype == 'A')
    {
        return(true);
    }
    return(false);
}


bool SamRecord::isStringType(char vtype) const
{
    if(vtype == 'Z')
    {
        return(true);
    }
    return(false);
}


void SamRecord::clearTags()
{
    if(extras.Entries() != 0)
    {
        extras.Clear();
    }
    strings.Clear();
    integers.Clear();
    doubles.Clear();
    myTagBufferSize = 0;
    resetTagIter();
}


bool SamRecord::rmTag(const char* tag, char type)
{
    // Check the length of tag.
    if(strlen(tag) != 2)
    {
        // Tag is the wrong length.
        myStatus.setStatus(SamStatus::INVALID, 
                           "rmTag called with tag that is not 2 characters\n");
        return(false);
    }

    myStatus = SamStatus::SUCCESS;
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            return(false);
        }
    }

    // Construct the key.
    int key = MAKEKEY(tag[0], tag[1], type);
    // Look to see if the key exsists in the hash.
    int offset = extras.Find(key);

    if(offset < 0)
    {
        // Not found, so return true, successfully removed since
        // it is not in tag.
        return(true);
    }

    int rmBuffSize = 0;
    
    // Offset is set, so recalculate the buffer size without this entry.
    // Do NOT remove from strings, integers, or doubles because then
    // extras would need to be updated for all entries with the new indexes
    // into those variables.
    switch(type)
    {
        case 'A':
        case 'c':
        case 'C':
            rmBuffSize = 4;
            break;
        case 's':
        case 'S':
            rmBuffSize = 5;
            break;
        case 'i':
        case 'I':
            rmBuffSize = 7;
            break;
        case 'f':
            rmBuffSize = 7;
            break;
        case 'Z':
            rmBuffSize = 4 + getString(offset).Length();
            break;
        default:
            myStatus.setStatus(SamStatus::INVALID, 
                               "rmTag called with unknown type.\n");
            return(false);
            break;
    };

    // The buffer tags are now out of sync.
    myNeedToSetTagsInBuffer = true;
    myIsTagsBufferValid = false;
    myIsBufferSynced = false;
    myTagBufferSize -= rmBuffSize;

    // Remove from the hash.
    extras.Delete(offset);
    return(true);
}


// Return the error after a failed SamRecord call.
const SamStatus& SamRecord::getStatus()
{
    return(myStatus);
}


String & SamRecord::getString(const char * tag)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            // TODO - what do we want to do on failure?
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], 'Z');
    int offset = extras.Find(key);

    int value;
    if (offset < 0)
    {
        // TODO - what do we want to do on failure?
        return(NOT_FOUND_TAG_STRING);
    }
    else
        value = extras[offset];

    return strings[value];
}

int & SamRecord::getInteger(const char * tag)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            // TODO - what do we want to do on failure?
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], 'i');
    int offset = extras.Find(key);

    int value;
    if (offset < 0)
    {
        // TODO - what do we want to do on failure?
        return NOT_FOUND_TAG_INT;
    }
    else
        value = extras[offset];

    return integers[value];
}

double & SamRecord::getDouble(const char * tag)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            // TODO - what do we want to do on failure?
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], 'f');
    int offset = extras.Find(key);

    int value;
    if (offset < 0)
    {
        // TODO - what do we want to do on failure?
        return NOT_FOUND_TAG_DOUBLE;
    }
    else
        value = extras[offset];

    return doubles[value];
}


bool SamRecord::checkTag(const char * tag, char type)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            return("");
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], type);

    return (extras.Find(key) != LH_NOTFOUND);
}


// Return the number of bases in this read that overlap the passed in
// region.  (start & end are 0-based)
uint32_t SamRecord::getNumOverlaps(int32_t start, int32_t end)
{
    // Determine whether or not the cigar has been parsed, which sets up
    // the cigar roller.  This is determined by checking the alignment length.
    if(myAlignmentLength == -1)
    {
        parseCigar();
    }
    return(myCigarRoller.getNumOverlaps(start, end, get0BasedPosition()));
}


// // Get the Extra field: <tag>:<vtype>:<value> for the key.
// void SamRecord::getSamExtraFieldFromKey(int key, String& extraField) 
// {
//     myStatus = SamStatus::SUCCESS;
//     // Extract the tag from the key.
//     char tag[3];
//     getTag(key, tag);

//     char vtype;
//     getVtype(key, vtype);

//     // Look up the key to get the offset into the appropriate array.
//     int offset = extras.Find(key);

//     // Get the value associated with the key based on the vtype.
//     switch (vtype)
//     {
//         case 'A' :
//             extraField = tag;
//             extraField += ":A:";
//             extraField += (char)(getInteger(offset));
//             break;
//         case 'f' :
//             extraField = tag;
//             extraField += ":f:";
//             extraField += getDouble(offset);
//             break;
//         case 'c' :
//         case 'C' :
//         case 's' :
//         case 'S' :
//         case 'i' :
//         case 'I' :
//             extraField = tag;
//             extraField += ":i:";
//             extraField += getInteger(offset);
//             break;
//         case 'Z' :
//             //         String valueString = getString(offset).c_str();
//             extraField = tag;
//             extraField += ":Z:";
//             extraField += getString(offset);
//             break;
//         default:
//             myStatus.setStatus(SamStatus::FAIL_PARSE,
//                                "Unknown tag type.");
//             break;
//     }
// }


// Allocate space for the record - does a realloc.  
// The passed in size is the size of the entire record including the
// block size field.
bool SamRecord::allocateRecordStructure(int size)
{
    if (allocatedSize < size)
    {
        bamRecordStruct* tmpRecordPtr = 
            (bamRecordStruct *)realloc(myRecordPtr, size);
        if(tmpRecordPtr == NULL)
        {
            // FAILED to allocate memory
            fprintf(stderr, "FAILED TO ALLOCATE MEMORY!!!");
            myStatus.addError(SamStatus::FAIL_MEM, "Failed Memory Allocation.");
            return(false);
        }
        // Successfully allocated memory, so set myRecordPtr.
        myRecordPtr = tmpRecordPtr;
        allocatedSize = size;
    }
    return(true);
}


// Index is the index into the strings array.
void* SamRecord::getStringPtr(int index)
{
    int value = extras[index];

    return &(strings[value]);
}

void* SamRecord::getIntegerPtr(int offset)
{
    int value = extras[offset];

    return &(integers[value]);
}

void* SamRecord::getDoublePtr(int offset)
{
    int value = extras[offset];

    return &(doubles[value]);
}


// Fixes the buffer to match the variable length fields if they are set.
bool SamRecord::fixBuffer(SequenceTranslation translation)
{
    // Check to see if the buffer is already synced.
    if(myIsBufferSynced &&
       (myBufferSequenceTranslation == translation))
    {
        // Already synced, nothing to do.
        return(true);
    }
   
    // Set the bin if necessary.
    if(!myIsBinValid)
    {
        // The bin that is set in the record is not valid, so
        // reset it.
        myRecordPtr->myBin = 
            bam_reg2bin(myRecordPtr->myPosition, get1BasedAlignmentEnd());      
        myIsBinValid = true;
    }

    // Not synced.
    bool status = true;

    // First determine the size the buffer needs to be.
    uint8_t newReadNameLen = getReadNameLength();
    uint16_t newCigarLen = getCigarLength();
    int32_t newReadLen = getReadLength();
    uint32_t newTagLen = getTagLength();
    uint32_t bamSequenceLen = (newReadLen+1)/2;

    // The buffer size extends from the start of the record to data
    // plus the length of the variable fields,
    // Multiply the cigar length by 4 since it is the number of uint32_t fields.
    int newBufferSize = 
        ((unsigned char*)(&(myRecordPtr->myData)) - 
         (unsigned char*)myRecordPtr) +
        newReadNameLen + ((newCigarLen)*4) +
        newReadLen + bamSequenceLen + newTagLen;
   
    if(!allocateRecordStructure(newBufferSize))
    {
        // Failed to allocate space.
        return(false);
    }

    // Now that space has been added to the buffer, check to see what if
    // any fields need to be extracted from the buffer prior to starting to
    // overwrite it.  Fields need to be extracted from the buffer if the 
    // buffer is valid for the field and a previous variable length field has
    // changed length.
    bool readNameLenChange = (newReadNameLen != myRecordPtr->myReadNameLength);
    bool cigarLenChange = (newCigarLen != myRecordPtr->myCigarLength);
    bool readLenChange = (newReadLen != myRecordPtr->myReadLength);

    // If the tags are still stored in the buffer and any other fields changed
    // lengths, they need to be extracted.
    if(myIsTagsBufferValid &&
       (readNameLenChange | cigarLenChange | readLenChange))
    {
        status &= setTagsFromBuffer();
        // The tag buffer will not be valid once the other fields
        // are written, so set it to not valid.
        myIsTagsBufferValid = false;
    }

    // If the sequence or quality strings are still stored in the buffer, and
    // any of the previous fields have changed length, extract it from the
    // current buffer.
    if((myIsQualityBufferValid | myIsSequenceBufferValid) && 
       (readNameLenChange | cigarLenChange | readLenChange))
    {
        setSequenceAndQualityFromBuffer();
        // The quality and sequence buffers will not be valid once the other
        // fields are written, so set them to not valid.
        myIsQualityBufferValid = false;
        myIsSequenceBufferValid = false;
    }

    // If the cigar is still stored in the buffer, and any of the
    // previous fields have changed length, extract it from the current buffer.
    if((myIsCigarBufferValid) && 
       (readNameLenChange))
    {
        status &= parseCigarBinary();
        myIsCigarBufferValid = false;
    }

    // Set each value in the buffer if it is not already valid.
    if(!myIsReadNameBufferValid)
    {
        memcpy(&(myRecordPtr->myData), myReadName.c_str(), 
               newReadNameLen);
   
        // Set the new ReadNameLength.
        myRecordPtr->myReadNameLength = newReadNameLen;
        myIsReadNameBufferValid = true;
    }

    unsigned char * readNameEnds = (unsigned char*)(&(myRecordPtr->myData)) +
        myRecordPtr->myReadNameLength;
   
    // Set the Cigar.  Need to reformat from the string to 
    unsigned int * packedCigar = (unsigned int *) (void *) readNameEnds;
      
    if(!myIsCigarBufferValid)
    {
        // The cigar was already parsed when it was set, so just copy
        // data from the temporary buffer.
        myRecordPtr->myCigarLength = newCigarLen;
        memcpy(packedCigar, myCigarTempBuffer, 
               myRecordPtr->myCigarLength * sizeof(uint32_t));
      
        myIsCigarBufferValid = true;
    }

    unsigned char * packedSequence = readNameEnds + 
        myRecordPtr->myCigarLength * sizeof(int);
    unsigned char * packedQuality = packedSequence + bamSequenceLen;
   
    if(!myIsSequenceBufferValid || !myIsQualityBufferValid || 
       (myBufferSequenceTranslation != translation))
    {
        myRecordPtr->myReadLength = newReadLen;
        // Determine if the quality needs to be set and is just a * and needs to
        // be set to 0xFF.
        bool noQuality = false;
        if((myQuality.Length() == 1) && (myQuality[0] == '*'))
        {
            noQuality = true;
        }
      
        const char* translatedSeq = NULL;
        // If the sequence is not valid in the buffer or it is not
        // properly translated, get the properly translated sequence
        // that needs to be put into the buffer.
        if((!myIsSequenceBufferValid) ||
           (translation != myBufferSequenceTranslation))
        {
            translatedSeq = getSequence(translation);
        }

        for (int i = 0; i < myRecordPtr->myReadLength; i++) 
        {
            if((!myIsSequenceBufferValid) ||
               (translation != myBufferSequenceTranslation))
            {
                // Sequence buffer is not valid, so set the sequence.
                int seqVal = 0;
                switch(translatedSeq[i])
                {
                    case '=':
                        seqVal = 0;
                        break;
                    case 'A':
                    case 'a':
                        seqVal = 1;
                        break;
                    case 'C':
                    case 'c':
                        seqVal = 2;
                        break;
                    case 'G':
                    case 'g':
                        seqVal = 4;
                        break;
                    case 'T':
                    case 't':
                        seqVal = 8;
                        break;
                    case 'N':
                    case 'n':
                    case '.':
                        seqVal = 15;
                        break;
                    default:
                        myStatus.addError(SamStatus::FAIL_PARSE,
                                          "Unknown Sequence character found.");
                        status = false;
                        break;
                };
            
                if(i & 1)
                {
                    // Odd number i's go in the lower 4 bits, so OR in the
                    // lower bits
                    packedSequence[i/2] |= seqVal;
                }
                else
                {
                    // Even i's go in the upper 4 bits and are always set first.
                    packedSequence[i/2] = seqVal << 4;
                }
            }

            if(!myIsQualityBufferValid)
            {
                // Set the quality.
                if((noQuality) || (myQuality.Length() <= i))
                {
                    // No quality or the quality is smaller than the sequence,
                    // so set it to 0xFF
                    packedQuality[i] = 0xFF;
                }
                else
                {
                    // Copy the quality string.
                    packedQuality[i] = myQuality[i] - 33;
                }
            }
        }
        myIsQualityBufferValid = true;
        myIsSequenceBufferValid = true;
        myBufferSequenceTranslation = translation;
    }

    if(!myIsTagsBufferValid)
    {
        status &= setTagsInBuffer();
    }

    // Set the lengths in the buffer.
    myRecordPtr->myReadNameLength = newReadNameLen;
    myRecordPtr->myCigarLength = newCigarLen;
    myRecordPtr->myReadLength = newReadLen;

    // Set the buffer block size to the size of the buffer minus the
    // first field.
    myRecordPtr->myBlockSize = newBufferSize - sizeof(int32_t);

    if(status)
    {
        myIsBufferSynced = true;
    }

    return(status);
}


// Sets the Sequence and Quality strings from the buffer.
// They are done together in one method because they require the same
// loop, so might as well be done at the same time.
void SamRecord::setSequenceAndQualityFromBuffer()
{
    // NOTE: If the sequence buffer is not valid, do not set the sequence
    // string from the buffer.
    // NOTE: If the quality buffer is not valid, do not set the quality string
    // from the buffer.

    // Extract the sequence if the buffer is valid and the string's length is 0.
    bool extractSequence = false;
    if(myIsSequenceBufferValid && (mySequence.Length() == 0))
    {
        extractSequence = true;
    }

    // Extract the quality if the buffer is valid and the string's length is 0.
    bool extractQuality = false;
    if(myIsQualityBufferValid && (myQuality.Length() == 0))
    {
        extractQuality = true;
    }

    // If neither the quality nor the sequence need to be extracted,
    // just return.
    if(!extractSequence && !extractQuality)
    {
        return;
    }

    // Set the sequence and quality strings..
    if(extractSequence)
    {
        mySequence.SetLength(myRecordPtr->myReadLength);
    }
    if(extractQuality)
    {
        myQuality.SetLength(myRecordPtr->myReadLength);
    }
   
    unsigned char * readNameEnds = 
        (unsigned char *)myRecordPtr->myData + myRecordPtr->myReadNameLength;
    unsigned char * packedSequence = 
        readNameEnds + myRecordPtr->myCigarLength * sizeof(int);
    unsigned char * packedQuality = 
        packedSequence + (myRecordPtr->myReadLength + 1) / 2;

    const char * asciiBases = "=AC.G...T......N";

    // Flag to see if the quality is specified - the quality contains a value
    // other than 0xFF.  If all values are 0xFF, then there is no quality.
    bool qualitySpecified = false;

    for (int i = 0; i < myRecordPtr->myReadLength; i++)
    {
        if(extractSequence)
        {
            mySequence[i] = i & 1 ?
                asciiBases[packedSequence[i / 2] & 0xF] :
                asciiBases[packedSequence[i / 2] >> 4];
        }

        if(extractQuality)
        {
            if(packedQuality[i] != 0xFF)
            {
                // Quality is specified, so mark the flag.
                qualitySpecified = true;
            }

            myQuality[i] = packedQuality[i] + 33;
        }
    }

    // If the read length is 0, then set the sequence and quality to '*'
    if(myRecordPtr->myReadLength == 0)
    {
        if(extractSequence)
        {
            mySequence = "*";
        }
        if(extractQuality)
        {
            myQuality = "*";
        }
    }
    else if(extractQuality && !qualitySpecified)
    {
        // No quality was specified, so set it to "*"
        myQuality = "*";
    }
}


// Parse the cigar to calculate the alignment/unclipped end.
bool SamRecord::parseCigar()
{
    // Determine if the cigar string or cigar binary needs to be parsed.
    if(myCigar.Length() == 0)
    {
        // The cigar string is not yet set, so parse the binary.
        return(parseCigarBinary());
    }
    return(parseCigarString());
}

// Parse the cigar to calculate the alignment/unclipped end.
bool SamRecord::parseCigarBinary()
{
    // Only need to parse if the string is not already set.
    // The length of the cigar string is set to zero when the 
    // record is read from a file into the buffer.
    if(myCigar.Length() != 0)
    {
        // Already parsed.
        return(true);
    }

    unsigned char * readNameEnds = 
        (unsigned char *)myRecordPtr->myData + myRecordPtr->myReadNameLength;
   
    unsigned int * packedCigar = (unsigned int *) (void *) readNameEnds;

    myCigarRoller.Set(packedCigar, myRecordPtr->myCigarLength);
    
    myCigarRoller.getCigarString(myCigar);

    myAlignmentLength = myCigarRoller.getExpectedReferenceBaseCount();

    myUnclippedStartOffset = myCigarRoller.getNumBeginClips();
    myUnclippedEndOffset = myCigarRoller.getNumEndClips();

    // if the cigar length is 0, then set the cigar string to "*"
    if(myRecordPtr->myCigarLength == 0)
    {
        myCigar = "*";
        return(true);
    }

    // Copy the cigar into a temporary buffer.
    int newBufferSize = myRecordPtr->myCigarLength * sizeof(uint32_t);
    if(newBufferSize > myCigarTempBufferAllocatedSize)
    {
        uint32_t* tempBufferPtr = 
            (uint32_t*)realloc(myCigarTempBuffer, newBufferSize);
        if(tempBufferPtr == NULL)
        {
            // Failed to allocate memory.
            // Do not parse, just return.
            fprintf(stderr, "FAILED TO ALLOCATE MEMORY!!!");
            myStatus.addError(SamStatus::FAIL_MEM,
                              "Failed to Allocate Memory.");
            return(false);
        }
        myCigarTempBuffer = tempBufferPtr;
        myCigarTempBufferAllocatedSize = newBufferSize;
    }

    memcpy(myCigarTempBuffer, packedCigar, 
           myRecordPtr->myCigarLength * sizeof(uint32_t));

    // Set the length of the temp buffer.
    myCigarTempBufferLength = myRecordPtr->myCigarLength;

    return(true);
}

// Parse the cigar string to calculate the cigar length and alignment end.
bool SamRecord::parseCigarString()
{
    myCigarTempBufferLength = 0;
    if(myCigar == "*")
    {
        // Cigar is empty, so initialize the variables.
        myAlignmentLength = 0;
        myUnclippedStartOffset = 0;
        myUnclippedEndOffset = 0;
        myCigarRoller.clear();
        return(true);
    }

    myCigarRoller.Set(myCigar);
    
    myAlignmentLength = myCigarRoller.getExpectedReferenceBaseCount();
    
    myUnclippedStartOffset = myCigarRoller.getNumBeginClips();
    myUnclippedEndOffset = myCigarRoller.getNumEndClips();

    // Check to see if the Temporary Cigar Buffer is large enough to contain
    // this cigar.  If we make it the size of the length of the cigar string,
    // it will be more than large enough.
    int newBufferSize = myCigar.Length() * sizeof(uint32_t);
    if(newBufferSize > myCigarTempBufferAllocatedSize)
    {
        uint32_t* tempBufferPtr = 
            (uint32_t*)realloc(myCigarTempBuffer, newBufferSize);
        if(tempBufferPtr == NULL)
        {
            // Failed to allocate memory.
            // Do not parse, just return.
            fprintf(stderr, "FAILED TO ALLOCATE MEMORY!!!");
            myStatus.addError(SamStatus::FAIL_MEM,
                              "Failed to Allocate Memory.");
            return(false);
        }
        myCigarTempBuffer = tempBufferPtr;
        myCigarTempBufferAllocatedSize = newBufferSize;
    }

    // Track if there were any errors.
    bool status = true;

    // Track the index into the cigar string that is being parsed.
    char *cigarOp;
    const char* cigarEntryStart = myCigar.c_str();
    int opLen = 0;
    int op = 0;

    unsigned int * packedCigar = myCigarTempBuffer;
    // TODO - maybe one day make a cigar list... or maybe make a 
    // reference cigar string for ease of lookup....
    const char* endCigarString = cigarEntryStart + myCigar.Length();
    while(cigarEntryStart < endCigarString)
    {
        bool validCigarEntry = true;
        // Get the opLen from the string.  cigarOp will then point to 
        // the operation.
        opLen = strtol(cigarEntryStart, &cigarOp, 10);
        // Switch on the type of operation.
        switch(*cigarOp)
        {
            case('M'):
                op = 0;
                break;
            case('I'):
                // Insert into the reference position, so do not increment the
                // reference end position.
                op = 1;
                break;
            case('D'):
                op = 2;
                break;
            case('N'):
                op = 3;
                break;
            case('S'):
                op = 4;
                break;
            case('H'):
                op = 5;
                break;
            case('P'):
                op = 6;
                break;
            default:
                fprintf(stderr, "ERROR parsing cigar\n");
                validCigarEntry = false;
                status = false;
                myStatus.addError(SamStatus::FAIL_PARSE,
                                  "Unknown operation found when parsing the Cigar.");
                break;
        };
        if(validCigarEntry)
        {
            // Increment the cigar length.
            ++myCigarTempBufferLength;
            *packedCigar = (opLen << 4) | op;
            packedCigar++;
        }
        // The next Entry starts at one past the cigar op, so set the start.
        cigarEntryStart = ++cigarOp;
    }

    // Check clipLength to adjust the end position.
    return(status);
}


bool SamRecord::setTagsFromBuffer()
{
    // If the tags do not need to be set from the buffer, return true.
    if(myNeedToSetTagsFromBuffer == false)
    {
        // Already been set from the buffer.
        return(true);
    }

    // Mark false, as they are being set now.
    myNeedToSetTagsFromBuffer = false;

    unsigned char * readNameEnds = 
        (unsigned char *)myRecordPtr->myData + myRecordPtr->myReadNameLength;
    unsigned char * packedSequence = 
        readNameEnds + myRecordPtr->myCigarLength * sizeof(int);
    unsigned char * packedQuality = 
        packedSequence + (myRecordPtr->myReadLength + 1) / 2;

    unsigned char * extraPtr = packedQuality + myRecordPtr->myReadLength;

    // Default to success, will be changed to false on failure.
    bool status = true;

    // Clear any previously set tags.
    clearTags();
    while (myRecordPtr->myBlockSize + 4 - 
           (extraPtr - (unsigned char *)myRecordPtr) > 0)
    {
        int key;
        int value;
        void * content = extraPtr + 3;
        int tagBufferSize = 0;
      
        switch (extraPtr[2])
        {
            case 'A' :
                key = MAKEKEY(extraPtr[0], extraPtr[1], 'A');
                value = integers.Length();
                integers.Push(* (char *) content);
                extraPtr += 4;
                tagBufferSize += 4;
                break;
            case 'c' :
                key = MAKEKEY(extraPtr[0], extraPtr[1], 'c');
                value = integers.Length();
                integers.Push(* (char *) content);
                extraPtr += 4;
                tagBufferSize += 4;
                break;
            case 'C' :
                key = MAKEKEY(extraPtr[0], extraPtr[1], 'C');
                value = integers.Length();
                integers.Push(* (unsigned char *) content);
                extraPtr += 4;
                tagBufferSize += 4;
                break;
            case 's' :
                key = MAKEKEY(extraPtr[0], extraPtr[1], 's');
                value = integers.Length();
                integers.Push(* (short *) content);
                extraPtr += 5;
                tagBufferSize += 5;
                break;
            case 'S' :
                key = MAKEKEY(extraPtr[0], extraPtr[1], 'S');
                value = integers.Length();
                integers.Push(* (unsigned short *) content);
                extraPtr += 5;
                tagBufferSize += 5;
                break;
            case 'i' :
                key = MAKEKEY(extraPtr[0], extraPtr[1], 'i');
                value = integers.Length();
                integers.Push(* (int *) content);
                extraPtr += 7;
                tagBufferSize += 7;
                break;
            case 'I' :
                key = MAKEKEY(extraPtr[0], extraPtr[1], 'I');
                value = integers.Length();
                integers.Push((int) * (unsigned int *) content);
                extraPtr += 7;
                tagBufferSize += 7;
                break;
            case 'Z' :
                key = MAKEKEY(extraPtr[0], extraPtr[1], 'Z');
                value = strings.Length();
                strings.Push((const char *) content);
                extraPtr += 4 + strings.Last().Length();
                tagBufferSize += 4 + strings.Last().Length();
                break;
            case 'f' :
                key = MAKEKEY(extraPtr[0], extraPtr[1], 'f');
                value = doubles.Length();
                doubles.Push(* (float *) content);
                fprintf(stderr, "\n\nFLOAT!!!\n\n");
                extraPtr += 7;
                tagBufferSize += 7;
                break;
            default :
                fprintf(stderr, 
                        "parsing BAM - Unknown custom field of type %c%c:%c\n",
                        extraPtr[0], extraPtr[1], extraPtr[2]);
                // Failed on read.
                // Increment extraPtr just by the size of the 3 known fields
                extraPtr += 3;
                myStatus.addError(SamStatus::FAIL_PARSE,
                                  "Unknown tag type.");
                status = false;
        }

        // Only add the tag if it has so far been successfully processed.
        if(status)
        {
            extras.Add(key, value);
            myTagBufferSize += tagBufferSize;
        }
    }
    return(status);
}


bool SamRecord::setTagsInBuffer()
{
    // The buffer size extends from the start of the record to data
    // plus the length of the variable fields,
    // Multiply the cigar length by 4 since it is the number of uint32_t fields.
    int bamSequenceLength = (myRecordPtr->myReadLength+1)/2;
    int newBufferSize = ((unsigned char*)(&(myRecordPtr->myData)) - 
                         (unsigned char*)myRecordPtr) +  
        myRecordPtr->myReadNameLength + ((myRecordPtr->myCigarLength)*4) +
        myRecordPtr->myReadLength + bamSequenceLength + myTagBufferSize;

    // Make sure the buffer is big enough.
    if(!allocateRecordStructure(newBufferSize))
    {
        // Failed to allocate space.
        return(false);
    }

    unsigned char * readNameEnds = (unsigned char*)(&(myRecordPtr->myData)) +
        myRecordPtr->myReadNameLength;
    unsigned char * packedSequence = readNameEnds + 
        myRecordPtr->myCigarLength * sizeof(int);
    unsigned char * packedQuality = 
        packedSequence + bamSequenceLength;
   
    char * extraPtr = (char*)packedQuality + myRecordPtr->myReadLength;

    bool status = true;

    // Set the tags in the buffer.
    if (extras.Entries())
    {
        for (int i = 0; i < extras.Capacity(); i++)
        {
            if (extras.SlotInUse(i))
            {
                int key = extras.GetKey(i);
                getTag(key, extraPtr);
                extraPtr += 2;
                getVtype(key, extraPtr[0]);
                char vtype = extraPtr[0];

                // increment the pointer to where the value is.
                extraPtr += 1;

                switch (vtype)
                {
                    case 'A' :
                        *(char*)extraPtr = (char)getInteger(i);
                        // sprintf(extraPtr, "%d", getInteger(i));
                        extraPtr += 1;
                        break;
                    case 'c' :
                        *(int8_t*)extraPtr = (int8_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 1;
                        break;
                    case 'C' :
                        *(uint8_t*)extraPtr = (uint8_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 1;
                        break;
                    case 's' :
                        *(int16_t*)extraPtr = (int16_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 2;
                        break;
                    case 'S' :
                        *(uint16_t*)extraPtr = (uint16_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 2;
                        break;
                    case 'i' :
                        *(int32_t*)extraPtr = (int32_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 4;
                        break;
                    case 'I' :
                        *(uint32_t*)extraPtr = (uint32_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 4;
                        break;
                    case 'Z' :
                        sprintf(extraPtr, "%s", getString(i).c_str());
                        extraPtr += getString(i).Length() + 1;
                        break;
                    case 'f' :
                        // TODO, figure out if %f is good enough.
                        sprintf(extraPtr, "%f", getDouble(i));
                        extraPtr += 4;
                        break;
                    default :
                        myStatus.addError(SamStatus::FAIL_PARSE,
                                          "Unknown tag type.");
                        status = false;
                        break;
                }
            }
        }
    }

    // Validate that the extra pointer is at the end of the allocated buffer.
    // If not then there was a problem.
    if(extraPtr != (char*)myRecordPtr + newBufferSize)
    {
        fprintf(stderr, "ERROR updating the buffer.  Incorrect size.");
        myStatus.addError(SamStatus::FAIL_PARSE,
                          "ERROR updating the buffer.  Incorrect size.");
        status = false;
    }


    // The buffer tags are now in sync.
    myNeedToSetTagsInBuffer = false;
    myIsTagsBufferValid = true;
    return(status);
}


// Reset the variables for a newly set buffer.  The buffer must be set first
// since this looks up the reference ids in the buffer to set the reference
// names.
void SamRecord::setVariablesForNewBuffer(SamFileHeader& header)
{
    // Lookup the reference name & mate reference name associated with this
    // record.
    myReferenceName = 
        header.getReferenceLabel(myRecordPtr->myReferenceID);
    myMateReferenceName = 
        header.getReferenceLabel(myRecordPtr->myMateReferenceID);      

    // Clear the SAM Strings that are now not in-sync with the buffer.
    myReadName.SetLength(0);
    myCigar.SetLength(0);
    mySequence.SetLength(0);
    mySeqWithEq.clear();
    mySeqWithoutEq.clear();
    myQuality.SetLength(0);
    myNeedToSetTagsFromBuffer = true;
    myNeedToSetTagsInBuffer = false;

    //Set that the buffer is valid.
    myIsBufferSynced = true;
    // Set that the variable length buffer fields are valid.
    myIsReadNameBufferValid = true;
    myIsCigarBufferValid = true;
    myIsSequenceBufferValid = true;
    myBufferSequenceTranslation = NONE;
    myIsQualityBufferValid = true;
    myIsTagsBufferValid = true;
    myIsBinValid = true;
}


// Extract the vtype from the key.
void SamRecord::getVtype(int key, char& vtype) const
{
    // Extract the vtype from the key.
    vtype = (key >> 16) & 0xFF;
}

// Extract the tag from the key.
void SamRecord::getTag(int key, char* tag) const
{
    // Extract the tag from the key.
    tag[0] = key & 0xFF;
    tag[1] = (key >> 8) & 0xFF;
    tag[2] = 0;
}


// Index is the index into the strings array.
String & SamRecord::getString(int index)
{
    int value = extras[index];

    return strings[value];
}

int & SamRecord::getInteger(int offset)
{
    int value = extras[offset];

    return integers[value];
}

double & SamRecord::getDouble(int offset)
{
    int value = extras[offset];

    return doubles[value];
}


