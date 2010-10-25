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

#include "BamInterface.h"
#include "CharBuffer.h"

BamInterface::BamInterface()
{
}


BamInterface::~BamInterface()
{
}


// Read a BAM file's header.
SamStatus::Status BamInterface::readHeader(IFILE filePtr, SamFileHeader& header)
{
    if(filePtr == NULL)
    {
        // File is not open, return false.
        return(SamStatus::FAIL_ORDER);
    }

    // Clear the passed in header.
    header.resetHeader();

    SamStatus::Status status = header.setHeaderFromBamFile(filePtr);
    if(status != SamStatus::SUCCESS)
    {
        return(status);
    }

    int referenceCount;
    // Read the number of references sequences.
    ifread(filePtr, &referenceCount, sizeof(int));

//     header.referenceContigs.Dimension(referenceCount);
//     header.referenceLengths.Dimension(referenceCount);
    CharBuffer refName;

    // Read each reference sequence
    for (int i = 0; i < referenceCount; i++)
    {
        int nameLength;
        // Read the length of the reference name.
        ifread(filePtr, &nameLength, sizeof(int));
      
        // Read the name.
        refName.readFromFile(filePtr, nameLength);

        // Read the length of the reference sequence.
        int32_t refLen;
        ifread(filePtr, &refLen, sizeof(int));

        header.addReferenceInfo(refName.c_str(), refLen);
    }

    // Successfully read the file.
    return(SamStatus::SUCCESS);
}


SamStatus::Status BamInterface::writeHeader(IFILE filePtr, 
                                            SamFileHeader& header)
{
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return failure.
        return(SamStatus::FAIL_ORDER);
    }

    char magic[4];
    magic[0] = 'B';
    magic[1] = 'A';
    magic[2] = 'M';
    magic[3] = 1;

    // Write magic to the file.
    ifwrite(filePtr, magic, 4);

    ////////////////////////////////
    // Write the header to the file.
    ////////////////////////////////
    // Construct a string containing the entire header.
    std::string headerString = "";
    header.getHeaderString(headerString);

    int32_t headerLen = headerString.length();
    int numWrite = 0;
    
    // Write the header length.
    numWrite = ifwrite(filePtr, &headerLen, sizeof(int32_t));
    if(numWrite != sizeof(int32_t))
    {
        return(SamStatus::FAIL_IO);
    }
   
    // Write the header to the file.
    numWrite = ifwrite(filePtr, headerString.c_str(), headerLen);
    if(numWrite != headerLen)
    {
        return(SamStatus::FAIL_IO);
    }
    
    ////////////////////////////////////////////////////////
    // Write the Reference Information.
    const SamReferenceInfo* refInfo = header.getReferenceInfo();
    if(refInfo == NULL)
    {
        // Failed to get the reference info.
        return(SamStatus::INVALID);
    }

    // Get the number of sequences.    
    int32_t numSeq = refInfo->getNumEntries();

    // Check to see if the reference info has not yet been set.
    // If not, check for the SQ header.
    if(numSeq == 0)
    {
        header.generateReferenceInfo();
        numSeq = refInfo->getNumEntries();
    }
    ifwrite(filePtr, &numSeq, sizeof(int32_t));

    // Write each reference sequence
    for (int i = 0; i < numSeq; i++)
    {
        const char* refName = refInfo->getReferenceName(i);
        // Add one for the null value.
        int32_t nameLength = strlen(refName) + 1;
        // Write the length of the reference name.
        ifwrite(filePtr, &nameLength, sizeof(int32_t));
      
        // Write the name.
        ifwrite(filePtr, refName, nameLength);
        // Write the length of the reference sequence.
        int32_t refLen = refInfo->getReferenceLength(i);
        ifwrite(filePtr, &refLen, sizeof(int32_t));
    }

    return(SamStatus::SUCCESS);
}


void BamInterface::readRecord(IFILE filePtr, SamFileHeader& header,
                              SamRecord& record, 
                              SamStatus& samStatus)
{
    // TODO - need to validate there are @SQ lines in both sam/bam - MAYBE!

    // Reset the record prior to reading a new one.
    record.resetRecord();

    if(record.setBufferFromFile(filePtr, header) != SamStatus::SUCCESS)
    {
        // Failed, so add the error message.
        samStatus.addError(record.getStatus());
    }
}

SamStatus::Status BamInterface::writeRecord(IFILE filePtr, 
                                            SamFileHeader& header,
                                            SamRecord& record)
{
    // Write the file, returning the status.
    return(record.writeRecordBuffer(filePtr));
}


