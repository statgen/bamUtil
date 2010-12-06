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

#include "SamInterface.h"

#include <limits>
#include <stdint.h>

SamInterface::SamInterface()
{
}


SamInterface::~SamInterface()
{
}


// Read a SAM file's header.
SamStatus::Status SamInterface::readHeader(IFILE filePtr, SamFileHeader& header)
{
    if(filePtr == NULL)
    {
        // File is not open.
        return(SamStatus::FAIL_ORDER);
    }

    // Clear the passed in header.
    header.resetHeader();

    do {
        StringIntHash tags;
        StringArray   values;
        buffer.ReadLine(filePtr);
      
        // Stop reading header lines if at the end of the file or
        // if the line is not blank and does not start with an @.
        if ( ifeof(filePtr) || 
             ((buffer.Length() != 0) && (buffer[0] != '@')) )
        {
            break;
        }
      
        // This is a header line, so add it to header.
        header.addHeaderLine(buffer.c_str());

        // Continue to the next line if this line is less than 3 characters
        // or is not an SQ line.
        if ((buffer.Length() < 3) || (buffer[1] != 'S') || (buffer[2] != 'Q'))
            continue;
      
        ParseHeaderLine(tags, values);
      
        int name = tags.Integer("SN");
        int length = tags.Integer("LN");
      
        if (name < 0 || length < 0) continue;

        header.addReferenceInfo(values[name], 
                                values[length].AsInteger());
      
    } while (1);
   
    // Store the first record since it was read.
    myFirstRecord = buffer;

    // Successfully read.
    return(SamStatus::SUCCESS);
}

SamStatus::Status SamInterface::writeHeader(IFILE filePtr,
                                            SamFileHeader& header)
{
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return failure.
        return(SamStatus::FAIL_ORDER);
    }

    ////////////////////////////////
    // Write the header to the file.
    ////////////////////////////////
    // Construct a string containing the entire header.
    std::string headerString = "";
    header.getHeaderString(headerString);
    
    int32_t headerLen = headerString.length();
    int numWrite = 0;
    
    // Write the header to the file.
    numWrite = ifwrite(filePtr, headerString.c_str(), headerLen);
    if(numWrite != headerLen)
    {
        return(SamStatus::FAIL_IO);
    }
    return(SamStatus::SUCCESS);
}


void SamInterface::readRecord(IFILE filePtr, SamFileHeader& header,
                              SamRecord& record, 
                              SamStatus& samStatus)
{
    // Initialize the status to success - will be set to false on failure.
    samStatus = SamStatus::SUCCESS;

    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open.
        samStatus.addError(SamStatus::FAIL_ORDER, 
                           "filePtr does not point to an open file.");
        return;
    }
    
    // If the first record has been set, use that and clear it,
    // otherwise read the record from the file.
    if(myFirstRecord.Length() != 0)
    {
        buffer = myFirstRecord;
        myFirstRecord.Clear();
    }
    else
    {
        // Read the next record.
        buffer.Clear();
        buffer.ReadLine(filePtr);
        // If the end of the file and nothing was read, return false.
        if ((ifeof(filePtr)) && (buffer.Length() == 0))
        {
            // end of the file and nothing to process.
            samStatus.addError(SamStatus::NO_MORE_RECS, 
                               "No more records in the file.");
            return;
        }
    }
    
    tokens.ReplaceColumns(buffer, '\t');
    
    
    // Error string for reporting a parsing failure.
    String errorString = "";
    
    if (tokens.Length() < 11)
    {
        errorString = "Too few columns (";
        errorString += tokens.Length();
        errorString += ") in the Record, expected at least 11.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
        return;
    }
        
    // Reset the record before setting any fields.
    record.resetRecord();

    if(!record.setReadName(tokens[0]))
    {
        samStatus.addError(record.getStatus());
    }
    
    long flagInt = 0;
    if(!tokens[1].AsInteger(flagInt))
    {
        errorString = "flag, ";
        errorString += tokens[1].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if((flagInt < 0) || (flagInt > UINT16_MAX))
    {
        errorString = "flag, ";
        errorString += tokens[1].c_str();
        errorString += ", is not between 0 and (2^16)-1 = 65535.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.setFlag(flagInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setReferenceName(header, tokens[2]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    long posInt = 0;
    if(!tokens[3].AsInteger(posInt))
    {
        errorString = "position, ";
        errorString += tokens[3].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if((posInt < INT32_MIN) || (posInt > INT32_MAX))
    {
        // If it is not in this range, it cannot fit into a 32 bit int.
        errorString = "position, ";
        errorString += tokens[3].c_str();
        errorString += ", does not fit in a 32 bit signed int.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.set1BasedPosition(posInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    long mapInt = 0;
    if(!tokens[4].AsInteger(mapInt))
    {
        errorString = "map quality, ";
        errorString += tokens[4].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if((mapInt < 0) || (mapInt > UINT8_MAX))
    {
        errorString = "map quality, ";
        errorString += tokens[4].c_str();
        errorString += ", is not between 0 and (2^8)-1 = 255.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.setMapQuality(mapInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setCigar(tokens[5]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setMateReferenceName(header, tokens[6]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    long matePosInt = 0;
    if(!tokens[7].AsInteger(matePosInt))
    {
        errorString = "mate position, ";
        errorString += tokens[7].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.set1BasedMatePosition(matePosInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    long insertInt = 0;
    if(!tokens[8].AsInteger(insertInt))
    {
        errorString = "insert size, ";
        errorString += tokens[8].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.setInsertSize(insertInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setSequence(tokens[9]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setQuality(tokens[10]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }
    
    // Clear the tag fields.
    record.clearTags();
    
    // Add the tags to the record.
    for (int i = 11; i < tokens.Length(); i++)
    {
        String & nugget = tokens[i];
        
        if (nugget.Length() < 6 || nugget[2] != ':' || nugget[4] != ':')
        {
            // invalid tag format.
            errorString = "Invalid Tag Format: ";
            errorString += nugget.c_str();
            errorString += ", should be cc:c:x*.";
            samStatus.addError(SamStatus::FAIL_PARSE,
                               errorString.c_str());
            continue;
        }
        
        // Valid tag format.
        // Add the tag.
        if(!record.addTag((const char *)nugget, nugget[3],
                          (const char *)nugget + 5))
        {
            samStatus.addError(record.getStatus().getStatus(),
                               record.getStatus().getStatusMessage());
        }
    }

    return;
}


SamStatus::Status SamInterface::writeRecord(IFILE filePtr,
                                            SamFileHeader& header, 
                                            SamRecord& record,
                                            SamRecord::SequenceTranslation translation)
{
    // Store all the fields into a string, then write the string.
    String recordString = record.getReadName();
    recordString += "\t";
    recordString += record.getFlag();
    recordString += "\t";
    recordString += record.getReferenceName();
    recordString += "\t";
    recordString += record.get1BasedPosition();
    recordString += "\t";
    recordString += record.getMapQuality();
    recordString += "\t";
    recordString += record.getCigar();
    recordString += "\t";
    recordString += record.getMateReferenceNameOrEqual();
    recordString += "\t";
    recordString += record.get1BasedMatePosition();
    recordString += "\t";
    recordString += record.getInsertSize();
    recordString += "\t";
    recordString += record.getSequence(translation);
    recordString += "\t";
    recordString += record.getQuality();
   
    char tag[3];
    char vtype;
    void* value;

    // Reset the tag iterator to ensure that all the tags are written.
    record.resetTagIter();

    // While there are more tags, write them to the recordString.
    while(record.getNextSamTag(tag, vtype, &value) != false)
    {
        recordString += "\t";
        recordString += tag;
        recordString += ":"; 
        recordString += vtype;
        recordString += ":";
        if(record.isIntegerType(vtype))
        {
            recordString += (int)*(int*)value;
        }
        else if(record.isDoubleType(vtype))
        {
            recordString += (double)*(double*)value;
        }
        else if(record.isCharType(vtype))
        {
            recordString += (char)*(char*)value;
        }
        else
        {
            // String type.
            recordString += (String)*(String*)value;
        }
    }

    recordString += "\n";
   
   
    // Write the record.
    ifwrite(filePtr, recordString.c_str(), recordString.Length());
    return(SamStatus::SUCCESS);
}


void SamInterface::ParseHeaderLine(StringIntHash & tags, StringArray & values)
{
    tags.Clear();
    values.Clear();

    tokens.AddColumns(buffer, '\t');

    for (int i = 1; i < tokens.Length(); i++)
    {
        tags.Add(tokens[i].Left(2), i - 1);
        values.Push(tokens[i].SubStr(3));
    }
}

