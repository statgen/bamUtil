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

#include <stdio.h>
#include <string.h>
#include "SamFile.h"
#include "SamFlag.h"

// print Usage
void printUsage() {
    fprintf(stderr,"Usage: squeezeBam [inFile] [outFile]\n");
    fprintf(stderr,"squeezeBam will modify the file dropping OQ fields, duplicates and shortening Read Names.\n");
}

// main function
int main(int argc, char ** argv)
{
    SamFile samIn;
    SamFile samOut;

    if ( argc != 3 )
    {
        printUsage();
        fprintf(stderr, "ERROR: Failed to properly specify the arguments\n");
        return(-1);
    }

    if ( ! samIn.OpenForRead(argv[1]) )
    {

        fprintf(stderr, "%s\n", samOut.GetStatusMessage());
        return(samOut.GetStatus());
    }
  
    if(!samOut.OpenForWrite(argv[2]))
    {
        fprintf(stderr, "%s\n", samOut.GetStatusMessage());
        return(samOut.GetStatus());
    }
  
    fprintf(stderr,"Arguments in effect: \n");
    fprintf(stderr,"\tInput file : %s\n",argv[1]);
    fprintf(stderr,"\tOutput file : %s\n",argv[2]);
  
    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);
    // Write the sam header.
    samOut.WriteHeader(samHeader);
  
    // Set returnStatus to success.  It will be changed
    // to the failure reason on failure.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;

    SamRecord samRecord;
  
    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Successfully read a record from the file.

        // Remove the record if it is a duplicate.
        if(SamFlag::isDuplicate(samRecord.getFlag()))
        {
            // Duplicate, so do not write it to the output
            // file and just continue to the next record.
            continue;
        }

        // Remove the OQ tag.
        if(!samRecord.rmTag("OQ", 'Z'))
        {
            // Failed to remove a tag.
            SamStatus errorStatus = samRecord.getStatus();
            fprintf(stderr, "%s\n", errorStatus.getStatusMessage());
            returnStatus = errorStatus.getStatus();
        }

        if(!samOut.WriteRecord(samHeader, samRecord))
        {
            // Failed to write a record.
            fprintf(stderr, "Failure in writing record %s\n", samOut.GetStatusMessage());
            returnStatus = samOut.GetStatus();
        }
    }
   
    if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        // Failed to read a record.
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        returnStatus = samOut.GetStatus();
    }   
   
    std::cerr << std::endl << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;
    std::cerr << "Number of records written = " << 
        samOut.GetCurrentRecordCount() << std::endl;

    // Since the reads were successful, return the status based
    // on the status of the reads/writes.  If any failed, return
    // their failure status.
    samIn.Close();
    samOut.Close();
    return returnStatus;
}
