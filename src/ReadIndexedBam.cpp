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

//////////////////////////////////////////////////////////////////////////
// This file contains the processing for the executable option "readIndexedBam"
// which reads an indexed BAM file by chromosome and writes it into a new
// file sorted from reference id -1 to maxRefID.

#include "ReadIndexedBam.h"
#include "SamFile.h"
#include "Parameters.h"
#include "SamValidation.h"
#include "PhoneHome.h"

void ReadIndexedBam::readIndexedBamDescription()
{
    std::cerr << " readIndexedBam - Read Indexed BAM By Reference and write it from reference id -1 to maxRefId" << std::endl;
}


void ReadIndexedBam::description()
{
    readIndexedBamDescription();
}


void ReadIndexedBam::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam readIndexedBam <inputFilename> <outputFile.sam/bam> <bamIndexFile>" << std::endl;
}


int ReadIndexedBam::execute(int argc, char ** argv)
{
    if(argc != 5)
    {
        String noPhArg = "--noPhoneHome";
        if((argc != 6) || (noPhArg.SlowCompareToStem(argv[5]) != 0))
        {
            usage();
            exit(-1);
        }
    }
    else
    {
        PhoneHome::checkVersion(getProgramName(), VERSION);
    }
    return(readIndexedBam(argv[2], argv[3], argv[4]));
}

int ReadIndexedBam::readIndexedBam(const char* inputFilename,
                                   const char* outputFilename,
                                   const char* indexFilename)
{
    // Open the input file for reading.
    SamFile samIn;
    samIn.OpenForRead(inputFilename);

    // Open the bam index file for reading.
    samIn.ReadBamIndex(indexFilename);

    // Open the output file for writing.
    SamFile samOut;
    samOut.OpenForWrite(outputFilename);

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    // Write the sam header.
    samOut.WriteHeader(samHeader);

    SamRecord samRecord;
    SamValidationErrors samInvalidErrors;

    // Get the number of references.
    int numReferences = samHeader.getReferenceInfo().getNumEntries();

    // Loop through each Reference.
    for(int i = -1; i < numReferences; i++)
    {
        int numSectionRecords = 0;
        samIn.SetReadSection(i);
        // Keep reading records until they aren't more.
        while(samIn.ReadRecord(samHeader, samRecord))
        {
            numSectionRecords++;
            // Successfully read a record from the file, so write it.
            samOut.WriteRecord(samHeader, samRecord);
        }

        std::cerr << "Reference ID " << i << " has " << numSectionRecords 
                  << " records" << std::endl;
    }
   
    std::cerr << "Number of records = " << samIn.GetCurrentRecordCount() 
              << std::endl;
   
    return(0);
}


