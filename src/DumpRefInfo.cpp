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
// This file contains the processing for the executable option "dumpRefInfo"
// which prints the SAM/BAM Reference Name information to the screen.

#include "DumpRefInfo.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"

void DumpRefInfo::dumpRefInfoDescription()
{
    std::cerr << " dumpRefInfo - Print SAM/BAM Reference Name Information" 
              << std::endl;
}


void DumpRefInfo::description()
{
    dumpRefInfoDescription();
}


void DumpRefInfo::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam dumpRefInfo --in <inputFilename> [--noeof] [--printRecordRefs] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in               : the SAM/BAM file to be read" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--noeof            : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--printRecordRefs  : print the reference information for the records in the file (grouped by reference)." << std::endl;
    std::cerr << "\t\t--params           : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


// Dump the reference information from specified SAM/BAM file.
int DumpRefInfo::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    bool noeof = false;
    bool printRecordRefs = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("printRecordRefs", &printRecordRefs)
        LONG_PARAMETER("params", &params)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);

    // If no eof block is required for a bgzf file, set the bgzf file type to 
    // not look for it.
    if(noeof)
    {
        // Set that the eof block is not required.
        BgzfFileType::setRequireEofBlock(false);
    }

    // Check to see if the in file was specified, if not, report an error.
    if(inFile == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--in is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Open the input file for reading.
    SamFile samIn;
    samIn.OpenForRead(inFile);

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    const SamReferenceInfo& refInfo = samHeader.getReferenceInfo();
    int numReferences = refInfo.getNumEntries();
    
    for(int i = 0; i < numReferences; i++)
    {
        std::cout << "Reference Index " << i;
        std::cout << "; Name: " << refInfo.getReferenceName(i)
                  << std::endl;
    }
    if(numReferences == 0)
    {
        // There is no reference info.
        std::cerr << "The header contains no reference information.\n";
    }

    // If we are to print the references as found in the records, loop
    // through reading the records.
    if(printRecordRefs)
    {
        SamRecord samRecord;

        // Track the prev name/id.
        std::string prevName = "";
        int prevID = -2;
        int recCount = 0; // track the num records in a ref.
        // Keep reading records until ReadRecord returns false.
        while(samIn.ReadRecord(samHeader, samRecord))
        {
            const char* name = samRecord.getReferenceName();
            int id = samRecord.getReferenceID();
            if((strcmp(name, prevName.c_str()) != 0) || (id != prevID))
            {
                if(prevID != -2)
                {
                    std::cout << "\tRef ID: " << prevID
                              << "\tRef Name: " << prevName 
                              << "\tNumRecs: " << recCount
                              << std::endl;
                }
                recCount = 0;
                prevID = id;
                prevName = name;
            }
            ++recCount;
        }
        // Print the last index.
        if(prevID != -2)
        {
            std::cout << "\tRef ID: " << prevID
                      << "\tRef Name: " << prevName 
                      << "\tNumRecs: " << recCount
                      << std::endl;
        }
    }
    return(SamStatus::SUCCESS);
}
