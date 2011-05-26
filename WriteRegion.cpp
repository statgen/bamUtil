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
// This file contains the processing for the executable option "writeRegion"
// which writes a file with the reads in the specified region.

#include "WriteRegion.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"

void WriteRegion::writeRegionDescription()
{
    std::cerr << " writeRegion - Write a file with reads in the specified region" << std::endl;
}

void WriteRegion::description()
{
    writeRegionDescription();
}

void WriteRegion::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam writeRegion --in <inputFilename>  --out <outputFilename> [--bamIndex <bamIndexFile>] [--noeof] [--refName <reference Name> | --refID <reference ID>] [--start <0-based start pos>] [--end <0-based end psoition>] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in       : the BAM file to be read" << std::endl;
    std::cerr << "\t\t--out      : the SAM/BAM file to write to" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--noeof  : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--bamIndex : the path/name of the bam index file" << std::endl;
    std::cerr << "\t\t             (if not specified, uses the --in value + \".bai\")" << std::endl;
    std::cerr << "\t\t--refName  : the BAM reference Name to read" << std::endl;
    std::cerr << "\t\t             Either this or refID can be specified." << std::endl;
    std::cerr << "\t\t             Defaults to all references." << std::endl;
    std::cerr << "\t\t--refID    : the BAM reference ID to read." << std::endl;
    std::cerr << "\t\t             Either this or refName can be specified." << std::endl;
    std::cerr << "\t\t             Defaults to all references." << std::endl;
    std::cerr << "\t\t             Specify -1 for unmapped" << std::endl;
    std::cerr << "\t\t--start    : inclusive 0-based start position." << std::endl;
    std::cerr << "\t\t             Defaults to -1: meaning from the start of the reference." << std::endl;
    std::cerr << "\t\t             Only applicable if refName/refID is set." << std::endl;
    std::cerr << "\t\t--end      : exclusive 0-based end position." << std::endl;
    std::cerr << "\t\t             Defaults to -1: meaning til the end of the reference." << std::endl;
    std::cerr << "\t\t             Only applicable if refName/refID is set." << std::endl;
    std::cerr << "\t\t--readName : print only reads with this name." << std::endl;
    std::cerr << "\t\t--params   : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int WriteRegion::execute(int argc, char **argv)
{
    // Extract command line arguments.
    static const int UNSPECIFIED_INT = -1;
    static const int UNSET_REF = -2;
    String inFile = "";
    String outFile = "";
    String indexFile = "";
    String refName = "";
    String readName = "";
    int refID = UNSET_REF;
    int start = UNSPECIFIED_INT;
    int end = UNSPECIFIED_INT;
    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
        LONG_STRINGPARAMETER("refName", &refName)
        LONG_INTPARAMETER("refID", &refID)
        LONG_INTPARAMETER("start", &start)
        LONG_INTPARAMETER("end", &end)
        LONG_STRINGPARAMETER("readName", &readName)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    inputParameters.Read(argc-1, &(argv[1]));

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
        // mandatory argument was not specified.
        inputParameters.Status();
        std::cerr << "Missing mandatory argument: --in" << std::endl;
        return(-1);
    }
    if(outFile == "")
    {
        usage();
        // mandatory argument was not specified.
        inputParameters.Status();
        std::cerr << "Missing mandatory argument: --out" << std::endl;
        return(-1);
    }
    
    if(indexFile == "")
    {
        // In file was not specified, so set it to the in file
        // + ".bai"
        indexFile = inFile + ".bai";
    }

    if(refID != UNSET_REF && refName.Length() != 0)
    {
        std::cerr << "Can't specify both refID and refName" << std::endl;
        inputParameters.Status();
        return(-1);
    }

    if(params)
    {
        inputParameters.Status();
    }

    SamFile samIn;
    // Open the file for reading.   
    samIn.OpenForRead(inFile);

    // If refName is set, use that.
    if(refName.Length() != 0)
    {
        // Use Reference Name.
        samIn.SetReadSection(refName.c_str(), start, end);
    }
    else if(refID != UNSET_REF)
    {
        // Use Reference ID
        samIn.SetReadSection(refID, start, end);
    }

    // Open the output file for writing.
    SamFile samOut;
    samOut.OpenForWrite(outFile);

    // Open the bam index file for reading.
    samIn.ReadBamIndex(indexFile);

    // Read & write the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);
    samOut.WriteHeader(samHeader);

    // Read the sam records.
    SamRecord samRecord;
    // Track the status.
    int numSectionRecords = 0;

    // Set returnStatus to success.  It will be changed
    // to the failure reason if any of the writes fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;

    // Keep reading records until they aren't anymore.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        if(!readName.IsEmpty())
        {
            // Check for readname.
            if(strcmp(samRecord.getReadName(), readName.c_str()) != 0)
            {
                // not a matching read name, so continue to the next record.
                continue;
            }
        }
        
        // Successfully read a record from the file, so write it.
        samOut.WriteRecord(samHeader, samRecord);
        ++numSectionRecords;
    }

    std::cerr << "Wrote " << outFile << " with " << numSectionRecords
              << " records.\n";
    return(returnStatus);
}
