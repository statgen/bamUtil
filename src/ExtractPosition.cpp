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
// This file contains the processing for the executable option "extractPosition"
// which writes a file with the reads in the specified region.

#include "ExtractPosition.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"
#include "BaseUtilities.h"

ExtractPosition::ExtractPosition()
    : BamExecutable(),
      myWithinReg(false),
      myWroteReg(false),
      myStart(UNSPECIFIED_INT),
      myEnd(UNSPECIFIED_INT),
      myPrevStart(UNSPECIFIED_INT),
      myPrevEnd(UNSPECIFIED_INT),
      myRefID(UNSPECIFIED_INT),
      myRefName(),
      myPrevRefName(),
      myBedRefID(SamReferenceInfo::NO_REF_ID),
      myBedFile(NULL)
{
    
}

void ExtractPosition::extractPositionDescription()
{
    std::cerr << " extractPosition - Print position specific data" << std::endl;
}

void ExtractPosition::description()
{
    extractPositionDescription();
}

void ExtractPosition::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam writeRegion --in <inputFilename> [--bamIndex <bamIndexFile>] "
              << "[--refName <reference Name>] [--start <0-based start pos>] "
              << "[--params] [--noeof]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in        : the BAM file to be read" << std::endl;
    std::cerr << "\tOptional Parameters for Specifying a Region:" << std::endl;
    std::cerr << "\t\t--bamIndex  : the path/name of the bam index file" << std::endl;
    std::cerr << "\t\t              (if not specified, uses the --in value + \".bai\")" << std::endl;
    std::cerr << "\t\t--refName   : the BAM reference Name to read" << std::endl;
    std::cerr << "\t\t              Either this or refID can be specified." << std::endl;
    std::cerr << "\t\t              Defaults to all references." << std::endl;
    std::cerr << "\t\t--start     : inclusive 0-based start position." << std::endl;
    std::cerr << "\t\t              Defaults to -1: meaning from the start of the reference." << std::endl;
    std::cerr << "\t\t              Only applicable if refName/refID is set." << std::endl;
    std::cerr << "\tOptional Parameters For Other Operations:\n";
    std::cerr << "\t\t--params    : print the parameter settings" << std::endl;
    std::cerr << "\t\t--noeof     : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << std::endl;
}


int ExtractPosition::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String indexFile = "";
    String readName = "";
    myStart = UNSPECIFIED_INT;
    myEnd = UNSPECIFIED_INT;
    myPrevStart = UNSPECIFIED_INT;
    myPrevEnd = UNSPECIFIED_INT;
    myRefID = UNSET_REF;
    myRefName.Clear();
    myPrevRefName.Clear();
    myBedRefID = SamReferenceInfo::NO_REF_ID;
    bool noeof = false;
    bool params = false;
    myWithinReg = false;
    myWroteReg = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_PARAMETER_GROUP("Optional Region Parameters")        
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
        LONG_STRINGPARAMETER("refName", &myRefName)
        LONG_INTPARAMETER("start", &myStart)
        LONG_PARAMETER_GROUP("Optional Other Parameters")
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
    
    if(indexFile == "")
    {
        // In file was not specified, so set it to the in file
        // + ".bai"
        indexFile = inFile + ".bai";
    }

    if(params)
    {
        inputParameters.Status();
    }

    if((myRefName.Length() == 0) || (myStart == -1))
    {
        usage();
        // mandatory argument was not specified.
        inputParameters.Status();
        std::cerr << "Missing mandatory argument: --refName/--start" << std::endl;
        return(-1);
    }

    // Open the file for reading.   
    mySamIn.OpenForRead(inFile);

    // Open the bam index file for reading if a region was specified.
    mySamIn.ReadBamIndex(indexFile);

    // Read & write the sam header.
    mySamIn.ReadHeader(mySamHeader);

    // Read the sam records.
    SamRecord samRecord;

    // Set returnStatus to success.  It will be changed
    // to the failure reason if any of the writes fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;
        
    mySamIn.SetReadSection(myRefName.c_str(), myStart, myStart+1);
    // Keep reading records until they aren't anymore.
    while(mySamIn.ReadRecord(mySamHeader, samRecord))
    {
        // Successfully read a record from the file, so write it.
        // Write position.
        Cigar* cigar = samRecord.getCigarInfo();
        int32_t index = cigar->getQueryIndex(myStart, samRecord.get0BasedPosition());
        if(index != Cigar::INDEX_NA)
        {
            std::cerr << samRecord.getReadName() << "\t"
                      << samRecord.getSequence(index) << "\t"
                      << (int)BaseUtilities::getPhredBaseQuality(samRecord.getQuality(index)) << "\n";
        }
    }
    return(returnStatus);
}
