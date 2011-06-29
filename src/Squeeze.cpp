/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan
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

#include "Squeeze.h"
#include "BgzfFileType.h"
// #include <stdio.h>
// #include <string.h>
// #include "SamFile.h"
#include "SamFlag.h"

Squeeze::Squeeze()
{
}

Squeeze::~Squeeze()
{
}

void Squeeze::squeezeDescription()
{
    std::cerr << " squeeze -  reduces files size by dropping OQ fields, duplicates, specified tags, using '=' when a base matches the reference." << std::endl;
}


void Squeeze::description()
{
    squeezeDescription();
}


// print Usage
void Squeeze::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam squeeze --in <inputFile> --out <outputFile.sam/bam/ubam (ubam is uncompressed bam)> [--cigar] [--qual] [--keepTags] [--rmBQ] [--rmTags <Tag:Type[;Tag:Type]*>] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in         : the SAM/BAM file to be read" << std::endl;
    std::cerr << "\t\t--out        : the SAM/BAM file to be written" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--refFile    : reference file name used to convert any bases that match the reference to '-'" << std::endl;
    std::cerr << "\t\t--keepOQ     : keep the OQ tag rather than removing it.  Default is to remove it." << std::endl;
    std::cerr << "\t\t--keepDups   : keep duplicates rather than removing records marked duplicate.  Default is to remove them." << std::endl;
    std::cerr << "\t\t--rmTags     : Remove the specified Tags formatted as Tag:Type;Tag:Type;Tag:Type..." << std::endl;
    std::cerr << "\t\t--noeof      : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params     : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}

// main function
int Squeeze::execute(int argc, char ** argv)
{
    String inFile = "";
    String outFile = "";
    String refFile = "";
    bool noeof = false;
    bool params = false;
    bool keepOQ = false;
    bool keepDups =  false;
    String rmTags = "";

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_PARAMETER("keepOQ", &keepOQ)
        LONG_PARAMETER("keepDups", &keepDups)
        LONG_STRINGPARAMETER("rmTags", &rmTags)
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
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--in is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(outFile == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--out is a mandatory argument, "
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

    // Open the output file for writing.
    SamFile samOut;
    samOut.OpenForWrite(outFile);
    // Check to see if the ref file was specified.
    // Open the reference.
    GenomeSequence* refPtr = NULL;
    if(refFile != "")
    {
        refPtr = new GenomeSequence(refFile);
        // Since a reference was specified, convert matching bases to '='.
        samOut.SetWriteSequenceTranslation(SamRecord::EQUAL);
        // Set the reference for the output file.
        samOut.SetReference(refPtr);
    }

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    // Write the sam header.
    samOut.WriteHeader(samHeader);

    SamRecord samRecord;

    // Set returnStatus to success.  It will be changed to the
    // failure reason if any of the writes or updates fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;
  
    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Successfully read a record from the file.

        // Remove the record if it is a duplicate and we are not
        // supposed to keep duplicates.
        if(!keepDups && SamFlag::isDuplicate(samRecord.getFlag()))
        {
            // Duplicate, so do not write it to the output
            // file and just continue to the next record.
            continue;
        }

        // Remove the OQ tag if we are not supposed to remove OQ tags.
        if(!keepOQ)
        {
            if(!samRecord.rmTag("OQ", 'Z'))
            {
                // Failed to remove a tag.
                SamStatus errorStatus = samRecord.getStatus();
                fprintf(stderr, "%s\n", errorStatus.getStatusMessage());
                returnStatus = errorStatus.getStatus();
            }
        }

        // Remove any specified tags.
        if(!rmTags.IsEmpty())
        {
            if(!samRecord.rmTags(rmTags.c_str()))
            {
                // Failed to remove the specified tags.
                fprintf(stderr, "%s\n", samIn.GetStatusMessage());
                returnStatus = samIn.GetStatus();
            }
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
