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
// This file contains the processing for the executable option "convert"
// which reads an SAM/BAM file and writes a SAM/BAM file (it can convert 
// between SAM and BAM formats).

#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"
#include "SamValidation.h"


void convertDescription()
{
    std::cerr << " convert - Convert SAM/BAM to SAM/BAM:" << std::endl;
    std::cerr << "\t./bam convert --in <inputFile> --out <outputFile.sam/bam/ubam (ubam is uncompressed bam)> [--refFile] [--seqBases|--seqEquals|--seqOrig] [--noeof] [--params]" << std::endl;
}


void convertUsage()
{
    convertDescription();
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in       : the SAM/BAM file to be read" << std::endl;
    std::cerr << "\t\t--out      : the SAM/BAM file to be written" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--refFile  : reference file name" << std::endl;
    std::cerr << "\t\t--noeof    : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params   : print the parameter settings" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\tOptional Sequence Parameters (only specify one):" << std::endl;
    std::cerr << "\t\t--seqOrig  : Leave the sequence as is (default & used if reference is not specified)." << std::endl;
    std::cerr << "\t\t--seqBases : Convert any '=' in the sequence to the appropriate base using the reference (requires --ref)." << std::endl;
    std::cerr << "\t\t--seqEquals : Convert any bases that match the reference to '=' (requires --ref)." << std::endl;
}


int convert(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    String refFile = "";
    bool noeof = false;
    bool params = false;

    bool seqBases = false;
    bool seqEquals = false;
    bool seqOrig = false;
    
    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PARAMETER_GROUP("SequenceConversion")
            EXCLUSIVE_PARAMETER("seqBases", &seqBases)
            EXCLUSIVE_PARAMETER("seqEquals", &seqEquals)
            EXCLUSIVE_PARAMETER("seqOrig", &seqOrig)
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
        convertUsage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--in is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(outFile == "")
    {
        convertUsage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--out is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    // Check to see if the ref file was specified.
    // Open the reference.
    GenomeSequence* refPtr = NULL;
    if(refFile != "")
    {
        refPtr = new GenomeSequence(refFile);
    }

    SamRecord::SequenceTranslation translation;
    if((seqBases) && (refPtr != NULL))
    {
        translation = SamRecord::BASES;
    }
    else if((seqEquals) && (refPtr != NULL))
    {
        translation = SamRecord::EQUAL;
    }
    else
    {
        seqOrig = true;
        translation = SamRecord::NONE;
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
    samOut.SetWriteSequenceTranslation(translation);
    samOut.SetReference(refPtr);

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    // Write the sam header.
    samOut.WriteHeader(samHeader);

    SamRecord samRecord;
    SamValidationErrors samInvalidErrors;

    // Set returnStatus to success.  It will be changed
    // to the failure reason if any of the writes fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;

    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Successfully read a record from the file, so write it.
        if(!samOut.WriteRecord(samHeader, samRecord))
        {
            // Failed to write a record.
            fprintf(stderr, "%s\n", samOut.GetStatusMessage());
            returnStatus = samOut.GetStatus();
        }
    }

    std::cerr << std::endl << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;
    std::cerr << "Number of records written = " << 
        samOut.GetCurrentRecordCount() << std::endl;

    if(refPtr != NULL)
    {
        delete(refPtr);
    }

    // Since the reads were successful, return the status based
    // on the status of the writes.  If any failed, return
    // their failure status.
    return(returnStatus);
}


