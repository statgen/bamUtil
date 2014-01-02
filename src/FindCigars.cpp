/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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
// This file contains the processing for the executable option "findCigars"
// which reads a SAM/BAM file and writes a SAM/BAM file with just the 
// reads that contain cigars that have any of the specified CIGAR operations.

#include <bitset>
#include "FindCigars.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"

FindCigars::FindCigars()
{
}

void FindCigars::findCigarsDescription()
{
    std::cerr << " findCigars - Output just the reads that contain any of the specified CIGAR operations." << std::endl;
}


void FindCigars::description()
{
    findCigarsDescription();
}


void FindCigars::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam findCigars --in <inputFile> --out <outputFile.sam/bam/ubam (ubam is uncompressed bam)> [--cinsert] [--cdel] [--cpad] [--cskip] [--chardClip] [--csoftClip] [--nonM] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in         : the SAM/BAM file to be read" << std::endl;
    std::cerr << "\t\t--out        : the SAM/BAM file to be written" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--cinsert     : output reads that contain insertions ('I')." << std::endl;
    std::cerr << "\t\t--cdel        : output reads that contain deletions ('D')." << std::endl;
    std::cerr << "\t\t--cpad        : output reads that contain pads ('P')." << std::endl;
    std::cerr << "\t\t--cskip       : output reads that contain skips ('N')." << std::endl;
    std::cerr << "\t\t--chardClip   : output reads that contain hard clips ('H')." << std::endl;
    std::cerr << "\t\t--csoftClip   : output reads that contain soft clips ('S')." << std::endl;
    std::cerr << "\t\t--nonM       : output reads that contain any non match/mismatch (anything other than 'M', '=', 'X')" << std::endl;
    std::cerr << "\t\t--noeof      : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params     : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int FindCigars::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    bool insert = false;
    bool del = false;
    bool pad = false;
    bool skip = false;
    bool hardClip = false;
    bool softClip = false;
    bool nonM = false;
    bool noeof = false;
    bool params = false;

    std::bitset<Cigar::MAX_OP_VALUE+1> desiredOps;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_PARAMETER("cinsert", &insert)
        LONG_PARAMETER("cdel", &del)
        LONG_PARAMETER("cpad", &pad)
        LONG_PARAMETER("cskip", &skip)
        LONG_PARAMETER("chardClip", &hardClip)
        LONG_PARAMETER("csoftClip", &softClip)
        LONG_PARAMETER("nonM", &nonM)
        LONG_PARAMETER("noeof", &noeof)
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

    if(outFile == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--out is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(nonM)
    {
        insert = true;
        del = true;
        pad = true;
        skip = true;
        hardClip = true;
        softClip = true;
    }

    if(insert)
    {
        desiredOps[Cigar::insert] = true;
    }
    if(del)
    {
        desiredOps[Cigar::del] = true;
    }
    if(pad)
    {
        desiredOps[Cigar::pad] = true;
    }
    if(skip)
    {
        desiredOps[Cigar::skip] = true;
    }
    if(hardClip)
    {
        desiredOps[Cigar::hardClip] = true;
    }
    if(softClip)
    {
        desiredOps[Cigar::softClip] = true;
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Check that there is at least one cigar character to search for.
    if(!insert && !del && !pad && !skip && !hardClip && !softClip)
    {
        std::cerr << "no CIGAR operations were specified, so no reads to output." 
                  << std::endl;
        return(-1);
    }

    // Open the input file for reading.
    SamFile samIn;
    samIn.OpenForRead(inFile);

    // Open the output file for writing.
    SamFile samOut;
    samOut.OpenForWrite(outFile);

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
        // Check the cigar.
        Cigar* cigarPtr = samRecord.getCigarInfo();
        if(cigarPtr == NULL)
        {
            std::cerr << "Failed to get a cigar.";
            return(-1);
        }
        for(int i = 0; i < cigarPtr->size(); i++)
        {
            Cigar::Operation op = cigarPtr->getOperator(i).operation;

            if(desiredOps[op])
            {
                // This cigar has a desired operation, so write the record.
                if(!samOut.WriteRecord(samHeader, samRecord))
                {
                    // Failed to write a record.
                    fprintf(stderr, "%s\n", samOut.GetStatusMessage());
                    returnStatus = samOut.GetStatus();
                }
                // Found a desired operation, so skip to the next string.
                break;
            }
        }
    }
    std::cerr << std::endl << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;
    std::cerr << "Number of records written = " << 
        samOut.GetCurrentRecordCount() << std::endl;

    // Since the reads were successful, return the status based
    // on the status of the writes.  If any failed, return
    // their failure status.
    return(returnStatus);
}
