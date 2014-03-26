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
// This file contains the processing for the executable option "revert"
// which reads an SAM/BAM file and writes a SAM/BAM file with the 
// specified previous values restored if the values are known.

#include "Revert.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"
#include "SamTags.h"

Revert::Revert()
    : myKeepTags(false)
{
}

void Revert::revertDescription()
{
    std::cerr << " revert - Revert SAM/BAM replacing the specified fields with their previous values (if known) and removes specified tags" << std::endl;
}


void Revert::description()
{
    revertDescription();
}


void Revert::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam revert --in <inputFile> --out <outputFile.sam/bam/ubam (ubam is uncompressed bam)> [--cigar] [--qual] [--keepTags] [--rmBQ] [--rmTags <Tag:Type[,Tag:Type]*>] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in         : the SAM/BAM file to be read" << std::endl;
    std::cerr << "\t\t--out        : the SAM/BAM file to be written" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--cigar      : update the cigar and the position based on the OC & OP tags." << std::endl;
    std::cerr << "\t\t--qual       : update the quality based on the OQ tag." << std::endl;
    std::cerr << "\t\t--keepTags   : keep the tags that are used to update the record.  Default is to remove them." << std::endl;
    std::cerr << "\t\t--rmBQ       : Remove the BQ Tag." << std::endl;
    std::cerr << "\t\t--rmTags     : Remove the specified Tags formatted as Tag:Type,Tag:Type,Tag:Type..." << std::endl;
    std::cerr << "\t\t--noeof      : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params     : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int Revert::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    bool cigar = false;
    bool qual = false;
    bool noeof = false;
    bool params = false;
    bool rmBQ = false;
    String rmTags = "";
    myKeepTags = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_PARAMETER("cigar", &cigar)
        LONG_PARAMETER("qual", &qual)
        LONG_PARAMETER("keepTags", &myKeepTags)
        LONG_PARAMETER("rmBQ", &rmBQ)
        LONG_STRINGPARAMETER("rmTags", &rmTags)
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
        // Update the cigar & position.
        if(cigar)
        {
            if(!updateCigar(samRecord))
            {
                // Failed to update the cigar & position.
                fprintf(stderr, "%s\n", samIn.GetStatusMessage());
                returnStatus = samIn.GetStatus();
            }
        }
        if(qual)
        {
            if(!updateQual(samRecord))
            {
                // Failed to update the quality.
                fprintf(stderr, "%s\n", samIn.GetStatusMessage());
                returnStatus = samIn.GetStatus();
            }
        }

        if(rmBQ)
        {
            if(!removeBQ(samRecord))
            {
                // Failed to remove BQ.
                fprintf(stderr, "%s\n", samIn.GetStatusMessage());
                returnStatus = samIn.GetStatus();
            }
        }

        if(rmTags != "")
        {
            if(!samRecord.rmTags(rmTags.c_str()))
            {
                // Failed to remove the specified tags.
                fprintf(stderr, "%s\n", samIn.GetStatusMessage());
                returnStatus = samIn.GetStatus();
            }
        }

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

    // Since the reads were successful, return the status based
    // on the status of the writes.  If any failed, return
    // their failure status.
    return(returnStatus);
}

bool Revert::updateCigar(SamRecord& samRecord)
{
    // Get the OC tag, which is a string.
    const String* oldCigar = samRecord.getStringTag(SamTags::ORIG_CIGAR_TAG);
    // Get the OP tag, which is an integer.

    bool status = true;
    if(oldCigar != NULL)
    {
        // The old cigar was found, so set it in the record.
        status &= samRecord.setCigar((*oldCigar).c_str());

        if(!myKeepTags)
        {
            // Remove the tag.
            status &= samRecord.rmTag(SamTags::ORIG_CIGAR_TAG, SamTags::ORIG_CIGAR_TAG_TYPE);
        }
    }

    int oldPos = 0;
    if(samRecord.getIntegerTag(SamTags::ORIG_POS_TAG, oldPos))
    {
        // The old position was found, so set it in the record.
        status &= samRecord.set1BasedPosition(oldPos);

        if(!myKeepTags)
        {
            // Remove the tag.
            status &= samRecord.rmTag(SamTags::ORIG_POS_TAG, SamTags::ORIG_POS_TAG_TYPE);
        }
    }

    return(status);
}


bool Revert::updateQual(SamRecord& samRecord)
{
        // Get the OQ tag, which is a string.
    const String* oldQual = samRecord.getStringTag(SamTags::ORIG_QUAL_TAG);

    bool status = true;
    if(oldQual != NULL)
    {
        // The old quality was found, so set it in the record.
        status &= samRecord.setQuality((*oldQual).c_str());

        if(!myKeepTags)
        {
            // Remove the tag.
            samRecord.rmTag(SamTags::ORIG_QUAL_TAG, SamTags::ORIG_QUAL_TAG_TYPE);
        }
    }
    return(status);
}


bool Revert::removeBQ(SamRecord& samRecord)
{
    // Remove the tag.
    return(samRecord.rmTag(SamTags::BQ_TAG, SamTags::BQ_TAG_TYPE));
}

