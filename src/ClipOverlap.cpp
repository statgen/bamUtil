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
// This file contains the processing for the executable option "clipOverlap"
// which clips overlapping read pairs.
#include "ClipOverlap.h"
#include "SamFile.h"
#include "BgzfFileType.h"
#include "CigarHelper.h"


ClipOverlap::ClipOverlap()
    : BamExecutable()
{
}


void ClipOverlap::clipOverlapDescription()
{
    std::cerr << " clipOverlap - Clip overlapping read pairs in a SAM/BAM File already sorted by ReadName" << std::endl;
}


void ClipOverlap::description()
{
    clipOverlapDescription();
}


void ClipOverlap::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam clipOverlap --in <inputFile> --out <outputFile> [--storeOrig <tag>] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in : the SAM/BAM file to clip overlaping read pairs for" << std::endl;
    std::cerr << "\t\t--out        : the SAM/BAM file to be written" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--storeOrig   : Store the original cigar in the specified tag." << std::endl;
    std::cerr << "\t\t--noeof       : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params      : Print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int ClipOverlap::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    String storeOrig = "";

    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_STRINGPARAMETER("storeOrig", &storeOrig)
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

    // Check to see if the out file was specified, if not, report an error.
    if(outFile == "")
    {
        usage();
        inputParameters.Status();
        // Out file was not specified but it is mandatory.
        std::cerr << "--out is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if((storeOrig.Length() != 0) && (storeOrig.Length() != 2))
    {
        usage();
        inputParameters.Status();
        std::cerr << "--storeOrig tag name must be 2 characters.\n";
        return(-1);
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Open the files for reading/writing.
    SamFile samIn;
    samIn.OpenForRead(inFile);
    SamFile samOut;
    samOut.OpenForWrite(outFile);

    // Read/write the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);
    samOut.WriteHeader(samHeader);


    // Set returnStatus to success.  It will be changed
    // to the failure reason if any of the writes fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;

    // Read the sam records.
    SamRecord* prevSamRecord = NULL;
    SamRecord* samRecord = new SamRecord;
    SamRecord* tmpRecord = new SamRecord;
    if((samRecord == NULL) || (tmpRecord == NULL))
    {
        std::cerr << "Failed to allocate a SamRecord, so exit.\n";
        return(-1);
    }

    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, *samRecord))
    {
        if(prevSamRecord == NULL)
        {
            // Nothing to compare this record to, so set this record to the
            // previous, and the next record.
            prevSamRecord = samRecord;
            samRecord = tmpRecord;
            tmpRecord = NULL;
            continue;
        }

        // Check if the read name matches the previous read name.
        if(samRecord->getReadName() == prevSamRecord->getReadName())
        {
            // Read name match, so check for overlap, this one
            // starts between the previous records start & end.
            uint32_t prevAlignmentStart = prevSamRecord->get0BasedPosition();
            uint32_t prevAlignmentEnd = prevSamRecord->get0BasedAlignmentEnd();
            uint32_t alignmentStart = samRecord->get0BasedPosition();
            if((alignmentStart >= prevAlignmentStart) && 
               (alignmentStart <= prevAlignmentEnd))
            {
                // overlap, determine how much needs to be clipped.
                // Clip from the 
                if(alignmentStart < prevAlignmentStart)
                {
                    // This clip starts first, so clip from the end of it.
                    clipFromRecord(*samRecord, prevAlignmentStart, storeOrig);
                }
                else
                {
                    // The previous record starts first, so clip from the end of it.
                    clipFromRecord(*prevSamRecord, alignmentStart, storeOrig);
                }
            }

            // Found a read pair, so write both the previous record and this record.
            if(!samOut.WriteRecord(samHeader, *prevSamRecord))
            {
                // Failed to write a record.
                fprintf(stderr, "%s\n", samOut.GetStatusMessage());
                returnStatus = samOut.GetStatus();
            }
            if(!samOut.WriteRecord(samHeader, *samRecord))
            {
                // Failed to write a record.
                fprintf(stderr, "%s\n", samOut.GetStatusMessage());
                returnStatus = samOut.GetStatus();
            }
            // Setup for the next read with no previous.
            tmpRecord = prevSamRecord;
            prevSamRecord = NULL;
            
        }
        else
        {
            // Read name does not match, so write the previous record.
            if(!samOut.WriteRecord(samHeader, *prevSamRecord))
            {
                // Failed to write a record.
                fprintf(stderr, "%s\n", samOut.GetStatusMessage());
                returnStatus = samOut.GetStatus();
            }
            // Store this record as the previous.
            tmpRecord = prevSamRecord;
            prevSamRecord = samRecord;
            samRecord = tmpRecord;
            tmpRecord = NULL;
        }
    }

    // Write the previous record if there is one.
    if(prevSamRecord != NULL)
    {
        if(!samOut.WriteRecord(samHeader, *prevSamRecord))
        {
            // Failed to write a record.
            fprintf(stderr, "%s\n", samOut.GetStatusMessage());
            returnStatus = samOut.GetStatus();
        }
        delete prevSamRecord;
    }

    if(samRecord != NULL)
    {
        delete samRecord;
    }
    if(tmpRecord != NULL)
    {
        delete tmpRecord;
    }

    return(returnStatus);
}


void ClipOverlap::clipFromRecord(SamRecord& record, int32_t refPos, 
                                 const String& storeOrig)
{
    static CigarRoller newCigar; // holds updated cigar.
    // track clip position.
    int32_t clipPos = 
        CigarHelper::softClipEndByRefPos(record, refPos, newCigar);
    if(clipPos != CigarHelper::NO_CLIP)
    { 
        // Write the current cigar into a tag.
        if(!storeOrig.IsEmpty())
        {
            record.addTag(storeOrig, 'Z', record.getCigar());
        }
        // Update the cigar.
        record.setCigar(newCigar);
    }
}
