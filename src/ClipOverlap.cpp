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
#include "SamFlag.h"

ClipOverlap::ClipOverlap()
    : BamExecutable(),
      myStoreOrig("")
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
    myStoreOrig = "";

    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_STRINGPARAMETER("storeOrig", &myStoreOrig)
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

    if((myStoreOrig.Length() != 0) && (myStoreOrig.Length() != 2))
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
    uint16_t flag = 0;
    if((samRecord == NULL) || (tmpRecord == NULL))
    {
        std::cerr << "Failed to allocate a SamRecord, so exit.\n";
        return(-1);
    }

    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, *samRecord))
    {
        // Determine if the read is mapped.
        flag = samRecord->getFlag();
        if(!SamFlag::isMapped(flag))
        {
            // This record is not mapped, no clipping will be done, so
            // just write this record and move onto the next record without
            // storing this one.
            if(!samOut.WriteRecord(samHeader, *samRecord))
            {
                // Failed to write a record.
                fprintf(stderr, "%s\n", samOut.GetStatusMessage());
                returnStatus = samOut.GetStatus();
            }
            continue;
        }

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
        if(strcmp(samRecord->getReadName(), prevSamRecord->getReadName()) == 0)
        {
            // Read name match, so check for overlap, this one
            // starts between the previous records start & end.
            if(prevSamRecord->getReferenceID() == samRecord->getReferenceID())
            {
                // Determine which read starts first.
                if(prevSamRecord->get0BasedPosition() <= 
                   samRecord->get0BasedPosition())
                {
                    // The previous read starts at or before the current one.
                    clip(*prevSamRecord, *samRecord);
                }
                else
                {
                    // The current read starts before the previous one.
                    clip(*samRecord, *prevSamRecord);
                }
            }

            // Found a read pair, so write both records.
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


void ClipOverlap::clip(SamRecord& firstRecord, SamRecord& secondRecord)
{
    static CigarRoller newFirstCigar; // holds updated cigar.
    static CigarRoller newSecondCigar; // holds updated cigar.

    // Used for checking forward/reverse logic.
    uint16_t firstFlag = firstRecord.getFlag();
    uint16_t secondFlag = secondRecord.getFlag();

    // Check for overlap.
    // We already know that the first record starts at or before the 2nd,
    // so there is overlap if the first record ends at or after the 2nd
    // one starts.
    int32_t firstEnd = firstRecord.get0BasedAlignmentEnd();
    int32_t secondStart = secondRecord.get0BasedPosition();
    if(firstEnd >= secondStart)
    {
        // overlap, determine which record will get clipped by determining
        // which record has a lower base quality in the overlapping region.
        // First check clipping the first region.

        // Determine the clipping on the 1st record at the start of the 2nd.
        int32_t firstClipPos = 
            CigarHelper::softClipEndByRefPos(firstRecord, secondStart, 
                                             newFirstCigar);
        // Loop through counting the quality of the clipped bases.
        // They run from the firstClip to the length of the read.
        int32_t firstQualSum = getSumQual(firstRecord, firstClipPos, 
                                  firstRecord.getReadLength()-1);

        int32_t newPos = 0;
        int32_t secondClipPos =
            CigarHelper::softClipBeginByRefPos(secondRecord, firstEnd,
                                               newSecondCigar, newPos);
        // Loop through counting the quality of the clipped bases.
        // They run from the beginning until the secondClip(included).
        int32_t secondQualSum = getSumQual(secondRecord, 0, secondClipPos);

        // Check to see whether the 1st record or the 2nd one should be clipped
        // based on which has the lower quality.
        if(firstQualSum <= secondQualSum)
        {
            // First clip has lower or equal quality, so clip that.
            // Check to see if the entire read should be clipped by
            // checking forward/reverse.
            if(SamFlag::isReverse(firstFlag) && !SamFlag::isReverse(secondFlag))
            {
                // first record is reverse and 2nd is forward, so clip
                // those extending ends which means clipping all of the
                // reverse(first) strand since it's overlap has lower quality.
                clipEntire(firstRecord);

                // Soft clip the end of the forward(second) strand that
                // extends past the reverse strand (firstEnd + 1).
                if(CigarHelper::softClipEndByRefPos(secondRecord, firstEnd+1,
                                                    newSecondCigar)
                   != CigarHelper::NO_CLIP)
                {
                    // Write the original cigar into the specified tag.
                    if(!myStoreOrig.IsEmpty())
                    {
                        // Write original cigar.
                        secondRecord.addTag(myStoreOrig, 'Z', 
                                            secondRecord.getCigar());
                    }
                    secondRecord.setCigar(newSecondCigar);
                }
            }
            else
            {
                // No strand specific extended ends clipping is required,
                // so just clip the overlap.
                // Write the original cigar into the specified tag.
                if(!myStoreOrig.IsEmpty())
                {
                    // Write original cigar.
                    firstRecord.addTag(myStoreOrig, 'Z', 
                                        firstRecord.getCigar());
                }
                firstRecord.setCigar(newFirstCigar);
            }
        }
        else
        {
            // The 2nd clip has lower quality, so clip that.
            // Check to see if the entire read should be clipped by
            // checking forward/reverse.
            if(SamFlag::isReverse(firstFlag) && !SamFlag::isReverse(secondFlag))
            {
                // first record is reverse and 2nd is forward, so clip
                // those extending ends which means clipping all of the
                // forward(second) strand since it's overlap has lower quality.
                clipEntire(secondRecord);

                // Soft clip the front of the reverse(first) strand that
                // extends past the forward strand (secondStart - 1).
                if(CigarHelper::softClipBeginByRefPos(firstRecord, 
                                                      secondStart-1,
                                                      newFirstCigar,
                                                      newPos)
                   != CigarHelper::NO_CLIP)
                {
                    // Write the original cigar into the specified tag.
                    if(!myStoreOrig.IsEmpty())
                    {
                        // Write original cigar.
                        firstRecord.addTag(myStoreOrig, 'Z', 
                                            firstRecord.getCigar());
                    }
                    firstRecord.setCigar(newFirstCigar);
                    firstRecord.set0BasedPosition(newPos);
                 }
            }
            else
            {
                // No strand specific extended ends clipping is required,
                // so just clip the overlap.
                // Write the original cigar into the specified tag.
                if(!myStoreOrig.IsEmpty())
                {
                    // Write original cigar.
                    secondRecord.addTag(myStoreOrig, 'Z', 
                                         secondRecord.getCigar());
                }
                secondRecord.set0BasedPosition(newPos);
                secondRecord.setCigar(newSecondCigar);
            }
        }
    }
    else
    {
        // No overlap, but check to verify that the strands are not in
        // the wrong order, because in that case, they should both be clipped.
        if(SamFlag::isReverse(firstFlag) && !SamFlag::isReverse(secondFlag))
        {
            // first record is reverse and 2nd is forward, which means
            // the entire reverse record is before the forward one, so
            // clip them both.
            clipEntire(firstRecord);
            clipEntire(secondRecord);
        }
    }
}


void ClipOverlap::clipEntire(SamRecord& record)
{
    static CigarRoller newCigar; // holds updated cigar.

    // Clip the entire record.
    if(CigarHelper::softClipEndByRefPos(record, record.get0BasedPosition(),
                                        newCigar) != CigarHelper::NO_CLIP)
    { 
        // Write the original cigar into the specified tag.
        if(!myStoreOrig.IsEmpty())
        {
            // Write original cigar.
            record.addTag(myStoreOrig, 'Z', record.getCigar());
        }
        // Update the cigar.
        record.setCigar(newCigar);
    }
}


int32_t ClipOverlap::getSumQual(SamRecord& record, 
                                int32_t startPos, int32_t endPos)
{
    int32_t qualSum = 0;
    const char* quality = record.getQuality();
    for(int i = startPos; i <= endPos; i++)
    {
        qualSum += quality[i];
        // Check for the null terminator at the end of the quality string.
        if(quality[i] == 0)
        {
            // Qual is shorter than the read, so break.
            break;
        }
    }
    return(qualSum);
}    
