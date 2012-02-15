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
#include "SamHelper.h"


ClipOverlap::ClipOverlap()
    : BamExecutable(),
      myPool(),
      myClipsOnly(false),
      myStoreOrig(""),
      myStats(false),
      myOverlaps(),
      myNumForwardClips(0),
      myNumReverseClips(0),
      myNumOrientationClips(0),
      myNumMateFailures(0),
      myNumPoolFail(0),
      myPoolSkipClip(false)
{
}


void ClipOverlap::clipOverlapDescription()
{
    std::cerr << " clipOverlap - Clip overlapping read pairs in a SAM/BAM File already sorted by Coordinate or ReadName" << std::endl;
}


void ClipOverlap::description()
{
    clipOverlapDescription();
}


void ClipOverlap::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam clipOverlap --in <inputFile> --out <outputFile> [--storeOrig <tag>] [--readName] [--poolSize <numRecords allowed to allocate>] [--poolSkipClip] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in          : the SAM/BAM file to clip overlaping read pairs for" << std::endl;
    std::cerr << "\t\t--out         : the SAM/BAM file to be written" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--storeOrig   : Store the original cigar in the specified tag." << std::endl;
    std::cerr << "\t\t--readName    : Original file is sorted by Read Name instead of coordinate." << std::endl;
    std::cerr << "\t\t--stats       : Print some statistics on the overlaps." << std::endl;
    std::cerr << "\t\t--clipsOnly   : Only output the clipped reads" << std::endl;
    std::cerr << "\t\t--noeof       : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params      : Print the parameter settings to stderr" << std::endl;
    std::cerr << "\tClipping By Coordinate Optional Parameters:" << std::endl;
    std::cerr << "\t\t--poolSize    : Maximum number of records the program is allowed to allocate" << std::endl;
    std::cerr << "\t\t                for clipping on Coordinate sorted files. (Default: " << DEFAULT_POOL_SIZE << ")" << std::endl;
    std::cerr << "\t\t--poolSkipClip: Skip clipping reads to free of usable records when the" << std::endl;
    std::cerr << "\t\t                poolSize is hit. The default action is to just clip the" << std::endl;
    std::cerr << "\t\t                first read in a pair to free up the record." << std::endl;
    std::cerr << std::endl;
}


int ClipOverlap::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    myStoreOrig = "";

    bool readName = false;
    int poolSize = DEFAULT_POOL_SIZE;
    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_STRINGPARAMETER("storeOrig", &myStoreOrig)
        LONG_PARAMETER("readName", &readName)
        LONG_PARAMETER ("stats", &myStats)
        LONG_PARAMETER ("clipsOnly", &myClipsOnly)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PARAMETER_GROUP("Coordinate Clipping Optional Parameters")
        LONG_INTPARAMETER("poolSize", &poolSize)
        LONG_PARAMETER("poolSkipClip", &myPoolSkipClip)
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

    if(readName)
    {
        return(clipSortedByReadName(samIn, samOut));
    }

    myPool.setMaxAllocatedRecs(poolSize);

    return(clipSortedByCoord(samIn, samOut));
}


int ClipOverlap::clipSortedByReadName(SamFile& samIn, SamFile& samOut)
{
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
    uint16_t prevFlag = 0;
    if((samRecord == NULL) || (tmpRecord == NULL))
    {
        std::cerr << "Failed to allocate a SamRecord, so exit.\n";
        return(SamStatus::FAIL_MEM);
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
        if(strcmp(samRecord->getReadName(), prevSamRecord->getReadName()) == 0)
        {
            // Determine if the reads are mapped.
            flag = samRecord->getFlag();
            prevFlag = prevSamRecord->getFlag();
            
            // Read name match, so check if both reads are mapped and there
            // is an overlap, this one starts between the previous record's
            // start & end.
            bool clipped = false;
            if((SamFlag::isMapped(flag) && SamFlag::isMapped(prevFlag)) && 
               (prevSamRecord->getReferenceID() == samRecord->getReferenceID()) && 
               (prevSamRecord->get0BasedPosition() != -1) && (samRecord->get0BasedPosition() != -1))
            {
                // Determine which read starts first.
                if(prevSamRecord->get0BasedPosition() <= 
                   samRecord->get0BasedPosition())
                {
                    // The previous read starts at or before the current one.
                    clipped = clip(*prevSamRecord, *samRecord);
                }
                else
                {
                    // The current read starts before the previous one.
                    clipped = clip(*samRecord, *prevSamRecord);
                }
            }

            // Found a read pair, so write both records.
            if(!myClipsOnly || clipped)
            {
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
            }
            // Setup for the next read with no previous.
            tmpRecord = prevSamRecord;
            prevSamRecord = NULL;
            
        }
        else
        {
            // Read name does not match, so write the previous record
            // if we are writing all records.
            if(!myClipsOnly)
            {
                if(!samOut.WriteRecord(samHeader, *prevSamRecord))
                {
                    // Failed to write a record.
                    fprintf(stderr, "%s\n", samOut.GetStatusMessage());
                    returnStatus = samOut.GetStatus();
                }
            }
            // Store this record as the previous.
            tmpRecord = prevSamRecord;
            prevSamRecord = samRecord;
            samRecord = tmpRecord;
            tmpRecord = NULL;
        }
    }

    // Write the previous record if there is one.
    if((prevSamRecord != NULL) && !myClipsOnly)
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

    // Print Stats if applicable
    if(myStats)
    {
        std::cerr << "Overlap Statistics:" << std::endl;
        std::cerr << "Number of overlapping pairs: "
                  << myOverlaps.NumDataValues() << std::endl
                  << "Average # Reference Bases Overlapped: "
                  << myOverlaps.Mean() << std::endl
                  << "Variance of Reference Bases overlapped: "
                  << myOverlaps.Variance() << std::endl
                  << "Number of times the forward strand was clipped: "
                  << myNumForwardClips << std::endl
                  << "Number of times the reverse strand was clipped: "
                  << myNumReverseClips << std::endl
                  << "Number of times orientation causes additional clipping: "
                  << myNumOrientationClips << std::endl;
    }

    if(returnStatus == SamStatus::SUCCESS)
    {
        std::cerr << "Completed ClipOverlap Successfully.\n";
    }
    else
    {
        std::cerr << "Failed to complete ClipOverlap.\n";
    }
    return(returnStatus);
}


int ClipOverlap::clipSortedByCoord(SamFile& samIn, SamFile& samOut)
{
    SamCoordOutput outputBuffer(myPool);
    SamFileHeader samHeader;
    SamRecord* recordPtr;
    MateMapByCoord mateMap;
    int32_t pos = 0;

    myNumMateFailures = 0;
    myNumPoolFail = 0;

    // Read/write the sam header.
    samIn.ReadHeader(samHeader);
    samOut.WriteHeader(samHeader);

    // Set the output file in the output buffer.
    outputBuffer.setOutputFile(&samOut, &samHeader);

    // Set returnStatus to success.  It will be changed
    // to the failure reason if any of the writes fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;

    // Read the sam records.
    while(returnStatus == SamStatus::SUCCESS)
    {
        recordPtr = myPool.getRecord();
        if(recordPtr == NULL)
        {
            // Failed to allocate a new record, so try to free up records.
            if(!forceRecordFlush(mateMap, outputBuffer))
            {
                std::cerr << "Failed to flush the output buffer.\n";
                returnStatus = SamStatus::FAIL_IO;
                continue;
            }
            // Try to get a new record, one should have been cleared.
            recordPtr = myPool.getRecord();
            if(recordPtr == NULL)
            {
                std::cerr << "Failed to allocate any records.\n";
                returnStatus = SamStatus::FAIL_MEM;
                continue;
            }
        }

        if(!samIn.ReadRecord(samHeader, *recordPtr))
        {
            // Nothing to process, so continue.
            returnStatus = samIn.GetStatus();
            continue;
        }

        // Process this record, get the position before any potential clipping.
        // Chromosome ID will not be affected by clipping.
        pos = recordPtr->get0BasedPosition();
        handleCoordRead(*recordPtr, mateMap, outputBuffer);

        // Read a new record, cleanup/flush based on this read.
        if(!coordFlush(recordPtr->getReferenceID(),
                       pos,
                       mateMap, outputBuffer))
        {
            returnStatus = SamStatus::FAIL_IO;
        }
    }


    // Flush the rest of the unpaired reads and the output buffer.
    if(!coordFlush(-1, -1, mateMap, outputBuffer))
    {
        returnStatus = SamStatus::FAIL_IO;
    }

    // Output any mate errors.
    if(myNumMateFailures != 0)
    {
        std::cerr << "Failed to find expected overlapping mates for " 
                  << myNumMateFailures << " records." << std::endl;
    }
    if(returnStatus == SamStatus::NO_MORE_RECS)
    {
        returnStatus = SamStatus::SUCCESS;
    }
    
    if(myNumPoolFail != 0)
    {
        // Had to skip clipping some records due to running out of
        // memory and not being able to wait for the mate.
        if(myPoolSkipClip)
        {
            std::cerr << "Due to hitting the max record poolSize, had to skip clipping " 
                      << myNumPoolFail << " records." << std::endl;
        }
        else
        {
            std::cerr << "Due to hitting the max record poolSize, had to default clip " 
                      << myNumPoolFail << " records." << std::endl;
        }
        returnStatus = SamStatus::NO_MORE_RECS;
    }

    // Print Stats if applicable
    if(myStats)
    {
        std::cerr << "Overlap Statistics:" << std::endl;
        std::cerr << "Number of overlapping pairs: "
                  << myOverlaps.NumDataValues() << std::endl
                  << "Average # Reference Bases Overlapped: "
                  << myOverlaps.Mean() << std::endl
                  << "Variance of Reference Bases overlapped: "
                  << myOverlaps.Variance() << std::endl
                  << "Number of times the forward strand was clipped: "
                  << myNumForwardClips << std::endl
                  << "Number of times the reverse strand was clipped: "
                  << myNumReverseClips << std::endl
                  << "Number of times orientation causes additional clipping: "
                  << myNumOrientationClips << std::endl;
    }
    
    if(returnStatus == SamStatus::SUCCESS)
    {
        std::cerr << "Completed ClipOverlap Successfully.\n";
    }
    else if(returnStatus == SamStatus::NO_MORE_RECS)
    {
        std::cerr << "Completed ClipOverlap.\n";
    }
    else
    {
        std::cerr << "Failed to complete ClipOverlap.\n";
    }
    return(returnStatus);
}



void ClipOverlap::handleCoordRead(SamRecord& record,
                                  MateMapByCoord& mateMap,
                                  SamCoordOutput& outputBuffer)
{
    // Determine whether or not the reads overlap.
    int16_t flag = record.getFlag();
    // Do not clip if:
    //  1) the read is not paired.
    //  2) read and its mate are on different chromosome ids
    //  3) read is unmapped
    //  4) mate is unmapped.
    if(!SamFlag::isPaired(flag) || 
       (record.getMateReferenceID() != record.getReferenceID()) ||
       !SamFlag::isMapped(flag) || !SamFlag::isMateMapped(flag))
    {
        // No clipping is necessary
        if(!myClipsOnly)
        {
            // Just write it to the output buffer.
            outputBuffer.add(&record);
        }
        else
        {
            // Not writing it, so release the record.
            myPool.releaseRecord(&record);
        }
    }
    else
    {
        // Same chromosome and both reads are mapped
        // Check which read starts first.
        int32_t readStart = record.get0BasedPosition();
        int32_t mateStart = record.get0BasedMatePosition();

        // If either position is unknown (-1), then just add
        // the record to the output buffer because no clipping
        // needs to be done.
        if((readStart == -1) || (mateStart == -1))
        {
            if(!myClipsOnly)
            {
                outputBuffer.add(&record);
            }
            else
            {
                // Not writing it, so release the record.
                myPool.releaseRecord(&record);
            }
            return;
        }

        if(readStart < mateStart)
        {
            // This is the first read in the pair.
            // Check to see if there is an overlap.
            int32_t readEnd = record.get0BasedAlignmentEnd();
            if(readEnd < mateStart)
            {
                // This read finishes before the mate starts so there is no overlap.
                // If this read is the reverse and the other read is the forward
                // strand, then they completely passed each other and both should
                // be clipped.
                if(SamFlag::isReverse(flag) && !SamFlag::isMateReverse(flag))
                {
                    // Clip both.
                    clipEntire(record);
                    // Increment orientation clip counter on the first read
                    // in the pair.
                    if(myStats)
                    {
                        ++myNumOrientationClips;
                    }
                    // The whole read was clipped, so just write it to the output buffer.
                    outputBuffer.add(&record);
                }
                else if(!myClipsOnly)
                {
                    // No clipping is necessary, so just write it to the output buffer.
                    outputBuffer.add(&record);
                }
                else
                {
                    // Not writing it, so release the record.
                    myPool.releaseRecord(&record);
                }
            }
            else
            {
                // The reads overlap, so store this read so the overlap can be
                // clipped when the mate is found.
                mateMap.add(record);
            }
        }
        else
        {
            // This is the 2nd read in the pair or the reads have the
            // same start position.
            // Check the map for the mate.
            SamRecord* mate = mateMap.getMate(record);
            if(mate == NULL)
            {
                // Did not find the mate. 
                // If the start positions are the same, then just insert this
                // read to the mate map.
                if(readStart == mateStart)
                {
                    // Same start position, but the mate has not yet been read,
                    // so store this read.
                    mateMap.add(record);
                }
                else
                {
                    // The mate for this read has already been written to the
                    // output buffer, so write this read.
                    // Before writing, check to see if this entire read needs to be 
                    // clipped (this read is forward and the mate is reverse).
                    if(!SamFlag::isReverse(flag) && SamFlag::isMateReverse(flag))
                    {
                        // Clip both.
                        // Do not increment the orientation clips since
                        // it was incremented when the mate was clipped.
                        clipEntire(record);
                
                        // Did not find the mate but the whole read was already clipped, 
                        // so just write it to the output buffer.
                        outputBuffer.add(&record);
                    }
                    else if(!myClipsOnly)
                    {
                        // Did not find the mate and no clipping is necessary,
                        // so just write it to the output buffer.
                        outputBuffer.add(&record);
                    }
                    else
                    {
                        // Not writing it, so release the record.
                        myPool.releaseRecord(&record);
                    }
                }
            }
            else
            {
                // Found the mate, so clip the 2 reads and write then to the buffer.
                bool clipped = clip(*mate, record);
                if(clipped || !myClipsOnly)
                {
                    // Write both reads to the output buffer.
                    outputBuffer.add(mate);
                    outputBuffer.add(&record);
                }
                else
                {
                    // Not writing them, so release the records.
                    myPool.releaseRecord(&record);
                    myPool.releaseRecord(mate);
                }
            }
        }
    }
}


bool ClipOverlap::coordFlush(int32_t chromID, int32_t position,
                             MateMapByCoord& mateMap,
                             SamCoordOutput& outputBuffer)
{
    // We will flush the output buffer up to the first record left in the
    // mateMap.  If there are no records left in the mate map, then we
    // flush up to this record's position.

    // Track position to flush up to.
    int32_t flushChrom = chromID;
    int32_t flushPosition = position;

    // The current record's chromosome/position.  Used to determine
    // which records to cleanup from the mateMap.
    uint64_t chromPos = SamHelper::combineChromPos(chromID, position);

    // Cleanup any strangling reads at the beginning of the mate map
    // whose mate was not found at the position specified.
    // Stop after the first read is found whose mate has not yet been reached.
    SamRecord* firstRec = mateMap.first();
    while(firstRec != NULL)
    {
        uint64_t firstMateChromPos = 
            SamHelper::combineChromPos(firstRec->getMateReferenceID(),
                                       firstRec->get0BasedMatePosition());
        if((firstMateChromPos < chromPos) || (chromID == -1))
        {
            // Already past the mate's position, so note this read and
            // write it.
            ++myNumMateFailures;
            if(!myClipsOnly)
            {
                outputBuffer.add(firstRec);
            }
            else
            {
                myPool.releaseRecord(firstRec);
            }
            // Remove this record.
            mateMap.popFirst();
            firstRec = mateMap.first();
        }
        else
        {
            // The first record's mate position has not yet been passed, so
            // stop cleaning up the buffer.
            // We will flush up to the start of this record.
            flushChrom = firstRec->getReferenceID();
            flushPosition = firstRec->get0BasedPosition();
            break;
        }
    }

    ////////////////////////////
    // Flush the output buffer prior to this position.
    return(outputBuffer.flush(flushChrom, flushPosition));
}


bool ClipOverlap::forceRecordFlush(MateMapByCoord& mateMap,
                                   SamCoordOutput& outputBuffer)
{
    // Out of records, so ensure records get flushed.
    static CigarRoller newFirstCigar; // holds updated cigar.

    // The previous standard flush did not free up any records, so pop
    // the first record off of the mate map if there is one and flush.
    SamRecord* firstRec = mateMap.first();
    if(firstRec != NULL)
    {
        // Increment the number of records specially handled due to a pool failure.
        ++myNumPoolFail;

        bool clipped = false;

        // Ran out of records, so can't wait until the mate's position.
        if(!myPoolSkipClip)
        {
            // Rather than skipping this pair, just clip this read so it 
            // can be written and the pair will be clipped.
            clipped = true;
            // If this read is the reverse and the other read is the forward
            // strand, then the entire read should be clipped.
            int16_t flag = firstRec->getFlag();
            if(SamFlag::isReverse(flag) && !SamFlag::isMateReverse(flag))
            {
                if(myStats)
                {
                    ++myNumOrientationClips;

                    // Also increment the times that the reverse is clipped.
                    ++myNumReverseClips;
                    myOverlaps.Push(firstRec->get0BasedAlignmentEnd() -
                                    firstRec->get0BasedMatePosition() + 1);
                }
                clipEntire(*firstRec);
            }
            else
            {
                if(CigarHelper::softClipEndByRefPos(*firstRec, 
                                                    firstRec->get0BasedMatePosition(), 
                                                    newFirstCigar) != CigarHelper::NO_CLIP)
                {
                    // Update the number of clips.
                    if(myStats)
                    {
                        // Clipping the first record's strand.
                        if(SamFlag::isReverse(flag))
                        {
                            ++myNumReverseClips;
                        }
                        else
                        {
                            ++myNumForwardClips;
                        }
                        // Update the overlap length - the difference
                        // between the mate start and this read's end.
                        myOverlaps.Push(firstRec->get0BasedAlignmentEnd() -
                                        firstRec->get0BasedMatePosition() + 1);
                    }
                    // Write the original cigar into the specified tag.
                    if(!myStoreOrig.IsEmpty())
                    {
                        // Write original cigar.
                        firstRec->addTag(myStoreOrig, 'Z', firstRec->getCigar());
                    }
                    firstRec->setCigar(newFirstCigar);
                }
            }
        }

        // Add the record to the output buffer.
        if(clipped || !myClipsOnly)
        {
            outputBuffer.add(firstRec);
        }
        else
        {
            myPool.releaseRecord(firstRec);
        }
        // Remove this record.
        mateMap.popFirst();
        
        // Flush up to & including this record's position.
        return(outputBuffer.flush(firstRec->getReferenceID(), 
                                  firstRec->get0BasedPosition()));
    }
    else
    {
        // There was nothing in the mateMap, so flush everything from the output buffer.
        return(outputBuffer.flush(-1, -1));
    }
}


bool ClipOverlap::clip(SamRecord& firstRecord, SamRecord& secondRecord)
{
    static CigarRoller newFirstCigar; // holds updated cigar.
    static CigarRoller newSecondCigar; // holds updated cigar.

    // Track whether or not a clip occurs.
    bool clipped = false;
    
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
        // overlap
        clipped = true;
        // If we are keeping stats, update the stats with the overlap length.
        if(myStats)
        {
            // Both firstEnd & secondStart are 0-based inclusive.
            // Check to see which one ends last.
            int32_t secondEnd = secondRecord.get0BasedAlignmentEnd();
            int32_t end = firstEnd;
            if(firstEnd > secondEnd)
            {
                // The first read ends after the 2nd one, so the
                // overlap ends at the end of the 2nd read.
                end = secondEnd;
            }
            // If firstEnd = 5 & secondStart = 3, they overlap at
            // positions 3, 4, & 5.  So 5-3+1=3 positions.
            myOverlaps.Push(end - secondStart + 1);
        }

        // Determine which record will get clipped by determining
        // which record has a lower base quality in the overlapping region.
        // First check clipping the first region.
        // Determine the clipping on the 1st record at the start of the 2nd.
        int32_t firstClipPos = 
            CigarHelper::softClipEndByRefPos(firstRecord, secondStart, 
                                             newFirstCigar);
        // Loop through counting the quality of the clipped bases.
        // They run from the firstClip to the length of the read.
        double firstQualAvg = getAvgQual(firstRecord, firstClipPos, 
                                         firstRecord.getReadLength()-1);

        int32_t newPos = 0;
        int32_t secondClipPos =
            CigarHelper::softClipBeginByRefPos(secondRecord, firstEnd,
                                               newSecondCigar, newPos);
        // Loop through counting the quality of the clipped bases.
        // They run from the beginning until the secondClip(included).
        double secondQualAvg = getAvgQual(secondRecord, 0, secondClipPos);
        
        // Check to see whether the 1st record or the 2nd one should be clipped
        // based on which has the lower quality.
        if(firstQualAvg <= secondQualAvg)
        {
            // First read has lower or equal quality, so clip that.
            if(myStats)
            {
                if(SamFlag::isReverse(firstFlag))
                {
                    ++myNumReverseClips;
                }
                else
                {
                    ++myNumForwardClips;
                }
            }


            // Check to see if the entire read should be clipped by
            // checking forward/reverse.
            if(SamFlag::isReverse(firstFlag) && !SamFlag::isReverse(secondFlag))
            {
                // first record is reverse and 2nd is forward, so clip
                // those extending ends which means clipping all of the
                // reverse(first) strand since it's overlap has lower quality.
                if(myStats)
                {
                    // Increment that we are doing an orientation clip.
                    ++myNumOrientationClips;
                }

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
            if(myStats)
            {
                if(SamFlag::isReverse(secondFlag))
                {
                    ++myNumReverseClips;
                }
                else
                {
                    ++myNumForwardClips;
                }
            }
            // Check to see if the entire read should be clipped by
            // checking forward/reverse.
            if(SamFlag::isReverse(firstFlag) && !SamFlag::isReverse(secondFlag))
            {
                // first record is reverse and 2nd is forward, so clip
                // those extending ends which means clipping all of the
                // forward(second) strand since it's overlap has lower quality.
                if(myStats)
                {
                    // Increment that we are doing an orientation clip.
                    ++myNumOrientationClips;
                }
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
                    secondRecord.set0BasedMatePosition(newPos);
                    firstRecord.set0BasedPosition(newPos);
                    secondRecord.set0BasedMatePosition(newPos);
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
                firstRecord.set0BasedMatePosition(newPos);
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
            clipped = true;
            // Complete Overlap
            if(myStats)
            {
                ++myNumOrientationClips;
            }
        }
    }
    return(clipped);
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


double ClipOverlap::getAvgQual(SamRecord& record, 
                               int32_t startPos, int32_t endPos)
{
    int32_t qualSum = 0;
    const char* quality = record.getQuality();
    int numVals = 0;

    // Check for invalid start position.
    if((startPos < 0) || ((uint32_t)startPos > strlen(quality)))
    {
        // Invalid start position, just return 0.
        return(0);
    }

    for(int i = startPos; i <= endPos; i++)
    {
        // Check for the null terminator at the end of the quality string.
        if(quality[i] == 0)
        {
            // Qual is shorter than the read, so break.
            break;
        }
        qualSum += quality[i];
        // increment the number of values.
        ++numVals;
    }
    if(numVals != 0)
    {
        return(qualSum/(double)numVals);
    }
    else
    {
        return(0);
    }
}    
