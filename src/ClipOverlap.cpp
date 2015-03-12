/*
 *  Copyright (C) 2011-2013  Regents of the University of Michigan,
 *                           Yancy Lo
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
#include "OverlapClipLowerBaseQual.h"

ClipOverlap::ClipOverlap()
    : BamExecutable(),
      myOverlapHandler(NULL),
      myPool(),
      myOverlapsOnly(false),
      myIntExcludeFlags(0),
      myNumMateFailures(0),
      myNumPoolFail(0),
      myNumPoolFailNoHandle(0),
      myNumPoolFailHandled(0),
      myNumOutOfOrder(0),
      myPoolSkipOverlap(false)
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
    std::cerr << "\t./bam clipOverlap --in <inputFile> --out <outputFile> [--storeOrig <tag>] [--readName] [--stats] [--overlapsOnly] [--excludeFlags <flag>] [--poolSize <numRecords allowed to allocate>] [--poolSkipOverlap] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in           : the SAM/BAM file to clip overlaping read pairs for" << std::endl;
    std::cerr << "\t\t--out          : the SAM/BAM file to be written" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--storeOrig    : Store the original cigar in the specified tag." << std::endl;
    std::cerr << "\t\t--readName     : Original file is sorted by Read Name instead of coordinate." << std::endl;
    std::cerr << "\t\t--stats        : Print some statistics on the overlaps." << std::endl;
    std::cerr << "\t\t--overlapsOnly : Only output overlapping read pairs" << std::endl;
    std::cerr << "\t\t--excludeFlags : Skip records with any of the specified flags set, default 0x70C" << std::endl;
    std::cerr << "\t\t--unmapped     : Mark records that would be completely clipped as unmapped" << std::endl;
    std::cerr << "\t\t--noeof        : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params       : Print the parameter settings to stderr" << std::endl;
    std::cerr << "\tClipping By Coordinate Optional Parameters:" << std::endl;
    std::cerr << "\t\t--poolSize     : Maximum number of records the program is allowed to allocate" << std::endl;
    std::cerr << "\t\t                 for clipping on Coordinate sorted files. (Default: " << DEFAULT_POOL_SIZE << ")" << std::endl;
    std::cerr << "\t\t--poolSkipClip : Skip clipping reads to free of usable records when the" << std::endl;
    std::cerr << "\t\t                 poolSize is hit. The default action is to just clip the" << std::endl;
    std::cerr << "\t\t                 first read in a pair to free up the record." << std::endl;
    std::cerr << std::endl;
}


int ClipOverlap::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    String storeOrig = "";
    bool readName = false;
    bool stats = false;
    int poolSize = DEFAULT_POOL_SIZE;
    bool unmapped = false;
    bool noeof = false;
    bool params = false;
    String excludeFlags = "0x70C";

    // TODO, cleanup legacy parameters
    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_STRINGPARAMETER("storeOrig", &storeOrig)
        LONG_PARAMETER("readName", &readName)
        LONG_PARAMETER ("stats", &stats)
        LONG_PARAMETER ("overlapsOnly", &myOverlapsOnly)
        LONG_STRINGPARAMETER ("excludeFlags", &excludeFlags)
        LONG_PARAMETER("unmapped", &unmapped)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PARAMETER_GROUP("Coordinate Processing Optional Parameters")
        LONG_INTPARAMETER("poolSize", &poolSize)
        LONG_PARAMETER("poolSkipOverlap", &myPoolSkipOverlap)
        LONG_PHONEHOME(VERSION)
        BEGIN_LEGACY_PARAMETERS()
        LONG_PARAMETER ("clipsOnly", &myOverlapsOnly)
        LONG_PARAMETER("poolSkipClip", &myPoolSkipOverlap)
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

    myOverlapHandler = new OverlapClipLowerBaseQual();
    if(myOverlapHandler == NULL)
    {
        usage();
        inputParameters.Status();
        std::cerr << "Failed to allocate the overlap handler\n";
        return(-1);
    }

    if(unmapped)
    {
        myOverlapHandler->markAsUnmapped();
    }

    // Setup the overlap handler.
    myOverlapHandler->keepStats(stats);
    if(storeOrig.Length() != 0)
    {
        myOverlapHandler->storeOrigCigar(storeOrig);
    }

    myIntExcludeFlags = excludeFlags.AsInteger();

    if(params)
    {
        inputParameters.Status();
    }

    // For each step process the file.
    // Open the files & read/write the sam header.
    SamStatus::Status runStatus = SamStatus::SUCCESS;
    for(int i = 1; i <= myOverlapHandler->numSteps(); i++)
    {
        // Open the file for reading.
        mySamHeader.resetHeader();
        SamFile samIn(inFile, SamFile::READ, &mySamHeader);
        SamFile* samOutPtr = NULL;
        // Check if writing, if so, open the output file.
        if(i == myOverlapHandler->numSteps())
        {
            samOutPtr = new SamFile(outFile, SamFile::WRITE, &mySamHeader);
        }

        if(readName)
        {
            samIn.setSortedValidation(SamFile::QUERY_NAME);
            runStatus = handleSortedByReadName(samIn, samOutPtr);
        }
        else
        {
            // Coordinate sorted, so work with the pools.
            samIn.setSortedValidation(SamFile::COORDINATE);
            myPool.setMaxAllocatedRecs(poolSize);

            // Reset the number of failures
            myNumMateFailures = 0;
            myNumPoolFail = 0;
            myNumPoolFailNoHandle = 0;
            myNumPoolFailHandled = 0;
            myNumOutOfOrder = 0;

            // Run by coordinate
            if(samOutPtr != NULL)
            {
                // Setup the output buffer for writing.
                SamCoordOutput outputBuffer(myPool);
                outputBuffer.setOutputFile(samOutPtr, &mySamHeader);
                runStatus = handleSortedByCoord(samIn, &outputBuffer);

                // Cleanup the output buffer.
                if(!outputBuffer.flushAll())
                {
                    std::cerr << "ERROR: Failed to flush the output buffer\n";
                    runStatus = SamStatus::FAIL_IO;
                }
            }
            else
            {
                runStatus = handleSortedByCoord(samIn, NULL);
            }
        }

        if(runStatus != SamStatus::SUCCESS)
        {
            break;
        }
        // Close the input file, it will be reopened if there are 
        // multiple steps.
        samIn.Close();
        if(samOutPtr != NULL)
        {
            samOutPtr->Close();
            delete samOutPtr;
            samOutPtr = NULL;
        }
    }

    // Done processing.
    // Print Stats
    myOverlapHandler->printStats();

    if(myNumMateFailures != 0)
    {
        std::cerr << "WARNING: did not find expected overlapping mates for "
                  << myNumMateFailures << " records." << std::endl;
    }
    if(myNumPoolFail != 0)
    {
        // Had to skip clipping some records due to running out of
        // memory and not being able to wait for the mate.
        std::cerr << "WARNING: " << myNumPoolFail 
                  << " record pool failures\n";
        if(myNumPoolFailNoHandle != 0)
        {
            std::cerr << "Due to hitting the max record poolSize, skipped handling " 
                      << myNumPoolFailNoHandle << " records." << std::endl;
        }
        if(myNumPoolFailHandled != 0)
        {
            std::cerr << "Due to hitting the max record poolSize, default handled " 
                      << myNumPoolFailHandled << " records." << std::endl;
        }
        if(myNumOutOfOrder != 0)
        {
            std::cerr << "WARNING: Resulting File out of Order by " 
                      << myNumOutOfOrder << " records.\n";
        }
    }

    if(runStatus == SamStatus::SUCCESS)
    {
        if(myNumPoolFail == 0)
        {
            std::cerr << "Completed ClipOverlap Successfully.\n";
        }
        else
        {
            runStatus = SamStatus::NO_MORE_RECS;
            std::cerr << "Completed ClipOverlap with WARNINGS.\n";
        }
    }
    else
    {
        std::cerr << "Failed to complete ClipOverlap.\n";
    }
    return(runStatus);
}


SamStatus::Status ClipOverlap::handleSortedByReadName(SamFile& samIn, 
                                                      SamFile* samOutPtr)
{
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
        return(SamStatus::FAIL_MEM);
    }

    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(mySamHeader, *samRecord))
    {
        int16_t flag = samRecord->getFlag();
        if((flag & myIntExcludeFlags) != 0)
        {
            // This read should not be checked for overlaps.

            // Check if there is a previous SamRecord.
            if(prevSamRecord != NULL)
            {
                // There is a previous record.
                // If it has a different read name, write it.
                if(strcmp(samRecord->getReadName(), 
                          prevSamRecord->getReadName()) != 0)
                {
                    // Different read name, so write the previous record.
                    if((samOutPtr != NULL) && !myOverlapsOnly)
                    {
                        if(!samOutPtr->WriteRecord(mySamHeader, *prevSamRecord))
                        {
                            // Failed to write a record.
                            fprintf(stderr, "%s\n", samOutPtr->GetStatusMessage());
                            returnStatus = samOutPtr->GetStatus();
                        }
                    }
                    // Clear the previous record info.
                    tmpRecord = prevSamRecord;
                    prevSamRecord = NULL;
                } 
                // If it has the same read name, leave it in case there is another read with that name
            }
            // This record is not being checked for overlaps, so just write it and continue
            if((samOutPtr != NULL) && !myOverlapsOnly)
            {
                if(!samOutPtr->WriteRecord(mySamHeader, *samRecord))
                {
                    // Failed to write a record.
                    fprintf(stderr, "%s\n", samOutPtr->GetStatusMessage());
                    returnStatus = samOutPtr->GetStatus();
                }
            }
            continue;
        }

        if(prevSamRecord == NULL)
        {
            // Nothing to compare this record to, so set this
            // record to the previous, and the next record.
            prevSamRecord = samRecord;
            samRecord = tmpRecord;
            tmpRecord = NULL;
            continue;
        }

        // Check if the read name matches the previous read name.
        if(strcmp(samRecord->getReadName(), 
                  prevSamRecord->getReadName()) == 0)
        {
            bool overlap = false;
            // Same Read Name, so check clipping.
            OverlapHandler::OverlapInfo prevClipInfo = 
                myOverlapHandler->getOverlapInfo(*prevSamRecord);
            OverlapHandler::OverlapInfo curClipInfo = 
                myOverlapHandler->getOverlapInfo(*samRecord);
            
            // If either indicate a complete clipping, clip both.
            if((prevClipInfo == OverlapHandler::NO_OVERLAP_WRONG_ORIENT) ||
               (curClipInfo == OverlapHandler::NO_OVERLAP_WRONG_ORIENT))
            {
                overlap = true;
                myOverlapHandler->handleNoOverlapWrongOrientation(*prevSamRecord);
                // Don't update stats since this is the 2nd in the pair
                myOverlapHandler->handleNoOverlapWrongOrientation(*samRecord, 
                                                                  false);
            }
            else if((prevClipInfo == OverlapHandler::OVERLAP) ||
                    (prevClipInfo == OverlapHandler::SAME_START))
            {
                // The previous read starts at or before the current one.
                overlap = true;
                myOverlapHandler->handleOverlapPair(*prevSamRecord,
                                                    *samRecord);
            }
            else if(curClipInfo == OverlapHandler::OVERLAP)
            {
                // The current read starts before the previous one.
                overlap = true;
                myOverlapHandler->handleOverlapPair(*samRecord,
                                                    *prevSamRecord);
            }
            
            // Found a read pair, so write both records if: 
            //   1) output file is specified
            //   AND
            //     2a) all records should be written
            //     OR
            //     2b) the pair overlaps
            if((samOutPtr != NULL) && (!myOverlapsOnly || overlap))
            {
                if(!samOutPtr->WriteRecord(mySamHeader, *prevSamRecord))
                {
                    // Failed to write a record.
                    fprintf(stderr, "%s\n", samOutPtr->GetStatusMessage());
                    returnStatus = samOutPtr->GetStatus();
                }
                if(!samOutPtr->WriteRecord(mySamHeader, *samRecord))
                {
                    // Failed to write a record.
                    fprintf(stderr, "%s\n", samOutPtr->GetStatusMessage());
                    returnStatus = samOutPtr->GetStatus();
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
            if((samOutPtr != NULL) && !myOverlapsOnly)
            {
                if(!samOutPtr->WriteRecord(mySamHeader, *prevSamRecord))
                {
                    // Failed to write a record.
                    fprintf(stderr, "%s\n", samOutPtr->GetStatusMessage());
                    returnStatus = samOutPtr->GetStatus();
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
    if((samOutPtr != NULL) && (prevSamRecord != NULL) && !myOverlapsOnly)
    {
        if(!samOutPtr->WriteRecord(mySamHeader, *prevSamRecord))
        {
            // Failed to write a record.
            fprintf(stderr, "%s\n", samOutPtr->GetStatusMessage());
            returnStatus = samOutPtr->GetStatus();
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

    if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        return(samIn.GetStatus());
    }
    return(returnStatus);
}


SamStatus::Status ClipOverlap::handleSortedByCoord(SamFile& samIn, 
                                                   SamCoordOutput* outputBufferPtr)
{
    MateMapByCoord mateMap;

    // Get record & track success/fail
    SamStatus::Status returnStatus = SamStatus::SUCCESS;
    
    OverlapHandler::OverlapInfo overlapInfo = OverlapHandler::UNKNOWN_OVERLAP;
    SamRecord* matePtr = NULL;

    // Get the first read.
    SamRecord* recordPtr = myPool.getRecord();
    if(recordPtr == NULL)
    {
        // No records in the pool, can't process.
        std::cerr <<
            "ERROR: No records in the pool, can't process by coordinate.\n";
        return(SamStatus::FAIL_MEM);
    }

    // Track the original chrom/position of the record being processed
    // for flushing.
    int32_t chrom = -1;
    int32_t position = -1;

    bool overlap = false;
    while(returnStatus == SamStatus::SUCCESS)
    {
        // Get the next Record.
        returnStatus = readCoordRecord(samIn, &recordPtr, 
                                       mateMap, outputBufferPtr);
        if(returnStatus != SamStatus::SUCCESS)
        {
            break;
        }

        chrom = recordPtr->getReferenceID();
        position = recordPtr->get0BasedPosition();

        // Cleanup the mate map based on the newly read record.
        cleanupMateMap(mateMap, outputBufferPtr, chrom, position);

        // Check the read for overlaps.
        overlapInfo = myOverlapHandler->getOverlapInfo(*recordPtr, myIntExcludeFlags);
        // Handle the types of overlaps.
        switch(overlapInfo)
        {
            case OverlapHandler::OVERLAP:
                overlap = true;
                // 1st read, so store it in the mate map.
                mateMap.add(*recordPtr);
                // Clear the pointer so a new one is used next time.
                recordPtr = NULL;
                break;
            case OverlapHandler::NO_OVERLAP_WRONG_ORIENT:
                overlap = true;
                myOverlapHandler->handleNoOverlapWrongOrientation(*recordPtr);
                break;
            case OverlapHandler::SAME_START:
                overlap = true;
                // First check the mate map for the mate.
                matePtr = mateMap.getMate(*recordPtr);
                if(matePtr != NULL)
                {
                    // Mate was found, so handle the overlap.
                    myOverlapHandler->handleOverlapPair(*matePtr, *recordPtr);
                }
                else
                {
                    // Mate not found, so store this one.
                    mateMap.add(*recordPtr);
                    // Clear the pointer so a new one is used next time.
                    recordPtr = NULL;
                }
                break;
            case OverlapHandler::UNKNOWN_OVERLAP:
                matePtr = mateMap.getMate(*recordPtr);
                if(matePtr != NULL)
                {
                    // Mate was found, there is an overlap.
                    overlap = true;
                    myOverlapHandler->handleOverlapPair(*matePtr, *recordPtr);
                }
                else
                {
                    // No overlap if mate not found.
                    overlap = false;
                }
                break;
            case OverlapHandler::UNKNOWN_OVERLAP_WRONG_ORIENT:
                overlap = true;
                matePtr = mateMap.getMate(*recordPtr);
                if(matePtr != NULL)
                {
                    // Mate was found, there is an overlap..
                    myOverlapHandler->handleOverlapPair(*matePtr, *recordPtr);
                }
                else
                {
                    // Mate not found, so handle wrong orientation.
                    // Don't update stats since this is the 2nd in the pair
                    myOverlapHandler->handleNoOverlapWrongOrientation(*recordPtr,
                                                                      false);
                }
                break;
            case OverlapHandler::NO_OVERLAP:
            default:
                // No overlap
                overlap = false;
                break;
        }

        // Handle writing the record if necessary.
        if(outputBufferPtr != NULL)
        {
            // Writing step.
            if((matePtr != NULL) && (!myOverlapsOnly || overlap))
            {
                // Add the mate to the output buffer.
                outputBufferPtr->add(matePtr);
                matePtr = NULL;
            }
            if((recordPtr != NULL) && (!myOverlapsOnly || overlap))
            {
                // Add this record to the output buffer.
                outputBufferPtr->add(recordPtr);
                recordPtr = NULL;
            }

            // Flush the output buffer
            if(!flushOutputBuffer(mateMap, *outputBufferPtr,
                                  chrom, position))
            {
                std::cerr << "ERROR: Failed to flush the output buffer\n";
                returnStatus = SamStatus::FAIL_IO;
            }
        }
        else
        {
            // Not writing.
            // Release the mate if it is set.
            if(matePtr != NULL)
            {
                myPool.releaseRecord(matePtr);
            }
        }
    }

    // Done with the file, cleanup the mate map.
    // The calling method will cleanup the output buffer.
    cleanupMateMap(mateMap, outputBufferPtr);

    if(returnStatus != SamStatus::NO_MORE_RECS)
    {
        // Failure.
        std::cerr << "ERROR reading file, exiting.\n";
        return(returnStatus);
    }
    return(SamStatus::SUCCESS);
}


///////////////////////////////////////////////////////////////////
// Methods to handle Coordinate Specific Clipping Operations.

SamStatus::Status ClipOverlap::readCoordRecord(SamFile& samIn,
                                               SamRecord** recordPtr, 
                                               MateMapByCoord& mateMap,
                                               SamCoordOutput* outputBufferPtr)
{
    // Null pointer, so get a new pointer.
    if(*recordPtr == NULL)
    {
        *recordPtr = myPool.getRecord();
        if(*recordPtr == NULL)
        {
            // Failed to allocate a new record.
            // Try to free up records from the mate map
            if(!forceRecordFlush(mateMap, outputBufferPtr))
            {
                std::cerr << "Failed to flush the output buffer.\n";
                return(SamStatus::FAIL_IO);
            }
            // Try to get a new record, one should have been cleared.
            *recordPtr = myPool.getRecord();
            if(*recordPtr == NULL)
            {
                std::cerr << "Failed to allocate any records.\n";
                return(SamStatus::FAIL_MEM);
            }
        }
    }

    // RecordPtr is set.
    if(!samIn.ReadRecord(mySamHeader, **recordPtr))
    {
        // Nothing to process, so return.
        return(samIn.GetStatus());
    }
    return(SamStatus::SUCCESS);
}


///////////////////////////////////////////////////////////////////
// Methods to handle flushing records from the mate map and/or
// the output buffer.

///////////////////////////////////////////////////////////////////
// Methods to handle flushing records from the mate map and/or
// the output buffer.
bool ClipOverlap::forceRecordFlush(MateMapByCoord& mateMap,
                                   SamCoordOutput* outputBufferPtr)
{
    // The previous standard flush did not free up any records, so pop
    // the first record off of the mate map if there is one and process
    // it without waiting for its mate.
    SamRecord* firstRec = mateMap.first();
    int32_t flushChrom = -1;
    int32_t flushPos = -1;
    bool updated = false;

    // Increment number of pool failures.
    ++myNumPoolFail;

    if(firstRec != NULL)
    {
        // Flush up to & including this record's position.
        flushChrom = firstRec->getReferenceID();
        flushPos = firstRec->get0BasedPosition();

        // Remove this record.
        mateMap.popFirst();

        // Ran out of records, so can't wait until the mate's position.
        if(!myPoolSkipOverlap)
        {
            // Handle the single entry of the pair.
            updated = 
                myOverlapHandler->handleOverlapWithoutMate(*firstRec);
        }

        if(!updated)
        {
            ++myNumPoolFailNoHandle;
        }
        else
        {
            ++myNumPoolFailHandled;
        }

        // Add the record to the output buffer.
        if((outputBufferPtr != NULL) && (updated || !myOverlapsOnly))
        {
            outputBufferPtr->add(firstRec);
        }
        else
        {
            myPool.releaseRecord(firstRec);
        }
    }
    else
    {
        ++myNumOutOfOrder;
    }

    if(outputBufferPtr != NULL)
    {
        return(outputBufferPtr->flush(flushChrom, flushPos));
    }
    // No output buffer, return true.
    return(true);
}


bool ClipOverlap::flushOutputBuffer(MateMapByCoord& mateMap,
                                    SamCoordOutput& outputBuffer,
                                    int32_t prevChrom,
                                    int32_t prevPos)
{
    // We will flush the output buffer up to the first record left in the
    // mateMap.  If there are no records left in the mate map, then we
    // flush everything up to the previous chrom/pos that was processed since
    // any new records will have a higher coordinate.
    SamRecord* firstRec = mateMap.first();
    if(firstRec != NULL)
    {
        return(outputBuffer.flush(firstRec->getReferenceID(), 
                                  firstRec->get0BasedPosition()));
    }
    // Otherwise, flush based on the previous 
    return(outputBuffer.flush(prevChrom, prevPos));
}

void ClipOverlap::cleanupMateMap(MateMapByCoord& mateMap,
                                 SamCoordOutput* outputBufferPtr,
                                 int32_t chrom, int32_t position)
{
    // Cleanup any reads in the mateMap whose mates are prior to the position
    // currently being processed in the file.  It means the mate was not found 
    // as expected.  Stop cleaning up once one is found that is not passed.
    uint64_t chromPos = 0;
    if((chrom != -1) && (position != -1))

    {
        chromPos = SamHelper::combineChromPos(chrom, position);
    }
    else
    {
        chrom = -1;
    }
    
    // Stop after the first read is found whose mate has not yet been reached.
    SamRecord* firstRec = mateMap.first();
    while(firstRec != NULL)
    {
        uint64_t firstMateChromPos = 
            SamHelper::combineChromPos(firstRec->getMateReferenceID(),
                                       firstRec->get0BasedMatePosition());
        if((firstMateChromPos < chromPos) || (chrom == -1))
        {
            // Already past the mate's position, so note this read and
            // write it.
            ++myNumMateFailures;
            if((outputBufferPtr != NULL) && !myOverlapsOnly)
            {
                outputBufferPtr->add(firstRec);
            }
            else
            {
                myPool.releaseRecord(firstRec);
            }
            // Remove this record.
            mateMap.popFirst();
            // Get the next record to check.
            firstRec = mateMap.first();
        }
        else
        {
            // The first record's mate position has not yet been passed, so
            // stop cleaning up the buffer.
            break;
        }
    }
}
