/*
 *  Copyright (C) 2010-2015  Regents of the University of Michigan
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

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include "SamFile.h"
#include "Dedup_LowMem.h"
#include "Logger.h"
#include "SamHelper.h"
#include "SamStatus.h"
#include "BgzfFileType.h"

const int Dedup_LowMem::DEFAULT_MIN_QUAL = 15;
const uint32_t Dedup_LowMem::CLIP_OFFSET = 1000;

Dedup_LowMem::~Dedup_LowMem()
{
    // clean up the maps.
    // First free any fragment records.
    myFragmentMap.clear();

    // Free any paired records.
    myPairedMap.clear();

    // Free any entries in the mate map.
    myMateMap.clear();
}

void Dedup_LowMem::dedup_LowMemDescription()
{
    std::cerr << " dedup_LowMem - Mark Duplicates using only a little memory\n";
}


void Dedup_LowMem::description()
{
    dedup_LowMemDescription();
}


void Dedup_LowMem::usage()
{
    std::cerr << "Usage: ./bam dedup_LowMem --in <InputBamFile> --out <OutputBamFile> [--minQual <minPhred>] [--log <logFile>] [--oneChrom] [--rmDups] [--force] [--excludeFlags <flag>] [--verbose] [--noeof] [--params] [--recab] ";
    myRecab.recabSpecificUsageLine();
    std::cerr << std::endl << std::endl;
    std::cerr << "Required parameters :" << std::endl;
    std::cerr << "\t--in <infile>   : Input BAM file name (must be sorted)" << std::endl;
    std::cerr << "\t--out <outfile> : Output BAM file name (same order with original file)" << std::endl;
    std::cerr << "Optional parameters : " << std::endl;
    std::cerr << "\t--minQual <int> : Only add scores over this phred quality when determining a read's quality (default: "
              << DEFAULT_MIN_QUAL << ")" << std::endl;
    std::cerr << "\t--log <logfile> : Log and summary statistics (default: [outfile].log, or stderr if --out starts with '-')" << std::endl;
    std::cerr << "\t--oneChrom      : Treat reads with mates on different chromosomes as single-ended." << std::endl;
    std::cerr << "\t--rmDups        : Remove duplicates (default is to mark duplicates)" << std::endl;
    std::cerr << "\t--force         : Allow an already mark-duplicated BAM file, unmarking any previously marked " << std::endl;
    std::cerr << "\t                  duplicates and apply this duplicate marking logic.  Default is to throw errors" << std::endl;
    std::cerr << "\t                  and exit when trying to run on an already mark-duplicated BAM" << std::endl;
    std::cerr << "\t--excludeFlags <flag>    : exclude reads with any of these flags set when determining or marking duplicates" << std::endl;
    std::cerr << "\t                           by default (0x304): exclude unmapped, secondary reads, and QC failures" << std::endl;
    std::cerr << "\t--verbose       : Turn on verbose mode" << std::endl;
    std::cerr << "\t--noeof         : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t--params        : Print the parameter settings" << std::endl;
    std::cerr << "\t--recab         : Recalibrate in addition to dedup_LowMem" << std::endl;
    myRecab.recabSpecificUsage();
    std::cerr << "\n" << std::endl;
}

int Dedup_LowMem::execute(int argc, char** argv) 
{
    /* --------------------------------
     * process the arguments
     * -------------------------------*/
    String inFile, outFile, logFile;
    myDoRecab = false;
    bool removeFlag = false;
    bool verboseFlag = false;
    myForceFlag = false;
    myNumMissingMate = 0;
    myMinQual = DEFAULT_MIN_QUAL;
    String excludeFlags = "0x304";
    uint16_t intExcludeFlags = 0;
    bool noeof = false;
    bool params = false;

    LongParamContainer parameters;
    parameters.addGroup("Required Parameters");
    parameters.addString("in", &inFile);
    parameters.addString("out", &outFile);
    parameters.addGroup("Optional Parameters");
    parameters.addInt("minQual", & myMinQual);
    parameters.addString("log", &logFile);
    parameters.addBool("oneChrom", &myOneChrom);
    parameters.addBool("recab", &myDoRecab);
    parameters.addBool("rmDups", &removeFlag);
    parameters.addBool("force", &myForceFlag);
    parameters.addString("excludeFlags", &excludeFlags);
    parameters.addBool("verbose", &verboseFlag);
    parameters.addBool("noeof", &noeof);
    parameters.addBool("params", &params);
    parameters.addPhoneHome(VERSION);
    myRecab.addRecabSpecificParameters(parameters);

    ParameterList inputParameters;
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            parameters.getLongParameterList()));
    
    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);
    
    // If no eof block is required for a bgzf file, set the bgzf file type to 
    // not look for it.
    if(noeof)
    {
        // Set that the eof block is not required.
        BgzfFileType::setRequireEofBlock(false);
    }

    if(inFile.IsEmpty())
    {
        usage();
        inputParameters.Status();
        std::cerr << "Specify an input file" << std::endl;
        return EXIT_FAILURE;
    }

    if(outFile.IsEmpty())
    {
        usage();
        inputParameters.Status();
        std::cerr << "Specify an output file" << std::endl;
        return EXIT_FAILURE;
    }

    intExcludeFlags = excludeFlags.AsInteger();

    if(myForceFlag && SamFlag::isDuplicate(intExcludeFlags))
    {
        usage();
        inputParameters.Status();
        std::cerr << "Cannot specify --force and Duplicate in the excludeFlags.  Since --force indicates to override"
                  << " previous duplicate setting and the excludeFlags says to skip those, you can't do both.\n";
        return EXIT_FAILURE;
    }

    if(logFile.IsEmpty())
    {
        logFile = outFile + ".log";
    }
    
    if(myDoRecab)
    {
        int status = myRecab.processRecabParam();
        if(status != 0)
        {
            inputParameters.Status();
            return(status);
        }
    }

    if(params)
    {
        inputParameters.Status();
    }
    
    Logger::gLogger = new Logger(logFile.c_str(), verboseFlag);

    /* -------------------------------------------------------------------
     * The arguments are processed.  Prepare the input BAM file,
     * instantiate dedup_LowMem, and construct the read group library map
     * ------------------------------------------------------------------*/

    SamFile samIn;

    samIn.OpenForRead(inFile.c_str());
    // If the file isn't sorted it will throw an exception.
    samIn.setSortedValidation(SamFile::COORDINATE);

    SamFileHeader header;
    samIn.ReadHeader(header);

    buildReadGroupLibraryMap(header);

    lastReference = -1;
    lastCoordinate = -1;

    // for keeping some basic statistics
    uint32_t recordCount = 0;
    uint32_t pairedCount = 0;
    uint32_t properPairCount = 0;
    uint32_t unmappedCount = 0;
    uint32_t reverseCount = 0;
    uint32_t qualCheckFailCount = 0;
    uint32_t secondaryCount = 0;
    uint32_t excludedCount = 0;

    // Now we start reading records
    SamRecord* recordPtr;
    SamStatus::Status returnStatus = SamStatus::SUCCESS;
    while(returnStatus == SamStatus::SUCCESS)
    {
        recordPtr = mySamPool.getRecord();
        if(recordPtr == NULL)
        {
            std::cerr << "Failed to allocate enough records\n";
            return(-1);
        }
        if(!samIn.ReadRecord(header, *recordPtr))
        {
            returnStatus = samIn.GetStatus();
            continue;
        }
        // Take note of properties of this record
        int flag = recordPtr->getFlag();
        if(SamFlag::isPaired(flag))     ++pairedCount;
        if(SamFlag::isProperPair(flag)) ++properPairCount;
        if(SamFlag::isReverse(flag))    ++reverseCount;
        if(SamFlag::isQCFailure(flag))  ++qualCheckFailCount;
        if(SamFlag::isSecondary(flag))  ++secondaryCount;
        if(!SamFlag::isMapped(flag))    ++unmappedCount;

        // put the record in the appropriate maps:
        //   single reads go in myFragmentMap
        //   paired reads go in myPairedMap
        recordCount = samIn.GetCurrentRecordCount();

        // if we have moved to a new position, look back at previous reads for duplicates
        if (hasPositionChanged(*recordPtr))
        {
            cleanupPriorReads(recordPtr);
        }

        // Determine if this read should be checked for duplicates.
        if((!SamFlag::isMapped(flag)) || ((flag & intExcludeFlags) != 0))
        {
            ++excludedCount;

            // No deduping done on this record, but still build the recab table.
            if(myDoRecab)
            {
                myRecab.processReadBuildTable(*recordPtr);
            }
            // Nothing more to do with this record, so
            // release the pointer.
            mySamPool.releaseRecord(recordPtr);
        }
        else
        {
            if(SamFlag::isDuplicate(flag) && !myForceFlag)
            {
                // Error: Marked duplicates, and duplicates aren't excluded.
                Logger::gLogger->error("There are records already duplicate marked.");
                Logger::gLogger->error("Use -f to clear the duplicate flag and start the dedup_LowMem procedure over");
            }

            checkDups(*recordPtr, recordCount);
            mySamPool.releaseRecord(recordPtr);
        }
        // let the user know we're not napping
        if (verboseFlag && (recordCount % 100000 == 0))
        {
            Logger::gLogger->writeLog("recordCount=%u singleKeyMap=%u pairedKeyMap=%u, dictSize=%u", 
                                      recordCount, myFragmentMap.size(), 
                                      myPairedMap.size(), 
                                      myMateMap.size());
        }
    }

    // we're finished reading record so clean up the duplicate search and
    //  close the input file
    cleanupPriorReads(NULL);
    samIn.Close();

    // print some statistics
    Logger::gLogger->writeLog("--------------------------------------------------------------------------");
    Logger::gLogger->writeLog("SUMMARY STATISTICS OF THE READS");
    Logger::gLogger->writeLog("Total number of reads: %u",recordCount);
    Logger::gLogger->writeLog("Total number of paired-end reads: %u",
                              pairedCount);
    Logger::gLogger->writeLog("Total number of properly paired reads: %u",
                              properPairCount);
    Logger::gLogger->writeLog("Total number of unmapped reads: %u",
                              unmappedCount);
    Logger::gLogger->writeLog("Total number of reverse strand mapped reads: %u",
                              reverseCount);
    Logger::gLogger->writeLog("Total number of QC-failed reads: %u",
                              qualCheckFailCount);
    Logger::gLogger->writeLog("Total number of secondary reads: %u",
                              secondaryCount);
    Logger::gLogger->writeLog("Size of singleKeyMap (must be zero): %u",
                              myFragmentMap.size());
    Logger::gLogger->writeLog("Size of pairedKeyMap (must be zero): %u",
                              myPairedMap.size());
    Logger::gLogger->writeLog("Total number of missing mates: %u",
                              myNumMissingMate);
    Logger::gLogger->writeLog("Total number of reads excluded from duplicate checking: %u",
                              excludedCount);
    Logger::gLogger->writeLog("--------------------------------------------------------------------------");
    Logger::gLogger->writeLog("Sorting the indices of %d duplicated records",
                              myDupList.size());

    // sort the indices of duplicate records
    std::sort(myDupList.begin(), myDupList.end(),
              std::less<uint32_t> ());

    // get ready to write the output file by making a second pass
    // through the input file
    samIn.OpenForRead(inFile.c_str());
    samIn.ReadHeader(header);

    SamFile samOut;
    samOut.OpenForWrite(outFile.c_str());
    samOut.WriteHeader(header);

    // If we are recalibrating, output the model information.
    if(myDoRecab)
    {
        myRecab.modelFitPrediction(outFile);
    }

    // an iterator to run through the duplicate indices
    int currentDupIndex = 0;
    bool moreDups = !myDupList.empty();

    // let the user know what we're doing
    Logger::gLogger->writeLog("\nWriting %s", outFile.c_str());

    // count the duplicate records as a check
    uint32_t singleDuplicates(0), pairedDuplicates(0);

    // start reading records and writing them out
    SamRecord record;
    while(samIn.ReadRecord(header, record))
    {
        uint32_t currentIndex = samIn.GetCurrentRecordCount();

        bool foundDup = moreDups &&
            (currentIndex == myDupList[currentDupIndex]);

        // modify the duplicate flag and write out the record,
        // if it's appropriate
        int flag = record.getFlag();
        if (foundDup)
        {   
            // this record is a duplicate, so mark it.
            record.setFlag( flag | 0x400 );
            currentDupIndex++;
            // increment duplicate counters to verify we found them all
            if ( ( ( flag & 0x0001 ) == 0 ) || ( flag & 0x0008 ) )
            { // unpaired or mate unmapped
                singleDuplicates++;
            }
            else
            {
                pairedDuplicates++;
            }
            // recalibrate if necessary.
            if(myDoRecab)
            {
                myRecab.processReadApplyTable(record);
            }

            // write the record if we are not removing duplicates
            if (!removeFlag ) samOut.WriteRecord(header, record);
        }
        else
        {
            if(myForceFlag)
            { 
                // this is not a duplicate we've identified but we want to
                // remove any duplicate marking
                record.setFlag( flag & 0xfffffbff ); // unmark duplicate
            }
            // Not a duplicate, so recalibrate if necessary.
            if(myDoRecab)
            {
                myRecab.processReadApplyTable(record);
            }
            samOut.WriteRecord(header, record);
        }
	
        // Let the user know we're still here
        if (verboseFlag && (currentIndex % 100000 == 0)) {
            Logger::gLogger->writeLog("recordCount=%u", currentIndex);
        }
    }

    // We're done.  Close the files and print triumphant messages.
    samIn.Close();
    samOut.Close();

    Logger::gLogger->writeLog("Successfully %s %u unpaired and %u paired duplicate reads", 
                              removeFlag ? "removed" : "marked" ,
                              singleDuplicates,
                              pairedDuplicates/2);
    Logger::gLogger->writeLog("\nDedup_LowMem complete!");
    return 0;
}

// Now that we've reached coordinate on chromosome reference, look back and
// clean up any previous positions from being tracked.
void Dedup_LowMem::cleanupPriorReads(SamRecord* record)
{
    static DupKey emptyKey;
    static DupKey tempKey2;

    // Set where to stop cleaning out the structures.
    // Initialize to the end of the structures.
    FragmentMap::iterator fragmentFinish = myFragmentMap.end();
    PairedMap::iterator pairedFinish = myPairedMap.end();
    uint64_t mateStopPos = 0;

    // If a record was specified, stop before this record.
    if(record != NULL)
    {
        int32_t reference = record->getReferenceID();
        int32_t coordinate = record->get0BasedPosition();
        tempKey2.cleanupKey(reference, coordinate);
        fragmentFinish = myFragmentMap.lower_bound(tempKey2);
        // Now do the same thing with the paired reads
        PairedKey pairedKey(emptyKey, tempKey2);
        pairedFinish = myPairedMap.lower_bound(pairedKey);
        mateStopPos =
            SamHelper::combineChromPos(reference, 
                                       coordinate);
    }

    // For each key k < fragmentFinish, release the record since we are
    // done with that position and it is not a duplicate.
    for(FragmentMap::iterator iter = myFragmentMap.begin(); 
        iter != fragmentFinish; iter++)
    {
        // If it is not paired, we are done with this record.
        // If it is paired, it will be handled separately.
        if(!iter->second.paired)
        {
            // Unpaired, non-duplicate, so perform any additional handling.
            handleNonDuplicate();
        }
    }
    // Erase the entries from the map.
    if(fragmentFinish != myFragmentMap.begin())
    {
        myFragmentMap.erase(myFragmentMap.begin(), fragmentFinish);
    }

    // Now do the same thing with the paired reads
    for(PairedMap::iterator iter = myPairedMap.begin(); 
        iter != pairedFinish; iter++)
    {
        // These are not duplicates, but we are done with them, 
        // so perform any additional handling.
        handleNonDuplicatePair();
    }
    // Erase the entries.
    if (pairedFinish != myPairedMap.begin())
    {
        myPairedMap.erase(myPairedMap.begin(), pairedFinish);
    }

    // Clean up the Mate map from any reads whose mates were not found.
    // Loop through the mate map and release records prior to this position.
    MateMap::iterator mateIter;
    for(mateIter = myMateMap.begin(); mateIter != myMateMap.end(); mateIter++)
    {
        // stop if a record was specified and we have gone past the mate
        // stop position.  If no record was specified, we want to clean
        // it all out.
        if((record != NULL) && (mateIter->first >= mateStopPos))
        {
            break;
        }
        // Passed the mate, but it was not found.
        handleMissingMate(mateIter->second.key.reference, 
                          mateIter->first >> 32);
    }
    // Erase the entries.
    if(mateIter != myMateMap.begin())
    {
        myMateMap.erase(myMateMap.begin(), mateIter);
    }
    return;
}


// determine whether the record's position is different from the previous record
bool Dedup_LowMem::hasPositionChanged(SamRecord& record)
{
    if (lastReference != record.getReferenceID() || 
        lastCoordinate < record.get0BasedPosition())
    {
        if (lastReference != record.getReferenceID())
        {
            lastReference = record.getReferenceID();
            Logger::gLogger->writeLog("Reading ReferenceID %d\n", lastReference);
        }
        lastCoordinate = record.get0BasedPosition();
        return true;
    }
    return false;
}

// When a record is read, check if it is a duplicate or
// store for future checking.
void Dedup_LowMem::checkDups(SamRecord& record, uint32_t recordCount)
{
    // Only inside this method if the record is mapped.

    // Get the key for this record.
    static DupKey key;
    key.initKey(record, getLibraryID(record));

    int flag = record.getFlag(); 
    bool recordPaired = SamFlag::isPaired(flag) && SamFlag::isMateMapped(flag);
    int sumBaseQual = getBaseQuality(record);

    int32_t chromID = record.getReferenceID();
    int32_t mateChromID = record.getMateReferenceID();

    // If we are one-chrom and the mate is not on the same chromosome, 
    // mark it as not paired.
    if(myOneChrom && (chromID != mateChromID))
    {
        recordPaired = false;
    }
    
    // Look in the fragment map to see if an entry for this key exists.
    FragmentMapInsertReturn ireturn = 
        myFragmentMap.insert(std::make_pair(key, FragData()));

    FragData* fragData = &(ireturn.first->second);

    // Enter the new record in the fragData if (any of the below):
    // 1) there is no previous entry for this key (ireturn.second == true)
    // or
    // 2) the previous entry is not paired 
    //    AND
    //     a) the new record is paired
    //     or
    //     b) the new record has higher quality
    if((ireturn.second == true) ||
       ((fragData->paired == false) && 
        (recordPaired || (sumBaseQual > fragData->sumBaseQual))))
    {
        // Check if this is a new key.
        if(ireturn.second == true)
        {
            // New entry, so build the recalibration table now.
            if(myDoRecab)
            {
                myRecab.processReadBuildTable(record);
            }
        }
        else if(fragData->paired == false)
        {
            // There was a previous record and it is not paired,
            // so mark it as a duplicate. 
            // Duplicate checking/marking for pairs is handled below.
            handleDuplicate(fragData->recordIndex);
        }

        // Store this record for later duplicate checking.
        fragData->sumBaseQual = sumBaseQual;
        fragData->recordIndex = recordCount;
        fragData->paired = recordPaired;
    }
    else
    {
        // Leave the old record in fragData.  
        // If the new record is not paired, handle it as a duplicate.
        if(recordPaired == false)
        {
            // This record is a duplicate, so mark it and release it.
            handleDuplicate(recordCount);
        }
    }

    // Only paired processing is left, so return if not paired.
    if(recordPaired == false)
    {
        // Not paired, no more operations required, so return.
        return;
    }
    
    // This is a paired record, so check for its mate.
    uint64_t readPos = 
        SamHelper::combineChromPos(chromID,
                                   record.get0BasedPosition());
    uint64_t matePos =
        SamHelper::combineChromPos(mateChromID, 
                                   record.get0BasedMatePosition());
    int mateIndex = -1;
    MateData* mateData = NULL;

    // Check to see if the mate is prior to this record.
    if(matePos <= readPos)
    {
        // The mate map is stored by the mate position, so look for this 
        // record's position.
        // The mate should be in the mate map, so find it.
        std::pair<MateMap::iterator,MateMap::iterator> matches =
            myMateMap.equal_range(readPos);
        // Loop through the elements that matched the pos looking for the mate.
        for(MateMap::iterator iter = matches.first; 
            iter != matches.second; iter++)
        {
            if(strcmp((*iter).second.readName.c_str(), 
                      record.getReadName()) == 0)
            {
                // Found the match.
                mateData = &((*iter).second);
                // Update the quality and track the mate record and index.
                sumBaseQual += mateData->sumBaseQual;
                mateIndex = mateData->recordIndex;
                // Remove the entry from the map.
                myMateMap.erase(iter);
                break;
            }
        }
    }
    if(mateData == NULL)
    {
        if(matePos >= readPos)
        {
            // Haven't gotten to the mate yet, so store this record.
            MateMap::iterator mateIter = 
                myMateMap.insert(std::make_pair(matePos, MateData()));
            mateIter->second.sumBaseQual = sumBaseQual;
            mateIter->second.recordIndex = recordCount;
            mateIter->second.key.copy(key);
            mateIter->second.readName = record.getReadName();
        }
        else
        {
            // Passed the mate, but it was not found.
            handleMissingMate(record.getReferenceID(), record.getMateReferenceID());
        }
        return;
    }

    // Make the paired key.
    PairedKey pkey(key, mateData->key);

    // Check to see if this pair is a duplicate.
    PairedMapInsertReturn pairedReturn = 
        myPairedMap.insert(std::make_pair(pkey,PairedData()));
    PairedData* storedPair = &(pairedReturn.first->second);

    // Get the index for "record 1" - the one with the earlier coordinate.
    int record1Index = getFirstIndex(key, recordCount,
                                     mateData->key, mateIndex);

    // Check if we have already found a duplicate pair.
    // If there is no duplicate found, there is nothing more to do.
    if(pairedReturn.second == false)
    {
        // Duplicate found.
        bool keepStored = true;
        if(pairedReturn.first->second.sumBaseQual < sumBaseQual)
        {
            // The new pair has higher quality, so keep that.
            keepStored = false;
        }
        else if(pairedReturn.first->second.sumBaseQual == sumBaseQual)
        {
            // Same quality, so keep the one with the earlier record1Index.
            if(record1Index < storedPair->record1Index)
            {
                // The new pair has an earlier lower coordinate read,
                // so keep that.
                keepStored = false;
            }
        }
        // Check to see which one should be kept by checking qualities.
        if(keepStored)
        {
            // The old pair had higher quality so mark the new pair as a
            // duplicate and release them.
            handleDuplicate(mateIndex);
            handleDuplicate(recordCount);
        }
        else
        {
            // The new pair has higher quality, so keep that.
            // First mark the previous one as duplicates and release them.
            handleDuplicate(storedPair->record1Index);
            handleDuplicate(storedPair->record2Index);
            // Store this pair's information.
            if(record1Index == mateIndex)
            {
                // Mate has a lower coordinate, so make mate
                // record1.
                storedPair->sumBaseQual = sumBaseQual;
                storedPair->record1Index = mateIndex;
                storedPair->record2Index = recordCount;
            }
            else
            {
                // This record has a lower coordinate, so make it
                // record1.
                storedPair->sumBaseQual = sumBaseQual;
                storedPair->record1Index = recordCount;
                storedPair->record2Index = mateIndex;
            }
        }
    }
    else
    {
        // Store this pair's information.
        storedPair->sumBaseQual = sumBaseQual;

        if(record1Index == mateIndex)
        {
            // Mate has a lower coordinate, so make mate
            // record1.
            storedPair->record1Index = mateIndex;
            storedPair->record2Index = recordCount;
        }
        else
        {
            // This record has a lower coordinate, so make it
            // record1.
            storedPair->record1Index = recordCount;
            storedPair->record2Index = mateIndex;
        }
    }
}


// Finds the total base quality of a read 
int Dedup_LowMem::getBaseQuality(SamRecord & record) {
    const char* baseQualities = record.getQuality();
    int readLength = record.getReadLength();
    int quality = 0.;
    if(strcmp(baseQualities, "*") == 0)
    {
        return(0);
    }
    for(int i=0; i < readLength; ++i) {
        int q = static_cast<int>(baseQualities[i])-33;
        if ( q >= myMinQual ) quality += q;
    }
    return quality;
}


// build the read group library map
void Dedup_LowMem::buildReadGroupLibraryMap(SamFileHeader& header) {
    rgidLibMap.clear();
    numLibraries = 0;
    std::map<std::string,uint32_t> libNameMap;
    
    SamHeaderRecord * headerRecord = header.getNextRGRecord();
    while(headerRecord != NULL) {
        std::string ID = headerRecord->getTagValue("ID");
        std::string LB = headerRecord->getTagValue("LB");

        if ( ID.empty() ) {
            std::string headerRecordString;
            headerRecord->appendString(headerRecordString);
            Logger::gLogger->error("Cannot find readGroup ID information in the header line %s", 
                                   headerRecordString.c_str());
        }
        if ( rgidLibMap.find(ID) != rgidLibMap.end() ) {
            Logger::gLogger->error("The readGroup ID %s is not a unique identifier",ID.c_str());
        }

        if ( LB.empty() ) {
            std::string headerRecordString;
            headerRecord->appendString(headerRecordString);
            Logger::gLogger->warning("Cannot find library information in the header line %s. Using empty string for library name",
                                     headerRecordString.c_str());
        }

        if ( libNameMap.find( LB ) != libNameMap.end() ) {
            rgidLibMap[ID] = libNameMap[LB];
        }
        else {
            numLibraries = libNameMap.size()+1;
            libNameMap[LB] = numLibraries;
            rgidLibMap[ID] = numLibraries;
        }
        headerRecord = header.getNextRGRecord();
    }

    if (numLibraries > 0xff) {
        Logger::gLogger->error("More than 255 library names are identified. Dedup_LowMem currently only allows up to 255 library names");
    }
}    

// get the libraryID of a record
uint32_t Dedup_LowMem::getLibraryID(SamRecord& record, bool checkTags) {
    if ( ( checkTags == false ) && ( numLibraries <= 1 ) ) {
        return 0; 
    } else {
        char tag[3];
        char vtype;
        void* value;
        std::string rgID;
        record.resetTagIter();
        while( record.getNextSamTag(tag,vtype,&value) != false ) {
            if ( ( tag[0] == 'R' ) && ( tag[1] == 'G' ) && ( vtype == 'Z' ) ) {
                if ( !rgID.empty() ) {
                    Logger::gLogger->error("Multiple RG tag found in one record. ReadName is %s",record.getReadName());
                }
                else if ( record.isStringType(vtype) ) {
                    String s = (String)*(String*)value;
                    rgID = s.c_str();
                }
                else {
                    Logger::gLogger->error("vtype is not string (Z) for RG tag");
                }
            }
        }
        if ( rgID.empty() ) {
            Logger::gLogger->error("No RG tag is found in read %s",record.getReadName());
            return 0;
        }
        else {
            std::map<std::string,uint32_t>::iterator it = rgidLibMap.find(rgID);
            if ( it != rgidLibMap.end() ) {
                return it->second;
            }
            else {
                Logger::gLogger->warning("RG tag %s does not exist in the header",rgID.c_str());
                return 0; // cannot be reached
            }
        }
    }
}


void Dedup_LowMem::handleMissingMate(int32_t refID, int32_t mateRefID)
{
    static bool firstDifferChrom = true;
    static bool firstSameChrom = true;

    ++myNumMissingMate;

    // Passed the mate, but it was not found.
    if(refID != mateRefID)
    {
        if(firstDifferChrom)
        {
            std::cerr << "Mate on different chromosome was not found.\n"
                      << "If you are running single chromosome, consider "
                      << "using --oneChrom to treat reads with mates on "
                      << "different chromosomes as single-ended.\n";
            firstDifferChrom = false;
        }
    }
    else if(firstSameChrom)
    {
        std::cerr << "ERROR: Records with missing mate can't be checked for "
                  << "duplicates.\n";
        firstSameChrom = false;
    }

    // Don't consider this record to be a duplicate.
    handleNonDuplicate();
}


void Dedup_LowMem::handleDuplicate(uint32_t index)
{
    // Add the index to the duplicate list.
    myDupList.push_back(index);
}
