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

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include "SamFile.h"
#include "SuperDeDuper.h"

Logger* Logger::gLogger = NULL;

const uint32_t SuperDeDuper::CLIP_OFFSET = 1000L;
const uint64_t SuperDeDuper::UNMAPPED_SINGLE_KEY = 0xffffffffffffffffULL;
const uint64_t SuperDeDuper::UNMAPPED_PAIRED_KEY = 0xfffffffffffffffeULL;
const uint32_t SuperDeDuper::EMPTY_RECORD_COUNT =  0xffffffffUL;
const int SuperDeDuper::PAIRED_QUALITY_OFFSET = 10000;
const int SuperDeDuper::MAX_REF_ID = 0xffff;
const int SuperDeDuper::LOOK_BACK = 1000;

// ::printUsage(std::ostream&) - print the usage of the program in the case of argument error
void printUsage(std::ostream& os) {
     os << "Usage: SuperDeDupper (options) --in=<InputBamFile> --out=<OutputBamFile>\n" << std::endl;
     os << "Required parameters :" << std::endl;
     os << "-i/--in [infile]  : input BAM file name (must be sorted)" << std::endl;
     os << "-o/--out [outfile] : output BAM file name (same order with original file)" << std::endl;
     os << "Optional parameters : (see SAM format specification for details)" << std::endl;
     os << "-l/--log [logfile] : log and summary statistics (default: [outfile].log)" << std::endl;
     os << "-r/--rm  : Remove duplicates (default is to mark duplicates)" << std::endl;
     os << "-f/--force-unmark  : Allow mark-duplicated BAM file and force unmarking the duplicates" << std::endl;
     os << "                     Default is to throw errors when trying to run a mark-duplicated BAM" << std::endl;
     os << "-c/--clip : Soft clip lower quality segment of overlapping paired end reads" << std::endl;
     os << "-s/--swap : Soft clip lower quality segment of overlapping paired end reads" << std::endl;
     os << "            Higher quality bases are swapped into the overlap of the unclipped end" << std::endl;
     os << "-v/--verbose  : Turn on verbose mode" << std::endl;
     os << "\n" << std::endl;
}

int main(int argc, char** argv) {
    /* --------------------------------
    / process the arguments
    / -------------------------------*/
    std::string inFile, outFile, logFile;
    bool removeFlag = false;
    bool verboseFlag = false;
    bool forceFlag = false;
    bool overlapFlag = false;
    bool overlapReplacementFlag = false;
    char c;
    int optionIndex = 0;
    
    static struct option getoptLongOptions[] = 
	{
	    { "in", required_argument, NULL, 'i' },
	    { "out", required_argument, NULL, 'o' },
	    { "rm", no_argument, NULL, 'r' },
	    { "log", no_argument, NULL, 'l' },
	    { "clip", no_argument, NULL, 'c' },
	    { "swap", no_argument, NULL, 'c' },
	    { "force-unmark", no_argument, NULL, 'c' },
	    { "verbose", no_argument, NULL, 'v' },
	    { NULL, 0, NULL, 0 },
	};
  
    while ( ( c = getopt_long(argc, argv, "i:o:rvl:csf", getoptLongOptions, &optionIndex) ) != -1 ) {
	switch(c) {
	case 'i':
	    inFile = optarg;
	    break;
	case 'o':
	    outFile = optarg;
	    break;
	case 'r':
	    removeFlag = true;
	    break;
	case 'v':
	    verboseFlag = true;
	    break;
	case 'l':
	    logFile = optarg;
	    break;
	case 's':
	    overlapReplacementFlag = true;
	case 'c':
	    overlapFlag = true;
	    break;
	case 'f':
	    forceFlag = true;
	    break;
	default:
	    printUsage(std::cerr);
	    std::cerr << "Illegal switch: " << c << std::endl;
	    abort();
	}
    }

    if (inFile.empty()) {
	printUsage(std::cerr);
	std::cerr << "Specify an input file" << std::endl;
	abort();
    }

    if (outFile.empty()) {
	printUsage(std::cerr);
	std::cerr << "Specify an output file" << std::endl;
	abort();
    }

    if (logFile.empty() ) {
	logFile = outFile + ".log";
    }
    Logger::gLogger = new Logger(logFile.c_str(), verboseFlag);

    /* -------------------------------------------------------------------
    / The arguments are processed.  Prepare the input BAM file,
    / instantiate superDeDuper, and construct the read group library map
    / ------------------------------------------------------------------*/

    SamFile samIn;

    samIn.OpenForRead(inFile.c_str());
    samIn.setSortedValidation(SamFile::COORDINATE);

    SamFileHeader header;
    samIn.ReadHeader(header);

    SuperDeDuper superDeDuper;
    superDeDuper.buildReadGroupLibraryMap(header);

    int lastReferenceID = -1;
    int lastPosition = -1;
    superDeDuper.lastReference = -1;
    superDeDuper.lastCoordinate = -1;

    // for keeping some basic statistics
    uint32_t recordCount = 0;
    uint32_t pairedCount = 0;
    uint32_t properPairCount = 0;
    uint32_t unmappedCount = 0;
    uint32_t reverseCount = 0;
    uint32_t qualCheckFailCount = 0;

    // improve this later
    superDeDuper.overlapFlag = overlapFlag;
    superDeDuper.overlapReplacementFlag = overlapReplacementFlag;

    // Now we start reading records
    SamRecord record; 
    while(samIn.ReadRecord(header, record)) {

	// Take note of properties of this record
	int flag = record.getFlag();
	if ( flag & 0x0001 ) ++pairedCount;
	if ( flag & 0x0002 ) ++properPairCount;
	if ( flag & 0x0004 ) ++unmappedCount;
	if ( flag & 0x0010 ) ++reverseCount;
	if ( flag & 0x0200 ) ++qualCheckFailCount;
	if ( ( ( flag & 0x0400 ) > 0 ) && ( ! forceFlag ) ) {
	    Logger::gLogger->error("There are records already duplicate marked.");
	    Logger::gLogger->error("Use -f to clear the duplicate flag and start the deduping procedure over");
	}

	// put the record in the appropriate maps:
	//   single reads go in fragmentMap
	//   paired reads go in pairedMap
	recordCount = samIn.GetCurrentRecordCount();
	superDeDuper.placeRecordInMaps(record, recordCount);

	// Let us know if we're moving to a new chromosome
	if (lastReferenceID != record.getReferenceID()) {
	    lastReferenceID = record.getReferenceID();
	    Logger::gLogger->writeLog("Reading ReferenceID %d\n", lastReferenceID);
	} else {
	    int position = record.get0BasedPosition();
	    // We're assuming the BAM file is sorted.  Otherwise, there is trouble.
            if ( ( position >= 0 ) && (position < lastPosition ) ) {
		Logger::gLogger->error("The BAM file is not sorted in the readName %s, reference sequence %s, position %d", 
				       record.getReadName(), 
				       header.getReferenceLabel(record.getReferenceID()).c_str(), 
				       position+1);
	    }
	    lastPosition = position;
        }

	// if we have moved to a new position, look back at previous reads for duplicates
	if (superDeDuper.hasPositionChanged(record)) {
	    superDeDuper.markDuplicatesBefore(record);
	}

	// let the user know we're not napping
	if (verboseFlag && (recordCount % 100000 == 0)) {
	    Logger::gLogger->writeLog("recordCount=%u singleKeyMap=%u pairedKeyMap=%u, dictSize=%u", 
				      recordCount, superDeDuper.fragmentMap.size(), 
				      superDeDuper.pairedMap.size(), 
				      superDeDuper.readDataMap.size());
	}
    }

    // we're finished reading record so clean up the duplicate search and close the input file
    superDeDuper.markDuplicatesBefore(SuperDeDuper::MAX_REF_ID, 0);
    samIn.Close();

    // print some statistics
    uint32_t totalDuplicates = superDeDuper.singleDuplicates + 2*superDeDuper.pairedDuplicates;
    Logger::gLogger->writeLog("--------------------------------------------------------------------------");
    Logger::gLogger->writeLog("SUMMARY STATISTICS OF THE READS");
    Logger::gLogger->writeLog("Total number of reads: %u",recordCount);
    Logger::gLogger->writeLog("Total number of paired-end reads: %u",pairedCount);
    Logger::gLogger->writeLog("Total number of properly paired reads: %u",properPairCount);
    Logger::gLogger->writeLog("Total number of unmapped reads : %u",unmappedCount);
    Logger::gLogger->writeLog("Total number of reverse strand mapped reads: %u",reverseCount);
    Logger::gLogger->writeLog("Total number of QC-failed reads : %u",qualCheckFailCount);
    Logger::gLogger->writeLog("Size of singleKeyMap (must be zero) : %u",superDeDuper.fragmentMap.size());
    Logger::gLogger->writeLog("Size of pairedKeyMap (must be zero) : %u",superDeDuper.pairedMap.size());
    Logger::gLogger->writeLog("Total number of duplicate single-ended or mate-unpaired reads : %u",superDeDuper.singleDuplicates);
    Logger::gLogger->writeLog("Total number of duplicate paired-end reads (both ends mapped) : %u",superDeDuper.pairedDuplicates);
    Logger::gLogger->writeLog("Total number of duplicated records: %u",totalDuplicates);
    if (overlapFlag) {
	Logger::gLogger->writeLog("Total number of overlapping pairs clipped: %u",superDeDuper.overlappingPairs);
    }
    Logger::gLogger->writeLog("--------------------------------------------------------------------------");
    Logger::gLogger->writeLog("Sorting the indices of %d duplicated records",superDeDuper.duplicateIndices.size());

    // sort the indices of duplicate records
    std::sort(superDeDuper.duplicateIndices.begin(), superDeDuper.duplicateIndices.end(), std::less<uint32_t> ());

    // get ready to write the output file by making a second pass through the input file
    samIn.OpenForRead(inFile.c_str());
    samIn.ReadHeader(header);

    SamFile samOut;
    samOut.OpenForWrite(outFile.c_str());
    samOut.WriteHeader(header);

    // an iterator to run through the duplicate indices
    std::vector<uint32_t>::iterator currentDupIndex = superDeDuper.duplicateIndices.begin();

    // let the user know what we're doing
    Logger::gLogger->writeLog("\nWriting %s", outFile.c_str());

    // count the duplicate records as a check
    uint32_t singleDuplicates(0), pairedDuplicates(0);

    // start reading records and writing them out
    while(samIn.ReadRecord(header, record)) {
	uint32_t currentIndex = samIn.GetCurrentRecordCount();

	bool foundDup = (currentIndex == *currentDupIndex);

	// if we are looking for overlapping pairs, check to see if we need to update this record
	//    we may need to update the cigar, sequence, and quality array
	if (superDeDuper.overlapFlag) {
	    SuperDeDuper::UInt32ToUpdateDataPointerIterator recordIterator = 
		superDeDuper.recordUpdates.find(currentIndex);
	    if (recordIterator != superDeDuper.recordUpdates.end()) {
		record.setCigar(recordIterator->second->cigar);
		record.setSequence(recordIterator->second->sequence.c_str());
		record.setQuality(recordIterator->second->baseQualities.c_str());
		delete recordIterator->second; 
		superDeDuper.recordUpdates.erase(samIn.GetCurrentRecordCount());
	    }

	}

	// modify the duplicate flag and write out the record, if it's appropriate
	int flag = record.getFlag();
	if (foundDup) {                     // this record is a duplicate
	    record.setFlag( flag | 0x400 ); // mark duplicate
	    currentDupIndex++;
	    // increment duplicate counters to verify we found them all
	    if ( ( ( flag & 0x0001 ) == 0 ) || ( flag & 0x0008 ) ) { // unpaired or mate unmapped
		singleDuplicates++;
	    }
	    else {
		pairedDuplicates++;
	    }
	    // write the record if we are not removing duplicates
	    if (!removeFlag ) samOut.WriteRecord(header, record);
	}
	else if (forceFlag ) { 
	    // this is not a duplicate we've identified but we want to remove any duplicate marking
	    record.setFlag( flag & 0xfffffbff ); // unmark duplicate
	    samOut.WriteRecord(header, record);
	}
	else { // not a duplicate we've identified and we aren't worried about existing marking
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
    Logger::gLogger->writeLog("\nSuperDeDuper complete!");
    return 0;
}

// Now that we've reached coordinate on chromosome reference, look back to find any duplicates
void SuperDeDuper::markDuplicatesBefore(uint32_t reference, uint32_t coordinate) {
    // Find the key corresponding to the current position
    // We will first search through single reads up to this position
    uint64_t key = makeKey(reference, coordinate, false, 0);
    SingleKeyToReadDataPointerMapIterator fragmentFinish = fragmentMap.lower_bound(key);

    // For each key k < fragmentFinish, look through reads indexed by that key to find the
    //   one with the highest quality, mark the others as duplicates, and clean up the maps
    for (SingleKeyToReadDataPointerMapIterator iterator = fragmentMap.begin(); iterator != fragmentFinish; iterator++) {
	if (iterator->second.size() > 1) {
	    int maxQuality = -1;
	    ReadDataPointerVectorIterator maxRead = iterator->second.end();
	    for (ReadDataPointerVectorIterator currentRead = iterator->second.begin(); 
		 currentRead != iterator->second.end(); 
		 currentRead++) {
		if ( (*currentRead)->getPairedBaseQuality() > maxQuality) {
		    maxQuality = (*currentRead)->getPairedBaseQuality();
		    maxRead = currentRead;
		}
	    }
	    if (maxRead == iterator->second.end()) {
		Logger::gLogger->error("In markDuplicates:  no positive best quality found.");
	    }
	    for (ReadDataPointerVectorIterator currentRead = iterator->second.begin(); 
		 currentRead != iterator->second.end(); 
		 currentRead++) {
		if (currentRead != maxRead) {
		    if ( (*currentRead)->paired == false) {
			duplicateIndices.push_back((*currentRead)->recordCount1);
			singleDuplicates++;
		    }
		}

		StringToReadDataPointerMapIterator readDataMapIterator = readDataMap.find( (*currentRead)->readName);
		if (readDataMapIterator == readDataMap.end()) {
		    Logger::gLogger->error("In markDuplicatesBefore:  Cannot find %s in readDataMap", ((*currentRead)->readName).c_str());
		}
		if (readDataMapIterator->second->paired == false) {
		    delete readDataMapIterator->second; 
		    readDataMap.erase(readDataMapIterator);
		}
	    }	    
	} else {  // in this case, there are no duplicates
	    StringToReadDataPointerMapIterator readDataMapIterator = readDataMap.find(iterator->second.front()->readName);
	    if (readDataMapIterator->second->paired == false) {
		delete readDataMapIterator->second; 
		readDataMap.erase(readDataMapIterator);
	    }
	}
			
    }
    if (fragmentFinish != fragmentMap.begin()) {
	fragmentMap.erase(fragmentMap.begin(), fragmentFinish);
    }

    // Now do the same thing with the paired reads
    PairedKey pairedKey(0, key);
    PairedKeyToReadDataPointerMapIterator pairedFinish = pairedMap.lower_bound(pairedKey);
    for (PairedKeyToReadDataPointerMapIterator iterator = pairedMap.begin(); iterator != pairedFinish; iterator++) {
	if (iterator->second.size() > 1) {
	    int maxQuality = -1;
	    ReadDataPointerVectorIterator maxRead = iterator->second.end();	
	    for (ReadDataPointerVectorIterator currentRead = iterator->second.begin(); 
		 currentRead != iterator->second.end(); 
		 currentRead++) {
		if ( (*currentRead)->getPairedBaseQuality() > maxQuality) {
		    maxQuality = (*currentRead)->getPairedBaseQuality();
		    maxRead = currentRead;
		}
	    }
	    if (maxRead == iterator->second.end()) {
		Logger::gLogger->error("In markDuplicates:  no positive best quality found.");
	    }

	    for (ReadDataPointerVectorIterator currentRead = iterator->second.begin(); 
		 currentRead != iterator->second.end(); 
		 currentRead++) {
		if (currentRead != maxRead) {
		    duplicateIndices.push_back((*currentRead)->recordCount1);
		    duplicateIndices.push_back((*currentRead)->recordCount2);
		    pairedDuplicates++;
		}

		StringToReadDataPointerMapIterator readDataMapIterator = readDataMap.find( (*currentRead)->readName);
		if (readDataMapIterator == readDataMap.end()) {
		    Logger::gLogger->error("In markDuplicatesBefore:  Cannot find %s in readDataMap", ((*currentRead)->readName).c_str());
		}
		delete readDataMapIterator->second;  
		readDataMap.erase(readDataMapIterator);
	    }
	} else { // there are no duplicates
	    StringToReadDataPointerMapIterator readDataMapIterator = readDataMap.find(iterator->second.front()->readName);
	    delete readDataMapIterator->second; 
	    readDataMap.erase(readDataMapIterator);
	}
    }
    if (pairedFinish != pairedMap.begin()) {
	pairedMap.erase(pairedMap.begin(), pairedFinish);
    }
    return;
}

// Look at reads before this record and determine duplicates
void SuperDeDuper::markDuplicatesBefore(SamRecord& record) {
    uint32_t coordinate = record.get0BasedPosition();
    uint32_t lookBackCoordinate = coordinate - LOOK_BACK;
    if (lookBackCoordinate < 0) {
	lookBackCoordinate = 0;
    }
    return markDuplicatesBefore(record.getReferenceID(), lookBackCoordinate);
}

// determine whether the position of record is different from the previous record
bool SuperDeDuper::hasPositionChanged(SamRecord& record) {
    if (lastReference < record.getReferenceID() || lastCoordinate < record.get0BasedPosition()) {
	lastReference = record.getReferenceID();
	lastCoordinate = record.get0BasedPosition();
	return true;
    }
    return false;
}

// when record is read, we will put it into the appropriate maps
void SuperDeDuper::placeRecordInMaps(SamRecord& record, uint32_t recordCount) {
    int flag = record.getFlag();

    if (flag & 0x0002) properPair++;
    if (flag & 0x0004) {  // fragment is unmapped so ignore it
	unmapped++;
	return; 
    }

    // save some basic information about the record in a ReadData structure
    uint64_t key = makeKeyFromRecord(record);
    ReadData* readDataPointer = new ReadData();
    readDataPointer -> key1 = key;
    readDataPointer -> key2 = UNMAPPED_SINGLE_KEY;
    readDataPointer -> recordCount1 = recordCount;
    readDataPointer -> recordCount2 = EMPTY_RECORD_COUNT;
    readDataPointer -> readName = record.getReadName();
    readDataPointer -> baseQuality = getBaseQuality(record);

    if ( ( (flag & 0x0001) == 0) || (flag & 0x0008) ) {  
	// if it's a single read or its mate is unmapped, 
	//   label it as not paired and put it into fragmentMap
	if ((flag & 0x0001) == 0) singleRead++;
	readDataPointer -> paired = false;
	if(readDataMap.insert( std::pair <std::string, ReadData*> (record.getReadName(), readDataPointer)).second == false) {
	    Logger::gLogger -> error("The fragment %s has already been read.", record.getReadName());
	    return;
	}
	fragmentMap[key].push_back(readDataPointer);
	return;
    }

    // We have a paired end read with its mate mapped
    readDataPointer -> paired = true;

    // Let's determine if we have already seen the mate
    StringToReadDataPointerMapIterator earlierRead = readDataMap.find(record.getReadName());

    if(earlierRead == readDataMap.end()) { // We haven't seen this read before
	firstPair++;
	if (overlapFlag) {
	    // store additional information to process overlaps
	    OverlapData* overlapDataPointer = new OverlapData();
	    overlapDataPointer->clippedEnd = record.get0BasedAlignmentEnd();
	    overlapDataPointer->cigar = *record.getCigarInfo();
	    overlapDataPointer->start = record.get0BasedPosition();
	    overlapDataPointer->baseQualities = std::string(record.getQuality());
	    overlapDataPointer->sequence = std::string(record.getSequence());
	    readDataPointer->overlapData = overlapDataPointer;
	}
	// put this in fragmentMap, record that we've seen it, and move on to the next record
	fragmentMap[key].push_back(readDataPointer);
	readDataMap.insert( std::pair <std::string, ReadData*> (record.getReadName(), readDataPointer));
	return;
    }

    // We have already seen this read's mate
    foundPair++;
    ReadData* readData = earlierRead->second;

    if (overlapFlag && record.get0BasedPosition() <= readData->overlapData->clippedEnd) {
	// We're in here if the two ends overlap
	if (readData->overlapData->clippedEnd >= record.get0BasedAlignmentEnd()) {
	    // something more interesting will happen here in the next version
	} else {
	    overlappingPairs++;
	    processOverlap(readData, record, recordCount);
	}
    }

    // Update the information in the ReadData structure
    fragmentMap[key].push_back(readDataPointer);
    readData -> baseQuality += getBaseQuality(record);

    if (key < readData -> key1) {
	readData -> key2 = readData -> key1;
	readData -> key1 = key;
	readData -> recordCount2 = readData -> recordCount1;
	readData -> recordCount1 = recordCount;
    } else {
	readData -> key2 = key;
	readData -> recordCount2 = recordCount;
    }

    // Put the ReadData structure in the paired key map and go to the next record
    PairedKey pkey(readData -> key1, readData -> key2);
    pairedMap[pkey].push_back(readData);

}

// When the read summarized in readData and its mate in record overlap, we'll do some work here
void SuperDeDuper::processOverlap(ReadData * readData, SamRecord & record, uint32_t recordCount) {

    // The first read has subscript 0, the second subscript 1
    // First, we record some positions in the reference coordinate
    int32_t start0 = readData->overlapData->start;
    int32_t end0 = readData->overlapData->clippedEnd;
    int32_t start1 = record.get0BasedPosition();
    
    // get the two cigars
    Cigar * cigar0 = &(readData->overlapData->cigar);
    Cigar * cigar1 = record.getCigarInfo();

    // For each read, find the range of query indices describing the overlap
    // First for read 0
    int32_t queryIndex0End   = cigar0->getQueryIndex(end0 - start0);
    int32_t queryIndex0Start = queryIndex0End;
    int32_t targetRefOffset0 = start1 - start0;
    for (int32_t i = queryIndex0End; i >= 0; i--) {
	int32_t refOffset = cigar0->getRefOffset(i);
	if (refOffset == Cigar::INDEX_NA) continue;
	if (refOffset < targetRefOffset0) break;
	queryIndex0Start = i;
    }

    // Then for read 1
    int32_t queryIndex1Start = cigar1->getQueryIndex(0);
    int32_t queryIndex1End   = queryIndex1Start;
    int32_t targetRefOffset1 = end0 - start1;
    for (int32_t i = queryIndex1Start; i < record.getReadLength(); i++) {
	int32_t refOffset = cigar1->getRefOffset(i);
	if (refOffset == Cigar::INDEX_NA) continue;
	if (refOffset > targetRefOffset1) break;
	queryIndex1End = i;
    }

    // Find the average base quality in the overlap for read 0
    float baseQuality0 = 0;
    for (int i = queryIndex0Start; i <= queryIndex0End; i++) {
	baseQuality0 += static_cast<float>(readData->overlapData->baseQualities[i]);
    }
    baseQuality0 /= (queryIndex0End - queryIndex0Start + 1);

    // Find the average base quality in the overlap for read 1
    float baseQuality1 = 0;
    std::string baseQualities1 = std::string(record.getQuality());
    for (int i = queryIndex1Start; i <= queryIndex1End; i++) {
	baseQuality1 += static_cast<float>(baseQualities1[i]);
    }
    baseQuality1 /= (queryIndex1End - queryIndex1Start + 1);

    if (baseQuality0 < baseQuality1) {
	// read 0 has a lower average base quality so it will be soft clipped in the overlap
	if (overlapReplacementFlag) {
	    // if this flag is set, we will replace lower quality bases in read 1 with bases from read 0
	    std::string expandedCigar0, expandedCigar1;
	    cigar0->getExpandedString(expandedCigar0);
	    cigar1->getExpandedString(expandedCigar1);
	    std::string sequence1 = std::string(record.getSequence());

	    // walk through corresponding bases in the two reads 
	    for (int32_t reference = start1; reference <= end0; reference++) {
		int32_t qi0 = cigar0->getQueryIndex(reference-start0);
		if (qi0 == Cigar::INDEX_NA) continue;
		int32_t qi1 = cigar1->getQueryIndex(reference-start1);
		if (qi1 == Cigar::INDEX_NA) continue;
		if (readData->overlapData->baseQualities[qi0] > baseQualities1[qi1]) {
		    // if the base is read 0 has a higher quality, we substitute it into read 1
		    baseQualities1.replace(qi1, 1, 1, readData->overlapData->baseQualities.at(qi0));
		    sequence1.replace(qi1, 1, 1, readData->overlapData->sequence.at(qi0));
		    expandedCigar1.replace(qi1, 1, 1, expandedCigar0.at(qi0));
		}
	    }

	    // Store the updated information about read 1
	    UpdateData* updateDataPointer = new UpdateData();
	    updateDataPointer->sequence = sequence1;
	    updateDataPointer->baseQualities = baseQualities1;
	    updateDataPointer->cigar = *rollupCigar(expandedCigar1);
	    recordUpdates.insert(std::pair<uint32_t, UpdateData*> (recordCount, updateDataPointer));
	}

	// clip read 0 and store the updated information
	cigar0 = insertClipsIntoCigar(cigar0, queryIndex0Start, queryIndex0End);
	
	UpdateData* updateDataPointer = new UpdateData();
	updateDataPointer->sequence = readData->overlapData->sequence;
	updateDataPointer->baseQualities = readData->overlapData->baseQualities;
	updateDataPointer->cigar = *cigar0;

	recordUpdates.insert(std::pair<uint32_t, UpdateData*> (readData->recordCount1, updateDataPointer));
	    
    } else {
	// In this case, read 1 is the lower quality read,
	//    so we do the same thing with the role of the reads reversed
	if (overlapReplacementFlag) {
	    std::string expandedCigar0, expandedCigar1;
	    cigar0->getExpandedString(expandedCigar0);
	    cigar1->getExpandedString(expandedCigar1);
	    std::string sequence1 = std::string(record.getSequence());

	    std::string baseQualities0 = readData->overlapData->baseQualities;
	    std::string sequence0 = readData->overlapData->sequence;

	    for (int32_t reference = start1; reference <= end0; reference++) {
		int32_t qi0 = cigar0->getQueryIndex(reference-start0);
		if (qi0 == Cigar::INDEX_NA) continue;
		int32_t qi1 = cigar1->getQueryIndex(reference-start1);
		if (qi1 == Cigar::INDEX_NA) continue;
		if (baseQualities0[qi0] < baseQualities1[qi1]) {
		    baseQualities0.replace(qi0, 1, 1, baseQualities1[qi1]);
		    sequence0.replace(qi0, 1, 1, sequence1.at(qi1));
		    expandedCigar0.replace(qi0, 1, 1, expandedCigar1.at(qi1));
		}
	    }
	    
	    UpdateData* updateDataPointer = new UpdateData();
	    updateDataPointer->sequence = sequence0;
	    updateDataPointer->baseQualities = baseQualities0;
	    updateDataPointer->cigar = *rollupCigar(expandedCigar0);
	    recordUpdates.insert(std::pair<uint32_t, UpdateData*> (readData->recordCount1, updateDataPointer));
	}

	UpdateData* updateDataPointer = new UpdateData();
	updateDataPointer->sequence = record.getSequence();
	updateDataPointer->baseQualities = baseQualities1;
	updateDataPointer->cigar = *insertClipsIntoCigar(cigar1, queryIndex1Start, queryIndex1End);
	recordUpdates.insert(std::pair<uint32_t, UpdateData*> (recordCount, updateDataPointer));

    }

    delete readData->overlapData; 
}

// turns an expanded cigar string in a rolled up cigar string
Cigar * SuperDeDuper::rollupCigar(std::string cigarString) {
    CigarRoller * cigarRoller = new CigarRoller();
    for (uint32_t i = 0; i < cigarString.length(); i++) {
	cigarRoller->Add(cigarString.at(i), 1);
    }
    return cigarRoller;
}

// replaces portions of a cigar with soft clips
Cigar * SuperDeDuper::insertClipsIntoCigar(Cigar * cigar, int32_t begin, int32_t end) {
    std::string cigarString;
    cigar->getExpandedString(cigarString);
    int32_t insertLength = end-begin+1;
    cigarString.replace(begin, insertLength, insertLength, 'S');

    return rollupCigar(cigarString);
}

// Finds the total base quality of a read 
int SuperDeDuper::getBaseQuality(SamRecord & record) {
    const char* baseQualities = record.getQuality();
    int readLength = record.getReadLength();
    int quality = 0.;
    for(int i=0; i < readLength; ++i) {
	int q = static_cast<int>(baseQualities[i])-33;
	if ( q >= 15 ) quality += q;
    }
    return quality;
}

// makes a key from the chromosome number, coordinate, orientation, and libraryID
//   single reads with equal keys are duplicates
uint64_t SuperDeDuper::makeKey(uint32_t referenceID, uint32_t coordinate, bool orientation, uint32_t libraryID) {
    return ( (0xffff000000000000 & ( static_cast<uint64_t>(referenceID) << 48))  |
	     (0x0000ffffffff0000 & ( static_cast<uint64_t>(coordinate + CLIP_OFFSET) << 16)) |
	     (0x000000000000ff00 & ( static_cast<uint64_t>(orientation) << 8)) |
	     (0x00000000000000ff & ( static_cast<uint64_t>(libraryID))));
}

// extract the relevant information from the key and make the record
uint64_t SuperDeDuper::makeKeyFromRecord(SamRecord& record) {
    int32_t referenceID = record.getReferenceID();
    bool orientation = (record.getFlag() & 0x0010) > 0;
    int32_t coordinate = orientation ? record.get0BasedUnclippedEnd() : record.get0BasedUnclippedStart();
    if ( ( referenceID < 0 ) || ( coordinate + CLIP_OFFSET < 0 ) ) {
	Logger::gLogger->error("SuperDeDuper::makeKeyFromRecord(record) - refID or coordinate is negative. refID = %d, coordinate = %d",
			       referenceID, coordinate + CLIP_OFFSET);
	return UNMAPPED_SINGLE_KEY; // this represents an unmapped read end
    }
    return makeKey(referenceID, coordinate, orientation, getLibraryID(record));
}

// build the read group library map
void SuperDeDuper::buildReadGroupLibraryMap(SamFileHeader& header) {
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
	Logger::gLogger->error("More than 255 library names are identified. SuperDeDuper currently only allows up to 255 library names");
    }

}    

// get the libraryID of a record
uint32_t SuperDeDuper::getLibraryID(SamRecord& record, bool checkTags) {
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

Logger::Logger(const char* filename, bool verbose) {
    b_verbose = verbose;
    fp_log = fopen(filename, "w");
    if ( fp_log == NULL ) {
	fprintf(stderr,"ERROR: Cannot open the log file %s. Check if the directory exists and you have the permission to create a file", filename);
	abort();
    }
    fp_err = stderr;
}

void Logger::writeLog(const char* format, ... ) {
    va_list args;
    va_start (args, format);
    vfprintf(fp_log, format, args);
    va_end (args);
    fprintf(fp_log, "\n");
    
    if ( b_verbose ) {
	va_start (args, format);
	vfprintf(fp_err, format, args);
	va_end (args);
	fprintf(fp_err, "\n");
    }
}

void Logger::error(const char* format, ... ) {
    va_list args;
    va_start (args, format);
    fprintf(fp_log, "ERROR: ");
    vfprintf(fp_log, format, args);
    va_end (args);
    fprintf(fp_log, "\n");
    
    va_start (args, format);
    fprintf(fp_err, "ERROR : ");
    vfprintf(fp_err, format, args);
    va_end (args);
    fprintf(fp_err, "\n");
    
    abort();
}

void Logger::warning(const char* format, ... ) {
    va_list args;
    va_start (args, format);
    fprintf(fp_log, "WARNING: ");
    vfprintf(fp_log, format, args);
    va_end (args);
    fprintf(fp_log, "\n");
    
    if ( b_verbose ) {
	va_start (args, format);
	fprintf(fp_err, "WARNING : ");
	vfprintf(fp_err, format, args);
	va_end (args);
	fprintf(fp_err, "\n");
    }
}
