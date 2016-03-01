/*
 * MateVectorByRN.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: fanzhang
 */

#include "MateVectorByRN.h"

#include <stdio.h>
#include "Bam2FastQ.h"
#include <algorithm>
#include <mutex>

std::mutex myLock2;

#define SHRINK_LIMIT_PARAMETER 100000
#define MAX_LIMIT_PARAMETER 500000
#define RECORD_POOL_PARAMETER 1000000

MateVectorByRN::MateVectorByRN() {
	// TODO Auto-generated constructor stub

}

MateVectorByRN::~MateVectorByRN() {
	// TODO Auto-generated destructor stub
	myMateBuffer[0].clear();
	myMateBuffer[1].clear();
	tmpFileNameList.clear();

}

MateVectorByRN::MateVectorByRN(int tmpVectorIndex, Bam2FastQ* conversionHost,
		std::string chromosome, long CoordinateStart, long CoordinateEnd) :
		bufferInUse(0), myNumPairs(0), count(0), SHRINK_LIMIT(
				SHRINK_LIMIT_PARAMETER), MAX_LIMIT(MAX_LIMIT_PARAMETER)/*,myPool(RECORD_POOL_PARAMETER)*/{

	prevRec = NULL;
	prevRN = "";
	//std::cerr<<conversionHost->myFirstFileNameExt<<"\t"<<conversionHost->mySecondFileNameExt<<std::endl;
	FirstFileNameExt = conversionHost->myFirstFileNameExt;
	SecondFileNameExt = conversionHost->mySecondFileNameExt;
	UnpairedFileNameExt = conversionHost->myUnpairedFileNameExt;
	host = conversionHost;
	chr = chromosome;
	begin = CoordinateStart;
	end = CoordinateEnd;
	vectorIndex = tmpVectorIndex;
	tmpFileIndex = 0;
	myPool = new SamRecordPool;

}

void MateVectorByRN::clearMyPool() {
	delete myPool;
	myPool = NULL;
}

int MateVectorByRN::HandlePairedRN(SamRecord& samRec, bool ready2Dump) { //all SamRecord* handled via this function will be released

	if (prevRec == NULL) {
		prevRec = &samRec;

	} else {
		if (strcmp(prevRec->getReadName(), samRec.getReadName()) != 0) {
			if (ready2Dump) {
				//++myNumMateFailures;
//	            if (SamFlag::isFirstFragment(samRec.getFlag()))
//	            	host->writeFastQ(samRec, tmpFile,  UnpairedFileNameExt, myPool, true, host->myFirstRNExt.c_str());
//	            else
//	            	host->writeFastQ(samRec, tmpFile,  UnpairedFileNameExt, myPool, true, host->mySecondRNExt.c_str());
				//ifwrite(tmpFile,record.second,record.first);
				tmpFile.WriteRecord(tmpHeader, *prevRec);
				myPool->releaseRecord(prevRec);
				prevRec = &samRec;
			} else {
				// Read Name does not match, error, did not find pair.
				myMateBuffer[1 - bufferInUse].push_back(
						std::make_pair(std::string(prevRec->getReadName()),
								prevRec));
				// Save this record to check against the next one.
				prevRec = &samRec;

			}
		} else {
			// Matching ReadNames.
			// Found the mate.

			// Check which is the first in the pair.
			if (SamFlag::isFirstFragment(samRec.getFlag())) {
//				if (SamFlag::isFirstFragment(prevRec->getFlag())) {
//					std::cerr << "Both reads of " << samRec.getReadName()
//							<< " are first fragment, so "
//							<< "splitting one to be in the 2nd fastq.\n";
//				}

				myLock2.lock();

				host->writeFastQ(samRec, host->myFirstFile, FirstFileNameExt,
						myPool, false, host->myFirstRNExt.c_str());
				host->writeFastQ(*prevRec, host->mySecondFile,
						SecondFileNameExt, myPool, false,
						host->mySecondRNExt.c_str());
				++host->myNumPairs;
				myLock2.unlock();

			} else {
//				if (!SamFlag::isFirstFragment(prevRec->getFlag())) {
//					std::cerr << "Neither read of " << samRec.getReadName()
//							<< " are first fragment, so "
//							<< "splitting one to be in the 2nd fastq.\n";
//				}

				myLock2.lock();

				host->writeFastQ(*prevRec, host->myFirstFile, FirstFileNameExt,
						myPool, false, host->myFirstRNExt.c_str());
				host->writeFastQ(samRec, host->mySecondFile, SecondFileNameExt,
						myPool, false, host->mySecondRNExt.c_str());
				++host->myNumPairs;
				myLock2.unlock();

			}
			// No previous record.
			prevRec = NULL;
		}
	}

	return 0;
}

bool MateVectorByRN::IsInTheRegion(SamRecord* record) {
	return (record->get1BasedAlignmentEnd() <= end);
}

int MateVectorByRN::Read(const std::string& bamFile) {		//end not included

	SamFile FIN(bamFile.c_str(), SamFile::READ, &tmpHeader);
	if (!FIN.IsOpen()) {
		std::cerr << "open " << bamFile << " failed!" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (chr == "*") {
		tmpFileNameBase = std::string(host->myTmpOutBase.c_str())
				+ std::string("_") + "unmap" + ":" + std::to_string(begin) + "_"
				+ std::to_string(end) + "_";
	} else
		tmpFileNameBase = std::string(host->myTmpOutBase.c_str())
				+ std::string("_") + chr + ":" + std::to_string(begin) + "_"
				+ std::to_string(end) + "_";
	FIN.ReadBamIndex();
	//std::cerr<<"Region:"<<chr<<":"<<begin<<"-"<<end<<"\t in process"<<std::endl;
	if (!FIN.SetReadSection(chr.c_str(), begin, end)) {
		std::cerr << "Region:" << chr << ":" << begin << "-" << end
				<< "\tdoes not exist" << std::endl;
		return 1;
	}
	SamRecord* recordPtr = myPool->getRecord();
	if (recordPtr == NULL) {
		std::cerr << "Get record from myPool failed!" << std::endl;
		exit(EXIT_FAILURE);
	}
	int alreadyRead = 0;
	while (FIN.ReadRecord(tmpHeader, *recordPtr)) {
		if (!IsInTheRegion(recordPtr))
			continue;
		Add(recordPtr);
		alreadyRead++;
		recordPtr = myPool->getRecord();

	}
	myPool->releaseRecord(recordPtr);
	DumpMateVector();
	ClearWorkingBuffer();
	FlipWorkingBuffer();
	assert(myMateBuffer[bufferInUse].size() == 0);		//ClearWorkingBuffer();

	return 0;
}

int MateVectorByRN::Add(SamRecord* record) {
	if (count >= MAX_LIMIT)			//ready to dump tmp files
			{
		DumpMateVector();
		SHRINK_LIMIT = SHRINK_LIMIT_PARAMETER;
		ClearWorkingBuffer();
		FlipWorkingBuffer();
		count = 0;
		//restart
		myMateBuffer[bufferInUse].push_back(
				std::make_pair(std::string(record->getReadName()), record));
		count++;
		return 2;			//finish job
	} else if (count >= SHRINK_LIMIT) {
		if (prevRec != NULL) {
			myMateBuffer[bufferInUse].push_back(
					std::make_pair(prevRec->getReadName(), prevRec));
			prevRec = NULL;
		}
		Sort();
		//if(count!=myMateBuffer[bufferInUse].size()) {std::cerr<<"Buffer vector size doesn't match!"<<std::endl;exit(EXIT_FAILURE);}
		for (uint i = 0; i < myMateBuffer[bufferInUse].size(); ++i) {
			HandlePairedRN(*myMateBuffer[bufferInUse][i].second, false);
		}
		SHRINK_LIMIT *= 2;
		ClearWorkingBuffer();
		FlipWorkingBuffer();
		count = myMateBuffer[bufferInUse].size();
		//restart
		myMateBuffer[bufferInUse].push_back(
				std::make_pair(std::string(record->getReadName()), record));
		count++;
		return 1;

	} else {
		myMateBuffer[bufferInUse].push_back(
				std::make_pair(std::string(record->getReadName()), record));
		count++;
	}

	return 0;
}

void MateVectorByRN::Sort() {
	//fprintf(stderr,"Sort %dth vector working buffer...\n",vectorIndex);
	std::sort(myMateBuffer[bufferInUse].begin(),
			myMateBuffer[bufferInUse].end(), RecordComparison);
}

int MateVectorByRN::DumpMateVector() {

	tmpFileName = tmpFileNameBase + std::to_string(vectorIndex) + "th.vector."
			+ std::to_string(tmpFileIndex) + std::string("th.dump.tmp");
	tmpFileNameList.push_back(tmpFileName);

	tmpFile.OpenForWrite(tmpFileName.c_str(), &tmpHeader);
	if (!tmpFile.IsOpen()) {
		std::cerr << "open " << tmpFileName
				<< " failed, check if there is enough space to create file!"
				<< std::endl;
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "Dumping %dth tmp file of %dth vector in %s section...\n",
			tmpFileIndex, vectorIndex, tmpFileNameBase.c_str());
	if (prevRec != NULL) {
		myMateBuffer[bufferInUse].push_back(
				std::make_pair(prevRec->getReadName(), prevRec));
		prevRec = NULL;
	}
	Sort();
	for (uint i = 0; i != myMateBuffer[bufferInUse].size(); ++i) {
		HandlePairedRN(*myMateBuffer[bufferInUse][i].second, true);
	}
	if (prevRec != NULL) {			//last one
		tmpFile.WriteRecord(tmpHeader, *prevRec);
		myPool->releaseRecord(prevRec);
		prevRec = NULL;
	}

	tmpFile.Close();
	tmpFileIndex++;

	return 0;

}

void MateVectorByRN::ClearWorkingBuffer() {
	myMateBuffer[bufferInUse].clear();
}

