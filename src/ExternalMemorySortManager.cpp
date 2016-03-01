/*
 * ExternalMemorySortManager.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: fanzhang
 */

#include "ExternalMemorySortManager.h"
#include "Bam2FastQ.h"
#include <thread>
#include <mutex>
#include <cstdio>

#define Section_Skip 10000000
std::mutex myLock;
ExternalMemorySortManager::ExternalMemorySortManager() {
	// TODO Auto-generated constructor stub

}

ExternalMemorySortManager::~ExternalMemorySortManager() {
	// TODO Auto-generated destructor stub
}

ExternalMemorySortManager::ExternalMemorySortManager(Bam2FastQ * tmpHost,
		std::string tmpBamFile, int nThread, int recordLimit,
		std::string tmpbedFile) {
	THREAD_NUM = nThread;
	RECORD_LIMIT = recordLimit;
	bamFile = tmpBamFile;
	bedFile = tmpbedFile;
	host = tmpHost;
}

static int parseRegion(std::string& line, std::string &chr, long &begin,
		long &end) {
	//chr2:100-200
	//chr2\t100\t200
	uint pos(0), pos2(0);
	if ((pos = line.find(":", 0)) != line.npos
			&& (pos2 = line.find("-", 0)) != line.npos)	//find ':'
					{
		chr = line.substr(0, pos);
		begin = atoi(line.substr(pos + 1, pos2 - pos + 1).c_str());
		end = atoi(line.substr(pos2, line.length() - pos2 + 1).c_str());
	} else if ((pos = line.find("\t", 0)) != line.npos)	//find '\t'
			{
		chr = line.substr(0, pos);
		pos2 = line.find_last_of("\t", pos + 1);
		begin = atoi(line.substr(pos + 1, pos2 - pos + 1).c_str());
		end = atoi(line.substr(pos2, line.length() - pos2 + 1).c_str());
	} else {
		std::cerr
				<< "Unrecognized region format! Please use format:\"chr2:100-200\" or \"chr2\t100\t200\""
				<< std::endl;
		exit(EXIT_FAILURE);
	}
	return 0;
}
int ExternalMemorySortManager::FindAllBeginPointer() {

	std::ifstream fin(bedFile.c_str());
	std::priority_queue<POPCORN, std::vector<POPCORN>, PopcornComparison> tmpHeap;
	if (bedFile == "") {	//create segments on the fly
		SamFileHeader tmpHeader;
		SamFile fin(bamFile.c_str(), SamFile::READ, &tmpHeader);
		int nSQ = tmpHeader.getNumSQs();
		std::string chrName;
		int seqLength(0);
		SamHeaderRecord* tmpHeaderRecord;
		int vectorIndex(0);
		int begin(1);
		for (int i = 0; i != nSQ; ++i) {
			tmpHeaderRecord = tmpHeader.getNextSQRecord();
			chrName = tmpHeaderRecord->getTagValue("SN");
			seqLength = atoi(tmpHeaderRecord->getTagValue("LN"));
			//std::cerr<<"Process:"<<chrName<<"\t"<<seqLength<<std::endl;
			begin = 1;
			//if(0)
			for (; begin + Section_Skip < seqLength; begin += Section_Skip)
			{
				MateVectorList.push_back(
						new MateVectorByRN(vectorIndex, host, chrName, begin,
								begin + Section_Skip - 1));
				masterTmpFileList.push_back(std::string(""));
				sortHeap.push_back(tmpHeap);
				vectorIndex++;
			}
			//else begin=(seqLength-Section_Skip)>0?(seqLength-Section_Skip):1;//for debug purpose
			MateVectorList.push_back(
					new MateVectorByRN(vectorIndex, host, chrName, begin,
							seqLength));
			masterTmpFileList.push_back(std::string(""));
			sortHeap.push_back(tmpHeap);
			vectorIndex++;
		}
		MateVectorList.push_back(
				new MateVectorByRN(vectorIndex, host, std::string("*"), begin,
						begin + 1));
		masterTmpFileList.push_back(std::string(""));
		sortHeap.push_back(tmpHeap);
		vectorIndex++;

	} else {	//divide task according to bed file

		if (!fin.is_open()) {
			std::cerr << "Open file " << bedFile << " failed!";
			exit(EXIT_FAILURE);
		}
		int vectorIndex(0);
		std::string line, chr;
		long begin(0), end(0);
		while (getline(fin, line)) {
			parseRegion(line, chr, begin, end);
			MateVectorList.push_back(
					new MateVectorByRN(vectorIndex, host, chr, begin, end));
			masterTmpFileList.push_back(std::string(""));
			sortHeap.push_back(tmpHeap);
			vectorIndex++;
		}
	}

	return 0;
}

int ExternalMemorySortManager::Dispatch() {

	if (THREAD_NUM > 1) {
		std::thread t[THREAD_NUM];

		int *workDone = new int[MateVectorList.size()];
		for (uint i = 0; i != MateVectorList.size(); ++i)
			workDone[i] = 0;

		for (int i = 0; i != THREAD_NUM; ++i)
			t[i] = std::thread(&ExternalMemorySortManager::ProcessVector, this,
					workDone);

		for (int i = 0; i != THREAD_NUM; ++i)
			t[i].join();

		delete[] workDone;
	} else {
		int *workDone = new int[MateVectorList.size()];
		for (uint i = 0; i != MateVectorList.size(); ++i)
			workDone[i] = 0;
		ProcessVector(workDone);
		delete[] workDone;
	}

	return 0;
}
int ExternalMemorySortManager::ProcessVector(int* workDone) {

	int prevIndex = -1;
	while (1) {
		uint nWorkDone = 0;
		int vectorIndex = 0;
		for (uint i = 0; i != MateVectorList.size(); ++i) {
			//mutex
			myLock.lock();
			if (workDone[i] == 0) {
				vectorIndex = i;
				workDone[i] = 1;
				myLock.unlock();
				break;
			} else {
				myLock.unlock();
				nWorkDone++;
			}
		}
		if (nWorkDone == MateVectorList.size())
			break;
		if (vectorIndex != prevIndex) {
			MateVectorList[vectorIndex]->Read(bamFile);
			miniMergeSort(MateVectorList[vectorIndex]);
			delete MateVectorList[vectorIndex];
			MateVectorList[vectorIndex] = NULL;
			prevIndex = vectorIndex;
		}
	}
	return 0;
}
int ExternalMemorySortManager::MergeSort(/*std::vector<std::string> tmpFileList*/) {//reading from signle tmp file(TODO:or heap), output to fianl file

	SamFileHeader tmpHeader;
	RECORD tmpRecord;
	tmpRecord.second = nullptr;
	std::vector<std::string> tmpFileList = masterTmpFileList;
	int nFile = tmpFileList.size();
	std::vector<SamFile*> tmpFileHandle;
	std::cerr << "Enter MergeSort Stage..." << std::endl;
	for (uint i = 0; i < tmpFileList.size(); ++i) {
		SamFile *fin = new SamFile(tmpFileList[i].c_str(), SamFile::READ,
				&tmpHeader);
		if (!fin->IsOpen()) {
			std::cerr << tmpFileList[i] << " file open failed...\n"
					<< std::endl;
			exit(EXIT_FAILURE);
		}

		if (fin->IsEOF()) //only contain header
		{
			fin->Close();
			delete fin;
			nFile--;
			continue;
		}

		tmpFileHandle.push_back(fin);
		tmpRecord.second = host->myPool.getRecord();
		fin->ReadRecord(tmpHeader, *tmpRecord.second);
		tmpRecord.first = tmpRecord.second->getReadName();
		masterHeap.push(std::make_pair(tmpRecord, fin));
	}

	if (nFile <= 0) {
		std::cerr << "Fatal error, no all.tmp file exists!\n" << std::endl;

		for (uint i = 0; i != tmpFileList.size(); ++i) {
			std::remove(tmpFileList[i].c_str());
		}
		exit(EXIT_FAILURE);
		//return 1; //no existing tmp files
	}
	int closedFile = 0;
	POPCORN tmpPop;
	SamFile* tmpFin = nullptr;
	int round = 0;
	while (closedFile < nFile) {
		round++;
		tmpPop = masterHeap.top();
		masterHeap.pop();
		tmpRecord = tmpPop.first;
		tmpFin = tmpPop.second;
		host->handlePairedRN(*tmpRecord.second); //ensure this is the only place call this function
		if (tmpFin->IsEOF()) {
			closedFile++;
			tmpFin->Close();
		} else {

			tmpRecord.second = host->myPool.getRecord();
			tmpFin->ReadRecord(tmpHeader, *tmpRecord.second);
			tmpRecord.first = tmpRecord.second->getReadName();
			masterHeap.push(std::make_pair(tmpRecord, tmpFin));
		}

	}
//	if(host->handlePairedRN::prevRec!=NULL)//we didnt pick out the record in local static variable in handlePairedRN functions
//	{
//		tmpVector.tmpFile.WriteRecord(tmpHeader,*tmpVector.prevRec); //writing bam1_t to tmp.all file
//		tmpVector.myPool->releaseRecord(tmpVector.prevRec);
//		tmpVector.prevRec=NULL;
//	}

	//std::cerr<<"How many pops in masterHeap?\t"<<masterHeap.size()<<std::endl;

	for (uint i = 0; i != tmpFileHandle.size(); ++i) {
		tmpFileHandle[i]->Close();
		delete tmpFileHandle[i];
	}

	for (uint i = 0; i != tmpFileList.size(); ++i) {
		std::remove(tmpFileList[i].c_str());
	}
	return 0;
}

int ExternalMemorySortManager::miniMergeSort(MateVectorByRN* tmpVector) { //reading from tmp file, ouput single to single tmp file (TODO:or heap)
	std::vector<std::string> tmpFileList = tmpVector->tmpFileNameList;

	RECORD tmpRecord;
	tmpRecord.second = nullptr;
	SamFileHeader tmpHeader;
	int vectorIndex = tmpVector->vectorIndex;
	std::cerr << "Enter " << vectorIndex << "th vector's miniMergeSort stage..."
			<< std::endl;
	masterTmpFileList[vectorIndex] = tmpVector->tmpFileNameBase
			+ std::to_string(vectorIndex) + "th.vector.all.tmp";
	tmpVector->tmpFile.OpenForWrite(masterTmpFileList[vectorIndex].c_str(),
			&tmpVector->tmpHeader);
	if (!tmpVector->tmpFile.IsOpen()) {
		std::cerr << "open " << masterTmpFileList[vectorIndex]
				<< " failed, check if there is enough space to create file!"
				<< std::endl;
		exit(EXIT_FAILURE);
	}

	int nFile = tmpFileList.size();
	std::vector<SamFile*> tmpFileHandle;
	for (uint i = 0; i < tmpFileList.size(); ++i) {

		SamFile * fin = new SamFile(tmpFileList[i].c_str(), SamFile::READ,
				&tmpHeader);

		if (!fin->IsOpen()) {
			std::cerr << tmpFileList[i] << " file open failed...\n"
					<< std::endl;
			exit(EXIT_FAILURE);
		}

		if (fin->IsEOF()) { //only contain header
			fin->Close();
			delete fin;
			nFile--;
			continue;
		}
		tmpFileHandle.push_back(fin);
		tmpRecord.second = tmpVector->myPool->getRecord();
		fin->ReadRecord(tmpHeader, *tmpRecord.second);
		tmpRecord.first = tmpRecord.second->getReadName();
		sortHeap[vectorIndex].push(POPCORN(tmpRecord, fin));
	}

	if (nFile <= 0) {

		for (uint i = 0; i != tmpFileList.size(); ++i) {
			std::remove(tmpFileList[i].c_str());
		}
		tmpVector->tmpFile.Close();
		tmpVector->clearMyPool();
		return 1; //no existing tmp files
	}

	int closedFile = 0;
	POPCORN tmpPop;
	SamFile* tmpFin = nullptr;
	while (closedFile < nFile) {

		tmpPop = sortHeap[vectorIndex].top();
		sortHeap[vectorIndex].pop();
		tmpRecord = tmpPop.first;
		if (tmpVector->prevRec == NULL) {
			tmpVector->prevRec = tmpRecord.second;
		} else {
			if (strcmp(tmpRecord.second->getReadName(),
					tmpVector->prevRec->getReadName()) == 0) {
				SamRecord* samRec = tmpRecord.second;
				if (SamFlag::isFirstFragment(samRec->getFlag())) {
//					if (SamFlag::isFirstFragment(
//							tmpVector.prevRec->getFlag())) {
//						std::cerr << "Both reads of " << samRec->getReadName()//<<"\tand\t"<<tmpVector.prevRec->getReadName()<<"\t"<<samRec<<"\t"<<tmpVector.prevRec
//								<< " are first fragment, so "
//								<< "splitting one to be in the 2nd fastq.\n";
//					}
					myLock.lock();

					host->writeFastQ(*samRec, host->myFirstFile,
							host->myFirstFileNameExt, tmpVector->myPool, false,
							host->myFirstRNExt.c_str());
					host->writeFastQ(*tmpVector->prevRec, host->mySecondFile,
							host->mySecondFileNameExt, tmpVector->myPool, false,
							host->mySecondRNExt.c_str());
					++host->myNumPairs;
					myLock.unlock();

				} else {
//					if (!SamFlag::isFirstFragment(
//							tmpVector.prevRec->getFlag())) {
//						std::cerr << "Neither read of " << samRec->getReadName()//<<"\tand\t"<<tmpVector.prevRec->getReadName()<<"\t"<<samRec<<"\t"<<tmpVector.prevRec
//								<< " are first fragment, so "
//								<< "splitting one to be in the 2nd fastq.\n";
//					}

					myLock.lock();

					host->writeFastQ(*tmpVector->prevRec, host->myFirstFile,
							host->myFirstFileNameExt, tmpVector->myPool, false,
							host->myFirstRNExt.c_str());
					host->writeFastQ(*samRec, host->mySecondFile,
							host->mySecondFileNameExt, tmpVector->myPool, false,
							host->mySecondRNExt.c_str());
					++host->myNumPairs;
					myLock.unlock();

				}
				// No previous record.
				tmpVector->prevRec = NULL;
			} else {
				tmpVector->tmpFile.WriteRecord(tmpHeader, *tmpVector->prevRec); //writing bam1_t to tmp.all file
				tmpVector->myPool->releaseRecord(tmpVector->prevRec);
				tmpVector->prevRec = tmpRecord.second;

			}
		}
		tmpFin = tmpPop.second;
		if (tmpFin->IsEOF()) {
			closedFile++;
			tmpFin->Close();
		} else {

			tmpRecord.second = tmpVector->myPool->getRecord();
			tmpFin->ReadRecord(tmpHeader, *tmpRecord.second);
			tmpRecord.first = tmpRecord.second->getReadName();
			sortHeap[vectorIndex].push(POPCORN(tmpRecord, tmpFin));
		}

	}
	if (tmpVector->prevRec != NULL) {
		tmpVector->tmpFile.WriteRecord(tmpHeader, *tmpVector->prevRec); //writing bam1_t to tmp.all file
		tmpVector->myPool->releaseRecord(tmpVector->prevRec);
		tmpVector->prevRec = NULL;
	}

	tmpVector->tmpFile.Close();
	tmpVector->clearMyPool();

	for (uint i = 0; i != tmpFileHandle.size(); ++i) {
		tmpFileHandle[i]->Close();
		delete tmpFileHandle[i];
	}
	for (uint i = 0; i != tmpFileList.size(); ++i) {
		std::remove(tmpFileList[i].c_str());
	}

	return 0;
}

