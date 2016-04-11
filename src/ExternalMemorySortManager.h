/*
 * ExternalMemorySortManager.h
 *
 *  Created on: Feb 17, 2016
 *      Author: fanzhang
 */

#ifndef SRC_BAMUTIL_SRC_EXTERNALMEMORYSORTMANAGER_H_
#define SRC_BAMUTIL_SRC_EXTERNALMEMORYSORTMANAGER_H_

#include "MateVectorByRN.h"

class Bam2FastQ;

typedef std::pair<RECORD,SamFile*> POPCORN;//elements in heap
class PopcornComparison
{
  bool reverse;
public:
  PopcornComparison(const bool& revparam=false)
    {reverse=revparam;}
  bool operator() ( POPCORN& lhs,  POPCORN&rhs)const
  {
    if (reverse) return (std::string(lhs.first.second->getReadName()) > std::string(rhs.first.second->getReadName()));
    else return (std::string(lhs.first.second->getReadName()) < std::string(rhs.first.second->getReadName()));
  }
};

typedef std::priority_queue<POPCORN,std::vector<POPCORN>,PopcornComparison> HEAP;

class ExternalMemorySortManager {
public:
	ExternalMemorySortManager();
	virtual ~ExternalMemorySortManager();
	ExternalMemorySortManager(Bam2FastQ * tmpHost,std::string tmpBamFile,int nThread, int recordLimit,std::string bedFile);

	int FindAllBeginPointer();

	int Dispatch();

	int miniMergeSort(MateVectorByRN* tmpVector);

	int MergeSort(/*std::vector<std::string> tmpFileList*/);

	int ProcessVector(int* workDone);

	bool readRecord(IFILE fin, int32_t* blockSize, char* tmpBam);

	//void ReadFastq(std::ifstream * fin, FASTQ * tmp, std::string& line);

private:

	int THREAD_NUM;

	int RECORD_LIMIT;//number of records in each section

	std::string bamFile;

	std::string bedFile;

	//std::queue<SamRecord*> FIFO;//interface to caller

	HEAP masterHeap;

	std::vector<HEAP> sortHeap;//size is equal to the number of tmp files

	std::vector<MateVectorByRN*> MateVectorList;

	std::vector<std::string> masterTmpFileList;

	//std::vector<std::ifstream> fileHandleList;

	Bam2FastQ * host;//interface to final output

};

#endif /* SRC_BAMUTIL_SRC_EXTERNALMEMORYSORTMANAGER_H_ */
