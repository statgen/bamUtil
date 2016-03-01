/*
 * MateVectorByRN.h
 *
 *  Created on: Feb 16, 2016
 *      Author: fanzhang
 */

#ifndef SRC_BAMUTIL_SRC_MATEVECTORBYRN_H_
#define SRC_BAMUTIL_SRC_MATEVECTORBYRN_H_

#include "SamFile.h"
#include <map>
#include <string>
#include <fstream>

#include "SamRecordPool.h"
#include "SamFlag.h"

typedef std::pair<std::string, SamRecord*> RECORD;
typedef std::vector<RECORD> MATE_VECTOR;

class Bam2FastQ;

struct {
	bool operator()(const RECORD& lhs, const RECORD&rhs) const {
		return (lhs.first > rhs.first); //sort it descendingly, because heap pop out the max
	}
} RecordComparison;

class MateVectorByRN {
public:
	MateVectorByRN(int vectorIndex, Bam2FastQ* host, std::string chr,
			long begin, long end);
	MateVectorByRN();
	virtual ~MateVectorByRN();

	/// Return the mate of the specified record, or NULL if mate is not found.
	/// The mate is removed from this map if it is found.
	//SamRecord* GetMate(SamRecord& record);

	int Read(const std::string& bamFile);

	bool IsInTheRegion(SamRecord* record);

	//SamRecord* SetRecordFromBamT(bam1_t * buffer);
	/// Add the specified record to this mate vector
	int Add(SamRecord* record);
	/// sort record for buffer in use, we can consider binary search plus insertion when speed is slow
	void Sort();

	int DumpMateVector();

	int HandlePairedRN(SamRecord& Record, bool ready2Dump);
	//int ReadNameReduction();//optional
	//clear buffer, reset shrink_limit
	void ClearWorkingBuffer();

	inline void FlipWorkingBuffer() {
		bufferInUse = 1 - bufferInUse;
	}

	void clearMyPool();

	MATE_VECTOR myMateBuffer[2];//use these two alternatively, to effectively free memory

	int bufferInUse;

	int myNumPairs;

	int myNumMateFailures;

	long count;

	long SHRINK_LIMIT;

	long MAX_LIMIT;

	std::string chr;

	long begin;

	long end;

	int vectorIndex;

	SamRecord* prevRec;

	std::string prevRN;
	SamRecordPool* myPool;

	std::string tmpFileName;
	std::string tmpFileNameBase;
	SamFileHeader tmpHeader;
	SamFile tmpFile;
	int tmpFileIndex;
	std::vector<std::string> tmpFileNameList;

	std::string FirstFileNameExt;
	std::string SecondFileNameExt;
	std::string UnpairedFileNameExt;

	Bam2FastQ* host;	//interface to final output
private:
};

#endif /* SRC_BAMUTIL_SRC_MATEVECTORBYRN_H_ */
