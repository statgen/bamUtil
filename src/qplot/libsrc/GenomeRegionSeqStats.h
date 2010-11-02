#ifndef __GENOMEREGIONSEQSTATS_H__
#define __GENOMEREGIONSEQSTATS_H__


#include "GenomeSequence.h"
#include  "StringArray.h"
#include "InputFile.h"
#include <map>
#include <cmath>

class Stats
{
 public:
  uint64_t nReads_all;
  uint64_t nReads_complete;
  uint64_t nBases;
  double meanDepth;
  int segCount;
  uint32_t totalLen;
  uint32_t atCnt;
  uint32_t gcCnt;
  double gcContent;
  
 public:
  Stats(){
    nReads_all = 0;
    nReads_complete = 0;
    nBases = 0;
    segCount = 0;
    totalLen = 0;
    atCnt=gcCnt=0;
    gcContent = 0.0;
  }
};

class GenomeRegionSeqStats
{
public:
  GenomeSequence referencegenome;
  std::vector<String> contigs;
  std::map<String, std::vector< std::pair<int, int> > > genomeRegions;
  std::map<String, std::vector<String> > genomeRegions_lines;
  std::map<String, std::vector<String> > genomeRegionGroups;
  std::map<String, Stats> genomeRegionStats;
  std::map<String, int>   genomeRegions_currentIndex;
  std::map<String, Stats> groupStats;
  std::map<int, uint64_t> regionReadCountInGCBins;
  std::map<int, uint64_t> groupReadCountInGCBins;
  std::vector<String> groups;
  std::vector<uint16_t> depth;
  std::vector<uint32_t> regionIndicator;
  int minOverlapLen;
  unsigned int contigFinishedCnt;
  uint64_t nReads;
  uint64_t nUnMapped;
  uint64_t nMapped2Targets;
  uint64_t nMappedOutTargets;
public:
  GenomeRegionSeqStats();
  GenomeRegionSeqStats(String &);
  ~GenomeRegionSeqStats();
  void Init();
  void SetMinOverlapLen(int overlapLen) { minOverlapLen = overlapLen; }
  void LoadRegionList(String &);
  void LoadReferenceGenome(String &);
  void SetRegionIndicator(GenomeSequence &);
  void CalcRegionGCContent();
  double CalcRegionGCContent(GenomeSequence &ref, uint32_t, uint32_t, int &, int &);
  void CalcRegionReadCountInGCBins();
  void CalcGroupReadCountInGCBins();
  void CalcClusters(String &, int);
  void CalcClusters(StringArray &, int);
  void Reset();
  bool IsInRegions(String &, int start, int end);
  bool PartialUpdate(String &chr, int start, int end);
  bool PartialUpdate_Unique(String &chr, int start, int end);
  bool UpdateRegionStats(String &, int start, int end);
  void UpdateDepth(int start, int end);
  void OutputRegionStats(String &outFile);
  void OutputGroupStats(String &outFile);
  void OutputRegionGCContent(String &outFile);
  void OutputClusters(String &, int, double, int);
  void CalcRegionStats(String &bamFiles);
  void CalcRegionStats(StringArray &bamFiles);
  int toGCBin(double);
};

#endif

