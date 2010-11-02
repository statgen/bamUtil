#ifndef __GCCONTENT_H__
#define __GCCONTENT_H__

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

#include <stdint.h>

#include "IntArray.h"
#include "StringHash.h"


#include "GenomeSequence.h"

class GCContent
{
 public:
  uint8_t *gcCount;
  uint8_t *atCount;
  uint64_t *statsVecByGC;
  double *statsVsGC;
  uint8_t *gcVec;
  uint32_t *gcContentVec;
  uint32_t length;
  int windowSize;
  GenomeSequence *genome;
  StringIntHash GCFreq;
  std::vector<bool> regionIndicator;
  //std::map<int, uint64_t> depthDist;
 public:
  GCContent();
  ~GCContent();
  
  void SetGenomeSequence(GenomeSequence *genomeSeq, int);
  void SetWindowSize(int);
  void ResetWindowSize(int);
  void SetZeroCount();
  void LoadRegions(String&, GenomeSequence &);
  void ReadGCContent(String);
  void CalcGCCount();
  void CalcGCFreq();
  void CalcStatsByGC(uint8_t *);
  void IncrementStatsByGC(uint32_t start, int len);
  void OutputFreq(String);
  void OutputGCContent(String & genome, int windowSize, String &gcContentFile, String &regionFile);
};


class ReadDepth
{
 public:
  uint8_t *depth;
  uint32_t length;
  //StringIntHash depthFreq;
  uint32_t *depthFreq;
  int LIM;
 public:
  ReadDepth();
  ReadDepth(uint32_t);
  ~ReadDepth();
  void AllocateMemory(uint32_t);
  void SetZeroCount();
  void SetDepthLimit(int);
  void IncrementCount(uint32_t);
  void IncrementCount(uint32_t, int);
  void IncrementCount(uint32_t, int,int);
  void CalcDepthFreq();
  void OutputFreq(String);
};
#endif
