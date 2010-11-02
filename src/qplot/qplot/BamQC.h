#ifndef __BAMQC_H__
#define __BAMQC_H__

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>

#include "StringArray.h"
#include "GenomeSequence.h"
#include "GenomeRegionSeqStats.h"
#include "QCStats.h"

#define MAXQ 50

class BamQC
{
 public:
  StringArray bamFiles;
  String GCInputFile;
  String label;
  String bamLabel;
  String lanes;
  std::map<int, int> lanes2Process;
  int size;
  GenomeRegionSeqStats regions;
  std::vector<bool> dbSNPIndicator;
  std::vector<bool> regionIndicator;
  std::vector<bool> genomePosCovered;
  GenomeSequence referencegenome;
  uint32_t refBaseNCount;
  ReadDepth depthVec;
  GCContent GC;
  QCStats *stats;
  int nRecords2Process;
  bool noDepth; // Indicator of whether to process depth dist
  bool noGC; // Indicator of whether to process GC content
  int page;

 public:
  BamQC();
  BamQC(StringArray &);
  BamQC(StringArray &, int);
  ~BamQC();
  void Init(StringArray &, int);
  void SetGCInputFile(String &in){ GCInputFile = in; };
  void SetQCStatsReferencePtr();
  void SetLabel(String &lb) { label = lb; }
  void SetBamLabels(String &lb) {bamLabel = lb; }
  void SetLanes2Process(String &);
  void LoadRegions(String &);
  void LoadGenomeSequence(String & refGenomeFile);
  void LoaddbSNP(String & dbSNPFile);
  void CalcNBaseCount();
  void SetNumRecords2Process(int n) { nRecords2Process=n; }
  void CalculateQCStats(QSamFlag &filter, double minMapQual);
  void OutputStats(String &);
  //void OutputPlotData(String &);
  void Plot(String &, FILE*);
  String GenRscript_EPSvsPhred_Plot();
  String GenRscript_EPSvsCycle_Plot();
  String GenRscript_DepthVsGC_Plot();
  String GenRscript_DepthDist_Plot();
  String GenRscript_BaseComp_Plot();
  String GenRscript_InsertSize_Plot();
  String GenRscript_DepthCoverage_Plot();
  String GenRscript_Q20_Plot();
  String GenRscript_GeneralStats_Plot();

  String GenRscript_EPSvsPhred_Data(int idx);
  String GenRscript_EPSvsCycle_Data(int idx);
  String GenRscript_PhredDist_Data(int idx);
  String GenRscript_DepthVsGC_Data(int idx);
  String GenRscript_DepthDist_Data(int idx);
  String GenRscript_BaseComp_Data(int idx);
  String GenRscript_InsertSize_Data(int idx);
  String GenRscript_CycleDist_Data(int idx);

/*  
  String GenPlotData_EPSvsPhred(int idx);
  String GenPlotData_EPSvsCycle(int idx);
  String GenPlotData_DepthVsGC(int idx);
  String GenPlotData_DepthDist(int idx);
  //String GenPlotData_BaseComp_Data(int idx);
*/
};

#endif
