#include "QCStats.h"
#include "BamQC.h"

#include "Parameters.h"

int main(int argc, char *argv[])
{
  bool unpaired = false;
  bool read1 = false;
  bool read2 = false;
  bool paired = false;
  bool keepDup = false;
  bool keepQCFail = false;
  double minMapQuality = 0;
  int nRecords = -1;

  String reference = "/data/local/ref/karma.ref/human.g1k.v37.umfa";
  String dbSNPFile = "/home/bingshan/data/db/dbSNP/dbSNP130.UCSC.coordinates.tbl";
  String gcContentFile = "/share/swg/bingshan/db/gccontent/human.g1k.w100.gc";
  String regions;
  String gcContentFile_create;
  int windowSize = 100;

  ParameterList pl;

  String statsFile;
  String plotFile;
  String RcodeFile;
  String label;
  String bamLabel;
  String lanes;
    
  bool noGC = false;
  bool noDepth = false;
  int page = 2;
    
  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("References")
      LONG_STRINGPARAMETER("reference",&reference)
      LONG_STRINGPARAMETER("dbsnp", &dbSNPFile)
      LONG_STRINGPARAMETER("gccontent", &gcContentFile)
     LONG_PARAMETER_GROUP("Create gcContent file")
       LONG_STRINGPARAMETER("create_gc",&gcContentFile_create)
       LONG_INTPARAMETER("winsize", &windowSize)
     LONG_PARAMETER_GROUP("Region list")
       LONG_STRINGPARAMETER("regions", &regions)
    LONG_PARAMETER_GROUP("Flag filters")
      LONG_PARAMETER("read1_skip", &read1)
      LONG_PARAMETER("read2_skip", &read2)
      LONG_PARAMETER("paired_skip", &paired)
      LONG_PARAMETER("unpaired_skip", &unpaired)
    LONG_PARAMETER_GROUP("Dup and QCFail")
      LONG_PARAMETER("dup_keep", &keepDup)
      LONG_PARAMETER("qcfail_keep", &keepQCFail)
    LONG_PARAMETER_GROUP("Mapping filters")
      LONG_DOUBLEPARAMETER("minMapQuality", &minMapQuality)
    LONG_PARAMETER_GROUP("Records to process") 
      LONG_INTPARAMETER("first_n_record", &nRecords)
    LONG_PARAMETER_GROUP("Lanes to process") 
      LONG_STRINGPARAMETER("lanes", &lanes)   
    LONG_PARAMETER_GROUP("Output files")
      LONG_STRINGPARAMETER("plot", &plotFile)
      LONG_STRINGPARAMETER("stats", &statsFile)
      LONG_STRINGPARAMETER("Rcode", &RcodeFile)
    LONG_PARAMETER_GROUP("Plot labels")
      LONG_STRINGPARAMETER("label", &label)
      LONG_STRINGPARAMETER("bamLabel", &bamLabel)
   END_LONG_PARAMETERS();
  

  pl.Add(new LongParameters("\n", longParameters));
  
  StringArray bamFiles;

  int in = pl.ReadWithTrailer(argc, argv, 1);
  for (int i=in+1; i<argc; i++){
    bamFiles.Push(argv[i]);
  }
  
  pl.Status();

  if(reference.Length()==0)
   error("Reference not provided!\n");

  if(gcContentFile.Length()==0 && gcContentFile_create.Length()==0)
  {
     error("No GC content file provided! You can create one based on window size of N by command option \n\t --create_gc gcContentFile --winsize N\n\n");
    }
         
  if(gcContentFile_create.Length()>0)
   {
     GCContent GC;
     fprintf(stderr, "Creating GC content file...\n");
     GC.OutputGCContent(reference, windowSize, gcContentFile_create, regions);
     return(0);
    }

  if(bamFiles.Length()==0)
    error("No SAM/BAM files provided!\n");


 fprintf(stderr, "The following files are to be processed...\n\n");  
  for(int i=0; i<bamFiles.Length();i++)
    fprintf(stderr, "%s\n", bamFiles[i].c_str());
 fprintf(stderr, "\n");
  
  if(plotFile.Length()==0)
    warning("No plot will be generated!\n");

 if(bamLabel.Length()>0) {
  StringArray bamLabelArray;
  bamLabelArray.ReplaceTokens(bamLabel, ",");
  if(bamLabelArray.Length()<bamFiles.Length())
  	error("BAM/SAM file number larger than lable number!\n");
  if(bamLabelArray.Length()>bamFiles.Length())
  	warning("BAM/SAM file number smaller than lable number and extra lables ignored!\n");
 }

  FILE *RCODE=NULL;
  FILE *pf = NULL;
  FILE *STATSFH = NULL;

  if(RcodeFile.Length()>0){
    RCODE = fopen(RcodeFile.c_str(), "w");
    if(RCODE==NULL)
      error("Open Rcode file for output failed!\n", RcodeFile.c_str());
  }

  if(plotFile.Length()>0)
  {
   pf = popen("Rscript --vanilla -", "w");
   if(pf==NULL)
    error("Open Rscript failed!\n", plotFile.c_str());
  }

  if(statsFile.Length()>0)
  {
   STATSFH = fopen(statsFile.c_str(), "w");
   if(STATSFH==NULL)
    error("Open stats file %s failed!\n", statsFile.c_str());
   fclose(STATSFH);
  }

  QSamFlag filter;

  if(paired && unpaired)
    warning("The filter --unpaired overrides --paired\n");
  
  if(unpaired) paired=true;
  filter.SetRead1(read1);
  filter.SetRead2(read2);
  filter.SetPaired(paired);
  filter.SetUnPaired(unpaired);
  filter.SetDuplicate(!keepDup);
  filter.SetQCFail(!keepQCFail);

  BamQC qc(bamFiles);

  qc.noGC = noGC;
  qc.noDepth = noDepth;
  qc.page = page;

  qc.SetLanes2Process(lanes);
  qc.SetNumRecords2Process(nRecords);
  qc.SetGCInputFile(gcContentFile);
  qc.SetLabel(label);
  qc.SetBamLabels(bamLabel);
  qc.LoadGenomeSequence(reference);
  qc.LoadRegions(regions);
  qc.LoaddbSNP(dbSNPFile);
  qc.CalculateQCStats(filter, minMapQuality);
  qc.OutputStats(statsFile);

  if(RcodeFile.Length()>0){
   qc.Plot(plotFile, RCODE);
  }

  if(plotFile.Length()>0) 
  {
   qc.Plot(plotFile, pf);
  }

  return(0);

} //END of main
