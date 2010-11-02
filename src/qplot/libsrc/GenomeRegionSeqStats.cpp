#include "StringArray.h"
#include "GenomeRegionSeqStats.h"
#include "Error.h"
#include "SamFile.h"
#include "SamFlag.h"

#define MAXDP 65535

GenomeRegionSeqStats::GenomeRegionSeqStats()
{
  Init();
}
GenomeRegionSeqStats::GenomeRegionSeqStats(String &input)
{
  Init();
  LoadRegionList(input);
}
GenomeRegionSeqStats::~GenomeRegionSeqStats() {}

void GenomeRegionSeqStats::Init()
{
  minOverlapLen = 1; 
  contigFinishedCnt = 0;
  nReads = 0;
  nUnMapped = 0;
  nMapped2Targets = 0;
  nMappedOutTargets = 0;
}

void GenomeRegionSeqStats::CalcRegionStats(StringArray &bamFiles)
{
 for(int i=0; i<bamFiles.Length(); i++) 
   CalcRegionStats(bamFiles[i]);
}

void GenomeRegionSeqStats::CalcRegionStats(String &bamFile)
{
  SamFile sam;
  SamRecord samRecord;
  SamFileHeader samHeader;

  if(!sam.OpenForRead(bamFile.c_str()))
    error("Open BAM file %s failed!\n", bamFile.c_str());

  if(!sam.ReadHeader(samHeader)) {
      error("Read BAM file header %s failed!\n", bamFile.c_str());
  }
  
  String contigLabel;
  int start, end;
  Reset();
  while(sam.ReadRecord(samHeader, samRecord))
    {
      nReads++;
      if(samRecord.getFlag() & SamFlag::UNMAPPED) nUnMapped++;

      if(contigFinishedCnt>=contigs.size()) continue;

      CigarRoller cigar(samRecord.getCigar());

      int nonClipSequence = 0;

      if(cigar.size()!=0 && cigar[0].operation==Cigar::softClip)
          nonClipSequence = cigar[0].count;

      contigLabel = samRecord.getReferenceName();
      start = nonClipSequence + samRecord.get0BasedPosition();  // start is 0-based
      end = start + samRecord.getReadLength() - 1;
      if(UpdateRegionStats(contigLabel, start, end)) nMapped2Targets++;
    }
    CalcRegionReadCountInGCBins();
    CalcGroupReadCountInGCBins();
    std::cout << "Total reads : " << nReads << std::endl;
}

void GenomeRegionSeqStats::CalcClusters(StringArray &bamFiles, int minMapQuality)
{
 for(int i=0; i<bamFiles.Length(); i++) 
   CalcClusters(bamFiles[i], minMapQuality);
}

void GenomeRegionSeqStats::CalcClusters(String &bamFile, int minMapQuality)
{
  SamFile sam;
  SamRecord samRecord;
  SamFileHeader samHeader;

  if(!sam.OpenForRead(bamFile.c_str()))
    error("Open BAM file %s failed!\n", bamFile.c_str());

  if(!sam.ReadHeader(samHeader)) {
      error("Read BAM file header %s failed!\n", bamFile.c_str());
  }
  
  if(depth.size()==0) depth.resize(referencegenome.sequenceLength());
  
  String contigLabel;
  uint32_t start;
  uint32_t gstart;
  Reset();
  while(sam.ReadRecord(samHeader, samRecord))
    {
      nReads++;
      if(samRecord.getFlag() & SamFlag::UNMAPPED) nUnMapped++;

      if(samRecord.getMapQuality() < minMapQuality) continue;

      CigarRoller cigar(samRecord.getCigar());

      int nonClipSequence = 0;

      if(cigar.size()!=0 && cigar[0].operation==Cigar::softClip)
          nonClipSequence = cigar[0].count;

      contigLabel = samRecord.getReferenceName();
      start = nonClipSequence + samRecord.get0BasedPosition();  // start is 0-based

      gstart = referencegenome.getGenomePosition(contigLabel.c_str(), start);

      if(IsInRegions(contigLabel, start, start+samRecord.getReadLength())) continue;

      for(uint32_t i=gstart; i<gstart+samRecord.getReadLength(); i++)
       if(depth[i]<MAXDP)
        depth[i]++;
      nMappedOutTargets++;
    }
}

void GenomeRegionSeqStats::LoadReferenceGenome(String &reference)
{
{
  bool memoryMap = false;
  fprintf(stderr,"Loading reference... ");  
  referencegenome.setReferenceName(reference.c_str());
  referencegenome.useMemoryMap(memoryMap);
  if (referencegenome.open())
  {
    fprintf(stderr, "Failed to open reference index and is creating one...\n");
    if(referencegenome.create(false))
    	error("Creating reference index failed!\n");
    if(referencegenome.open())
    	error("Open newly created reference failed!\n");
   }
  fprintf(stderr, "DONE! Total sequence length %u\n", referencegenome.sequenceLength());
 }
}

void GenomeRegionSeqStats::LoadRegionList(String &inputList)
{
  FILE *in = fopen(inputList.c_str(), "r");
  if(in==NULL) error("Open region input file %s failed!\n", inputList.c_str());
  StringArray tokens;
  String buffer;
  while(!feof(in))
    {
      buffer.ReadLine(in);
      if (buffer.IsEmpty() || buffer[0] == '#') continue;
      tokens.ReplaceTokens(buffer);
      if(tokens.Length()<3)
	error("Too few columns: %s\n", buffer.c_str());
      
      String CSE = tokens[0]+":"+tokens[1]+":"+tokens[2];
      std::pair<int, int> start_end;
      start_end.first = tokens[1].AsInteger();
      start_end.second = tokens[2].AsInteger();
      if(start_end.first>=start_end.second) // positions are 0-based. Otherwise == is valid
      	error("Region end is equal or smaller than the start: %s!\n", buffer.c_str());      
      genomeRegions_lines[tokens[0]].push_back(buffer);
      genomeRegions[tokens[0]].push_back(start_end);
      genomeRegions_currentIndex[tokens[0]] = 0; 

      if(tokens.Length()>3) {
	groupStats[tokens[3]].segCount++;
	groupStats[tokens[3]].totalLen += (start_end.second - start_end.first);
	genomeRegionGroups[CSE].push_back(tokens[3]);
      }
    }
  
  fclose(in);
  
  // Chromosome info
  contigs.clear();
  std::map<String, vector<std::pair<int, int> > >::iterator p;
  for(p=genomeRegions.begin(); p!=genomeRegions.end(); p++)
    {
      contigs.push_back(p->first);
      for(unsigned int i=1; i<genomeRegions[p->first].size(); i++)
	if(genomeRegions[p->first][i].first<genomeRegions[p->first][i-1].first)
	  error("Input coordinates are not in order: %s %d %d!\n", p->first.c_str(),genomeRegions[p->first][i].first,genomeRegions[p->first][i].second);
    }
  // Group info such as gene names
  groups.clear();
  std::map<String, Stats>::iterator p2;
  for(p2=groupStats.begin(); p2!=groupStats.end(); p2++)
    groups.push_back(p2->first);
}

void GenomeRegionSeqStats::Reset()
{
  contigFinishedCnt = 0;
  for(unsigned int i=0; i<contigs.size(); i++)
    genomeRegions_currentIndex[contigs[i]] = 0;
}

void GenomeRegionSeqStats::CalcRegionGCContent()
{
  double gc;
  int start, end;
  genomeIndex_t gstart, gend;
  String CSE; //Chr_Start_End
  for(unsigned int i=0; i<contigs.size(); i++)
    {
      for(unsigned int j=0; j<genomeRegions[contigs[i]].size(); j++)
	{
	  start = genomeRegions[contigs[i]][j].first;
	  end   = genomeRegions[contigs[i]][j].second;
	  CSE = contigs[i]+":"+start+":"+end;
	  gstart = referencegenome.getGenomePosition(contigs[i].c_str(), start);
	  gend = gstart + (end-start);
	  int atCnt, gcCnt;
	  gc = CalcRegionGCContent(referencegenome,gstart, gend, atCnt, gcCnt);
	  genomeRegionStats[CSE].gcContent = gc;
	  genomeRegionStats[CSE].atCnt = atCnt;
	  genomeRegionStats[CSE].gcCnt = gcCnt;
	  
	  if(groups.size()>0) {
	    for(unsigned int g=0; g<genomeRegionGroups[CSE].size(); g++)
	    {
	     groupStats[genomeRegionGroups[CSE][g]].atCnt+=atCnt;
	     groupStats[genomeRegionGroups[CSE][g]].gcCnt+=gcCnt;
	    }
	  }
	}
    }
}

double GenomeRegionSeqStats::CalcRegionGCContent(GenomeSequence &ref, uint32_t start, uint32_t end, int &atCnt, int &gcCnt)
{
  atCnt=0;
  gcCnt=0;
  for(uint32_t i=start; i<end; i++)
    {
      if(toupper(ref[i])=='A' || toupper(ref[i])=='T') atCnt++;
      if(toupper(ref[i])=='C' || toupper(ref[i])=='G') gcCnt++;
    }
  double gcContent = 100*double(gcCnt)/(atCnt+gcCnt);
  return(gcContent);
}

void GenomeRegionSeqStats::CalcRegionReadCountInGCBins()
{
  if(referencegenome.sequenceLength()==0) return;
  double gc;
  int start, end;

  String CSE; //Chr_Start_End
  for(unsigned int i=0; i<contigs.size(); i++)
    {
      for(unsigned int j=0; j<genomeRegions[contigs[i]].size(); j++)
	{
	  start = genomeRegions[contigs[i]][j].first;
	  end   = genomeRegions[contigs[i]][j].second;
	  CSE = contigs[i]+":"+start+":"+end;
	  gc = genomeRegionStats[CSE].gcContent;
	  regionReadCountInGCBins[toGCBin(gc)] += genomeRegionStats[CSE].nReads_all;
	}
    }
}

void GenomeRegionSeqStats::CalcGroupReadCountInGCBins()
{
  double gc;
  for(unsigned int i=0; i<groups.size(); i++)
    {
      gc =  double(groupStats[groups[i]].gcCnt)/(groupStats[groups[i]].atCnt+groupStats[groups[i]].gcCnt)*100;
      groupStats[groups[i]].gcContent = gc;
      groupReadCountInGCBins[toGCBin(gc)] += groupStats[groups[i]].nReads_all;
    }
}

int GenomeRegionSeqStats::toGCBin(double gc)
{
 return(floor(gc+0.5));
}

// start and end are 0-based
bool GenomeRegionSeqStats::IsInRegions(String &chr, int start, int end)
{
  if(start>end){ fprintf(stderr, "WARNING: Start position %d is larger than the end position %d of %s!\n", start, end, chr.c_str()); return(false); }

  int32_t &currentIndex = genomeRegions_currentIndex[chr];

  if(currentIndex<0) return(false);
  if(currentIndex>=(int)genomeRegions[chr].size()){
    currentIndex=-1;
    contigFinishedCnt++;
    return(false);
  }

  vector< std::pair<int, int> > &region = genomeRegions[chr];
  if(end<=region[currentIndex].first) return(false);

  if(start>=region[currentIndex].first && end<=region[currentIndex].second) return(true);  
  if(start<=region[currentIndex].first && end>region[currentIndex].first) return(true);
  if(start<region[currentIndex].second && end>=region[currentIndex].second) return(true);
  if(start<=region[currentIndex].first && end>=region[currentIndex].second) return(true);

  while(start>=region[currentIndex++].second)
    {
      if(currentIndex>=(int)genomeRegions[chr].size()) return(false);

      if(start>=region[currentIndex].first && end<=region[currentIndex].second) return(true);  
      if(start<=region[currentIndex].first && end>region[currentIndex].first) return(true);
      if(start<region[currentIndex].second && end>=region[currentIndex].second) return(true);
      if(start<=region[currentIndex].first && end>=region[currentIndex].second) return(true);

      if(end<=region[currentIndex].first) return(false);
    }
  return(false);
}


bool GenomeRegionSeqStats::UpdateRegionStats(String &chr, int start, int end)
{
  int &currentIndex = genomeRegions_currentIndex[chr];
  if(currentIndex<0) return(false);
  if(currentIndex>=(int)genomeRegions[chr].size()){
    currentIndex=-1;
    contigFinishedCnt++;
    return(false);
  }
  vector< std::pair<int, int> > &region = genomeRegions[chr];

  if(end<region[currentIndex].first) return(false);

  bool overlap = PartialUpdate(chr, start, end);

  if(overlap==true) return(true);

  // No overlap found and find the next overlapping region
  // If found, update the stats and adjust index
  // If not, adjust the index for next sequence
  while(start>=region[currentIndex++].second)
    {
      if(currentIndex>=(int)genomeRegions[chr].size()) return(false);
      if(PartialUpdate(chr, start, end))  return(true);
      if(end<=region[currentIndex].first) return(false);
    }
  return(false);
}

// If two or more regions overlap, a single read mapped to the overlapping regions
// will be counted to ALL of the overlapping regions
// Note: start is 0-based
bool GenomeRegionSeqStats::PartialUpdate(String &chr, int start, int end)
{
  int &currentIndex = genomeRegions_currentIndex[chr];
  vector< std::pair<int, int> > &region = genomeRegions[chr];
  
  String CSE;
  vector<String> group;
  vector<int> counter;
  bool overlap = false;
  for(unsigned int i=currentIndex; i<region.size(); i++)
    {
      //  The second condition handles the situation where a read spans >1 region.
      if(region[i].first>region[currentIndex].second && region[i].first>end) break;
      
      CSE = chr+":"+region[i].first+":"+region[i].second;
      group = genomeRegionGroups[CSE];
       counter.resize(group.size(), 0); 
       
       // completely mapped to a region 
       //    ****** read 
       // --------------- region 
       if(start>=region[i].first && end<=region[i].second)
	   {
	    genomeRegionStats[CSE].nReads_complete++; 
	    genomeRegionStats[CSE].nReads_all++; 
	    genomeRegionStats[CSE].nBases += (end-start); //positions are 0-based overlap=true; 
	    if(group.size()>0) {
	  	 for(unsigned int g=0; g<group.size(); g++) 
	  	 { 
	  	  if((++counter[g])>1) continue;
	     groupStats[group[g]].nReads_complete++;
	     groupStats[group[g]].nReads_all++;
	     groupStats[group[g]].nBases += (end-start); //positions are 0-based
	     }
	     }
	   continue;
	  }

      // completely cover the entire region
      // ******************* read
      //    ---------        region
      if(start<=region[i].first  && end>=region[i].second)
	{
	  genomeRegionStats[CSE].nReads_all++;
	  genomeRegionStats[CSE].nReads_complete++;
	  genomeRegionStats[CSE].nBases += (region[i].second-region[i].first); // positions are 0-based
	  overlap=true;
	  if(group.size()>0) {
	    for(unsigned int g=0; g<group.size(); g++)
	    {
	  	 if((++counter[g])>1) continue;
	    groupStats[group[g]].nReads_complete++;
	    groupStats[group[g]].nReads_all++;
	    groupStats[group[g]].nBases += (region[i].second-region[i].first); //positions are 0-based
	    }
	  }
	  continue;
	}
 
      // partial mapping left
      //  *************                          read
      //       --------------------------        region
      if(start<region[i].first && end>(region[i].first+minOverlapLen-1))
	{
	  genomeRegionStats[CSE].nReads_all++;
	  genomeRegionStats[CSE].nBases += (end-region[i].first); // positions are 0-based
	  overlap=true;
	  if(group.size()>0) {
	  for(unsigned int g=0; g<group.size(); g++)
	  {
	  	 if((++counter[g])>1) continue;
	    groupStats[group[g]].nReads_all++;
	    groupStats[group[g]].nBases += (end-region[i].first); //positions are 0-based
	  }
	  }
	  continue;
	}
	
	// printf("%d\t%d\n", region[i].second, minOverlapLen);
      // partial mapping right
      //              ****************      read
      //  -----------------------           region
      if(start<(region[i].second-minOverlapLen+1) && end>region[i].second) 
	{
	  genomeRegionStats[CSE].nReads_all++;
	  genomeRegionStats[CSE].nBases += (region[i].second-start); // positions are 0-based
	  overlap=true;
	  if(group.size()>0) {
	  for(unsigned int g=0; g<group.size(); g++)
	  {
	  	 if((++counter[g])>1) continue;
	    groupStats[group[g]].nReads_all++;
	    groupStats[group[g]].nBases += (region[i].second-start); //positions are 0-based
	  }
	  }
	  continue;
	}
    }
  return(overlap);
}

// One read is counted ONLY ONCE. If two regions overlap, reads mapped to both regions
// will be counted to the region with samller start position
bool GenomeRegionSeqStats::PartialUpdate_Unique(String &chr, int start, int end)
{
  int &currentIndex = genomeRegions_currentIndex[chr];
  vector< std::pair<int, int> > &region = genomeRegions[chr];
  String CSE = chr+":"+region[currentIndex].first+":"+region[currentIndex].second;
  vector<String> group = genomeRegionGroups[CSE];
  vector<int> counter(group.size(), 0);
  //  completely mapped to a region
  //    ******            read
  //  ---------------     region
  if(start>=region[currentIndex].first && end<=region[currentIndex].second) 
    {
      genomeRegionStats[CSE].nReads_complete++;
      genomeRegionStats[CSE].nReads_all++;
      genomeRegionStats[CSE].nBases += (end-start);
      if(group.size()>0) {
  for(unsigned int g=0; g<group.size();g++)
  {
	  	 if((++counter[g])>1) continue;
	groupStats[group[g]].nReads_complete++;
	groupStats[group[g]].nReads_all++;
	groupStats[group[g]].nBases += (end-start); //positions are 0-based
	}
      }
      return(true);
    }
  
  // completely cover the entire region
  // ******************* read
  //    ---------        region
  if(start<=region[currentIndex].first  && end>=region[currentIndex].second)
    {
      genomeRegionStats[CSE].nReads_all++;
      genomeRegionStats[CSE].nReads_complete++;
      genomeRegionStats[CSE].nBases += (region[currentIndex].second-region[currentIndex].first); // positions are 0-based
      if(group.size()>0) {
   for(unsigned int g=0; g<group.size(); g++)
   {
	  	 if((++counter[g])>1) continue;
	groupStats[group[g]].nReads_complete++;
	groupStats[group[g]].nReads_all++;
	groupStats[group[g]].nBases += (region[currentIndex].second-region[currentIndex].first); //positions are 0-based
	}
      }
      return(true);
    }
  // partial mapping left
  //  *************                          read
  //       --------------------------        region
  if(start<region[currentIndex].first && end>(region[currentIndex].first+minOverlapLen-1))
    {
      genomeRegionStats[CSE].nReads_all++;
      genomeRegionStats[CSE].nBases += (end-region[currentIndex].first);
      if(group.size()>0) {
         for(unsigned int g=0; g<group.size(); g++)
         {
	  	 if((++counter[g])>1) continue;
	groupStats[group[g]].nReads_all++;
	groupStats[group[g]].nBases += (end-region[currentIndex].first); //positions are 0-based	
  }
      }
      return(true);
    }
  // partial mapping right
  //              ****************      read
  //  -----------------------           region
  if(start<(region[currentIndex].second-minOverlapLen+1) && end>region[currentIndex].second) 
    {
      genomeRegionStats[CSE].nReads_all++;
      genomeRegionStats[CSE].nBases += (region[currentIndex].second-start); // positions are 0-based
      if(group.size()>0) {
         for(unsigned int g=0; g<group.size(); g++){
	  	 if((++counter[g])>1) continue;
	groupStats[group[g]].nReads_all++;
	groupStats[group[g]].nBases += (region[currentIndex].second-start); //positions are 0-based
	}
      }
      return(true);
    }
  return(false);
}

/*
void  GenomeRegionSeqStats::CalcReadsMapped2Targets()
{
 nReadsMapped2Targets = 0;
 for(int i=0; i<contigs.size(); i++)
  {
   for(int j=0; j<genomeRegions[contigs[i]].size(); j++)
    {
     start = genomeRegions[contigs[i]][j].first;
     end   = genomeRegions[contigs[i]][j].second;
     String CSE = contigs[i]+":"+start+":"+end;
     int regionLen = end - start; // poistions are 0-based
     nReadsMapped2Targets += genomeRegionStats[CSE].nReads_all;
    }
  }
}
*/

void GenomeRegionSeqStats::OutputRegionStats(String &outFile)
{
  FILE *fh;
  if(outFile.Length()==0) fh = stdout;
  else fh = fopen(outFile.c_str(), "w");
  if(fh==NULL)
    error("Open out file %s failed!\n", outFile.c_str());

  int start, end;
  fprintf(fh, "##Total reads: %lu\n##Mapped reads: %u\n##Mapped rate(%%): %.2f\n##Mapped2targets reads: %lu\n##Mapped2Targets rate(%%): %.2f\n", nReads, (unsigned int) (nReads-nUnMapped), double(nReads-nUnMapped)/nReads*100, nMapped2Targets, double(nMapped2Targets)/nReads*100);
  fprintf(fh,"#Chr\tStart\tEnd");
  if(groupStats.size()>0)
    fprintf(fh, "\tGene/Group");
  fprintf(fh, "\tLen\tReadCnt\tRPM\tRPKM");
  if(referencegenome.sequenceLength()>0)
    fprintf(fh, "\t%%GC\tnReadGCBin\n");

  for(unsigned int i=0; i<contigs.size(); i++)
    {
      for(unsigned int j=0; j<genomeRegions[contigs[i]].size(); j++)
	{
	  start = genomeRegions[contigs[i]][j].first;
	  end   = genomeRegions[contigs[i]][j].second;
	  String CSE = contigs[i]+":"+start+":"+end;
	  int regionLen = end - start; // poistions are 0-based
	  double RPM = double(genomeRegionStats[CSE].nReads_all)/(double(nMapped2Targets)/1000000);
	  double RPKM = RPM/regionLen*1000;
	  //fprintf(fh, "%s\t%d\t%u\t%u\t%.2f", genomeRegions_lines[contigs[i]][j].c_str(), end-start, genomeRegionStats[CSE].nReads_all, genomeRegionStats[CSE].nReads_complete, double(genomeRegionStats[CSE].nBases)/regionLen);
	  fprintf(fh, "%s\t%d\t%u\t%.3f\t%.3f\t%.3f", genomeRegions_lines[contigs[i]][j].c_str(), end-start, (unsigned int) genomeRegionStats[CSE].nReads_all, double(genomeRegionStats[CSE].nBases)/regionLen, RPM, RPKM);
	  if(referencegenome.sequenceLength()>0)
	  	fprintf(fh, "\t%.2f\t%u", genomeRegionStats[CSE].gcContent, (unsigned int) regionReadCountInGCBins[toGCBin(genomeRegionStats[CSE].gcContent)]);
	  /*
	  if(group.Length()>0)
	    fprintf(fh, "\t%u\t%u\t%.2f\t%d\t%d", groupStats[group].nReads_all, groupStats[group].nReads_complete, double(groupStats[group].nBases)/groupStats[group].totalLen, groupStats[group].segCount, groupStats[group].totalLen);
	  */
	   fprintf(fh, "\n");
	  }
    }
  fclose(fh);
}

void GenomeRegionSeqStats::OutputGroupStats(String &outFile)
{
  if(groups.size()==0) { fprintf(stderr, "WARNING: Gene/Group (column 4) is not in the input regions and no stats output for genes/groups!\n"); return; }
  
  FILE *fh = fopen(outFile.c_str(), "w");
  if(fh==NULL) error("Open group outfile %s failed!\n", outFile.c_str());
  
  fprintf(fh, "##Total reads: %lu\n##Mapped reads: %u\n##Mapped rate(%%): %.2f\n##Mapped2targets reads: %lu\n##Mapped2Targets rate(%%): %.2f\n", nReads, (unsigned int) (nReads-nUnMapped), double(nReads-nUnMapped)/nReads*100, nMapped2Targets, double(nMapped2Targets)/nReads*100);
  fprintf(fh, "#Name\tLen\tReadCnt\tRPM\tRPKM");
  if(referencegenome.sequenceLength()>0)
    fprintf(fh, "\tGC%%\tnReadGCBin");
  fprintf(fh, "\n");
  
  //if(nReadsMapped2Targets==0) CalcReadsMapped2Targets();

  double RPM, RPKM;
  for(unsigned int i=0; i<groups.size(); i++)
    {
      RPM = double(groupStats[groups[i]].nReads_all)/(double(nMapped2Targets)/1000000);
      RPKM = RPM/groupStats[groups[i]].totalLen*1000;
      double gcContent = groupStats[groups[i]].gcContent;
 
      fprintf(fh, "%s\t%d\t%u\t%.3f\t%.3f", groups[i].c_str(), groupStats[groups[i]].totalLen, (unsigned int) groupStats[groups[i]].nReads_all, RPM, RPKM);
      if(referencegenome.sequenceLength()>0)
 	     fprintf(fh, "\t%.2f\t%u", gcContent, (unsigned int) groupReadCountInGCBins[toGCBin(gcContent)]);
      fprintf(fh, "\n");
    }
  fclose(fh);
}

void GenomeRegionSeqStats::OutputClusters(String &outFile, int minDepth, double minAvgDepth, int minClusterSize)
{
 FILE *fh = fopen(outFile.c_str(), "w");
 if(fh==NULL) error("Open cluster output file %s failed!\n", outFile.c_str());
 
 bool inCluster = false;
 int winSize = 0;
 double avgDepth = 0;
 int totalDepth = 0;
 uint32_t cstart, cend; //cluster start and cluster end
 for(uint32_t i=0; i<depth.size(); i++)
 {
  if(depth[i]>=minDepth)
  {
   if(inCluster==false)
   {
    inCluster = true;
    cstart = i;
   }
   winSize++;
   totalDepth += depth[i]; 
  }
  else
  {
    if(inCluster==false) continue;

    //std::cout<<totalDepth<<" "<<winSize<<" "<<inCluster<<":"<<minClusterSize<<":"<<avgDepth<<":"<<minAvgDepth<<std::endl;
    avgDepth = double(totalDepth)/winSize;
    if(winSize<minClusterSize || avgDepth<minAvgDepth) 
    {
      winSize = 0;
      totalDepth = 0;
      inCluster = false;
      continue;
     }
    winSize = 0;
    totalDepth = 0;
    inCluster = false;
     cend = i-1;
    int atCnt, gcCnt;
    double gc = CalcRegionGCContent(referencegenome, cstart, cend, atCnt, gcCnt);

    String chr; int start, end;
    referencegenome.getChromosomeAndIndex(chr, cstart);

    StringArray tokens;
    tokens.ReplaceTokens(chr, ":");
	 chr = tokens[0];
	 start = tokens[1].AsInteger();
	 end = start + (cend-cstart);

    fprintf(fh, "%s\t%u\t%u\t%.2f\t%.2f\n", chr.c_str(), start, end, avgDepth, gc);
  }
 }
 fclose(fh);
}

void GenomeRegionSeqStats::OutputRegionGCContent(String &outFile)
{
  FILE *fh;
  if(outFile.Length()==0) fh = stdout;
  else fh = fopen(outFile.c_str(), "w");
  if(fh==NULL)
    error("Open out file %s failed!\n", outFile.c_str());
  
  fprintf(fh,"#Chr\tStart\tEnd");
  if(groupStats.size()>0)
    fprintf(fh, "\tGene");
  fprintf(fh, "\tLen\tAT_Cnt\tGC_Cnt\t%%GC\n");
  
  int start, end;
  for(unsigned int i=0; i<contigs.size(); i++)
    {
      for(unsigned int j=0; j<genomeRegions[contigs[i]].size(); j++)
	{
	  start = genomeRegions[contigs[i]][j].first;
	  end   = genomeRegions[contigs[i]][j].second;
	  String CSE = contigs[i]+":"+start+":"+end;
	  fprintf(fh, "%s\t%d\t%u\t\%u\t%.2f\n", genomeRegions_lines[contigs[i]][j].c_str(), end-start, genomeRegionStats[CSE].atCnt, genomeRegionStats[CSE].gcCnt, genomeRegionStats[CSE].gcContent);  
	}
    }
  fclose(fh);
}

void GenomeRegionSeqStats::SetRegionIndicator(GenomeSequence &referencegenome)
{
 if(referencegenome.sequenceLength()==0) return;
 regionIndicator.resize(referencegenome.sequenceLength());
 genomeIndex_t startGenomeIndex, endGenomeIndex;
 
 int start,end;
 for(unsigned int i=0; i<contigs.size(); i++)
   {
    for(unsigned int j=0; j<genomeRegions[contigs[i]].size(); j++)
       {
         start = genomeRegions[contigs[i]][j].first;
         end   = genomeRegions[contigs[i]][j].second;
         startGenomeIndex = referencegenome.getGenomePosition(contigs[i].c_str(), start);
         endGenomeIndex = startGenomeIndex+(end-start);
         if(startGenomeIndex >= regionIndicator.size() ) {
           continue;
         }
         for(uint32_t k=startGenomeIndex; k<startGenomeIndex+(end-start); k++)
            regionIndicator[k] = true;
       }
     }
}
