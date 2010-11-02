#include "GCContent.h"
#include "Error.h"
#include "InputFile.h"
#include <cmath>
#include <fstream>

GCContent::GCContent()
{
 length=0;
 gcCount=NULL;
 statsVecByGC = NULL;
}

void GCContent::SetZeroCount()
{
  for(uint32_t i=0; i<length; i++){
     gcCount[i] = 0;
    //atCount[i] = 0;
  for(int i=0; i<=windowSize; i++)
     statsVecByGC[i] = 0;
        
   }
}

void GCContent::SetGenomeSequence(GenomeSequence *genomeSeq, int ws) {
  if(ws > 255)
    error("Window size can not be larger than 255!\n");
  windowSize = ws;
  statsVecByGC = new uint64_t[windowSize+1];
  for(int i=0; i<=windowSize; i++)
   statsVecByGC[i] = 0;
   
  genome = genomeSeq;
  gcCount = new uint8_t[genome->sequenceLength()];
  //atCount = new uint8_t[len];
  if(gcCount==NULL)
      error("Allocating memory for gcCount failed!\n");
  length = genome->sequenceLength();

   gcVec = new uint8_t[101];
   gcContentVec = new uint32_t[101];
     
}

void GCContent::SetWindowSize(int newSize)
{

}

void GCContent::ResetWindowSize(int newSize)
{
  if(newSize > 255)
    error("Window size can not be larger than 255!\n");
  windowSize = newSize;
  if(statsVecByGC!=NULL)
	 delete [] statsVecByGC;
  statsVecByGC = new uint64_t[windowSize+1];
  for(int i=0; i<=windowSize; i++)
   statsVecByGC[i] = 0;
}

GCContent::~GCContent(){ 
	if(length>0) {
	delete [] gcCount;
	//delete [] atCount;
	 delete [] statsVecByGC;
	 delete [] gcVec;
	 delete [] gcContentVec;
	    
	}
};

void GCContent::LoadRegions(String & regionsFile, GenomeSequence &genome)
{
  if(regionsFile.Length()==0) return;
  if(genome.sequenceLength()==0) error("No reference genome loaded!\n");
    
  IFILE fhRegions;
  fhRegions = ifopen(regionsFile.c_str(),"r");
  if(fhRegions==NULL)
    error("Open regions file %s failed!\n", regionsFile.c_str());
  
  regionIndicator.resize(genome.sequenceLength());
  
  StringArray tokens;
  String buffer;
  int len;
  
  fprintf(stderr, "Loading region list...");
  
  while (!ifeof(fhRegions)){
    buffer.ReadLine(fhRegions);
    if (buffer.IsEmpty() || buffer[0] == '#') continue;
    
    tokens.AddTokens(buffer, WHITESPACE);
    if(tokens.Length() < 3) continue;
    
    genomeIndex_t startGenomeIndex = 0;   
    int chromosome = -1;
    int chromosomeIndex = tokens[1].AsInteger();
    
     startGenomeIndex = genome.getGenomePosition(tokens[0].c_str(), chromosomeIndex);
    
    if(startGenomeIndex >= regionIndicator.size() ) {
      //fprintf(stderr, "WARNING: region list section %s position %u is not found in the reference and skipped...\n", tokens[0].c_str(), chromosomeIndex);    
      continue;
    }
    
    len = tokens[2].AsInteger() - tokens[1].AsInteger() + 1;
    for(uint32_t i=startGenomeIndex; i<startGenomeIndex+len; i++)
      regionIndicator[i] = true;  
    
    tokens.Clear();
    buffer.Clear();            
  }
  ifclose(fhRegions);
  fprintf(stderr, "DONE!\n");
}


void GCContent::ReadGCContent(String infile)
{
  FILE *in = fopen(infile.c_str(), "rb");
  if(in==NULL) error("Open GC content file %s failed!\n", infile.c_str());
  uint32_t result = fread(gcCount, 1, genome->sequenceLength(), in);
  if(result!=genome->sequenceLength())
    error("Reading GC content from file %f failed!\n",infile.c_str());
  
  result = fread(gcContentVec,4,101, in);
  
  if(result!=101) error("Reading GC content error!\n");
  
  fclose(in);
  
}

void GCContent::CalcGCCount()
{
  int cnt = 0;
//  int cnt2 = 0;
  for(int i=0; i<windowSize; i++)
    {
      if(toupper((*genome)[i]) == 'G' || toupper((*genome)[i]) == 'C')
	cnt++;
      //if(toupper((*genome)[i]) == 'A' || toupper((*genome)[i]) == 'T')
      //	cnt2++;
      gcCount[i/2] = cnt;
      //atCount[i/2] = cnt2;
    }
  for(uint i=windowSize; i<length; i++)
    {
      if(toupper((*genome)[i])=='G' || toupper((*genome)[i]) == 'C')
	cnt++;
      //if(toupper((*genome)[i])=='A' || toupper((*genome)[i]) == 'T')
      //   cnt2++;
      if(toupper((*genome)[i-windowSize]) == 'G' ||toupper((*genome)[i-windowSize]) == 'C')
	cnt--;
      //if(toupper((*genome)[i-windowSize]) == 'A' ||toupper((*genome)[i-windowSize]) == 'T')
      //   cnt2--;
      
      gcCount[i-windowSize/2] = cnt;
      //atCount[i-windowSize/2] = cnt2;
    }
}

void GCContent::OutputGCContent(String & reference, int windowSize, String &gcContentFile, String &regionFile)
{
  
  GenomeSequence genome;
  bool memoryMap = false;
  
  genome.setReferenceName(reference.c_str());
  
  if (genome.open())
    {
      fprintf(stderr, "Failed to open reference index and is creating one...\n");
      if(genome.create())
	error("Failed to create reference index!\n");
    }
  
  genome.useMemoryMap(memoryMap);
  if(genome.open())
    error("Open  reference failed...!\n");
  

  LoadRegions(regionFile, genome);

  printf("%u %u\n", regionIndicator.size(), genome.sequenceLength());
  
  FILE * fh = fopen(gcContentFile.c_str(), "wb");
  
  int gccnt = 0;
  int atcnt = 0;
  uint8_t gc = 0;
  
  uint32_t gcContentDist[101];
  for(int i=0; i<101; i++) gcContentDist[i] = 0;
  
  for(int i=0; i<windowSize; i++)
    {
      if(toupper((genome)[i]) == 'G' || toupper((genome)[i]) == 'C')
	gccnt++;
      if(toupper((genome)[i]) == 'A' || toupper((genome)[i]) == 'T')
	atcnt++;
      
    }
  
  for(int i=0; i<windowSize/2; i++)
  {
    fwrite(&gc, 1, sizeof(uint8_t), fh);
    gcContentDist[gc]++;
  }
  
  for(uint32_t i=windowSize; i<genome.sequenceLength(); i++)
    {
      if(toupper((genome)[i])=='G' || toupper((genome)[i]) == 'C')
	gccnt++;
      if(toupper((genome)[i-windowSize]) == 'G' ||toupper((genome)[i-windowSize]) == 'C')
	gccnt--;
      if(toupper((genome)[i])=='A' || toupper((genome)[i]) == 'T')
	atcnt++;
      if(toupper((genome)[i-windowSize]) == 'A' ||toupper((genome)[i-windowSize]) == 'T')
	atcnt--;

      if(gccnt+atcnt<windowSize) gc = 0;
      else gc = floor(double(gccnt)/(gccnt+atcnt)*100+0.5);

      fwrite(&gc, 1, sizeof(uint8_t), fh);

      if(regionIndicator.size()>0 && regionIndicator[i]==false) gc = 0;
      gcContentDist[gc]++;
    }
  
  gc = 0;
 
  for(int i=0; i<windowSize/2; i++)
  {
    fwrite(&gc, 1, sizeof(uint8_t), fh);
    gcContentDist[gc]++;
  }
  
  fwrite(&gcContentDist, 1, sizeof(gcContentDist), fh);
  
  for(int i=0; i<101; i++) printf("%u\n", gcContentDist[i]);
   
  fclose(fh);
}

void GCContent::CalcGCFreq()
{
  String s;
  for(uint32_t i=0; i<length-windowSize/2; i++)
    {
      //if(gcCount[i]==0) { gcCount[i] = -1; continue; }
      //gcPerc = int(double(1000*gcCount[i])/double((gcCount[i]+atCount[i])));
      s.printf("%u", gcCount[i]);
      GCFreq.IncrementCount(s);
      //gcCount[i] = gcPerc;
    }
  //delete [] atCount;
}


void GCContent::CalcStatsByGC(uint8_t * v)
{
  for(int i=0; i<=windowSize; i++)
    statsVecByGC[gcCount[i]]=0;
  
  for(uint32_t i=0; i<length; i++)
    {
      //statsDist[v[i]]++;
      if(v[i]==0)    continue;
      //if(gcCount[i]==0) continue;
      //gcPerc = int(double(100*gcCount[i])/double((gcCount[i]+atCount[i])));
      //s.printf("%d", gcPerc);
      //GCFreq.IncrementCount(s);
      statsVecByGC[gcCount[i]] += (v[i]);
    }
}

void GCContent::IncrementStatsByGC(uint32_t start, int len)
{
  for(uint32_t i=start; i<start+len; i++)
    {
      //printf("%d\n", gcCount[i]);
      //if(gcCount[i]==-1) continue;
      statsVecByGC[gcCount[i]]++;
    }
}

void GCContent::OutputFreq(String file)
{
    // use fstream here so we don't have to worry about convoluted
    // printf compatibility with varying types below
  std::ofstream out;
  out.open(file.c_str(), std::ios_base::out | std::ios_base::trunc);
  if(!out.is_open())
    {
      std::cerr << "Unable to open file '" << file << "'.\n";
      return;
    }
  for(int i=0; i<GCFreq.Capacity(); i++)
    if(GCFreq.SlotInUse(i))
      {
        if(statsVecByGC!=NULL)
            out << GCFreq[i].c_str() << "\t" << (uint32_t)(GCFreq.Integer(i)) << "\t" << statsVecByGC[GCFreq[i].AsInteger()] << "\t" << (double)statsVecByGC[GCFreq[i].AsInteger()]/double(GCFreq.Integer(i)) << "\n";

//	fprintf(f, "%s\t%u\t%u\t%f\n", GCFreq[i].c_str(), (uint32_t)(GCFreq.Integer(i)), statsVecByGC[GCFreq[i].AsInteger()], (double)statsVecByGC[GCFreq[i].AsInteger()]/double(GCFreq.Integer(i)));
        else
          out << GCFreq[i].c_str() << "\t++\t" << (uint32_t)GCFreq.Integer(i) << "\n";
//          fprintf(f, "%s\t++\t%u\n", GCFreq[i].c_str(), (uint32_t)GCFreq.Integer(i));
      }
  out.close();
}


ReadDepth::ReadDepth()
{
  LIM = 255;
  depth = NULL;
  depthFreq = NULL;
}

void ReadDepth::SetDepthLimit(int lim) { LIM = lim; }

ReadDepth::ReadDepth(uint32_t len){AllocateMemory(len);}

void ReadDepth::AllocateMemory(uint32_t len)
{
  length = len;
  depth = new uint8_t[len];
  if(depth==NULL)
    error("Allocating memory for depth failed!\n");
  depthFreq = new uint32_t[LIM+1];
  //SetZeroCount();
}

ReadDepth::~ReadDepth() { 
  if(depth!=NULL)
    delete [] depth;
  if(depthFreq!=NULL) 
    delete [] depthFreq;
}

void ReadDepth::SetZeroCount()
{
  for(uint32_t i=0; i<length; i++)
    depth[i] = 0;
  for(int i=0; i<=LIM; i++)
    depthFreq[i] = 0;
}
void ReadDepth::IncrementCount(uint32_t idx)
{
  if(depth[idx]<LIM)
    ++depth[idx];
}
void ReadDepth::IncrementCount(uint32_t idx, int cnt)
{
  if((depth[idx]+cnt)<LIM)
    depth[idx] += cnt;
}
void ReadDepth::IncrementCount(uint32_t start, int len, int amount)
{
  for(uint32_t i=start; i<start+len; i++)
    if((depth[i]+amount)<LIM)
      depth[i] += amount;
}

void ReadDepth::CalcDepthFreq()
{
  String s;
  for(uint32_t i=0; i<length; i++)
    {
      //s.printf("%u", depth[i]);
      //depthFreq.IncrementCount(s);
      depthFreq[depth[i]]++;
    }
}

void ReadDepth::OutputFreq(String file)
{
  FILE *f = fopen(file.c_str(), "w");
  for(int i=0; i<=LIM; i++)
    //if(depthFreq.SlotInUse(i))
    fprintf(f, "%d\t%u\n", i, depthFreq[i]);
  
  fclose(f);
}
