#ifndef __GEN_MATRIX_H__
#define __GEN_MATRIX_H__

#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>

#include "VarMatrix.h"
#include "base/string_tokenizer.h"

///////////////////////////////////////////////////////////////
// VarMatrix class

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>

#include "base/logger.h"
#include "base/string_tokenizer.h"

#define DEFAULT_BUFFER_SIZE 1UL
#define DEFAULT_MAX_LINE_LENGTH 1000000UL

class VarMatrix;
class VcfMatrix;
class BedMatrix;

class VarIndividual;

class VarMarker;
class VcfMarker;
class BedMarker;

typedef struct {
  // 1 + 7 + 1 + 7 bit
  // 1 : 
  uint16_t genotype;
} BedData;

typedef struct {
  uint16_t genotype;
  float genotypeQuality;
  uint32_t depth;
} VcfData;

class Individual;
class VarMarker;

class VarMatrix {
  VarIndividual* vpIndividuals;
  VarMarker* vpMarkers;
  uint32_t nMarkerWindow;
  uint32_t nMarkerHead;
  uint32_t nMarkerTail;
};

class VarMatrix {
 public:
  std::vector<Individual*> vpIndividuals;
  std::vector<VarMarker*> vpVarMarkers;

  uint32_t nMarkerStartIndex; // start index of marker buffered
  uint32_t nMarkerEndIndex;   // end index of marker buffered
  uint32_t nBufferSize;
  bool bInfiniteBuffer;
  uint32_t nMaxLineLength;

  // genotype value has the following format for diploid
  // leftmost bit - 0 : unphased, 1 : phased
  // bit2-8   - allele1 representation (0-128)
  // bit9     - 0 : non-missing, 1 : missing
  // bit10-16 - allele2 representation (0-128)
  std::vector<uint16_t> vnGenotypes;

 VarMatrix() : nMarkerStartIndex(0), nMarkerEndIndex(0), nBufferSize(DEFAULT_BUFFER_SIZE), bInfiniteBuffer(false), nMaxLineLength(DEFAULT_MAX_LINE_LENGTH), pcGenotypes(NULL) {}  virtual ~VarMatrix();

  void setBufferSize(uint32_t newBufferSize) { nBufferSize = newBufferSize; }
  void setMaxLineLength(uint32_t newMaxLineLength) { nMaxLineLength = newMaxLineLength; }
  void setInfiniteBuffer() { bInfiniteBuffer = true; }

  virtual bool iterateMarker() = 0;
  char* currentGenotypes();
  char currentGenotype(uint32_t idx);
  inline VarMarker* currentMarker() { return vpVarMarkers[nMarkerEndIndex-1]; }
};

class BedMatrix : public VarMatrix {
 public:
  std::ifstream inFile;

  virtual ~BedMatrix();
 BedMatrix() {}

  void open(const char* famFile, const char* bimFile, const char* bedFile, const char* fastaFile = NULL);
  void openMarkerOnly(const char* bimFile, const char* fastaFile = NULL);
  virtual bool iterateMarker();
};

class VcfMatrix : public VarMatrix {
 public:
  std::ifstream inFile;
  double* pdGenotypeQualities;
  uint32_t* pnReadDepths;
  double* pdGenotypeLikelihoods;
  char* pcBuf;

  virtual ~VcfMatrix();
 VcfMatrix() : pdGenotypeQualities(NULL), pnReadDepths(NULL), pdGenotypeLikelihoods(NULL) {}
  void open(const char* vcfFile);
  virtual bool iterateMarker();
};

// Individual information, containing basic pedigree
class Individual {
 public:
  std::string sFamID;
  std::string sIndID;
  std::string sFatID;
  std::string sMotID;
  enum Gender { UnknownGender, Male, Female } gender;
};

// Maker Information (without actual genotype values)
// Contains all the fields compatible to VCF 4.0 (and also PLINK & MERLIN format)
class VarMarker {
 public:
  std::string sChrom;
  uint32_t nBasePosition;
  double dCenterMorgan;
  std::string sMarkerID;
  std::string sRef;
  std::string sAllele1;
  std::string sAllele2;
  double dAlleleFrequency;
  double dQual;
  std::vector<std::string> vsFilters;
  std::vector<std::string> vsInfoKeys;
  std::vector<std::string> vsInfoValues;

 VarMarker() : nBasePosition(0), dCenterMorgan(0), dAlleleFrequency(0), dQual(0) {}

  uint64_t getGenomicPosition();
  uint32_t getChromosomeIndex();
};

VarMatrix::~VarMatrix() {
  for(uint32_t i=0; i < vpIndividuals.size(); ++i) {
    delete vpIndividuals[i];
  }
  for(uint32_t i=0; i < vpVarMarkers.size(); ++i) {
    delete vpVarMarkers[i];
  }
}

BedMatrix::~BedMatrix() {
  if ( pcGenotypes != NULL ) {
    free(pcGenotypes);
  }
}

VcfMatrix::~VcfMatrix() {
  if ( pcGenotypes != NULL ) {
    free(pcGenotypes);
  }
  if ( pdGenotypeQualities != NULL ) {
    free(pdGenotypeQualities);
  }
  if ( pnReadDepths != NULL ) {
    free(pnReadDepths);
  }
  if ( pdGenotypeLikelihoods != NULL ) {
    free(pdGenotypeLikelihoods);
  }
  if ( pcBuf != NULL ) {
    delete [] pcBuf;
  }
}

// open BED files without reading actual genotypes
void BedMatrix::open(const char* famFile, const char* bimFile, const char* bedFile, const char* fastaFile) 
{
  char* buf = new char[nMaxLineLength];

  Logger::gLogger->writeLog("Opening binary PLINK files..\n\tFamFile = %s\n\tBimFile = %s\n\tBedFile = %s", famFile,bimFile,bedFile);

  // read .fam file and load into memory
  Logger::gLogger->writeLog("Reading %s", famFile);

  std::ifstream ifsFam(famFile);
  if ( ! ifsFam.is_open() ) {
    Logger::gLogger->error("openBED() : Error opening file %s", famFile);
  }

  while( !ifsFam.getline(buf, nMaxLineLength).eof() ) {
    int len = static_cast<int>(ifsFam.gcount());
    //Logger::gLogger->error("Line length = %d, %u", len, nMaxLineLength);

    CStringTokenizer t(buf, buf+len, "\t ");
    Individual* pInd = new Individual;
    int i;
    for(i=0; t.GetNext(); ++i) {
      switch(i) {
      case 0:
	pInd->sFamID = t.token().c_str();
	break;
      case 1:
	pInd->sIndID = t.token().c_str();
	break;
      case 2:
	pInd->sFatID = t.token().c_str();
	break;
      case 3:
	pInd->sMotID = t.token().c_str();
	break;
      case 4:
	{
	  const char* pToken = t.token().c_str();
	  if ( strlen(pToken) != 1 ) {
	    Logger::gLogger->error("BedMatrix::open() - Unrecognized gender %s in %s",pToken,famFile);
	  }
	  switch(pToken[0]) {
	  case '0':
	    pInd->gender = Individual::UnknownGender;
	    break;
	  case '1':
	    pInd->gender = Individual::Male;
	    break;
	  case '2':
	    pInd->gender = Individual::Female;
	    break;
	  default:
	    Logger::gLogger->error("BedMatrix::open() - Unrecognized gender %s in %s",pToken,famFile);
	    break;
	  }
	}
	break;
      case 5: // could be phenotype. do nothing now.
	break;
      default:
	Logger::gLogger->error("BedMatrix::open() - file %s contains more than 6 columns",famFile);
	break;
      }
    }

    if ( i < 5 ) {
      Logger::gLogger->error("BedMatrix::open() - file %s contains only %u columns",famFile, i);
    }

    vpIndividuals.push_back(pInd);
  }

  Logger::gLogger->writeLog("Finished Reading %s containing %u individuals' info", famFile, vpIndividuals.size());

  Logger::gLogger->writeLog("Reading %s", bimFile);
  // read .bim file and load into memory
  std::ifstream ifsBim(bimFile);
  while( !ifsBim.getline(buf, nMaxLineLength).eof() ) {
    int len = static_cast<int>(ifsBim.gcount());
    CStringTokenizer t(buf, buf+len, "\t ");
    VarMarker* pMarker = new VarMarker;
    int i;
    // we would expect six columns - CHROM, ID, CM, BP, AL1, AL2
    for(i=0; t.GetNext(); ++i) {
      switch(i) {
      case 0:
	pMarker->sChrom = t.token().c_str();
	break;
      case 1:
	pMarker->sMarkerID = t.token().c_str();
	break;
      case 2:
	pMarker->dCenterMorgan = atof(t.token().c_str());
	break;
      case 3:
	pMarker->nBasePosition = static_cast<uint32_t>(atoi(t.token().c_str()));
	break;
      case 4:
	pMarker->sAllele1 = t.token().c_str();
	break;
      case 5:
	pMarker->sAllele2 = t.token().c_str();
	break;
      case 6:
	pMarker->dAlleleFrequency = atof(t.token().c_str());
	break;
      default:
	Logger::gLogger->error("BedMatrix::open() - file %s contains more than 6 columns",bimFile);
	break;
      }
    }

    if ( ( i != 4 ) && ( i != 6 ) && ( i != 7 ) ) {
      Logger::gLogger->error("BedMatrix::open() - file %s contains less than 6 columns",bimFile);
    }

    vpVarMarkers.push_back(pMarker);
  }

  Logger::gLogger->writeLog("Finished Reading %s containing %u markers' info", bimFile, vpVarMarkers.size());

  Logger::gLogger->writeLog("Opening %s and checking the magic numbers", bedFile);
  // open the BED file and check the magic numbers
  inFile.open(bedFile, std::ifstream::in);
  int magicNumbers[3] = { 0x6c, 0x1b, 0x01 };
  for(int i=0; i < 3; ++i) {
    char c;
    inFile.read(&c,1);
    if ( magicNumbers[i] != c ) {
      Logger::gLogger->error("BedMatrix::open() - failed to validate the magic number of %s and direction (SNPmajor)",bedFile);
    }
  }
  Logger::gLogger->writeLog("BED magic numbers / SNP-major orientation were verified");

  // allocate buffers for genotypes only
  if ( nBufferSize == 0 ) {
    Logger::gLogger->error("BedMatrix::open() - buffer size is zero");
  }
  Logger::gLogger->writeLog("Allocating initial genotype buffer of %u lines", nBufferSize);
  pcGenotypes = (char*) malloc(sizeof(char) * nBufferSize * vpIndividuals.size());
  delete [] buf;
}

// open BED files without reading actual genotypes
void BedMatrix::openMarkerOnly(const char* bimFile, const char* fastaFile) 
{
  char* buf = new char[nMaxLineLength];

  Logger::gLogger->writeLog("Opening PLINK bimFile = %s",bimFile);

  Logger::gLogger->writeLog("Reading %s", bimFile);
  // read .bim file and load into memory
  std::ifstream ifsBim(bimFile);
  while( !ifsBim.getline(buf, nMaxLineLength).eof() ) {
    int len = static_cast<int>(ifsBim.gcount());
    CStringTokenizer t(buf, buf+len, "\t ");
    VarMarker* pMarker = new VarMarker;
    int i;
    // we would expect six columns - CHROM, ID, CM, BP, AL1, AL2
    for(i=0; t.GetNext(); ++i) {
      switch(i) {
      case 0:
	pMarker->sChrom = t.token().c_str();
	break;
      case 1:
	pMarker->sMarkerID = t.token().c_str();
	break;
      case 2:
	pMarker->dCenterMorgan = atof(t.token().c_str());
	break;
      case 3:
	pMarker->nBasePosition = static_cast<uint32_t>(atoi(t.token().c_str()));
	break;
      case 4:
	pMarker->sAllele1 = t.token().c_str();
	break;
      case 5:
	pMarker->sAllele2 = t.token().c_str();
	break;
      case 6:
	pMarker->dAlleleFrequency = atof(t.token().c_str());
	break;
      default:
	Logger::gLogger->error("BedMatrix::open() - file %s contains more than 6 columns",bimFile);
	break;
      }
    }

    if ( ( i != 6 ) && ( i != 4 ) && ( i != 7 ) ) {
      Logger::gLogger->error("BedMatrix::open() - file %s should contain 4 or 6 columns",bimFile);
    }

    vpVarMarkers.push_back(pMarker);
  }

  Logger::gLogger->writeLog("Finished Reading %s containing %u markers' info", bimFile, vpVarMarkers.size());

  delete [] buf;
}


// read a marker from BED file
bool BedMatrix::iterateMarker() {
  uint32_t nInds = vpIndividuals.size();
  uint32_t nBytes = static_cast<uint32_t>(ceil(nInds/4.));
  char* buf = new char[nBytes+1];
  inFile.read(buf,nBytes); // 0 - HOMREF, 1 - HET, 2 - MISSING, 3-HOMALT

  if ( !inFile.good() ) {
    delete [] buf;
    return false;
  }

  //for(uint32_t i=0; i < nBytes; ++i) {
  //  Logger::gLogger->writeLog("Byte %u/%u - %u",i,nBytes,static_cast<uint8_t>(buf[i]));
  //}

  // double the buffer size
  if ( ( bInfiniteBuffer ) && ( nMarkerEndIndex == nBufferSize ) ) {
    nBufferSize *= 2;
    pcGenotypes = (char*) realloc ( pcGenotypes, sizeof(char) * nBufferSize * nInds );
  }

  uint32_t curIdx = ( nMarkerEndIndex % nBufferSize) * nInds;
  //if ( curIdx != 0 ) Logger::gLogger->error("curIdx != 0");
  int c;
  for(uint32_t i=0; i < nInds; ++i) {
    switch ( c = ((buf[i/4] >> ( (i % 4) * 2 )) & ( 0x03 )) ) {
    case 0: // homozygote
      pcGenotypes[curIdx+i] = 1;
      break;
    case 1: // missing
      pcGenotypes[curIdx+i] = 0;
      //Logger::gLogger->error("Marker %u : missing genotype is observed at individual %u",nMarkerEndIndex, i);
      break;
    case 2: // heterozygote
      pcGenotypes[curIdx+i] = 2;
      break;
    case 3:
      pcGenotypes[curIdx+i] = 3;
      break;
    default:
      Logger::gLogger->error("BedMatrix::iterateMarker() - Unknown genotype value %u at marker %u, individual %u", c, nMarkerEndIndex, i);
    }
  }
  //Logger::gLogger->writeLog("%u %u %u %u %u %u %u %u %u",pcGenotypes[0],pcGenotypes[1],pcGenotypes[2],pcGenotypes[3],pcGenotypes[4],pcGenotypes[5],pcGenotypes[6],pcGenotypes[7],pcGenotypes[8]);
  ++nMarkerEndIndex;
  delete[] buf;
  return true;
}

char* VarMatrix::currentGenotypes() {
  return &pcGenotypes[ ((nMarkerEndIndex - 1) % nBufferSize) * vpIndividuals.size() ];
}

char VarMatrix::currentGenotype(uint32_t idx) {
  if ( idx < vpIndividuals.size() ) {
    return pcGenotypes[ ((nMarkerEndIndex - 1) % nBufferSize) * vpIndividuals.size() + idx ];
  }
  else {
    Logger::gLogger->error("VarMatrix::currentGenotype() - Out of bound in idx value %u",idx);
    return 0;
  }
}

uint64_t VarMarker::getGenomicPosition() {
  return (static_cast<uint64_t>(getChromosomeIndex()) << 32) & (static_cast<uint64_t>(nBasePosition));
}

uint32_t VarMarker::getChromosomeIndex() {
  uint32_t n = static_cast<uint32_t>(atoi(sChrom.c_str()));
  
  if ( n > 0 ) {
    return n;
  }
  else if ( sChrom.compare("X") == 0 ) {
    return 23;
  }
  else if ( sChrom.compare("Y") == 0 ) {
    return 24;
  }
  else if ( sChrom.compare("XY") == 0 ) {
    return 25;
  }
  else if ( ( sChrom.compare("MT") == 0 ) || ( sChrom.compare("Mt") == 0 ) || ( sChrom.compare("mt") == 0 ) ) {
    return 26;
  }
  else {
    Logger::gLogger->error("VarMarker::getChromosomeIndex() - Cannot recognize chromosome %s", sChrom.c_str());
    return 0;
  }
}


// open BED files without reading actual genotypes
void VcfMatrix::open(const char* vcfFile)
{
  pcBuf = new char[nMaxLineLength];

  Logger::gLogger->writeLog("Opening VCF file %s",vcfFile);

  inFile.open(vcfFile, std::ifstream::in);
  if ( ! inFile.is_open() ) {
    Logger::gLogger->error("VcfMatrix::open() : Error opening file %s", vcfFile);
  }

  while( !inFile.getline(pcBuf, nMaxLineLength).eof() ) {
    if ( ( pcBuf[0] == '#' ) && ( pcBuf[1] == '#' ) ) {  // comment line
      // ignore it
    }
    else if ( strncmp(pcBuf,"#CHROM\t",7) == 0 ) { // header comment line
      int len = static_cast<int>(inFile.gcount());
      CStringTokenizer t(pcBuf, pcBuf+len, "\t ");
      int i;
      for(i=0; t.GetNext(); ++i) {
	switch(i) {
	case 0:
	  if ( strcmp(t.token().c_str(),"#CHROM") != 0 ) {
	    Logger::gLogger->error("Unrecognized header column name %s at column %d", t.token().c_str(), i+1);
	  }
	  break;
	case 1:
	  if ( strcmp(t.token().c_str(),"POS") != 0 ) {
	    Logger::gLogger->error("Unrecognized header column name %s at column %d", t.token().c_str(), i+1);
	  }
	  break;
	case 2:
	  if ( strcmp(t.token().c_str(),"ID") != 0 ) {
	    Logger::gLogger->error("Unrecognized header column name %s at column %d", t.token().c_str(), i+1);
	  }
	  break;
	case 3:
	  if ( strcmp(t.token().c_str(),"REF") != 0 ) {
	    Logger::gLogger->error("Unrecognized header column name %s at column %d", t.token().c_str(), i+1);
	  }
	  break;
	case 4:
	  if ( strcmp(t.token().c_str(),"ALT") != 0 ) {
	    Logger::gLogger->error("Unrecognized header column name %s at column %d", t.token().c_str(), i+1);
	  }
	  break;
	case 5:
	  if ( strcmp(t.token().c_str(),"QUAL") != 0 ) {
	    Logger::gLogger->error("Unrecognized header column name %s at column %d", t.token().c_str(), i+1);
	  }
	  break;
	case 6:
	  if ( strcmp(t.token().c_str(),"FILTER") != 0 ) {
	    Logger::gLogger->error("Unrecognized header column name %s at column %d", t.token().c_str(), i+1);
	  }
	  break;
	case 7:
	  if ( strcmp(t.token().c_str(),"INFO") != 0 ) {
	    Logger::gLogger->error("Unrecognized header column name %s at column %d", t.token().c_str(), i+1);
	  }
	  break;
	case 8:
	  if ( strcmp(t.token().c_str(),"FORMAT") != 0 ) {
	    Logger::gLogger->error("Unrecognized header column name %s at column %d", t.token().c_str(), i+1);
	  }
	  break;
	default:
	  {
	    std::string tok = t.token().c_str();
	    tok.erase( tok.find_last_not_of(' ') + 1 );
	    if ( !tok.empty() ) {
	      Individual* pInd = new Individual;
	      pInd->sIndID = tok;
	      vpIndividuals.push_back(pInd);
	      Logger::gLogger->writeLog("%s Inserted - Total of %u individuals",tok.c_str(), vpIndividuals.size());
	    }
	    else {
	      Logger::gLogger->warning("Trailing space is found in header lines and ignored");
	    }
	  }
	  break;
	}
      }
      if ( i < 9 ) {
	Logger::gLogger->error("VcfMatrix::open() - too few # columns (%d) in the header line", i);
      }
      break;
    }
    else {
      Logger::gLogger->error("VcfMatrix::open() - no VCF header line was detected before reading the contents");
      // cannot be reached
    }
  }

  if ( ( vpIndividuals.size() == 0 ) || ( ! inFile.good() ) ) {
    Logger::gLogger->error("VcfMatrix::open() - no individual info has been loaded. ");
  }

  Logger::gLogger->writeLog("Finished reading header lines.. %u Individual IDs were recognized", vpIndividuals.size());

  // allocate buffers for genotypes only
  if ( nBufferSize == 0 ) {
    Logger::gLogger->error("VcfMatrix::open() - buffer size is zero");
  }
  Logger::gLogger->writeLog("Allocating initial genotype buffer of %u lines", nBufferSize);
  pcGenotypes = (char*) malloc(sizeof(char) * nBufferSize * vpIndividuals.size());
}

bool VcfMatrix::iterateMarker() {
  uint32_t nInds = vpIndividuals.size();
  if ( inFile.getline(pcBuf, nMaxLineLength).eof() ) {
    return false;
  }
  else {
    CStringTokenizer t(pcBuf, pcBuf+inFile.gcount(), "\t ");
    int i;
    VarMarker* pMarker = new VarMarker;
    int nGTidx = -1;

    // double the buffer size if needed
    if ( ( bInfiniteBuffer ) && ( nMarkerEndIndex == nBufferSize ) ) {
      nBufferSize *= 2;
      pcGenotypes = (char*) realloc ( pcGenotypes, sizeof(char) * nBufferSize * nInds );
    }
    uint32_t curIdx = (nMarkerEndIndex % nBufferSize) * nInds;

    for(i=0; t.GetNext(); ++i) {
      switch(i) {
      case 0: // CHROM
	pMarker->sChrom = t.token().c_str();
	break;
      case 1: // POS
	pMarker->nBasePosition = static_cast<uint32_t>(atoi(t.token().c_str()));
	break;
      case 2: // ID
	pMarker->sMarkerID = t.token().c_str();
	break;
      case 3: // REF
	pMarker->sAllele1 = t.token().c_str();
	pMarker->sRef = pMarker->sAllele1;
	break;
      case 4: // ALT
	pMarker->sAllele2 = t.token().c_str();
	break;
      case 5: // QUAL
	pMarker->dQual = atof(t.token().c_str());
	break;
      case 6: // FILTER
	pMarker->sFilter = t.token().c_str();
	break;
      case 7: // INFO 
	{
	  std::string tokenInfo = t.token().c_str();
	  StringTokenizer t2(tokenInfo,";");
	  while(t2.GetNext()) {
	    // copy the token
	    std::string tokenInfoItem = t2.token().c_str();
	    //Logger::gLogger->writeLog("Parsing INFO token %s",tokenInfoItem.c_str());
	    StringTokenizer t3(tokenInfoItem,"=");

	    if ( t3.GetNext() ) {
	      //Logger::gLogger->writeLog("INFO token %s",t3.token().c_str());
	      std::string key = t3.token();
	      if ( t3.GetNext() ) {
		std::string val = t3.token();
		//Logger::gLogger->writeLog("INFO token %s",t3.token().c_str());
		if ( t3.GetNext() ) {
		  //Logger::gLogger->writeLog("INFO token %s",t3.token().c_str());
		  Logger::gLogger->error("VcfMatrix::iterateMarker() : More than three columns observed in a INFO tag %s",t2.token().c_str());
		}
		vsInfoKeys.push_back(key);
		vsInfoValues.push_back(val);
	      }
	      else { // no value;
		vsInfoKeys.push_back(key);
		vsInfoValues.push_back("");
	      }
	    }
	    else {
	      Logger::gLogger->error("VcfMatrix::iterateMarker() : No token obtained in parsing a INFO tag %s - %s in marker %u",tokenInfo.c_str(),nMarkerEndIndex+1);
	    }
	  }
	}
	break;
      case 8: // FORMAT column
	{
	  StringTokenizer t2(t.token(),":");
	  for(int j=0;t2.GetNext();++j) {
	    if ( strcmp(t2.token().c_str(),"GT") == 0 ) {
	      nGTidx = j;
	    }
	  }
	  if ( nGTidx < 0 ) {
	    Logger::gLogger->error("VcfMatrix::open() : No GT field is format column %s",t.token().c_str());
	  }
	}
	break;
      default: // Actual values
	{
	  std::string tok = t.token().c_str();
	  tok.erase( tok.find_last_not_of(" ") + 1 );
	  if ( !tok.empty() ) {
	    StringTokenizer t2(tok,":");
	    int j;
	    for(j=0; t2.GetNext(); ++j) {
	      if ( j == nGTidx ) { // matches with GT
		const char* p = t2.token().c_str();
		if ( strlen(p) == 3 ) {
		  if ( ( p[1] == '/') || ( p[1] == '|' ) || ( p[1] == '\\' ) ) {
		    if ( ( p[0] == '.' ) && ( p[2] == '.' ) ) {
		      pcGenotypes[curIdx+i-8] = 0;
		    }
		    else if ( ( p[0] == '.' ) || ( p[2] == '.' ) ) {
		      Logger::gLogger->error("VcfMatrix::iterateMarker() : Unrecognized GT token %s at line %d column %d", p, nMarkerEndIndex, i+1);
		    }
		    else {
		      pcGenotypes[curIdx+i-8] = 1 + (p[0] > 0 ? 1 : 0) + (p[2] > 0 ? 1 : 0);  // supports only biallelic marker
		  }
		}
		else {
		  Logger::gLogger->error("VcfMatrix::iterateMarker() : Unrecognized GT token %s at line %d column %d", p, nMarkerEndIndex, i+1);
		}
	      }
	      else {
		Logger::gLogger->error("VcfMatrix::iterateMarker() : Unrecognized GT token %s at line %d column %d", p, nMarkerEndIndex, i+1);
	      }
	      }
	    }
	    if ( j < nGTidx ) { // not enough # columns
	      Logger::gLogger->error("VcfMatrix::iterateMarker() : Unrecognized Individual token %s at line %d column %d", t.token().c_str(), nMarkerEndIndex, i+1);
	    }
	  }
	  else {
	    Logger::gLogger->warning("Trailing space is found and ignored in iterating individual %d at marker %u",i,nMarkerEndIndex+1);
	    --i;
	  }
	}
	break;
      }
    }
    if ( i != static_cast<int>(nInds) + 9 ) {
      std::string tok = pcBuf;
      std::remove(tok.begin(), tok.end(), ' ');
      //tok.erase(std::remove_if(tok.begin(), tok.end(), std::isspace), tok.end());
      delete pMarker;

      if ( tok.empty() ) {
	Logger::gLogger->warning("VcfMatrix::iterateMarker() : Empty marker line detected and ignored");
      }
      else {
	Logger::gLogger->error("VcfMatrix::iterateMarker() : Number of tokens %d does not match to the expected number %d+9=%d at marker %d",i,nInds,nInds+9,nMarkerEndIndex+1);
      }
    }
    vpVarMarkers.push_back(pMarker);
    ++nMarkerEndIndex;
    return true;
  }
}

#endif // __GEN_MATRIX_H__
