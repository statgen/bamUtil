/*
 *  Copyright (C) 2010-2015  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include "CSG_MD5.h"
#include "SamFile.h"
#include "PolishBam.h"
#include "Logger.h"
#include "PhoneHome.h"

///////////////////////////////////////////////////////////////////////
//
//  FastaFile reading.
//
///////////////////////////////////////////////////////////////////////
class FastaFile {
public:
  std::string sFileName;
  std::vector<std::string> vsSequenceNames;
  std::vector<std::string> vsMD5sums;
  std::vector<uint32_t> vnSequenceLengths;
  std::ifstream ifsIn;
  char* pcLine;
  static const uint32_t nDefaultMaxLineLength = 100000;
  uint32_t nMaxLineLength;

  MD5_CTX curMD5Ctx;
  std::string curSeqName;
  uint32_t curSeqLength;
  int curSeqIndex;
  uint32_t nCurrentLine;

  FastaFile() : pcLine(NULL), nMaxLineLength(nDefaultMaxLineLength) {}
  ~FastaFile() {close();}
  bool open(const char* fileName);
  void close() {freeBuffer();}
  void allocateBuffer(int bufLength = nDefaultMaxLineLength);
  void freeBuffer();
  bool readLine();
  void readThru();
  void updateCurrentSequenceInfo();
  static uint32_t tokenizeString(const char* s, 
                                 std::vector<std::string>& tokens);
};

bool FastaFile::open(const char* fileName) {
  if ( sFileName.empty() && ( pcLine == NULL ) ) {
    sFileName = fileName;
    ifsIn.open(fileName, std::ios::in);
    if ( ifsIn.good() ) {
      nCurrentLine = 0;
      curSeqIndex = -1;
      curSeqLength = 0;
      allocateBuffer();
      return true;
    }
    else {
      return false;
    }
  }
  else {
    return false;
  }
}

void FastaFile::allocateBuffer(int bufLength) {
  freeBuffer();
  pcLine = new char[bufLength];
  nMaxLineLength = bufLength;
}

void FastaFile::freeBuffer() {
  if ( pcLine != NULL ) {
    delete[] pcLine;
    pcLine = NULL;
  }
}

bool FastaFile::readLine() {
  int n;

  ifsIn.getline(pcLine, nMaxLineLength);
  if ( pcLine[0] == '>' ) {
    if ( curSeqIndex >= 0 ) {
      updateCurrentSequenceInfo();
    }
    std::vector<std::string> tokens;
    tokenizeString(pcLine, tokens);
    curSeqName = tokens[0].substr(1); // make SeqName without leading '>'
    MD5Init(&curMD5Ctx);
    curSeqLength = 0;
    ++curSeqIndex;
  }
  else {
    n = strlen(pcLine);
    // convert to upper-case for calculating FASTA MD5
    for(int i = 0; i < n; ++i)
    {
        pcLine[i] = toupper(pcLine[i]);
    }
    MD5Update(&curMD5Ctx, (unsigned char*)pcLine, n);
    curSeqLength += n;
  }
  ++nCurrentLine;
  return ifsIn.good();
}

void FastaFile::readThru() {
    while(readLine());
    updateCurrentSequenceInfo();
}

void FastaFile::updateCurrentSequenceInfo() {
  unsigned char md5digest[1024];
  char md5string[1024];
  vsSequenceNames.push_back(curSeqName);
  vnSequenceLengths.push_back(curSeqLength);
  MD5Final(md5digest, &curMD5Ctx);
  for(int i=0; i < 16; ++i) {
    sprintf(&md5string[i*2],"%02x",static_cast<int>(md5digest[i]));
  }
  vsMD5sums.push_back(md5string);
}

uint32_t FastaFile::tokenizeString(const char* s, 
                                   std::vector<std::string>& tokens) {
    if ( !tokens.empty() ) {
        tokens.clear();
    }

    // Assuming delimiters are space of tab, 
    // use stringstream to tokenize
    // In order to accept other delimiters, 
    // this routine needs to be rewritten
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    std::string tok;
    
    ss << s;
    
    while( ss >> tok ) {
        tokens.push_back(tok); // each token is considerd as strings
    }
    
    return(tokens.size()); // return the number of tokenized elements
}


///////////////////////////////////////////////////////////////////////
//
//  polishBam
//
///////////////////////////////////////////////////////////////////////
/*
This is a C++ version of Tom Blackwell's awk script (detailed description
below) to add ReadGroup tag info. The C++ version takes each of the read group tag information as parameters instead of parsing the information from sequence.index file. Tom's awk script is modified to use as wrapper for this to be able to process BAM to BAM writing directly.
Hyun Min Kang, March 22, 2010

#  Add read group information to either a file or a stream (pipeline) 
#  in .sam format.  There are three possible arguments;  each can be 
#  set on the command line as  "var=value"  (no "-v" is needed).  
#  For most uses, default values for the arguments will suffice.
#  "indexfile" = path name to the 1000 Genomes "sequence.index" file.  
#  "keylen" = length of the substring to be used as a run identifier.  
#  "key" = identifier for current run, to be matched in column 3 of 
#  the sequence.index file, also used as the "RG" tag on every line.
#  "indexfile"  defaults to "sequence.index" if not set externally.  
#  "keylen"  defaults to 9 characters, in current 1000 Genomes usage.  
#  "key"  defaults to the first "keylen" characters from the name of 
#  the first sequence read in the .sam file.  
#  Tom Blackwell,  University of Michigan,   March 16, 2010
*/

void PolishBam::polishBamDescription()
{
    std::cerr << " polishBam - adds/updates header lines & adds the RG tag to each record" << std::endl;
}


void PolishBam::description()
{
    polishBamDescription();
}


void PolishBam::usage()
{
    BamExecutable::usage();
     std::cerr << "Usage: polishBam (options) --in <inBamFile> --out <outBamFile>\n" << std::endl;
     std::cerr << "Required parameters : " << std::endl;
     std::cerr << "-i/--in : input BAM file" << std::endl;
     std::cerr << "-o/--out : output BAM file" << std::endl;
     std::cerr << "Optional parameters :" << std::endl;
     std::cerr << "-v : turn on verbose mode" << std::endl;
     std::cerr << "-l/--log : writes logfile with specified name." << std::endl;
     std::cerr << "--HD : add @HD header line" << std::endl;
     std::cerr << "--RG : add @RG header line" << std::endl;
     std::cerr << "--PG : add @PG header line" << std::endl;
     std::cerr << "--CO : add @CO header line" << std::endl;
     std::cerr << "-f/--fasta : fasta reference file to compute MD5sums and update SQ tags" << std:: endl;
     std::cerr << "--AS : AS tag for genome assembly identifier" << std::endl;
     std::cerr << "--UR : UR tag for @SQ tag (if different from --fasta)" << std::endl;
     std::cerr << "--SP : SP tag for @SQ tag" << std:: endl;
     std::cerr << "--checkSQ : check the consistency of SQ tags (SN and LN) with existing header lines. Must be used with --fasta option" << std::endl;
     std::cerr << "\n" << std::endl;
}

void replaceTabChar(std::string& headerLine)
{
    // Convert "\t" to '\t'
    size_t tabPos = headerLine.find("\\t", 0);
        while(tabPos != std::string::npos)
        {
            headerLine.replace(tabPos, 2, "\t");
            tabPos = headerLine.find("\\t", tabPos);
        }

}

void checkHeaderStarts(std::vector<std::string>& headerLines, const char* type) {
    for(uint32_t i=0; i < headerLines.size(); ++i)
    {
        replaceTabChar(headerLines[i]);

        if ( headerLines[i].find(type) != 0 ) {
            Logger::gLogger->error("The following header line does not start with '%s'\n%s",type,headerLines[i].c_str());
        }
    }
}


void checkOrAddStarts(std::vector<std::string>& headerLines, const char* type) {
    for(uint32_t i=0; i < headerLines.size(); ++i)
    {
        replaceTabChar(headerLines[i]);

        size_t pos = headerLines[i].find(type);
        if(pos == std::string::npos) {
            // Does not start with prepend header;
            headerLines[i].insert(0, type);
        }else if(pos != 0 ) {
            Logger::gLogger->error("The following header line does not start with '%s'\n%s",type,headerLines[i].c_str());
        }
    }
}


int PolishBam::execute(int argc, char ** argv)
{
  static struct option getopt_long_options[] = 
    {
      // Input options
      { "fasta", required_argument, NULL, 'f'},
      { "in", required_argument, NULL, 'i'},
      { "out", required_argument, NULL, 'o'},
      { "verbose", no_argument, NULL, 'v'},
      { "log", required_argument, NULL, 'l'},
      { "clear", no_argument, NULL, 0},
      { "AS", required_argument, NULL, 0},
      { "UR", required_argument, NULL, 0},
      { "SP", required_argument, NULL, 0},
      { "HD", required_argument, NULL, 0},
      { "RG", required_argument, NULL, 0},
      { "PG", required_argument, NULL, 0},
      { "CO", required_argument, NULL, 0},
      { "checkSQ", no_argument, NULL, 0},
      { "noPhoneHome", no_argument, NULL, 'p'},
      { "nophonehome", no_argument, NULL, 'P'},
      { "phoneHomeThinning", required_argument, NULL, 't'},
      { "phonehomethinning", required_argument, NULL, 'T'},
      { NULL, 0, NULL, 0 },
    };

  // Shift arguments due to format being ./bam polishBam and then the args.
  ++argv;
  --argc;

  int n_option_index = 0, c;
  
  std::string sAS, sUR, sSP, sFasta, sInFile, sOutFile, sLogFile;
  bool bClear, bCheckSQ, bVerbose;
  std::vector<std::string> vsHDHeaders, vsRGHeaders, vsPGHeaders, vsCOHeaders;
  bool noPhoneHome = false;

  bCheckSQ = bVerbose = false;
  bClear = true;

  while ( (c = getopt_long(argc, argv, "vf:i:o:l:", getopt_long_options, &n_option_index)) != -1 ) {
      //    std::cout << getopt_long_options[n_option_index].name << "\t" << optarg << std::endl;
    if ( c == 'f' ) {
      sFasta = optarg;
    }
    else if ( c == 'i' ) {
      sInFile = optarg;
    }
    else if ( c == 'o' ) {
      sOutFile = optarg;
    }
    else if ( c == 'v' ) {
      bVerbose = true;
    }
    else if ( c == 'l' ) {
	sLogFile = optarg;
    }
    else if (( c == 'p' )||( c == 'P' )) {
        noPhoneHome = true;
    }
    else if (( c == 't' )||( c == 'T' )) {
        PhoneHome::allThinning = atoi(optarg);
    }
    else if ( strcmp(getopt_long_options[n_option_index].name,"AS") == 0 ) {
      sAS = optarg;
    }
    else if ( strcmp(getopt_long_options[n_option_index].name,"UR") == 0 ) {
      sUR = optarg;
    }
    else if ( strcmp(getopt_long_options[n_option_index].name,"SP") == 0 ) {
      sSP = optarg;
    }
    else if ( strcmp(getopt_long_options[n_option_index].name,"HD") == 0 ) {
      vsHDHeaders.push_back(optarg);
    }
    else if ( strcmp(getopt_long_options[n_option_index].name,"RG") == 0 ) {
      vsRGHeaders.push_back(optarg);
    }
    else if ( strcmp(getopt_long_options[n_option_index].name,"PG") == 0 ) {
      vsPGHeaders.push_back(optarg);
    }
    else if ( strcmp(getopt_long_options[n_option_index].name,"CO") == 0 ) {
      vsCOHeaders.push_back(optarg);
    }
    else if ( strcmp(getopt_long_options[n_option_index].name,"checkSQ") == 0 ) {
      bCheckSQ = true;
    }
    else {
      std::cerr << "Error: Unrecognized option " << getopt_long_options[n_option_index].name << std::endl;
      return(-1);
    }
  }

  if(!noPhoneHome)
  {
      PhoneHome::checkVersion(getProgramName(), VERSION);
  }

  if ((sLogFile.compare("__NONE__") == 0) ||  sLogFile.empty())
  {
      if(sOutFile.empty())
      {
          sLogFile = "-";
      }
      else
      {
          sLogFile = (sOutFile + ".log");
      }
  }

  Logger::gLogger = new Logger(sLogFile.c_str(), bVerbose);

  if ( optind < argc ) {
    usage();
    Logger::gLogger->error("non-option argument %s exist ",argv[optind]);
  }

  if ( sInFile.empty() || sOutFile.empty() ) {
    usage();
    Logger::gLogger->error("Input and output files are required");
  }

  if ( ( bCheckSQ ) && ( sFasta.empty() ) ) {
    usage();
    Logger::gLogger->error("--checkSQ option must be used with --fasta option");
  }

  // check whether each header line starts with a correct tag
  checkHeaderStarts(vsHDHeaders, "@HD\t");
  checkHeaderStarts(vsRGHeaders, "@RG\t");
  checkHeaderStarts(vsPGHeaders, "@PG\t");
  checkOrAddStarts(vsCOHeaders, "@CO\t");

  Logger::gLogger->writeLog("Arguments in effect:");
  Logger::gLogger->writeLog("\t--in [%s]",sInFile.c_str());
  Logger::gLogger->writeLog("\t--out [%s]",sOutFile.c_str());
  Logger::gLogger->writeLog("\t--log [%s]",sLogFile.c_str());
  Logger::gLogger->writeLog("\t--fasta [%s]",sFasta.c_str());
  Logger::gLogger->writeLog("\t--AS [%s]",sAS.c_str());
  Logger::gLogger->writeLog("\t--UR [%s]",sUR.c_str());
  Logger::gLogger->writeLog("\t--SP [%s]",sSP.c_str());
  Logger::gLogger->writeLog("\t--checkSQ [%s]",bClear ? "ON" : "OFF" );
  if ( vsHDHeaders.empty() ) {
    Logger::gLogger->writeLog("\t--HD []");
  }
  else {
    Logger::gLogger->writeLog("\t--HD [%s]",vsHDHeaders[0].c_str());
  }
  if ( vsRGHeaders.empty() ) {
    Logger::gLogger->writeLog("\t--RG []");
  }
  else {
    Logger::gLogger->writeLog("\t--RG [%s]",vsRGHeaders[0].c_str());
  }
  if ( vsPGHeaders.empty() ) {
    Logger::gLogger->writeLog("\t--PG []");
  }
  else {
    for(uint32_t i=0; i < vsPGHeaders.size(); ++i) {
      Logger::gLogger->writeLog("\t--PG [%s]",vsPGHeaders[i].c_str());
    }
  }
  if ( vsCOHeaders.empty() ) {
    Logger::gLogger->writeLog("\t--CO []");
  }
  else {
    for(uint32_t i=0; i < vsCOHeaders.size(); ++i) {
      Logger::gLogger->writeLog("\t--CO [%s]",vsCOHeaders[i].c_str());
    }
  }

  if ( (vsHDHeaders.empty() ) && ( vsRGHeaders.empty() ) && ( vsPGHeaders.empty() ) && ( vsCOHeaders.empty() ) && ( !bClear ) && ( sFasta.empty() ) ) {
    Logger::gLogger->warning("No option is in effect for modifying BAM files. The input and output files will be identical");
  }

  if ( ( vsHDHeaders.size() > 1 ) || ( vsRGHeaders.size() > 1 ) ) {
    Logger::gLogger->error("HD and RG headers cannot be multiple");
  }

  FastaFile fastaFile;
  if ( ! sFasta.empty() ) {
    if ( fastaFile.open(sFasta.c_str()) ) {
      Logger::gLogger->writeLog("Reading the reference file %s",sFasta.c_str());
      fastaFile.readThru();
      fastaFile.close();
      Logger::gLogger->writeLog("Finished reading the reference file %s",sFasta.c_str());      
    }
    else {
      Logger::gLogger->error("Failed to open reference file %s",sFasta.c_str());
    }
  }

  SamFile samIn;
  SamFile samOut;

  if ( ! samIn.OpenForRead(sInFile.c_str()) ) {
    Logger::gLogger->error("Cannot open BAM file %s for reading - %s",sInFile.c_str(), SamStatus::getStatusString(samIn.GetStatus()) );
  }
  if ( ! samOut.OpenForWrite(sOutFile.c_str()) ) {
    Logger::gLogger->error("Cannot open BAM file %s for writing - %s",sOutFile.c_str(), SamStatus::getStatusString(samOut.GetStatus()) );
  }

  SamFileHeader samHeader;
  SamHeaderRecord* pSamHeaderRecord;
  samIn.ReadHeader(samHeader);

  // check the sanity of SQ file
  // make sure the SN and LN matches, with the same order
  if ( bCheckSQ ) {
    unsigned int numSQ = 0;
    while( (pSamHeaderRecord = samHeader.getNextHeaderRecord()) != NULL ) {
      if ( pSamHeaderRecord->getType() == SamHeaderRecord::SQ ) {
	++numSQ;
      }
    }

    if ( numSQ != fastaFile.vsSequenceNames.size() ) {
      Logger::gLogger->error("# of @SQ tags are different from the original BAM and the reference file");
    }

    // iterator over all @SQ objects
    for(unsigned int i=0; i < numSQ; ++i) {
      pSamHeaderRecord = samHeader.getSQ(fastaFile.vsSequenceNames[i].c_str());
      if ( fastaFile.vsSequenceNames[i].compare(pSamHeaderRecord->getTagValue("SN")) != 0 ) {
	Logger::gLogger->error("SequenceName is not identical between fasta and input BAM file");
      }
      else if ( static_cast<int>(fastaFile.vnSequenceLengths[i]) != atoi(pSamHeaderRecord->getTagValue("LN")) ) {
	Logger::gLogger->error("SequenceLength is not identical between fasta and input BAM file");
      }
      else {
	if ( !sAS.empty() ) 
	  samHeader.setSQTag("AS",sAS.c_str(),fastaFile.vsSequenceNames[i].c_str());
	samHeader.setSQTag("M5",fastaFile.vsMD5sums[i].c_str(),fastaFile.vsSequenceNames[i].c_str());
	if ( !sUR.empty() ) 
	  samHeader.setSQTag("UR",sUR.c_str(),fastaFile.vsSequenceNames[i].c_str());
	if ( !sSP.empty() ) 
	  samHeader.setSQTag("SP",sSP.c_str(),fastaFile.vsSequenceNames[i].c_str());
      }
    }
    Logger::gLogger->writeLog("Finished checking the consistency of SQ tags");
  }
  else {
    Logger::gLogger->writeLog("Skipped checking the consistency of SQ tags");
  }

  // go over the headers again, 
  // assuming order of HD, SQ, RG, PG, CO, and put proper tags at the end of the original tags

  Logger::gLogger->writeLog("Creating the header of new output file");
  //SamFileHeader outHeader;
  samHeader.resetHeaderRecordIter();

  Logger::gLogger->writeLog("Adding %d HD, %d RG, %d PG, and %d CO headers",vsHDHeaders.size(), vsRGHeaders.size(), vsPGHeaders.size(), vsCOHeaders.size());
  int numHdSuccess = 0;
  for(unsigned int i=0; i < vsHDHeaders.size(); ++i) {
      if(samHeader.addHeaderLine(vsHDHeaders[i].c_str()))
      {
          ++numHdSuccess;
      }
  }

  /*
  for(int i=0; i < fastaFile.vsSequenceNames.size(); ++i) {
    std::string s("@SQ\tSN:");
    char buf[1024];
    s += fastaFile.vsSequenceNames[i];
    sprintf(buf,"\tLN:%d",fastaFile.vnSequenceLengths[i]);
    s += buf;
    if ( !sAS.empty() ) {
      sprintf(buf,"\tAS:%s",sAS.c_str());
      s += buf;
    }
    if ( !sUR.empty() ) {
      sprintf(buf,"\tUR:%s",sUR.c_str());
      s += buf;
    }
    sprintf(buf,"\tM5:%s",fastaFile.vsMD5sums[i].c_str());
    s += buf;
    if ( !sSP.empty() ) {
      sprintf(buf,"\tSP:%s",sSP.c_str());
      s += buf;
    }
    outHeader.addHeaderLine(s.c_str());
    }*/

  int numRgSuccess = 0;
  for(unsigned int i=0; i < vsRGHeaders.size(); ++i) {
      if(samHeader.addHeaderLine(vsRGHeaders[i].c_str()))
      {
          ++numRgSuccess;
      }
  }

  int numPgSuccess = 0;
  for(unsigned int i=0; i < vsPGHeaders.size(); ++i) {
      if(samHeader.addHeaderLine(vsPGHeaders[i].c_str()))
      {
          ++numPgSuccess;
      }
  }

  int numCoSuccess = 0;
  for(unsigned int i=0; i < vsCOHeaders.size(); ++i) {
      if(samHeader.addHeaderLine(vsCOHeaders[i].c_str()))
      {
          ++numCoSuccess;
      }
  }

  samOut.WriteHeader(samHeader);
  Logger::gLogger->writeLog("Successfully added %d HD, %d RG, %d PG, and %d CO headers",numHdSuccess, numRgSuccess, numPgSuccess, numCoSuccess);
  Logger::gLogger->writeLog("Finished writing output headers");

  // parse RG tag and get RG ID to append
  std::string sRGID;
  if ( ! vsRGHeaders.empty() ) {
    std::vector<std::string> tokens;
    FastaFile::tokenizeString( vsRGHeaders[0].c_str(), tokens );
    for(unsigned int i=0; i < tokens.size(); ++i) {
      if ( tokens[i].find("ID:") == 0 ) {
	sRGID = tokens[i].substr(3);
      }
    }
  }
  
  Logger::gLogger->writeLog("Writing output BAM file");
  SamRecord samRecord;
  while (samIn.ReadRecord(samHeader, samRecord) == true) {
    if ( !sRGID.empty() ) {
      if ( samRecord.addTag("RG",'Z',sRGID.c_str()) == false ) {
	Logger::gLogger->error("Failed to add a RG tag %s",sRGID.c_str());
      }
      // temporary code added
      if ( strncmp(samRecord.getReadName(),"seqcore_",8) == 0 ) {
	char buf[1024];
	sprintf(buf,"UM%s",samRecord.getReadName()+8);
	samRecord.setReadName(buf);
      }
    }
    samOut.WriteRecord(samHeader, samRecord);
    //if ( samIn.GetCurrentRecordCount() == 1000 ) break;
  }
  samOut.Close();
  Logger::gLogger->writeLog("Successfully written %d records",samIn.GetCurrentRecordCount());
  delete Logger::gLogger;
  return 0;
}
