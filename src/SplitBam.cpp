/*
 *  Copyright (C) 2010  Hyun Min Kang
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
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include "SplitBam.h"
#include "SamFile.h"

////////////////////////////////////////////////////////////////////////
// SplitBam : 
//   Split a BAM file into multiple BAM files based on ReadGroup
// 
// Copyright (c) 2010 Hyun Min Kang
// Last modified Jun 10, 2010
// Modified 1/16/2012 by Mary Kate Trost to put into bamUtil.
//
// SplitBam splits a BAM file into multiple BAM files based on
//   ReadGroup according to the following details.
// (1) Creates multiple output files named [outprefix].[RGID].bam, for
//     each ReadGroup ID (RGID) existing in the bam file
// (2) Headers are a copy of the original file, removing @RG and @PG
//     headers where IDs match with the other ReadGroup IDs.
// (3) Copy each of the original file's BAM record to one of the output
//     file where the ReadGroup ID matches
//
// Usage : bam splitBam [-v] [-i inputBAMFile] [-o outPrefix] [-L logFile]
// Required arguments:
//    -i inputBAMFile : Original BAM file containing readGroup info
//    -o outprefix  : prefix of output bam files of [outprefix].[RGID].bam
// Optional arguments:
//    -L logFile  : log file name. default is listFile.log
//    -v : turn on verbose mode
///////////////////////////////////////////////////////////////////////


// Logger class for logging/error/warning
class Logger {
 protected:
  FILE* fp_log;
  FILE* fp_err;
  bool b_verbose;

  Logger() {} // default constructor prohibited
 public:
  static Logger* gLogger;
  Logger(const char* filename, bool verbose);
  void writeLog(const char* format, ...);
  void error(const char* format, ...);
  void warning(const char* format, ...);
};

// ReadGroup w/ RG ID and full header line
class ReadGroup {
public:
  std::string s_id;
  std::string s_header_line;
};

// Global variables
Logger* Logger::gLogger = NULL;
uint64_t MAX_GENOMIC_COORDINATE = 0xffffffffffffffffULL;
uint64_t UNMAPPED_GENOMIC_COORDINATE = 0xfffffffffffffffeULL;

// function declarations
bool parseListFile(std::string& listFile, vector<std::string>& bamFiles, vector<ReadGroup>& readGroups, vector<ReadGroup>& uniqReadGroups);
bool equalHeaders(SamFileHeader& header1, SamFileHeader& header2);
uint64_t getGenomicCoordinate(SamRecord& r);
void addReadGroupToHeader(SamFileHeader& header, ReadGroup& rg);
void addReadGroupTag(SamRecord& record, ReadGroup& rg);
uint32_t addTokenizedStrings(const std::string& str, const std::string& delimiters, vector<std::string>& tokens);


void SplitBam::splitBamDescription()
{
    std::cerr << " splitBam - Split a BAM file into multiple BAM files based on ReadGroup" << std::endl;
}


void SplitBam::description()
{
    splitBamDescription();
}


// print Usage
void SplitBam::usage()
{
    BamExecutable::usage();
  std::cerr << "\t ./bam splitBam [-v] -i <inputBAMFile> -o <outPrefix> [-L logFile]" << std::endl;
  std::cerr << "splitBam splits a BAM file into multiple BAM files based on" << std::endl;
  std::cerr << "ReadGroup according to the following details." << std::endl;
  std::cerr << "\t(1) Creates multiple output files named [outprefix].[RGID].bam, for" << std::endl;
  std::cerr << "\teach ReadGroup ID (RGID) in the BAM record" << std::endl;
  std::cerr << "\t(2) Headers are a copy of the original file, removing @RG and @PG" << std::endl;
  std::cerr << "\theaders where IDs match with the other ReadGroup IDs." << std::endl;
  std::cerr << "\t(3) Copy each of the original file's BAM record to one of the output" << std::endl;
  std::cerr << "file where the ReadGroup ID matches" << std::endl;
  std::cerr << "Required arguments:" << std::endl;
  std::cerr << "-i/--in [inputBAMFile] : Original BAM file containing readGroup info" << std::endl;
  std::cerr << "-o/--out [outPrefix] : prefix of output bam files of [outprefix].[RGID].bam" << std::endl;
  std::cerr << "Optional arguments:" << std::endl;
  std::cerr << "-L/--log [logFile]  : log file name. default is listFile.log" << std::endl;
  std::cerr << "-v/--verbose : turn on verbose mode" << std::endl;
}

// main function
int SplitBam::execute(int argc, char ** argv)
{
  static struct option getopt_long_options[] = 
    {
      // Input options
      { "in", required_argument, NULL, 'i'},
      { "out", required_argument, NULL, 'o'},
      { "verbose", no_argument, NULL, 'v'},
      { "log", required_argument, NULL, 'L'},
      { NULL, 0, NULL, 0 },
    };

  int n_option_index = 0;
  char c;
  bool b_verbose = false;

  std::string s_in, s_out, s_logger;

  while ( ( c = getopt_long(argc-1, &(argv[1]), "i:o:vL:", getopt_long_options, &n_option_index) ) != -1 ) {
    switch(c) {
    case 'i':
      s_in = optarg;
      break;
    case 'o':
      s_out = optarg;
      break;
    case 'v':
      b_verbose = true;
      break;
    case 'L':
      s_logger = optarg;
      break;
    default:
      fprintf(stderr,"Unrecognized option %s",getopt_long_options[n_option_index].name);
      abort();
    }
  }

  if ( s_logger.empty() ) {
    s_logger = s_out + ".log";
  }

  // create a logger object, now possible to write logs/warnings/errors
  Logger::gLogger = new Logger(s_logger.c_str(), b_verbose);

  // every argument must correspond to an option
  if ( optind < (argc-1) ) {
    usage();
    Logger::gLogger->error("non-option argument exist");
  }

  // check the required arguments are nonempty
  if ( s_in.empty() || s_out.empty() ) {
    usage();
    Logger::gLogger->error("At least one of the required argument is missing");
  }

  Logger::gLogger->writeLog("Input BAM file : %s",s_in.c_str());
  Logger::gLogger->writeLog("Output BAM prefix : %s",s_out.c_str());
  Logger::gLogger->writeLog("Output log file : %s",s_logger.c_str());
  Logger::gLogger->writeLog("Verbose mode    : %s",b_verbose ? "On" : "Off");
  
  SamFile inBam;
  SamFileHeader inHeader;
  std::map<std::string,uint32_t> msRGidx;
  std::vector<std::string> vsRGIDs;
  std::vector<SamFile*> vpOutBams;
  std::vector<SamFileHeader*> vpOutHeaders;

  if ( ! (inBam.OpenForRead(s_in.c_str()))  ) {
    Logger::gLogger->error("Cannot open BAM file %s for reading - %s",s_in.c_str(), SamStatus::getStatusString(inBam.GetStatus()) );
  }
  inBam.ReadHeader(inHeader);

  SamHeaderRecord* pSamHeaderRecord;
  while( (pSamHeaderRecord = inHeader.getNextHeaderRecord()) != NULL ) {
      SamHeaderRecord::SamHeaderRecordType sHeaderRecordType(pSamHeaderRecord->getType());
    if ( sHeaderRecordType == SamHeaderRecord::RG) {
      std::string sRGID = pSamHeaderRecord->getTagValue("ID");
      if ( sRGID.empty() ) {
	Logger::gLogger->error("Readgroup ID is empty");
      }
      vsRGIDs.push_back(sRGID);
      uint32_t idx = msRGidx.size();
      msRGidx[sRGID] = idx;
      SamFile* pNewFile = new SamFile;
      vpOutBams.push_back(pNewFile);
      
      std::string outFileName = s_out + "." + sRGID + ".bam";
      if ( !pNewFile->OpenForWrite(outFileName.c_str()) ) {
	Logger::gLogger->error("Cannot open BAM file %s for writing",outFileName.c_str());
      }
      
      SamFileHeader* pNewHeader = new SamFileHeader(inHeader);
      vpOutHeaders.push_back(pNewHeader);
    }
  }

  Logger::gLogger->writeLog("The following ReadGroup IDs are identified. Splitting into %u BAM files..",vsRGIDs.size());
  for(uint32_t i=0; i < vsRGIDs.size(); ++i) {
    Logger::gLogger->writeLog("\t%d: %s",i+1,vsRGIDs[i].c_str());
  }
  

  if ( vsRGIDs.size() == 0 ) {
    Logger::gLogger->error("Only %u readGroups are observed",vsRGIDs.size());
  }
  else if ( vsRGIDs.size() ==1 ) {
    Logger::gLogger->warning("Only %u readGroups are observed",vsRGIDs.size());
  }

  // remove non-relevant readGroups
  for(uint32_t i=0; i < vsRGIDs.size(); ++i) {
    for(uint32_t j=0; j < vsRGIDs.size(); ++j) {
      if ( i != j ) {
	vpOutHeaders[i]->removeRG(vsRGIDs[j].c_str());
	vpOutHeaders[i]->removePG(vsRGIDs[j].c_str());
      }
    }
  }

  // write headers to the output file
  for(uint32_t i=0; i < vsRGIDs.size(); ++i) {
    vpOutBams[i]->WriteHeader(*vpOutHeaders[i]);
  }

  SamRecord record;
  while( inBam.ReadRecord(inHeader, record) == true ) {
    char tag[3];
    char vtype;
    void* value;
    bool found = false;
    while( record.getNextSamTag(tag, vtype, &value) != false ) {
      if ( strcmp(tag,"RG") == 0 ) {
	found = true;
	if ( vtype == 'Z' ) {
	  std::string sValue = ((String)*(String*)value).c_str();
	  if ( msRGidx.find(sValue) != msRGidx.end() ) {
	    uint32_t idx = msRGidx[sValue];
	    if ( (idx >= 0 ) && ( idx < vsRGIDs.size () ) ) {
	      vpOutBams[idx]->WriteRecord(inHeader, record);
	    }
	    else {
	      Logger::gLogger->error("ReadGroup Index Lookup Failure");
	    }
	  }
	  else {
	    Logger::gLogger->error("ReadGroup ID %s cannot be found",sValue.c_str());
	  }
	}
	else {
	  Logger::gLogger->error("vtype of RG tag must be 'Z'");
	}
	break;
      }
    }
    if ( found == false ) {
      Logger::gLogger->error("Cannot find RG tag for readName %s",record.getReadName());
    }
  }

  for(uint32_t i=0; i < vsRGIDs.size(); ++i) {
    Logger::gLogger->writeLog("Successfully wrote %d record for readGroup %s",vpOutBams[i]->GetCurrentRecordCount(), vsRGIDs[i].c_str());
    vpOutBams[i]->Close();
    delete vpOutBams[i];
    delete vpOutHeaders[i];
  }

  delete Logger::gLogger;
  return 0;
}

// Constructor of logger
Logger::Logger(const char* filename, bool verbose) 
{
  b_verbose = verbose;
  fp_log = fopen(filename, "w");
  if ( fp_log == NULL ) {
    fprintf(stderr,"ERROR: Cannot open the log file %s. Check if the directory exists and you have the permission to create a file", filename);
    abort();
  }
  fp_err = stderr;
}

// Write a log to output file
void Logger::writeLog(const char* format, ... ) {
  va_list args;
  va_start (args, format);
  vfprintf(fp_log, format, args);
  va_end (args);
  fprintf(fp_log, "\n");
  fflush(fp_log);
  
  if ( b_verbose ) {
    va_start (args, format);
    vfprintf(fp_err, format, args);
    va_end (args);
    fprintf(fp_err, "\n");
  }
}

// Write error messages and abort
void Logger::error(const char* format, ... ) {
  va_list args;
  va_start (args, format);
  fprintf(fp_log, "ERROR: ");
  vfprintf(fp_log, format, args);
  va_end (args);
  fprintf(fp_log, "\n");
  fflush(fp_log);

  va_start (args, format);
  fprintf(fp_err, "ERROR : ");
  vfprintf(fp_err, format, args);
  va_end (args);
  fprintf(fp_err, "\n");

  abort();
}

// Write warning messages
void Logger::warning(const char* format, ... ) {
  va_list args;
  va_start (args, format);
  fprintf(fp_log, "WARNING: ");
  vfprintf(fp_log, format, args);
  va_end (args);
  fprintf(fp_log, "\n");
  fflush(fp_log);

  if ( b_verbose ) {
    va_start (args, format);
    fprintf(fp_err, "WARNING : ");
    vfprintf(fp_err, format, args);
    va_end (args);
    fprintf(fp_err, "\n");
  }
}

// parse a list file containing BAM file and readGroup info
bool parseListFile(std::string& listFile, vector<std::string>& bamFiles, vector<ReadGroup>& readGroups, vector<ReadGroup>& uniqReadGroups) {
  // use getline and stringstream to tokenize the input
  std::ifstream fin;
  std::string line;
  std::map<std::string, uint32_t> columnDict;

  fin.open(listFile.c_str(), std::ios::in);
  if ( ! fin.good() ) {
    Logger::gLogger->error("Cannot read the input file %s", listFile.c_str());
  }

  // read header lines and generate a dictionary
  if ( ! std::getline(fin, line) ) {
    Logger::gLogger->error("Cannot read the first line of the input file %s",listFile.c_str());
  }
  vector<std::string> tokens;
  std::string delimiters = "\t\n";
  addTokenizedStrings(line, delimiters, tokens);
  for(uint32_t i=0; i < tokens.size(); ++i) {
    columnDict[tokens[i]] = i;
  }
  /*
  std::stringstream ss (std::stringstream::in | std::stringstream::out);
  std::string tok;
  ss << line;
  for(uint32_t i=0; ss >> tok; ++i) {
    columnDict[tok] = i;  // columnDict contains the indices of each column name
    }*/

  // make sure that the required fields exists
  const char *requiredColumnNames[3] = {"BAM","ID","SM"};
  for(uint32_t i=0; i < 3; ++i) {
    if ( columnDict.find( requiredColumnNames[i] ) == columnDict.end() ) {
      Logger::gLogger->error("Required column name '%s' does not exist in the header of the list file", requiredColumnNames[i]);
    }
  }

  // read the next line and create a readgroup using existing fields
  //vector<std::string> tokens;
  const uint32_t NUM_COLS = 9;
  const char* columnNames[NUM_COLS] = {"ID","SM","LB","DS","PU","PI","CN","DT","PL"};
  while ( std::getline(fin, line) ) {
    tokens.clear();
    if ( addTokenizedStrings(line, delimiters, tokens) != columnDict.size() ) {
      Logger::gLogger->error("The column sizes are different between the lines in the input file %s. The header size is %d, and the line size is %d",listFile.c_str(),columnDict.size(), tokens.size());
    }

    /*
    std::stringstream ss2 (std::stringstream::in | std::stringstream::out);
    tokens.clear();
    ss2 << line;
    while( ss2 >> tok ) { tokens.push_back(tok); }
    if ( columnDict.size() != tokens.size() ) {
      Logger::gLogger->error("The column sizes are different between the lines in the input file %s. The header size is %d, and the line size is %d",listFile.c_str(),columnDict.size(), tokens.size());
      }*/

    bamFiles.push_back(tokens[columnDict["BAM"]]);
    std::string readGroupID = tokens[columnDict["ID"]]; 
    std::string headerLine = "@RG";
    for(uint32_t i=0; i < NUM_COLS; ++i) {
      if ( columnDict.find(columnNames[i]) != columnDict.end() ) {
	if ( !tokens[columnDict[columnNames[i]]].empty() ) {
	  headerLine += "\t";
	  headerLine += columnNames[i];
	  headerLine += ":";
	  headerLine += tokens[columnDict[columnNames[i]]];
	}
      }
    }
    // scan existing readGroup to detect duplicates
    bool foundDuplicates = false;
    for(uint32_t i=0; i < readGroups.size(); ++i) {
      if ( readGroupID.compare(readGroups[i].s_id) == 0 ) {
	foundDuplicates = true;
	if ( headerLine.compare(readGroups[i].s_header_line) != 0 ) {
	  Logger::gLogger->error("Duplicated readGroup ID %s contains inconsistent column info between them\nLine 1 : %s\nLine 2 : %s",readGroupID.c_str(),headerLine.c_str(),readGroups[i].s_header_line.c_str());
	}
      }
    }

    ReadGroup rg;
    readGroups.push_back(rg);
    readGroups.back().s_id = readGroupID;
    readGroups.back().s_header_line = headerLine;
    if ( ! foundDuplicates ) {
      uniqReadGroups.push_back(rg);
      uniqReadGroups.back().s_id = readGroupID;
      uniqReadGroups.back().s_header_line = headerLine;

      Logger::gLogger->writeLog("Adding a header line - %s",headerLine.c_str());
    }
    else {
      Logger::gLogger->writeLog("Skipping a duplicated header line - %s",headerLine.c_str());
    }
  }
  return true;
}

// Compare the identity of two headers
bool equalHeaders(SamFileHeader& header1, SamFileHeader& header2) {
  std::string s1, s2;
  header1.getHeaderString(s1);
  header2.getHeaderString(s2);
  return ( s1.compare(s2) == 0 );
}

// make a 64-bit genomic coordinate [24bit-chr][32bit-pos][8bit-orientation]
uint64_t getGenomicCoordinate(SamRecord& r) {
  // 64bit string consisting of
  // 24bit refID, 32bit pos, 8bit orientation
  if ( ( r.getReferenceID() < 0 ) || ( r.get0BasedPosition() < 0 ) ) {
    return UNMAPPED_GENOMIC_COORDINATE;
  }
  else {
    return ( ( static_cast<uint64_t>(r.getReferenceID()) << 40 ) | ( static_cast<uint64_t>(r.get0BasedPosition()) << 8 ) | static_cast<uint64_t>( r.getFlag() & 0x0010 ) );
  }
}

// add readgroup header line to the SamFileHeader
void addReadGroupToHeader(SamFileHeader& header, ReadGroup& rg) {
  if ( !header.addHeaderLine(rg.s_header_line.c_str()) ) {
    Logger::gLogger->error("Failed to add ID = %s, header line %s",rg.s_id.c_str(),rg.s_header_line.c_str());
  }
}

// add a readgroup tag to ech SamRecord
void addReadGroupTag(SamRecord& record, ReadGroup& rg) {
  if ( record.addTag("RG",'Z',rg.s_id.c_str()) == false ) {
    Logger::gLogger->error("Failed to add readgroup tag %s",rg.s_id.c_str());
  }
}

uint32_t addTokenizedStrings(const std::string& str, const std::string& delimiters, vector<std::string>& tokens)
{
  std::string::size_type delimPos = 0, tokenPos = 0, pos = 0;
  uint32_t numAddedTokens = 0;
  
  if( str.length() < 1 )  return 0;

  while( true ){
    delimPos = str.find_first_of(delimiters, pos);
    tokenPos = str.find_first_not_of(delimiters, pos);
    
    if(std::string::npos != delimPos){
      if(std::string::npos != tokenPos){
	if( tokenPos < delimPos ) {
	  tokens.push_back(str.substr(pos,delimPos-pos));
	  ++numAddedTokens;
	}
	else {
	  tokens.push_back("");
	  ++numAddedTokens;
	}
      }
      else {
	tokens.push_back("");
	++numAddedTokens;
      }
      pos = delimPos+1;
    } 
    else {
      if( std::string::npos != tokenPos ){
	tokens.push_back(str.substr(pos));
	++numAddedTokens;
      } 
      else {
	tokens.push_back("");
	++numAddedTokens;
      }
      break;
    }
  }
  return numAddedTokens;
}
