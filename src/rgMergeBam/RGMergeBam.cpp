#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include "SamFile.h"

////////////////////////////////////////////////////////////////////////
// RGMergeBAM : Merge multiple BAM files appending ReadGroup IDs
// 
// Copyright (c) 2010 Hyun Min Kang
//
// rgMergeBam merges multiple sorted BAM files into one BAM file like
//   'samtools merge' command, but merges BAM headers.
// (1) check the HD and SQ tags are identical across the BAM files
// (2) 
// (1) Add @RG headers from a tabular input file containing the fields' info
// (2) Add RG:Z:[RGID] tag for each record based on the source BAM file
// (3) Ensure that the headers are identical across the input files
//     and that input/output BAM records are sorted
//
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

// print Usage
void printUsage(std::ostream& os) {
     os << "Usage: rgMergeBam (options) --list=<RGAListFile> --out=<outBamFile>\n" << std::endl;
     os << "RGAMerge merges multiple sorted BAM files while appending readgroup taggs\n";
     os << "Required parameters :" << std::endl;
     os << "--out/-o : Output BAM file (sorted)" << std::endl;
     os << "--list/-l : RGAList File. Tab-delimited list consisting of following columns (with headers):" << std::endl;
     os << "\tBAM* : Input BAM file name to be merged" << std::endl;
     os << "\tID* : Unique read group identifier" << std::endl;
     os << "\tSM* : Sample name" << std::endl;
     os << "\tLB : Library name" << std::endl;
     os << "\tDS : Description" << std::endl;
     os << "\tPU : Platform unit" << std::endl;
     os << "\tPI : Predicted median insert size" << std::endl;
     os << "\tCN : Name of sequencing center producing the read" << std::endl;
     os << "\tDT : Date the rn was produced" << std::endl;
     os << "\tPL : Platform/technology used to produce the read" << std::endl;
     os << "\t* (Required fields)" << std::endl;
     os << "Optional parameters : " << std::endl;
     os << "--log/-l : Log file" << std::endl;
     os << "--verbose/-v : Turn on verbose mode" << std::endl;
}

// main function
int main(int argc, char ** argv)
{
  static struct option getopt_long_options[] = 
    {
      // Input options
      { "list", required_argument, NULL, 'l'},
      { "out", required_argument, NULL, 'o'},
      { "verbose", no_argument, NULL, 'v'},
      { "log", required_argument, NULL, 'L'},
      { NULL, 0, NULL, 0 },
    };

  int n_option_index = 0;
  char c;
  bool b_verbose = false;

  std::string s_list, s_out, s_logger;

  while ( ( c = getopt_long(argc, argv, "l:o:vL:", getopt_long_options, &n_option_index) ) != -1 ) {
    switch(c) {
    case 'l':
      s_list = optarg;
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
  if ( optind < argc ) {
    printUsage(std::cerr);
    Logger::gLogger->error("non-option argument exist");
  }

  // check the required arguments are nonempty
  if ( s_list.empty() || s_out.empty() ) {
    printUsage(std::cerr);
    Logger::gLogger->error("At least one of the required argument is missing");
  }

  Logger::gLogger->writeLog("Input list file : %s",s_list.c_str());
  Logger::gLogger->writeLog("Output BAM file : %s",s_out.c_str());
  Logger::gLogger->writeLog("Output log file : %s",s_logger.c_str());
  Logger::gLogger->writeLog("Verbose mode    : %s",b_verbose ? "On" : "Off");
  
  vector<std::string> vs_in_bam_files; // input BAM files
  vector<ReadGroup> v_readgroups;      // readGroups corresponding to BAM file
  vector<ReadGroup> v_uniq_readgroups; // unique readGroups written to header

  // parse the list file and fill the vectors above
  if ( parseListFile(s_list, vs_in_bam_files, v_readgroups, v_uniq_readgroups) == false ) {
    Logger::gLogger->error("Error in parsing the list file %s",s_list.c_str());
  }

  // sanity check
  uint32_t n_bams = vs_in_bam_files.size();
  Logger::gLogger->writeLog("Total of %d BAM files are being merged",n_bams);

  if ( n_bams != v_readgroups.size() ) {
    Logger::gLogger->error("parseListFile gave different size for vs_in_bam_files, v_readgroups");
  }
  else if ( n_bams < 2 ) {
    Logger::gLogger->error("At least two BAM files should exist in the list to merge");
  }

  // create SamFile and SamFileHeader object for each BAM file
  SamFile *p_in_bams = new SamFile[n_bams];
  SamFileHeader *p_headers = new SamFileHeader[n_bams];

  // read each BAM file and its header, 
  // making sure that the headers are identical
  for(uint32_t i=0; i < n_bams; ++i) {
    if ( ! p_in_bams[i].OpenForRead(vs_in_bam_files[i].c_str()) ) {
      Logger::gLogger->error("Cannot open BAM file %s for reading",vs_in_bam_files[i].c_str());
    }
    p_in_bams[i].setSortedValidation(SamFile::COORDINATE);

    p_in_bams[i].ReadHeader(p_headers[i]);
    if ( i > 0 ) {
      if ( ! equalHeaders(p_headers[0], p_headers[i]) ) {
	Logger::gLogger->error("The headers are not identical at index %d",i);
      }
    }
  }

  // first header will be the new header to be written to output
  // adding all possible readGroups to the new header
  for(uint32_t i=0; i < v_uniq_readgroups.size(); ++i) {
    addReadGroupToHeader(p_headers[0], v_uniq_readgroups[i]);
  }

  // Write an output file with new headers
  SamFile bam_out;
  if ( !bam_out.OpenForWrite(s_out.c_str()) ) {
    Logger::gLogger->error("Cannot open BAM file %s for writing",s_out.c_str());
  }
  bam_out.setSortedValidation(SamFile::COORDINATE);
  bam_out.WriteHeader(p_headers[0]);

  // create SamRecords and GenomicCoordinates for each input BAM file
  SamRecord* p_records = new SamRecord[n_bams];
  uint64_t* p_gcoordinates = new uint64_t[n_bams];

  // read the first record for every input BAM file
  for(uint32_t i=0; i < n_bams; ++i) {
    if ( p_in_bams[i].ReadRecord(p_headers[i],p_records[i]) ) {
      if ( p_records[i].isValid(p_headers[i]) ) {
	p_gcoordinates[i] = getGenomicCoordinate(p_records[i]);
      }
      else {
	Logger::gLogger->error("Invalid record found at the first line of file %u. Failure code is %d", i, static_cast<int>(p_in_bams[i].GetFailure()));
      }
    }
    else {
      if ( p_in_bams[i].GetFailure() == SamStatus::NO_MORE_RECS ) {
	// the BAM file has no record
	p_gcoordinates[i] = MAX_GENOMIC_COORDINATE;
      }
      else {
	Logger::gLogger->error("Invalid record found at the first line of file %u. Failure code is %d", i, static_cast<int>(p_in_bams[i].GetFailure()));
      }
    }
  }

  // Routine for writing output BAM file
  uint32_t nWrittenRecords = 0; // number of written BAM records
  while(true) {
    // scan the minimum index of genomic coordinate
    int min_idx = -1;
    uint64_t min_gcoordinate = MAX_GENOMIC_COORDINATE;
    for(uint32_t i=0; i < n_bams; ++i) {
      if ( min_gcoordinate > p_gcoordinates[i] ) {
	min_gcoordinate = p_gcoordinates[i];
	min_idx = static_cast<int>(i);
      }
    }

    // If every file eached EOF, exit the loop
    if ( min_idx < 0 ) break;

    // add readGroup tag to the record to write and write to output BAM file
    //Logger::gLogger->writeLog("%d",min_idx);
    addReadGroupTag(p_records[min_idx], v_readgroups[min_idx]);
    bam_out.WriteRecord(p_headers[0], p_records[min_idx]);
    ++nWrittenRecords;
    if ( nWrittenRecords % 1000000 == 0 ) {
      Logger::gLogger->writeLog("Writing %u records to the output file",nWrittenRecords);
    }

    // Read a record from the input BAM file 
    if ( p_in_bams[min_idx].ReadRecord(p_headers[min_idx], p_records[min_idx]) ) {
      if ( p_records[min_idx].isValid(p_headers[min_idx]) ) {
	p_gcoordinates[min_idx] = getGenomicCoordinate(p_records[min_idx]);
      }
      else { // if invalid record found
	Logger::gLogger->error("Invalid record found at recordCount %d of file %d. Failure code is %d", p_in_bams[min_idx].GetCurrentRecordCount(), min_idx, static_cast<int>(p_in_bams[min_idx].GetFailure()));
      }
    }
    else {
      if ( p_in_bams[min_idx].GetFailure() == SamStatus::NO_MORE_RECS ) {
	p_gcoordinates[min_idx] = MAX_GENOMIC_COORDINATE; // Mark that all record has been read
      }
      else {
	Logger::gLogger->error("Cannot read record at recordCount %d of file %d. Failure code is %d", p_in_bams[min_idx].GetCurrentRecordCount(), min_idx, static_cast<int>(p_in_bams[min_idx].GetFailure()));
      }
    }
  }

  // close files and free allocated memory
  Logger::gLogger->writeLog("Finished writing %d records into the output BAM file",bam_out.GetCurrentRecordCount());
  bam_out.Close();
  for(uint32_t i=0; i < n_bams; ++i) {
    p_in_bams[i].Close();
  }
  delete[] p_records;
  delete[] p_in_bams;
  delete[] p_headers;
  delete[] p_gcoordinates;
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
  char *requiredColumnNames[3] = {"BAM","ID","SM"};
  for(uint32_t i=0; i < 3; ++i) {
    if ( columnDict.find( requiredColumnNames[i] ) == columnDict.end() ) {
      Logger::gLogger->error("Required column name '%s' does not exist in the header of the list file", requiredColumnNames[i]);
    }
  }

  // read the next line and create a readgroup using existing fields
  //vector<std::string> tokens;
  const uint32_t NUM_COLS = 9;
  char* columnNames[NUM_COLS] = {"ID","SM","LB","DS","PU","PI","CN","DT","PL"};
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
