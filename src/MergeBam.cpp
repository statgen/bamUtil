/*
 *  Copyright (C) 2010-2014  Hyun Min Kang,
 *                           Regents of the University of Michigan
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
#include "SamFile.h"
#include "MergeBam.h"
#include "Logger.h"
#include "PhoneHome.h"

////////////////////////////////////////////////////////////////////////
// MergeBam : Merge multiple BAM files appending ReadGroup IDs if necessary
//
// mergeBam merges multiple sorted BAM files into one BAM file like
//   'samtools merge' command, but merges BAM headers.
// (1) check that all non-RG header lines are identical across the BAM files
// (2) add all RG lines from all BAMs into the output BAM.
// (3) Ensure that input/output BAM records are sorted
// Optionally add Read Groups
// (1) Add @RG headers from a tabular input file containing the fields' info
// (2) Add RG:Z:[RGID] tag for each record based on the source BAM file
//
///////////////////////////////////////////////////////////////////////

// ReadGroup w/ RG ID and full header line
class ReadGroup {
public:
  std::string s_id;
  std::string s_header_line;
};

// Global variables
uint64_t MAX_GENOMIC_COORDINATE = 0xffffffffffffffffULL;
uint64_t UNMAPPED_GENOMIC_COORDINATE = 0xfffffffffffffffeULL;

// function declarations
bool parseListFile(std::string& listFile, vector<std::string>& bamFiles, vector<ReadGroup>& readGroups, vector<ReadGroup>& uniqReadGroups);
void parseOutRG(SamFileHeader& header, std::string& noRgPgString, SamFileHeader* newHeader, bool ignorePI);
uint64_t getGenomicCoordinate(SamRecord& r);
void addReadGroupToHeader(SamFileHeader& header, ReadGroup& rg);
void addReadGroupTag(SamRecord& record, ReadGroup& rg);
uint32_t addTokenizedStrings(const std::string& str, const std::string& delimiters, vector<std::string>& tokens);

MergeBam::MergeBam()
    : BamExecutable(),
      myRegionArray(),
      myRegionArrayIndex(-1),
      myRegionFile(NULL)
{
}

void MergeBam::mergeBamDescription()
{
    std::cerr << " mergeBam - merge multiple BAMs and headers appending ReadGroupIDs if necessary" << std::endl;
}


void MergeBam::description()
{
    mergeBamDescription();
}


void MergeBam::usage()
{
    BamExecutable::usage();
     std::cerr << "Usage: mergeBam [-v] [--log logFile] [--ignorePI] --list <listFile> --out <outFile>\n" << std::endl;
     std::cerr << "Required parameters :" << std::endl;
     std::cerr << "--out/-o : Output BAM file (sorted)" << std::endl;
     std::cerr << "--in/-i  : BAM file to be input, must be more than one of these options." << std::endl;
     std::cerr << "            cannot be used with --list/-l" << std::endl;
     std::cerr << "--list/-l : RGAList File. Tab-delimited list consisting of following columns (with headers):" << std::endl;
     std::cerr << "\tBAM* : Input BAM file name to be merged" << std::endl;
     std::cerr << "\tID* : Unique read group identifier" << std::endl;
     std::cerr << "\tSM* : Sample name" << std::endl;
     std::cerr << "\tLB : Library name" << std::endl;
     std::cerr << "\tDS : Description" << std::endl;
     std::cerr << "\tPU : Platform unit" << std::endl;
     std::cerr << "\tPI : Predicted median insert size" << std::endl;
     std::cerr << "\tCN : Name of sequencing center producing the read" << std::endl;
     std::cerr << "\tDT : Date the rn was produced" << std::endl;
     std::cerr << "\tPL : Platform/technology used to produce the read" << std::endl;
     std::cerr << "\t* (Required fields)" << std::endl;
     std::cerr << "Optional parameters : " << std::endl;
     std::cerr << "--regions/-r : list of intervals, '<chr>:<start>-<end>', to merge separated by commas, ','\n";
     std::cerr << "--regionFile/-R : file containing list of intervals, '<chr>:<start>-<end>', to merge, one per line\n";
     std::cerr << "--ignorePI/-I : Ignore the RG PI field when comparing headers\n";
     std::cerr << "--log/-L : Log file" << std::endl;
     std::cerr << "--verbose/-v : Turn on verbose mode" << std::endl;
}

// main function
int MergeBam::execute(int argc, char ** argv)
{
  static struct option getopt_long_options[] = 
    {
      // Input options
      { "list", required_argument, NULL, 'l'},
      { "in", required_argument, NULL, 'i'},
      { "out", required_argument, NULL, 'o'},
      { "verbose", no_argument, NULL, 'v'},
      { "log", required_argument, NULL, 'L'},
      { "ignorePI", no_argument, NULL, 'I'},
      { "regions", required_argument, NULL, 'r'},
      { "regionFile", required_argument, NULL, 'R'},
      { "noPhoneHome", no_argument, NULL, 'p'},
      { "nophonehome", no_argument, NULL, 'P'},
      { "phoneHomeThinning", required_argument, NULL, 't'},
      { "phonehomethinning", required_argument, NULL, 'T'},
      { NULL, 0, NULL, 0 },
    };

  // Adjust the arguments since it is called as ./bam mergeBam instead of
  // just mergeBam.
  --argc;
  ++argv;

  int n_option_index = 0;
  char c;
  bool b_verbose = false;
  bool ignorePI = false;
  bool noPhoneHome = false;
  String regions = "";
  std::string regionFile = "";
  vector<std::string> vs_in_bam_files; // input BAM files

  std::string s_list, s_out, s_logger;

  while ( ( c = getopt_long(argc, argv, "l:i:o:vIL:r:R:", getopt_long_options, &n_option_index) ) != -1 ) {
    switch(c) {
    case 'i':
      vs_in_bam_files.push_back(optarg);
      break;
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
    case 'I':
      ignorePI = true;
      break;
    case 'r':
      regions = optarg;
      break;
    case 'R':
      regionFile = optarg;
      break;
    case 'p':
    case 'P':
      noPhoneHome = true;
      break;
    case 't':
    case 'T':
        PhoneHome::allThinning = atoi(optarg);
      break;
    default:
      fprintf(stderr,"Unrecognized option %s\n",getopt_long_options[n_option_index].name);
      return(-1);
    }
  }

  if(!noPhoneHome)
  {
      PhoneHome::checkVersion(getProgramName(), VERSION);
  }
  
  if ( s_logger.empty() ) {
      if(s_out.empty())
      {
          s_logger = "-";
      }
      else
      {
          s_logger = s_out + ".log";
      }
  }

  // create a logger object, now possible to write logs/warnings/errors
  Logger::gLogger = new Logger(s_logger.c_str(), b_verbose);

  // every argument must correspond to an option
  if ( optind < argc ) {
    usage();
    Logger::gLogger->error("non-option argument exist");
  }

  // check the required arguments are nonempty
  if ( (vs_in_bam_files.empty() && s_list.empty()) || s_out.empty() ) {
    usage();
    Logger::gLogger->error("At least one of the required argument is missing");
  }

  if(!vs_in_bam_files.empty() && !s_list.empty())
  {
      Logger::gLogger->error("Cannot specify both --in/-i and --list/-l");
  }

  // Check that both region and regionfile are not specified.
  if(!regions.IsEmpty() && !regionFile.empty())
  {
      Logger::gLogger->error("Cannot specify both --regions/-r and --regionFile/-R");
  }

  if(!regions.IsEmpty())
  {
      // Setup the region array
      myRegionArray.ReplaceColumns(regions, ',');
      myRegionArrayIndex = 0;
  }
  else if(!regionFile.empty())
  {
      myRegionFile = ifopen(regionFile.c_str(), "r");
      if(myRegionFile == NULL)
      {
          Logger::gLogger->error("Failed to open region file, %s", regionFile.c_str());
      }
  }

  if(!s_list.empty())
  {
      Logger::gLogger->writeLog("Input list file : %s",s_list.c_str());
  }
  else
  {
      std::string bamList = "";
      for(unsigned int i = 0; i < vs_in_bam_files.size(); i++)
      {
          if(i != 0)
          {
              bamList += ", ";
          }
          bamList += vs_in_bam_files[i];
      }
      Logger::gLogger->writeLog("Input list file : %s", bamList.c_str());
  }
  Logger::gLogger->writeLog("Output BAM file : %s",s_out.c_str());
  Logger::gLogger->writeLog("Output log file : %s",s_logger.c_str());
  Logger::gLogger->writeLog("Verbose mode    : %s",b_verbose ? "On" : "Off");
  
  vector<ReadGroup> v_readgroups;      // readGroups corresponding to BAM file
  vector<ReadGroup> v_uniq_readgroups; // unique readGroups written to header

  // If the list file is being used instead of the individual bams, parse it.
  if(!s_list.empty())
  {
      // parse the list file and fill the vectors above
      if ( parseListFile(s_list, vs_in_bam_files, v_readgroups, v_uniq_readgroups) == false ) {
          Logger::gLogger->error("Error in parsing the list file %s",s_list.c_str());
      }
      if ( vs_in_bam_files.size() != v_readgroups.size() ) {
          Logger::gLogger->error("parseListFile gave different size for vs_in_bam_files, v_readgroups: %d, %d", vs_in_bam_files.size(), v_readgroups.size());
      }
  }

  // sanity check
  uint32_t n_bams = vs_in_bam_files.size();
  Logger::gLogger->writeLog("Total of %d BAM files are being merged",n_bams);

  if ( n_bams < 2 )
  {
      Logger::gLogger->error("At least two BAM files must be specified for merging");
  }

  // create SamFile and SamFileHeader object for each BAM file
  SamFile *p_in_bams = new SamFile[n_bams];
  SamFileHeader *p_headers = new SamFileHeader[n_bams];

  // read each BAM file and its header, 
  // making sure that the headers are identical

  std::string firstHeaderNoRGPG = "";
  std::string headerNoRGPG = "";
  SamFileHeader newHeader;

  std::string firstHeaderString = "";
  for(uint32_t i=0; i < n_bams; ++i)
  {
      if ( ! p_in_bams[i].OpenForRead(vs_in_bam_files[i].c_str()) )
      {
          Logger::gLogger->error("Cannot open BAM file %s for reading",vs_in_bam_files[i].c_str());
      }
      p_in_bams[i].setSortedValidation(SamFile::COORDINATE);
      
      p_in_bams[i].ReadHeader(p_headers[i]);

      if((myRegionFile != NULL) || (myRegionArrayIndex != -1))
      {
          // Reading just specific regions, so read each index file.
          p_in_bams[i].ReadBamIndex();
      }

      // Extract the RGs from this header.
      if(i == 0)
      {
          // First header, so store it as the first header
          newHeader = p_headers[i];
          // Determine the header without RG.
          parseOutRG(p_headers[i], firstHeaderNoRGPG, NULL, ignorePI);
      }
      else
      {
          parseOutRG(p_headers[i], headerNoRGPG, &newHeader, ignorePI);
          if(firstHeaderNoRGPG != headerNoRGPG)
          {
              Logger::gLogger->error("The headers are not identical at index %d",i);
          }
          if(newHeader.getReferenceInfo() != p_headers[i].getReferenceInfo())
          {
              Logger::gLogger->error("The headers are not identical at index %d",i);
          }
      }
  }

  // first header will be the new header to be written to output
  // adding all possible readGroups to the new header
  for(uint32_t i=0; i < v_uniq_readgroups.size(); ++i)
  {
    addReadGroupToHeader(newHeader, v_uniq_readgroups[i]);
  }

  // Write an output file with new headers
  SamFile bam_out;
  if ( !bam_out.OpenForWrite(s_out.c_str()) )
  {
    Logger::gLogger->error("Cannot open BAM file %s for writing",s_out.c_str());
  }
  bam_out.setSortedValidation(SamFile::COORDINATE);
  bam_out.WriteHeader(newHeader);

  // create SamRecords and GenomicCoordinates for each input BAM file
  SamRecord* p_records = new SamRecord[n_bams];
  uint64_t* p_gcoordinates = new uint64_t[n_bams];

  // Loop through, processing each section.
  while(getNextSection(p_in_bams, n_bams))
  {
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
          
          // If every file reached EOF, exit the loop
          if ( min_idx < 0 ) break;
          
          
          // If adding read groups, add the tag.
          if(!v_readgroups.empty())
          {
              // add readGroup tag to the record to write and write to output BAM file
              //Logger::gLogger->writeLog("%d",min_idx);
              addReadGroupTag(p_records[min_idx], v_readgroups[min_idx]);
          }
          bam_out.WriteRecord(newHeader, p_records[min_idx]);
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
  const char* columnNames[NUM_COLS] =
      {"ID","SM","LB","DS","PU","PI","CN","DT","PL"};
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


void parseOutRG(SamFileHeader& header, std::string& noRgPgString, SamFileHeader* newHeader, bool ignorePI)
{
    noRgPgString.clear();
    // strings for comparing if two RGs with same ID are the same.
    static std::string prevString = "";
    static std::string newString = "";

    SamHeaderRecord* rec = header.getNextHeaderRecord();
    while(rec != NULL)
    {
        if(rec->getType() == SamHeaderRecord::RG)
        {
            if(newHeader != NULL)
            {
                // This is an RG line.
                // First check if this RG is already included in the new header.
                SamHeaderRG* prevRG = newHeader->getRG(rec->getTagValue("ID"));
                
                if(prevRG != NULL)
                {
                    // This RG already exists, check that they are the same.
                    // If they are the same, there is nothing to do.
                    bool status = true;
                    prevString.clear();
                    newString.clear();
                    status &= prevRG->appendString(prevString);
                    status &= rec->appendString(newString);

                    if(prevString != newString)
                    {
                        if(!ignorePI)
                        {
                            Logger::gLogger->error("Failed to add readgroup to "
                                                   "header, duplicate, but "
                                                   "non-identical RG ID, %s\n"
                                                   "prev:\t(%s)\nnew:\t(%s)",
                                                   rec->getTagValue("ID"),
                                                   prevString.c_str(),
                                                   newString.c_str());
                        }
                        else
                        {
                            // Check for a PI string.
                            size_t prevPIStart = prevString.find("PI:");
                            size_t newPIStart = newString.find("PI:");

                            // If they are both npos, then PI was not found
                            // so fail.
                            if((prevPIStart == std::string::npos) &&
                               (newPIStart == std::string::npos))
                            {
                                // They are not identical, so report an error.
                                Logger::gLogger->error("Failed to add readgroup"
                                                       " to header, duplicate,"
                                                       " but non-identical RG"
                                                       " ID, %s\n"
                                                       "prev:\t(%s)\nnew:\t(%s)",
                                                       rec->getTagValue("ID"),
                                                       prevString.c_str(),
                                                       newString.c_str());
                            }
                            else
                            {
                                // PI found in one or both strings.
                                size_t prevPIEnd;
                                size_t newPIEnd;
                                if(prevPIStart == std::string::npos)
                                {
                                    // new string has PI, so compare to the start of that.
                                    prevPIStart = newPIStart;
                                    prevPIEnd = newPIStart;
                                }
                                else
                                {
                                    prevPIEnd = prevString.find('\t', prevPIStart) + 1;
                                }
                                if(newPIStart == std::string::npos)
                                {
                                    // new string has PI, so compare to the start of that.
                                    newPIStart = prevPIStart;
                                    newPIEnd = newPIStart;
                                }
                                else
                                {
                                    newPIEnd = newString.find('\t', newPIStart) + 1;
                                }
                                // Compare before PI.
                                if((newString.compare(0, newPIStart, prevString, 0, prevPIStart) != 0) ||
                                   (newString.compare(newPIEnd, std::string::npos, prevString,
                                                      prevPIEnd, std::string::npos) != 0))
                                {
                                    // They are not identical, so report an error.
                                    Logger::gLogger->error("Failed to add readgroup to header, "
                                                           "duplicate, but non-identical RG ID, %s, "
                                                           "even when ignoring PI\n"
                                                           "prev:\t(%s)\nnew:\t(%s)",
                                                           rec->getTagValue("ID"),
                                                           prevString.c_str(),
                                                           newString.c_str());
                                }
                                else
                                {
                                    Logger::gLogger->warning("Warning: ignoring non-identical PI field "
                                                             "for RG ID, %s",
                                                             rec->getTagValue("ID"));
                                }
                            }
                        }
                    }
                }
                else
                {
                    // This RG does not exist yet, so add it to the new header.
                    if(!newHeader->addRecordCopy((SamHeaderRG&)(*rec)))
                    {
                        // Failed to add the RG, exit.
                        Logger::gLogger->error("Failed to add readgroup to header, %s",
                                               newHeader->getErrorMessage());
                    }
                }
            }
        }
        else if(rec->getType() == SamHeaderRecord::PG)
        {
            if(newHeader != NULL)
            {
                // This is a PG line.
                // First check if this PG is already included in the new header.
                SamHeaderPG* prevPG = newHeader->getPG(rec->getTagValue("ID"));
                
                if(prevPG != NULL)
                {
                    // This PG already exists, check if they are the same.
                    // If they are the same, there is nothing to do.
                    bool status = true;
                    prevString.clear();
                    newString.clear();
                    status &= prevPG->appendString(prevString);
                    status &= rec->appendString(newString);
                    if(prevString != newString)
                    {
                        // They are not identical, ignore for now.
                        // TODO: change the ID, and add it.
                        Logger::gLogger->warning("Warning: dropping duplicate, "
                                                 "but non-identical PG ID, %s",
                                                 rec->getTagValue("ID"));
                    }
                }
                else
                {
                    // This PG does not exist yet, so add it to the new header.
                    if(!newHeader->addRecordCopy((SamHeaderPG&)(*rec)))
                    {
                        // Failed to add the PG, exit.
                        Logger::gLogger->error("Failed to add PG to header, %s",
                                               newHeader->getErrorMessage());
                    }
                }
            }
        }
        else
        {
            rec->appendString(noRgPgString);
        }
        rec = header.getNextHeaderRecord();
    }

    // Append the comments.
    header.appendCommentLines(noRgPgString);
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

// add a readgroup tag to each SamRecord
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


bool MergeBam::getNextSection(SamFile *in_bams, uint32_t numBams)
{
    static bool done = false;

    if(done == true)
    {
        return(false);
    }

    if((myRegionFile == NULL) && (myRegionArrayIndex == -1))
    {
        // Not reading by section, so don't set a section, just return true.
        // Set the done flag so on the next call it will return false -
        //  done = false means nothing more to process.
        done = true;
        return(true);
    }

    // Get the next region from the specified file/string
    String regionStr;

    if(myRegionArrayIndex != -1)
    {
        if(myRegionArrayIndex >= myRegionArray.Length())
        {
            // Done processing, so set done and return false - no more regions.
            done = true;
            return(false);
        }
        // Get the next region string, then increment the index
        // for the next loop.
        regionStr = myRegionArray[myRegionArrayIndex++];
    }
    else
    {
        // Read the string from the file.
        if(regionStr.ReadLine(myRegionFile) < 0)
        {
            // error/EOF.
            if(myRegionFile->ifeof())
            {
                // done processing, so set done and return false
                // no more regions.
                done = true;
                return(false);
            }
            else
            {
                // TODO - other error.
                Logger::gLogger->error("Failed while reading %s",
                                       myRegionFile->getFileName());
            }
        }
    }

    // Parse the region string
    String chr;
    int start = -1;
    int end = -1;
    int chrStrEnd = regionStr.FastFindChar(':',1);
    if(chrStrEnd < 1)
    {
        // This region is just a chromosome.
        chr = regionStr;
    }
    else
    {
        chr = regionStr.Left(chrStrEnd);
        // get the start position.
        int startStrEnd = regionStr.FastFindChar('-',chrStrEnd);
        String startStr;
        if(startStrEnd < chrStrEnd)
        {
            // No end region
            startStr = regionStr.SubStr(chrStrEnd+1);
        }
        else
        {
            startStr = regionStr.Mid(chrStrEnd+1, startStrEnd-1);
            String endStr = regionStr.SubStr(startStrEnd+1);
            // Convert to integer.
            if(!endStr.AsInteger(end) || end < 1)
            {
                String errorStr = "Error: Invalid region string, '";
                errorStr += regionStr;
                errorStr += "', the end position, '";
                errorStr += endStr;
                errorStr += "', is not an integer >= 1.";
                error(errorStr.c_str());
            }
        }
        // Convert the start string to an integer
        if(!startStr.AsInteger(start) || start < 1)
        {
            String errorStr = "Error: Invalid region string, '";
            errorStr += regionStr;
            errorStr += "', the start position, '";
            errorStr += startStr;
            errorStr += "', is not an integer >= 1.";
            error(errorStr.c_str());
        }
        // Subtract 1 from start since it is 1-based inclusive, but
        // SetReadSection expects 0-based inclusive
        --start;
        // Do not need to subtract from end since it is 1-based 
        // exclusive which is the same as 0-based exclusive that
        // SetReadSection expects.
    }
  
    // Set the next section in every bam.
    for(uint32_t i = 0; i < numBams; ++i)
    {
        // SetReadSection start parameter is 0-based inclusive, 
        // start is 1-based inclusive, so subtract 1 from start
        // SetReadSection end parameter is 0-based exclusive,
        // end is 1-based inclusive, which is the same as 0-based exclusive
        in_bams[i].SetReadSection(chr.c_str(), start, end);
    }
    return(true);
}
