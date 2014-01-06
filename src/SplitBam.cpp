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
#include "Logger.h"
#include "BgzfFileType.h"
#include "PhoneHome.h"

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
  std::cerr << "-n/--noeof : turn off the check for an EOF block at the end of a bam file" << std::endl;
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
      { "noeof", no_argument, NULL, 'n'},
      { "log", required_argument, NULL, 'L'},
      { "noPhoneHome", no_argument, NULL, 'p'},
      { "nophonehome", no_argument, NULL, 'P'},
      { "phoneHomeThinning", required_argument, NULL, 't'},
      { "phonehomethinning", required_argument, NULL, 'T'},
      { NULL, 0, NULL, 0 },
    };

  int n_option_index = 0;
  char c;
  bool b_verbose = false;
  bool noeof = false;
  bool noPhoneHome = false;

  std::string s_in, s_out, s_logger;

  while ( ( c = getopt_long(argc-1, &(argv[1]), "i:o:vLn:", getopt_long_options, &n_option_index) ) != -1 ) {
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
    case 'n':
      noeof = true;
      break;
    case 'L':
      s_logger = optarg;
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
      fprintf(stderr,"ERROR: Unrecognized option %s\n",getopt_long_options[n_option_index].name);
      return(-1);
    }
  }

  if ( s_logger.empty() )
  {
      if(s_out.empty())
      {
          s_logger = "-";
      }
      else
      {
          s_logger = s_out + ".log";
      }
  }
  
  if(!noPhoneHome)
  {
      PhoneHome::checkVersion(getProgramName(), VERSION);
  }
  
  if(noeof)
  {
      // Set that the eof block is not required.
      BgzfFileType::setRequireEofBlock(false);
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
  Logger::gLogger->writeLog("BGFZ EOF indicator : %s",noeof ? "Off" : "On");
  
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
