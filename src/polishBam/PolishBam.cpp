/*
 *  Copyright (C) 2010  Regents of the University of Michigan
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
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include "SamFile.h"
#include "FastaFile.h"
#include "Helper.h"

static Logger* gpLogger;

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

void printUsage(std::ostream& os) {
     os << "Usage: polishBam (options) --in <inBamFile> --out <outBamFile>\n" << std::endl;
     os << "Required parameters : " << std::endl;
     os << "-i/--in : input BAM file" << std::endl;
     os << "-o/--out : output BAM file" << std::endl;
     os << "Optional parameters :" << std::endl;
     os << "-v : turn on verbose mode" << std::endl;
     os << "-l/--log : writes logfile. <outBamFile>.log will be used if value is unspecified" << std::endl;
     os << "--HD : add @HD header line" << std::endl;
     os << "--RG : add @RG header line" << std::endl;
     os << "--PG : add @PG header line" << std::endl;
     os << "-f/--fasta : fasta reference file to compute MD5sums and update SQ tags" << std:: endl;
     os << "--AS : AS tag for genome assembly identifier" << std::endl;
     os << "--UR : UR tag for @SQ tag (if different from --fasta)" << std::endl;
     os << "--SP : SP tag for @SQ tag" << std:: endl;
     os << "--checkSQ : check the consistency of SQ tags (SN and LN) with existing header lines. Must be used with --fasta option" << std::endl;
     os << "\n" << std::endl;
}

void checkHeaderStarts(std::vector<std::string>& headerLines, const char* type) {
  for(uint32_t i=0; i < headerLines.size(); ++i) {
    if ( headerLines[i].find(type) != 0 ) {
      gpLogger->error("The following header line does not start with '%s'\n%s",type,headerLines[i].c_str());
    }
  }
}

int main(int argc, char ** argv)
{
  gpLogger = new Logger;

  static struct option getopt_long_options[] = 
    {
      // Input options
      { "fasta", required_argument, NULL, 'f'},
      { "in", required_argument, NULL, 'i'},
      { "out", required_argument, NULL, 'o'},
      { "verbose", no_argument, NULL, 'v'},
      { "log", optional_argument, NULL, 'l'},
      { "clear", no_argument, NULL, 0},
      { "AS", required_argument, NULL, 0},
      { "UR", required_argument, NULL, 0},
      { "SP", required_argument, NULL, 0},
      { "HD", required_argument, NULL, 0},
      { "RG", required_argument, NULL, 0},
      { "PG", required_argument, NULL, 0},
      { "checkSQ", no_argument, NULL, 0},
      { NULL, 0, NULL, 0 },
    };

  int n_option_index = 0, c;
  
  std::string sAS, sUR, sSP, sFasta, sInFile, sOutFile, sLogFile;
  bool bClear, bCheckSQ, bVerbose;
  std::vector<std::string> vsHDHeaders, vsRGHeaders, vsPGHeaders;

  bCheckSQ = bVerbose = false;
  bClear = true;

  while ( (c = getopt_long(argc, argv, "vf:i:o:l:", getopt_long_options, &n_option_index)) != -1 ) {
    std::cout << getopt_long_options[n_option_index].name << "\t" << optarg << std::endl;
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
      if ( getopt_long_options[n_option_index].has_arg ) {
	sLogFile = optarg;
      }
      else {
	sLogFile = "__NONE__";
      }
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
    else if ( strcmp(getopt_long_options[n_option_index].name,"checkSQ") == 0 ) {
      bCheckSQ = true;
    }
    else {
      std::cerr << "Error: Unrecognized option " << getopt_long_options[n_option_index].name << std::endl;
      abort();
    }
  }

  if ( optind < argc ) {
    printUsage(std::cerr);
    gpLogger->error("non-option argument %s exist ",argv[optind]);
  }

  if ( sInFile.empty() || sOutFile.empty() ) {
    printUsage(std::cerr);
    gpLogger->error("Input and output files are required");
  }

  if ( sLogFile.compare("__NONE__") == 0 ) {
    sLogFile = (sOutFile + ".log");
  }

  gpLogger->open(sLogFile.c_str(), bVerbose);

  if ( ( bCheckSQ ) && ( sFasta.empty() ) ) {
    printUsage(std::cerr);
    gpLogger->error("--checkSQ option must be used with --fasta option");
  }

  // check whether each header line starts with a correct tag
  checkHeaderStarts(vsHDHeaders, "@HD\t");
  checkHeaderStarts(vsRGHeaders, "@RG\t");
  checkHeaderStarts(vsPGHeaders, "@PG\t");

  gpLogger->write_log("Arguments in effect:");
  gpLogger->write_log("\t--in [%s]",sInFile.c_str());
  gpLogger->write_log("\t--out [%s]",sOutFile.c_str());
  gpLogger->write_log("\t--log [%s]",sLogFile.c_str());
  gpLogger->write_log("\t--fasta [%s]",sFasta.c_str());
  gpLogger->write_log("\t--AS [%s]",sAS.c_str());
  gpLogger->write_log("\t--UR [%s]",sUR.c_str());
  gpLogger->write_log("\t--SP [%s]",sSP.c_str());
  gpLogger->write_log("\t--checkSQ [%s]",bClear ? "ON" : "OFF" );
  if ( vsHDHeaders.empty() ) {
    gpLogger->write_log("\t--HD []");
  }
  else {
    gpLogger->write_log("\t--HD [%s]",vsHDHeaders[0].c_str());
  }
  if ( vsRGHeaders.empty() ) {
    gpLogger->write_log("\t--RG []");
  }
  else {
    gpLogger->write_log("\t--RG [%s]",vsRGHeaders[0].c_str());
  }
  if ( vsPGHeaders.empty() ) {
    gpLogger->write_log("\t--PG []");
  }
  else {
    for(uint32_t i=0; i < vsPGHeaders.size(); ++i) {
      gpLogger->write_log("\t--PG [%s]",vsPGHeaders[i].c_str());
    }
  }

  if ( (vsHDHeaders.empty() ) && ( vsRGHeaders.empty() ) && ( vsPGHeaders.empty() ) && ( !bClear ) && ( sFasta.empty() ) ) {
    gpLogger->warning("No option is in effect for modifying BAM files. The input and output files will be identical");
  }

  if ( ( vsHDHeaders.size() > 1 ) || ( vsRGHeaders.size() > 1 ) ) {
    gpLogger->error("HD and RG headers cannot be multiple");
  }

  FastaFile fastaFile;
  if ( ! sFasta.empty() ) {
    if ( fastaFile.open(sFasta.c_str()) ) {
      gpLogger->write_log("Reading the reference file %s",sFasta.c_str());
      fastaFile.readThru();
      fastaFile.close();
      gpLogger->write_log("Finished reading the reference file %s",sFasta.c_str());      
    }
    else {
      gpLogger->error("Failed to open reference file %s",sFasta.c_str());
    }
  }

  SamFile samIn;
  SamFile samOut;

  if ( ! samIn.OpenForRead(sInFile.c_str()) ) {
    gpLogger->error("Cannot open BAM file %s for reading - %s",sInFile.c_str(), SamStatus::getStatusString(samIn.GetStatus()) );
  }
  if ( ! samOut.OpenForWrite(sOutFile.c_str()) ) {
    gpLogger->error("Cannot open BAM file %s for writing - %s",sOutFile.c_str(), SamStatus::getStatusString(samOut.GetStatus()) );
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
      gpLogger->error("# of @SQ tags are different from the original BAM and the reference file");
    }

    // iterator over all @SQ objects
    for(unsigned int i=0; i < numSQ; ++i) {
      pSamHeaderRecord = samHeader.getSQ(fastaFile.vsSequenceNames[i].c_str());
      if ( fastaFile.vsSequenceNames[i].compare(pSamHeaderRecord->getTagValue("SN")) != 0 ) {
	gpLogger->error("SequenceName is not identical between fasta and input BAM file");
      }
      else if ( static_cast<int>(fastaFile.vnSequenceLengths[i]) != atoi(pSamHeaderRecord->getTagValue("LN")) ) {
	gpLogger->error("SequenceLength is not identical between fasta and input BAM file");
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
    gpLogger->write_log("Finished checking the consistency of SQ tags");
  }
  else {
    gpLogger->write_log("Skipped checking the consistency of SQ tags");
  }

  // go over the headers again, 
  // assuming order of HD, SQ, RG, PG, and put proper tags at the end of the original tags

  gpLogger->write_log("Creating the header of new output file");
  //SamFileHeader outHeader;
  samHeader.resetHeaderRecordIter();

  for(unsigned int i=0; i < vsHDHeaders.size(); ++i) {
    samHeader.addHeaderLine(vsHDHeaders[i].c_str());
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

  for(unsigned int i=0; i < vsRGHeaders.size(); ++i) {
    samHeader.addHeaderLine(vsRGHeaders[i].c_str());
  }

  for(unsigned int i=0; i < vsPGHeaders.size(); ++i) {
    samHeader.addHeaderLine(vsPGHeaders[i].c_str());
  }

  samOut.WriteHeader(samHeader);
  gpLogger->write_log("Adding %d HD, %d RG, and %d PG headers",vsHDHeaders.size(), vsRGHeaders.size(), vsPGHeaders.size());
  gpLogger->write_log("Finished writing output headers");

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
  
  gpLogger->write_log("Writing output BAM file");
  SamRecord samRecord;
  while (samIn.ReadRecord(samHeader, samRecord) == true) {
    if ( !sRGID.empty() ) {
      if ( samRecord.addTag("RG",'Z',sRGID.c_str()) == false ) {
	gpLogger->error("Failed to add a RG tag %s",sRGID.c_str());
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
  gpLogger->write_log("Successfully written %d records",samIn.GetCurrentRecordCount());
  delete gpLogger;
  return 0;
}
