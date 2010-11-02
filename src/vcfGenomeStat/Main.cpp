#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>

#include "GenomeSequence.h"
#include "base/logger.h"
#include "base/string_tokenizer.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

Logger* Logger::gLogger = NULL;

int main(int argc, char** argv) {
  std::string sRef, sInFile, sOutFile, sLogFile;
  bool bVerbose;
  uint32_t win;
  uint32_t maxLineLength;

  try {
    TCLAP::CmdLine cmd("Command description message", ' ', "1.0");
    TCLAP::ValueArg<std::string> argReference("r","reference","Karma's reference sequence",true,"","string");
    TCLAP::ValueArg<std::string> argIn("i","in","Input VCF file",true,"","string");
    TCLAP::ValueArg<std::string> argOut("o","out","Output file",true,"","string");
    TCLAP::ValueArg<std::string> argLog("l","log","Log file - [out].log is default",false,"","string");
    TCLAP::SwitchArg switchVerbose("v","verbose","Turn on verbose mode",cmd,false);
    TCLAP::ValueArg<uint32_t> argWindow("w","win","Window size for HomPolymerRun (default:10)",false,10,"integer");
    TCLAP::ValueArg<uint32_t> argMaxLineLength("L","max-line-length","Maximun Line Length of the VCF file (default:1048576)",false,1048576,"integer");

    cmd.add(argReference);
    cmd.add(argIn);
    cmd.add(argOut);
    cmd.add(argLog);
    cmd.add(argWindow);
    cmd.add(argMaxLineLength);

    cmd.parse(argc, argv);

    sRef = argReference.getValue();
    sInFile = argIn.getValue();
    sOutFile = argOut.getValue();
    sLogFile = argLog.getValue();
    bVerbose = switchVerbose.getValue();
    win = argWindow.getValue();
    maxLineLength = argMaxLineLength.getValue();

    if ( sLogFile.empty() ) {
      sLogFile = sOutFile + ".log";
    }
    Logger::gLogger = new Logger(sLogFile.c_str(), bVerbose);
    
    //    Logger::gLogger->writeLog("Arguments in effect: \n%s",cmd.toString().c_str());
    Logger::gLogger->writeLog("Arguments in effect: \n");
    std::ostringstream oss;
    std::list<TCLAP::Arg*> argList = cmd.getArgList();
    for(std::list<TCLAP::Arg*>::iterator i=argList.begin(); i != argList.end(); ++i) {
		oss << "\t" << (*i)->toString() << std::endl;
    }
    Logger::gLogger->writeLog("%s", oss.str().c_str());
  }
  catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    abort();
  }

  GenomeSequence genomeSequence( sRef );

  char* buf = new char[maxLineLength];
  char* leftSequence = new char[win+1];
  char* rightSequence = new char[win+1];

  std::vector<std::string> vChromosomes;
  std::vector<uint32_t> vCoordinates;
  std::vector<char> vRefBases;
  std::vector<std::string> vLeftSequences;
  std::vector<std::string> vRightSequences;

  std::map<std::string,uint32_t> mLeftCounts;
  std::map<std::string,uint32_t> mRightCounts;
  std::vector<uint32_t> vHomPolRuns;
  
  Logger::gLogger->writeLog("Reading input VCF file %s..",sInFile.c_str());

  std::ifstream ifsVcf(sInFile.c_str());
  if ( ! ifsVcf.is_open() ) {
    Logger::gLogger->error("Error opening file %s",sInFile.c_str());
  }

  while( !ifsVcf.getline(buf, maxLineLength).eof() ) {
    if ( buf[0] == '#' ) {
      // comment line
    }
    else {
      CStringTokenizer t(buf, buf + ifsVcf.gcount(), "\t ");

      t.GetNext();
      std::string sChrom = t.token().c_str();

      t.GetNext();
      uint32_t nPos = static_cast<uint32_t>(atoi(t.token().c_str()));

      genomeIndex_t markerIndex = genomeSequence.getGenomePosition(sChrom.c_str(), nPos);
      bool leftA = true, leftT = true, rightA = true, rightT = true;
      uint32_t leftRun = 0, rightRun = 0;

      if ( markerIndex < win ) { 
	Logger::gLogger->error("Marker %s:%u is way too close to the end of genome given window size %u",sChrom.c_str(),nPos,win);
      }

      for(uint32_t i=0; i < win; ++i) {
	rightSequence[i] = genomeSequence[markerIndex+i+1];
	leftSequence[win-i-1] = genomeSequence[markerIndex-i-1];

	switch(rightSequence[i]) {
	case 'A':
	  if ( rightA ) {
	    rightRun = i+1;
	  }
	  rightT = false;
	  break;
	case 'T':
	  if ( rightT ) {
	    rightRun = i+1;
	  }
	  rightA = false;
	  break;
	default:
	  rightT = rightA = false;
	}

	switch(leftSequence[win-i-1]) {
	case 'A':
	  if ( leftA ) {
	    leftRun = i+1;
	  }
	  leftT = false;
	  break;
	case 'T':
	  if ( leftT ) {
	    leftRun = i+1;
	  }
	  leftA = false;
	  break;
	default:
	  leftT = leftA = false;
	}
      }

      leftSequence[win] = rightSequence[win] = '\0';

      //fprintf(fpOut,"%s\t%u\t%d\t%s\t%c\t%s\n",sChrom.c_str(), nPos, (rightRun > leftRun) ? rightRun : leftRun, leftSequence, genomeSequence[markerIndex], rightSequence);

      vChromosomes.push_back(sChrom);
      vCoordinates.push_back(nPos);
      vRefBases.push_back(genomeSequence[markerIndex]);
      vLeftSequences.push_back(leftSequence);
      vRightSequences.push_back(rightSequence);
      ++(mLeftCounts[leftSequence]);
      ++(mRightCounts[rightSequence]);
      vHomPolRuns.push_back( (rightRun > leftRun) ? rightRun : leftRun );
    }
  }

  Logger::gLogger->writeLog("Writing output file %s",sOutFile.c_str());
  FILE* fpOut = fopen(sOutFile.c_str(),"w");

  for(uint32_t i=0; i < vChromosomes.size(); ++i) {
    uint32_t leftCount = mLeftCounts[vLeftSequences[i]];
    uint32_t rightCount = mRightCounts[vRightSequences[i]];
    fprintf(fpOut,"%s\t%u\t%u\t%u\t%.5lf\t%s\t%c\t%s\n",vChromosomes[i].c_str(),vCoordinates[i],vHomPolRuns[i], (leftCount > rightCount) ? leftCount : rightCount, static_cast<double>((leftCount > rightCount) ? leftCount : rightCount)/static_cast<double>(vChromosomes.size()), vLeftSequences[i].c_str(), vRefBases[i], vRightSequences[i].c_str());
  }
  
  delete[] buf;
  delete[] leftSequence;
  delete[] rightSequence;
}
