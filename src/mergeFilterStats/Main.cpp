////////////////////////////////////////////////////////////////////// 
// mergeFilterStats/Main.cpp 
// (c) 2010 Hyun Min Kang
// 
// This file is distributed as part of the vcfCooker source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile vcfCooker
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Sunday February 13, 2011

#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <vector>

#include "Parameters.h"
#include "InputFile.h"
#include "Error.h"
#include "Logger.h"
#include "VcfFile.h"
#include "FilterStat.h"

Logger* Logger::gLogger = NULL;

namespace
{
  bool appendStatVcf( FilterStat* pStat, const char* file )
   {
     return pStat->appendStatVcf(file);
   }
}


int main(int argc, char** argv) {
  printf("mergeFilterStats 1.0.0 -- Merge VCF INFO stats\n"
	 "(c) 2010 Hyun Min Kang\n\n");

  String sAnchorVcf;
  String sPrefix;
  String sSuffix;
  String sIndexf;
  String sOutVcf;
  int nThreads = 1;
  int indexSkip = 2;
  bool bVerbose = true;

  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Options")
    LONG_STRINGPARAMETER("anchor",&sAnchorVcf)
    LONG_STRINGPARAMETER("prefix",&sPrefix)
    LONG_STRINGPARAMETER("suffix",&sSuffix)
    LONG_STRINGPARAMETER("index",&sIndexf)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("outvcf",&sOutVcf)

    LONG_PARAMETER_GROUP("Multithreading Options")
    LONG_INTPARAMETER("nthreads",&nThreads)

    LONG_PARAMETER_GROUP("Input file formats")
    LONG_INTPARAMETER("skipIndex",&indexSkip)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // create objects for logging
  if ( sOutVcf.IsEmpty() ) {
    error("ERROR: --outvcf is empty\n");
    //abort();
  }

  Logger::gLogger = new Logger((sOutVcf+".log").c_str(), bVerbose);
  
  time_t t;
  time(&t);
  Logger::gLogger->writeLog("Analysis started on %s", ctime(&t));
  
  ////////////////////////////////////////////////////////////
  // check the compatibility of arguments
  ///////////////////////////////////////////////////////////
  // Check the sanity of input file arguments
  ///////////////////////////////////////////////////////////
  if ( sAnchorVcf.IsEmpty() || sPrefix.IsEmpty() || sSuffix.IsEmpty() || sIndexf.IsEmpty() ) {
    Logger::gLogger->error("All the --anchor, --prefix, --suffix, --index option must be present");
  }

  // Read index file and list the file names to read
  IFILE indexFile = ifopen( sIndexf.c_str(), "rb" );
  String line;
  std::vector<std::string> inputVcfs;
  while( line.ReadLine(indexFile) > 0 ) {
    //fprintf(stderr,"line = %s",line.c_str());
    StringArray tok;
    tok.ReplaceColumns(line,'\t');
    if ( tok.Length() < indexSkip + 1 ) {
      Logger::gLogger->error("Cannot recognize %s in the index file",line.c_str());
    }
    for(int i=indexSkip; i < tok.Length(); ++i) {
      std::string s(sPrefix.c_str());
      StringArray paths;
      paths.ReplaceColumns(tok[i],'/');
      s += paths[paths.Length()-1].c_str();
      s += sSuffix.c_str();
      inputVcfs.push_back(s);
    }
  }

  int nVcfs = inputVcfs.size();
  Logger::gLogger->writeLog("Merging statistics from the following %d files..",nVcfs);
  for(int i=0; i < (int)inputVcfs.size(); ++i) {
    Logger::gLogger->writeLog("%s",inputVcfs[i].c_str());
  }

  // Read anchor VCF files
  FilterStat fStat;
  fStat.loadAnchorVcf(sAnchorVcf.c_str());

  // Add VCF file to add in multi-threaded manner
  for(int i=0; i < nVcfs; i += nThreads) {
    boost::thread_group Tg;
    for(int j=i; (j < nVcfs) && ( j < i + nThreads ); ++j) {
      Tg.create_thread( boost::bind( appendStatVcf, &fStat, inputVcfs[j].c_str() ) );
    }
    Tg.join_all();
  }

  // Re-read anchor VCF and print out merged VCF
  fStat.writeMergedVcf(sOutVcf.c_str());
}
