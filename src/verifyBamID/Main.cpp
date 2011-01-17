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
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <values.h>
#include <math.h>
#include <string>

#include "CigarRoller.h"
#include "Generic.h"
#include "GenomeSequence.h"
#include "SamFile.h"
#include "genmat/GenMatrix.h"
#include "Logger.h"
#include "base/string_tokenizer.h"
#include "BgzfFileType.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

#define MAX_Q 100 // maximum baseQuality phred score

Logger* Logger::gLogger = NULL; // Message log 

int main(int argc, char** argv) {
  // configurable parameters
  std::string sRef, sInFile, sInPrefix, sInSuffix, sIndexFile, sOutFile, sLogFile, sBfile, sBimFile, sMixValues, sHomValues, sIbdValues;
  bool bVerbose, bNoEOF, bSelfOnly, bBimAF, bMemoryMap, bUCSC, bPrecise;
  uint32_t minQ, maxQ, minMapQ, maxDepth;
  double genoError, minAF, mixUnit, homUnit, ibdUnit;

  // parameters currently not set to be configuraable
  bool sameSMFlag = true;
  uint16_t includeSamFlag = 0x0000;
  uint16_t excludeSamFlag = 0x0704;
  
  // configure parameter settings
  try {
    TCLAP::CmdLine cmd("Command description message", ' ', "0.0.2");
    TCLAP::ValueArg<std::string> argReference("r","reference","Karma's reference sequence",true,"","string");
    TCLAP::ValueArg<std::string> argIn("i","in","Input BAM file. Must be sorted and indexed",false,"","string");
    TCLAP::ValueArg<std::string> argInprefix("p","inprefix","Prefix of input BAM file for multi-chromosome inputs",false,"","string");
    TCLAP::ValueArg<std::string> argInsuffix("s","insuffix","Suffix of input BAM file for multi-chromosome inputs",false,"","string");
    TCLAP::ValueArg<std::string> argIndex("I","index","Index of input BAM file - [inputBam].bai will be default value",false,"","string");
    TCLAP::ValueArg<std::string> argOut("o","out","Prefix of output files",true,"","string");
    TCLAP::ValueArg<std::string> argLog("l","log","Log file - default: [out].log",false,"","string");
    TCLAP::ValueArg<std::string> argBfile("b","bfile","Binary PLINK format genotype file prefix. Must be forward-stranded",false,"","string");
    TCLAP::ValueArg<std::string> argBimFile("B","bimpfile","PLINK format BIM file with allele frequency information at the last column",false,"","string");
    TCLAP::ValueArg<uint32_t> argMinQ("q","minQ","Minimum Phred-scale base quality value (default:20) - bases with lower quality will be ignored",false,20,"integer");
    TCLAP::ValueArg<uint32_t> argMaxQ("Q","maxQ","Maximum Phred-scale base quality value (default:40) - higher base quality will be enforced to be this value",false,40,"integer");
    TCLAP::ValueArg<uint32_t> argMinMapQ("m","minMapQ","Minimum mapping quality value (default:10) - reads with lower quality will be ignored",false,10,"integer");
    TCLAP::ValueArg<uint32_t> argMaxDepth("d","maxDepth","Maximum depth per site (default:20) - bases with higher depth will be ignored due to possible alignment artifacts. For deep coverage data, it is important to set this value to a sufficiently large number (e.g. 200)",false,20,"integer");
    TCLAP::ValueArg<double> argGenoError("g","genoError","Error rate in genotype data (default: 0.005)",false,5.0e-3,"double");
    TCLAP::ValueArg<double> argMinAF("f","minAF","Minimum allele frequency (default: 0.005). Markers with lower allele frequencies will be ignored",false,5.0e-3,"double");
    TCLAP::ValueArg<double> argMixUnit("","mixUnit","unit of % mixture (default:0.01)",false,0.01,"double");
    TCLAP::ValueArg<double> argHomUnit("","homUnit","unit of % excessive homozygosity (default:0.01)",false,0.01,"double");
    TCLAP::ValueArg<double> argIbdUnit("","ibdUnit","unit of IBD values (default: 0.01)",false,0.01,"double");
    TCLAP::ValueArg<std::string> argMixValues("","mixValues","comma-separated list of candidate values for testing % mixture",false,"","string");
    TCLAP::ValueArg<std::string> argHomValues("","homValues","comma-separated list of candidate values for testing % excessive homozygosity",false,"","string");
    TCLAP::ValueArg<std::string> argIbdValues("","ibdValues","comma-separated list of candidate values for testing % IBD",false,"","string");
    TCLAP::SwitchArg switchVerbose("v","verbose","Toggle verbose mode (default:off)",cmd,false);
    TCLAP::SwitchArg switchNoeof("n","noeof","Do not check EOF marker for BAM file (default:off)",cmd,false);
    TCLAP::SwitchArg switchSelfOnly("S","selfonly","compare the genotypes with SELF (annotated sample) only (default:off)",cmd,false);
    TCLAP::SwitchArg switchBimAF("F","bimAF","use the allele frequency information by loading .bimp file instead of .bim file",cmd,false);
    TCLAP::SwitchArg switchMemoryMap("M","mmap","toggle whether to use memory map (default:true)",cmd,true);
    TCLAP::SwitchArg switchUCSC("u","ucsc","use UCSC-type reference file, where the sequence name starts with a prefix 'chr'",cmd,false);
    TCLAP::SwitchArg switchPrecise("","precise","improve the precision of calculation for high-depth (>100x) data",cmd,false);
    
    cmd.add(argReference);
    cmd.add(argIn);
    cmd.add(argInprefix);
    cmd.add(argInsuffix);
    cmd.add(argIndex);
    cmd.add(argOut);
    cmd.add(argLog);
    cmd.add(argBfile);
    cmd.add(argBimFile);
    cmd.add(argMinQ);
    cmd.add(argMaxQ);
    cmd.add(argMinMapQ);
    cmd.add(argMaxDepth);
    cmd.add(argGenoError);
    cmd.add(argMinAF);
    cmd.add(argMixUnit);
    cmd.add(argHomUnit);
    cmd.add(argIbdUnit);
    cmd.add(argMixValues);
    cmd.add(argHomValues);
    cmd.add(argIbdValues);

    cmd.parse(argc, argv);

    sRef = argReference.getValue();
    sInFile = argIn.getValue();
    sInPrefix = argInprefix.getValue();
    sInSuffix = argInsuffix.getValue();
    sIndexFile = argIndex.getValue();
    sOutFile = argOut.getValue();
    sLogFile = argLog.getValue();
    sBfile = argBfile.getValue();
    sBimFile = argBimFile.getValue();
    minQ = argMinQ.getValue();
    maxQ = argMaxQ.getValue();
    minMapQ = argMinMapQ.getValue();
    maxDepth = argMaxDepth.getValue();
    genoError = argGenoError.getValue();
    minAF = argMinAF.getValue();
    mixUnit = argMixUnit.getValue();
    homUnit = argHomUnit.getValue();
    ibdUnit = argIbdUnit.getValue();
    sMixValues = argMixValues.getValue();
    sHomValues = argHomValues.getValue();
    sIbdValues = argIbdValues.getValue();
    bVerbose = switchVerbose.getValue();
    bNoEOF = switchNoeof.getValue();
    bSelfOnly = switchSelfOnly.getValue();
    bBimAF = switchBimAF.getValue();
    bMemoryMap = switchMemoryMap.getValue();
    bUCSC = switchUCSC.getValue();
    bPrecise = switchPrecise.getValue();

    // create log file
    if ( sLogFile.empty() ) {
      sLogFile = sOutFile + ".log";
    }
    Logger::gLogger = new Logger(sLogFile.c_str(), bVerbose);

    // print out effective arguments
    //    Logger::gLogger->writeLog("Arguments in effect: \n%s",cmd.toString().c_str());
    Logger::gLogger->writeLog("Arguments in effect: \n");
    std::ostringstream oss;
    std::list<TCLAP::Arg*> argList = cmd.getArgList();
    for(std::list<TCLAP::Arg*>::iterator i=argList.begin(); 
        i != argList.end(); ++i)
    {
        oss << "\t" << (*i)->toString() << std::endl;
    }
    Logger::gLogger->writeLog("%s", oss.str().c_str());
  }
  // abort if failed paring the arguments
  catch (TCLAP::ArgException &e)  { 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    abort();
  }

  // post-process arguments to get configurable parameters
  std::vector<double> vfIBDs;
  std::vector<double> vfMixFracs;
  std::vector<double> vfHomFracs;
  uint32_t numIBDs;
  uint32_t numHomFracs;
  uint32_t numMixFracs;
  
  // build the list of IBD values to test
  if ( sIbdValues.empty() ) {
    numIBDs = static_cast<uint32_t>(ceil(1./ibdUnit))+1;
    for(uint32_t i=0; i < numIBDs; ++i) {
      double fIBD = ibdUnit * static_cast<double>(i);
      if ( fIBD > 1. ) {
	fIBD = 1.;
      }
      vfIBDs.push_back(fIBD);
    }
  }
  else {
    StringTokenizer t(sIbdValues, ", ");
    for(numIBDs = 0; t.GetNext(); ++numIBDs) {
      vfIBDs.push_back(atof(t.token().c_str()));
    }
    numIBDs = vfIBDs.size();
  }

  // build the list of %MIX values to test
  if ( sMixValues.empty() ) {
    numMixFracs = static_cast<uint32_t>(ceil(.5/mixUnit))+1;
    for(uint32_t i=0; i < numMixFracs; ++i) {
      double fMix = mixUnit * static_cast<double>(i);
      if ( fMix > .5 ) {
	fMix = .5;
      }
      vfMixFracs.push_back(fMix);
    }
  }
  else {
    StringTokenizer t(sMixValues, ", ");
    for(numMixFracs = 0; t.GetNext(); ++numMixFracs) {
      vfMixFracs.push_back(atof(t.token().c_str()));
    }
  }

  // build the list of %HOM values to test
  if ( sHomValues.empty() ) {
    numHomFracs = static_cast<uint32_t>(ceil(1./homUnit))+1;
    for(uint32_t i=0; i < numHomFracs; ++i) {
      double fHom = homUnit * static_cast<double>(i);
      if ( fHom > 1. ) {
	fHom = 1.;
      }
      vfHomFracs.push_back(fHom);
    }
  }
  else {
    StringTokenizer t(sHomValues, ", ");
    for(numHomFracs = 0; t.GetNext(); ++numHomFracs) {
      vfHomFracs.push_back(atof(t.token().c_str()));
    }
  }
  Logger::gLogger->writeLog("# IBDs = %d, # MIX = %d, # HOM = %d",numIBDs, numMixFracs,numHomFracs);

  double fPhred2Err[MAX_Q+1];

  // precompute the conversion from Phred score to baseError
  for(int32_t i=0; i < MAX_Q+1; ++i) {
    // If Phred>=maxQ, assume that the base quality is
    // overestimated and apply an upper threshold.
    if ( i > static_cast<int32_t>(maxQ) ) {
      fPhred2Err[i] = fPhred2Err[maxQ]; 
    }
    else {
      fPhred2Err[i] = pow(10.,(0-i)/10.);
    }
  }

  // check the validity of input files
  if ( sBfile.empty() && sBimFile.empty() ) {
    Logger::gLogger->error("Either --bfile or --bimfile is required");
  }

  if ( sInFile.empty() && ( sInPrefix.empty() || sInSuffix.empty() ) ) {
    Logger::gLogger->error("--in or (--inprefix & --insuffix) are required arguments");
  }

  if ( ( !sInFile.empty() ) && ( !sInPrefix.empty() ) && ( !sInSuffix.empty() ) ) {
    Logger::gLogger->error("--in and (--inprefix & --insuffix) cannot be combined together");
  }


  std::map<std::string,uint32_t> msRGidx; // map from RGID to an index
  std::vector<std::string> vsRGIDs; // list of indexed RGID
  std::vector<std::string> vsSMIDs; // list of SMID
  std::string sBamSMID; // per-BAM sample ID, when all readGroups have consistent sample ID

  GenomeSequence genomeSequence; // genomeSequence object
  genomeSequence.setReferenceName(sRef.c_str());
  genomeSequence.useMemoryMap(bMemoryMap);
  // open and create one if the index does not exist
  if ( genomeSequence.open() ) {
    Logger::gLogger->warning("Failed to open karma-type index of reference file %s. Creating a new one..",sRef.c_str());
    if ( genomeSequence.create(false)) {
      Logger::gLogger->error("Cannot crate karma-type index of reference file of %s",sRef.c_str());
    }
    if ( genomeSequence.open() ) {
      Logger::gLogger->error("Failed to open indexed karma-type reference file %s",sRef.c_str());
    }
  }

  // open genotype file
  BedMatrix* pBedMatrix = new BedMatrix;
  if ( !sBfile.empty() ) {
    if ( bBimAF ) { 
      // genotype is available, but use estimated allele frequency
      // this is useful when genotype data has limited # of individuals
      // (e.g. with --selfOnly option)
      pBedMatrix->open((sBfile+".fam").c_str(),(sBfile+".bimp").c_str(),(sBfile+".bed").c_str());
    }
    else {
      pBedMatrix->open((sBfile+".fam").c_str(),(sBfile+".bim").c_str(),(sBfile+".bed").c_str());
    }
  }
  else if ( !sBimFile.empty() ) { 
    // in the case only allele frequency is available
    bBimAF = true;
    pBedMatrix->openMarkerOnly(sBimFile.c_str());
  }
  else {
    Logger::gLogger->error("Neither --bfile nor --bimfile was specified");
  }

  // sanity check if the .bimp file loaded contains at least one positive allele frequency
  if ( bBimAF ) {
    uint32_t nZeros = 0;
    for(uint32_t i=0; i < pBedMatrix->vpGenMarkers.size(); ++i) {
      if ( pBedMatrix->vpGenMarkers[i]->dAlleleFrequency == 0 ) {
	++nZeros;
      }
    }
    if (nZeros == pBedMatrix->vpGenMarkers.size()) {
      Logger::gLogger->error("Allele frequencies are all zero. BIM file should contain additional column containing allele frequencies, and at least some of them should be nonzero");
    }
    else {
      Logger::gLogger->writeLog("Identified %u markers with zero minor allele frequencies and they are set to genotype error rate %lf",nZeros,genoError);
    }
  }

  uint32_t nMarkers = pBedMatrix->vpGenMarkers.size();
  uint32_t nInds = pBedMatrix->vpIndividuals.size();

  Logger::gLogger->writeLog("Total of %u markers over %u individuals exist in the genotype file",nMarkers,nInds);

  if ( bNoEOF ) {
    BgzfFileType::setRequireEofBlock(false);
  }

  std::string curChrom; // variable needed for tracking multi-chromosome data
  if ( !sInPrefix.empty() ) { // set up the first file to open
    if ( !sIndexFile.empty() ) {
      Logger::gLogger->error("Cannot use --index option with --prefix and --suffix option");
    }

    // open current BAM file
    curChrom = pBedMatrix->vpGenMarkers[0]->sChrom;
    sInFile = sInPrefix + pBedMatrix->vpGenMarkers[0]->getNCBIChromosomeName() + sInSuffix;
    sIndexFile = sInFile + ".bai"; // index name is fixed with multi-chromosome case
  }

  // routines for reading headers of the BAM file
  SamFile inBam;
  SamFileHeader inHeader;
  if ( ! (inBam.OpenForRead(sInFile.c_str()))  ) {
    Logger::gLogger->error("Cannot open BAM file %s for reading - %s",sInFile.c_str(), SamStatus::getStatusString(inBam.GetStatus()) );
  }
  inBam.ReadHeader(inHeader);

  // set index name
  if ( sIndexFile.empty() ) {
    sIndexFile = sInFile + ".bai";
  }

  // read index file
  if ( ! inBam.ReadBamIndex( sIndexFile.c_str() ) ) {
    Logger::gLogger->error("Cannot open BAM file index %s for reading %s", sIndexFile.c_str(),sInFile.c_str());
  }

  Logger::gLogger->writeLog("Reading header Records");

  // read through header records to get the readGroup ID and corresponding sample ID. 
  // It also checkes whether the sequenceName matches with the specified convention
  SamHeaderRecord* pSamHeaderRecord;
  bool bObservedChrinSQ = false;
  while( (pSamHeaderRecord = inHeader.getNextHeaderRecord()) != NULL ) {
    if ( pSamHeaderRecord->getType() == SamHeaderRecord::RG ) {
      std::string sRGID = pSamHeaderRecord->getTagValue("ID");
      std::string sSMID = pSamHeaderRecord->getTagValue("SM");
      if ( sRGID.empty() ) {
	Logger::gLogger->error("Readgroup ID is empty");
      }
      if ( sSMID.empty() ) {
	Logger::gLogger->warning("SM tag is missing in read group %s",sRGID.c_str());
      }
      vsRGIDs.push_back(sRGID);
      vsSMIDs.push_back(sSMID);
      if ( sBamSMID.empty() ) {
	sBamSMID = sSMID;
      }
      else if ( sBamSMID.compare(sSMID) != 0 ) {
	Logger::gLogger->warning("SM is not identical across the readGroups. Ignoring .bestSM/.selfSM outputs");
	sameSMFlag = false;
      }

      uint32_t idx = msRGidx.size();
      msRGidx[sRGID] = idx;
    }
    else if ( pSamHeaderRecord->getType() == SamHeaderRecord::SQ ) {
      std::string sSN = pSamHeaderRecord->getTagValue("SN");
      if ( strncmp(sSN.c_str(), "chr", 3) == 0 ) {
	bObservedChrinSQ = true;
      }
    }
  }

  if ( bObservedChrinSQ != bUCSC ) {
    Logger::gLogger->warning("--ucsc option %s specificed but %s 'chr' in the sequence names. If the results contain no data, try to toggle the option", bUCSC ? "is" : "is not", bObservedChrinSQ ? "observed" : "not observed");
  }

  uint32_t nCompInds = 0; // number of individuals to compare
  int32_t nSelfIndex = -1; // set when --selfonly option is set and the corresponding individual is found
  if ( bSelfOnly ) {
    if ( sameSMFlag ) {
      // when --selfonly is set, copy the matching individual to the very beginning of the data
      for(uint32_t i=0; i < nInds; ++i) {
	if ( sBamSMID.compare(pBedMatrix->vpIndividuals[i]->sIndID) == 0 ) {
	  if ( nSelfIndex < 0 ) {
	    nSelfIndex = static_cast<int32_t>(i);
	    pBedMatrix->vpIndividuals[0]->sIndID = pBedMatrix->vpIndividuals[i]->sIndID;
	  }
	  else {
	    Logger::gLogger->error("Multiple individuals with matching individual IDs (%s) in the PLINK file",sBamSMID.c_str());
	  }
	}
      }
      // disregard the genotypes when matching individual was not found with --selfonly option
      if ( nSelfIndex < 0 ) {
	Logger::gLogger->warning("--selfOnly was used but the individual in the BAM file cannot be found in the genotype file, skipping genotype comparison");
	nCompInds = 0;
      }
      else {
	nCompInds = 1;
      }
    }
    else {
      Logger::gLogger->error("The input BAM file consist of multiple individuals and cannot be run with --selfOnly option");
    }
  }
  else {
    nCompInds = nInds;
  }

  Logger::gLogger->writeLog("The following %u ReadGroup IDs are identified.",vsRGIDs.size());
  for(uint32_t i=0; i < vsRGIDs.size(); ++i) {
    Logger::gLogger->writeLog("\t%d: %s",i+1,vsRGIDs[i].c_str());
  }
  
  SamRecord samRecord;
  std::vector<uint32_t> nRGIndices;
  std::vector<char> cBases;
  std::vector<char> cQuals;
  CigarRoller cigarRoller;
  std::string cigar;
  uint32_t nRGs = vsRGIDs.size();

  // 1-dimensional (#ind) * (#geno) matrix of marker count (row-major)
  uint32_t* nIndGenoMarkers = new uint32_t[nCompInds * 3]();
  // (#ind) * (#geno) matrix of base count given marker genotypes
  uint32_t* nRGIndGenoMarkerBases = new uint32_t[nRGs * nCompInds * 3]();
  // (#ind) * (#geno) matrix of reference count given marker genotype
  uint32_t* nRGIndGenoMarkerRefs = new uint32_t[nRGs * nCompInds * 3]();
  // (#ind) * (#geno) matrix of non-reference count given marker genotype  
  uint32_t* nRGIndGenoMarkerAlts = new uint32_t[nRGs * nCompInds * 3]();

  // (#ind) * (#geno) matrix of base count given marker genotypes
  uint32_t* nSMIndGenoMarkerBases = new uint32_t[nCompInds * 3]();
  // (#ind) * (#geno) matrix of reference count given marker genotype
  uint32_t* nSMIndGenoMarkerRefs = new uint32_t[nCompInds * 3]();
  // (#ind) * (#geno) matrix of non-reference count given marker genotype  
  uint32_t* nSMIndGenoMarkerAlts = new uint32_t[nCompInds * 3]();

  // genotype-sequence IBD likelioods per RG, and per SM
  double* fSumRGIndLLKs = new double[nRGs * nCompInds * numIBDs]();
  double* fRGGenoIBDLLKs = new double[nRGs * 3]();
  //double* fRGnonIBDLLKs = new double[nRGs]();

  double* fSumSMIndLLKs = new double[nCompInds * numIBDs]();
  double* fSMGenoIBDLLKs = new double[3]();
  //double fSMnonIBDLLK = (bPrecise ? 0 : 1.);

  // %MIX likelihood per RG, and per SM
  double* fRGHomMixLLKs = new double[nRGs * (numHomFracs + numMixFracs - 1)]();
  double* fSMHomMixLLKs = new double[numHomFracs + numMixFracs - 1]();
  double* fRGMixMarkerLLKs = new double[nRGs * 9];
  double* fSMMixMarkerLLKs = new double[9];

  // Depth distribution of the markers per RG
  //uint32_t* nRGMarkerDepth = new uint32_t[nRGs * (maxDepth+1)]();
  //uint32_t* nRGMarkerDepthRefOnly = new uint32_t[nRGs * (maxDepth+1)]();
  //uint32_t* nRGMarkerDepthAltOnly = new uint32_t[nRGs * (maxDepth+1)]();
  //uint32_t* nRGMarkerDepthRefAlt = new uint32_t[nRGs * (maxDepth+1)]();
  //double* fRGSumEHetDepth = new double[nRGs * (maxDepth+1)]();

  double* fRGSumPriorHets = new double[nRGs]();
  double* fRGSumPosteriorHets = new double[nRGs]();
  uint32_t* nRGMarkers = new uint32_t[nRGs]();
  uint32_t* nRGBases = new uint32_t[nRGs]();
  uint32_t* nRGMarkerBases = new uint32_t[nRGs]();

  // count of ref, alt, other bases per readgroup
  //uint32_t* nRGRefCnts = new uint32_t[nRGs]();
  //uint32_t* nRGAltCnts = new uint32_t[nRGs]();
  //uint32_t* nRGEtcCnts = new uint32_t[nRGs]();

  // Depth distribution of the markers per SM
  //uint32_t* nSMMarkerDepth = new uint32_t[(maxDepth+1)]();
  //uint32_t* nSMMarkerDepthRefOnly = new uint32_t[(maxDepth+1)]();
  //uint32_t* nSMMarkerDepthAltOnly = new uint32_t[(maxDepth+1)]();
  //uint32_t* nSMMarkerDepthRefAlt = new uint32_t[(maxDepth+1)]();
  //double* fSMSumEHetDepth = new double[(maxDepth+1)]();

  double fSMSumPriorHets = 0;
  double fSMSumPosteriorHets = 0;
  uint32_t nSMMarkers = 0;
  uint32_t nSMBases = 0;

  // count of ref, alt, other bases per SM
  //uint32_t nSMRefCnt = 0;
  //uint32_t nSMAltCnt = 0;
  //uint32_t nSMEtcCnt = 0;

  Logger::gLogger->writeLog("Iterating over the genotypes..");
  //for( uint32_t markerCnt = 0; pBedMatrix->iterateMarker() && (markerCnt < 1000 ); ++markerCnt ) { // left for debugging purpose

  // iterate each variable sites (marker)
  uint32_t nMarkerCnt;
  for( nMarkerCnt = 0; pBedMatrix->iterateMarker(); ++nMarkerCnt ) {
    if ( (nMarkerCnt+1) % 1000 == 0 )
      Logger::gLogger->writeLog("Processing %u / %u markers",nMarkerCnt+1,nMarkers);
    
    // In 1000G data, one BAM file is splitted into each chromosome, 
    // this is the routine to process the special cases
    // [inPrefix][ncbiChrName][outPrefix]
    if ( pBedMatrix->currentMarker()->sChrom.compare(curChrom) != 0 ) {
      curChrom = pBedMatrix->currentMarker()->sChrom;
      Logger::gLogger->writeLog("Opening chromosome %s",pBedMatrix->currentMarker()->getNCBIChromosomeName().c_str());
      
      if ( !sInPrefix.empty() )  {
	// if the current marker does not match to current BAM file,
	// close the current file, open new file with header
	inBam.Close();
	
	sInFile = sInPrefix + pBedMatrix->currentMarker()->getNCBIChromosomeName() + sInSuffix;
	if ( ! (inBam.OpenForRead(sInFile.c_str())) ) {
	  Logger::gLogger->error("Cannot open BAM file %s for reading - %s",sInFile.c_str(), SamStatus::getStatusString(inBam.GetStatus()) );	
	}
	inBam.ReadHeader(inHeader);
	inHeader.resetHeaderRecordIter();
	sIndexFile = sInFile + ".bai";

	if ( ! inBam.ReadBamIndex( sIndexFile.c_str() ) ) {
	  Logger::gLogger->error("Cannot open BAM file index %s for reading %s", sIndexFile.c_str(),sInFile.c_str());
	}
	
	for( uint32_t k=0, rgCnt = 0; (pSamHeaderRecord = inHeader.getNextHeaderRecord()) != NULL; ++k ) {
	  if ( pSamHeaderRecord->getType() == SamHeaderRecord::RG ) {
	    std::string sRGID = pSamHeaderRecord->getTagValue("ID");
	    std::string sSMID = pSamHeaderRecord->getTagValue("SM");
	    
	    if ( sRGID.empty() ) {
	      Logger::gLogger->error("Readgroup ID is empty");
	    }
	    if ( sSMID.empty() ) {
	      Logger::gLogger->warning("SM tag is missing in read group %s",sRGID.c_str());
	    }
	    
	    if ( vsRGIDs[rgCnt].compare(sRGID) != 0 ) {
	      Logger::gLogger->error("Nonidentical headers - ReadGroup ID is different at %u : %s vs %s",rgCnt,vsRGIDs[rgCnt].c_str(),sRGID.c_str());
	    }
	    if ( vsSMIDs[rgCnt].compare(sSMID) != 0 ) {
	      Logger::gLogger->error("Nonidentical headers - Readgroup SM is different at %u : %s vs %s",rgCnt,vsSMIDs[rgCnt].c_str(),sSMID.c_str());
	    }
	    ++rgCnt;
	  }
	}
      }
    }


    // Routine for processing genotype data at each marker

    // Read the genotype file
    char* bedGenos = pBedMatrix->currentGenotypes();
    char a1 = pBedMatrix->currentMarker()->sAllele1[0];
    char a2 = pBedMatrix->currentMarker()->sAllele2[0];
    genomeIndex_t markerIndex = genomeSequence.getGenomePosition( bUCSC ? pBedMatrix->currentMarker()->getUCSCChromosomeName().c_str() : pBedMatrix->currentMarker()->getNCBIChromosomeName().c_str(), pBedMatrix->currentMarker()->nBasePosition );
    char refBase = genomeSequence[markerIndex];
    bool flipFlag = false;

    // routines for a1/a2 flip
    // note that strand flip would not be automatically handled
    if ( refBase == a1 ) {
      // fine
    }
    else if ( refBase == a2 ) {
      // flip genotypes if reference base does not match
      char tmp = a2;
      a2 = a1;
      a1 = tmp;
      for(uint32_t i=0; i < nInds; ++i) {
	if ( bedGenos[i] == 1 ) 
	  bedGenos[i] = 3;
	else if ( bedGenos[i] == 3 ) 
	  bedGenos[i] = 1;
      }
      flipFlag = true;
    }
    else {
      Logger::gLogger->error("Marker alleles %s at %s:%u do not match to reference allele",pBedMatrix->currentMarker()->sMarkerID.c_str(), pBedMatrix->currentMarker()->getNCBIChromosomeName().c_str(), pBedMatrix->currentMarker()->nBasePosition);
    }

    // calculate non-reference allele frequency
    double alleleFreq = 0.0;
    uint32_t nMissing = 0;
    if ( !bBimAF ) {
      for(uint32_t i=0; i < nInds; ++i) {
	if ( bedGenos[i] == 0 ) {
	  ++nMissing;
	  alleleFreq *= static_cast<double>(nInds-nMissing+1)/static_cast<double>(nInds-nMissing);
	}
	else {
	  alleleFreq += static_cast<double>(bedGenos[i]-1)/static_cast<double>(2*(nInds-nMissing));
	}
      }
    }
    else {
      alleleFreq = pBedMatrix->currentMarker()->dAlleleFrequency;
      if ( !flipFlag ) {
	alleleFreq = 1. - alleleFreq;
      }
    }

    if ( ( alleleFreq < minAF ) || ( alleleFreq > 1.-minAF) ) {
      continue;
    }

    if ( alleleFreq < genoError ) alleleFreq = genoError;
    else if ( alleleFreq > 1.-genoError ) alleleFreq = 1.-genoError;

    //Logger::gLogger->writeLog("%c %c %c %lf",refBase,a1,a2,alleleFreq);

    // genotype-match likelihood
    for(uint32_t i=0; i < nRGs*3; ++i) {
      fRGGenoIBDLLKs[i] = (bPrecise ? 0. : 1.);
    }
    //for(uint32_t i=0; i < nRGs; ++i) {
    //  fRGnonIBDLLKs[i] = (bPrecise ? 0. : 1.);
    //}

    fSMGenoIBDLLKs[0] = fSMGenoIBDLLKs[1] = fSMGenoIBDLLKs[2] = (bPrecise ? 0. : 1.);


    // retrieve the bases aligned to the marker position
    //std::cerr << "foo" << std::endl;
    int refID = inHeader.getReferenceID( bUCSC ? pBedMatrix->currentMarker()->getUCSCChromosomeName().c_str() : pBedMatrix->currentMarker()->getNCBIChromosomeName().c_str() );
    //std::cerr << "bar" << std::endl;
    uint32_t pos = pBedMatrix->currentMarker()->nBasePosition; // 1-based position
    int numSectionRecords = 0;
    
    nRGIndices.clear();
    cBases.clear();
    cQuals.clear();

    // **** Rouitine for reading each read
    if ( refID >= 0 ) {
      inBam.SetReadSection( refID, pos-1, pos ); // set bam file to retrieve the reads overlapping with the particular genomic position chr(refID):bp(pos) 

      // Keep reading records until they aren't anymore.
      while(inBam.ReadRecord(inHeader, samRecord)) {
	++numSectionRecords;

	// filtering step - mapQ
	if ( samRecord.getMapQuality() < minMapQ ) 
	  continue;

	// skip flagged reads
	uint16_t samFlags = samRecord.getFlag();
	if ( includeSamFlag && ( ( samFlags & includeSamFlag ) != includeSamFlag ) )
	  continue;
	if ( excludeSamFlag && ( samFlags & excludeSamFlag ) )
	  continue;

	// obtain readGroup info and store to rgIdx
	char tag[3];
	char vtype;
	void* value;
	bool found = false;
	uint32_t rgIdx;
	while( samRecord.getNextSamTag(tag, vtype, &value) != false ) {
	  if ( strcmp(tag, "RG") == 0 ) {
	    found = true;
	    if ( vtype == 'Z' ) {
	      std::string sValue = ((String)*(String*)value).c_str();
	      if ( msRGidx.find(sValue) != msRGidx.end() ) {
		rgIdx = msRGidx[sValue];
	      }
	      else {
		Logger::gLogger->error("ReadGroup ID %s cannot be found",sValue.c_str());
	      }
	    }
	    else {
	      Logger::gLogger->error("vtype of RG tag must be 'Z'");
	    }
	  }
	}
	if ( found == false ) {
	  Logger::gLogger->error("Cannot find RG tag for readName %s",samRecord.getReadName());
	}

	// access the base calls and qualities
	uint32_t readStartPosition = samRecord.get1BasedPosition();
	int32_t offset = pos - readStartPosition;
	const char* readQuality = samRecord.getQuality();
	const char* readSequence = samRecord.getSequence();

	cigar = samRecord.getCigar();
	cigarRoller.Set(cigar.c_str());

	if ( offset >= 0 ) {
	  int32_t readIndex = cigarRoller.getQueryIndex(offset);
	  if ( readIndex != CigarRoller::INDEX_NA ) {
	    if ( ( static_cast<uint32_t>(readQuality[readIndex]) >= minQ + 33 ) && ( readSequence[readIndex] != 'N' ) ) {
	      nRGIndices.push_back(rgIdx);
	      cBases.push_back(readSequence[readIndex]);
	      cQuals.push_back(readQuality[readIndex]);
	    }
	  }
	}
      }

      // at the end, we will have
      // nRGIndices - vector of readGroup indices
      // cBases - vector of bases
      // cQual - vector of base qualities

      if ( bSelfOnly ) {  
	if ( nSelfIndex >= 0 ) {
	  bedGenos[0] = bedGenos[nSelfIndex];
	}
      }

      for(uint32_t i=0; i < nCompInds; ++i) {
	if ( bedGenos[i] > 0 ) {
	  ++nIndGenoMarkers[3*i + bedGenos[i]-1];
	}
      }

      /*
      std::string sIndices;
      for(uint32_t j=0; j < nRGIndices.size(); ++j) {
	char buf[255];
	if ( j > 0 ) sIndices += " ";
	sprintf(buf,"%u:%c:%d",nRGIndices[j],cBases[j],static_cast<int>(cQuals[j])-33);
	sIndices += buf;
      }

      if ( !sIndices.empty() ) {
	Logger::gLogger->writeLog("Reading %s:%u (%c,%c)... %d records are found - %s", pBedMatrix->currentMarker()->sChrom.c_str(), pos, a1, a2, numSectionRecords, sIndices.c_str());
      }*/

      memset(nRGMarkerBases, 0, sizeof(uint32_t)*nRGs);
      if ( nRGIndices.size() <= maxDepth ) {
	//double probReadIBD[4], probReadnonIBD;
	double probReadIBD[3]; // missing genotype is not allowed
	//for(uint32_t i=0; i < nRGs; ++i) {
	//  nRGRefCnts[i] = 0;
	//  nRGAltCnts[i] = 0;
	//  nRGEtcCnts[i] = 0;
	//}
	//nSMRefCnt = nSMAltCnt = nSMEtcCnt = 0;

	for(uint32_t i=0; i < nRGIndices.size(); ++i) {
	  //double baseError = pow(10.,(33.-static_cast<double>(cQuals[i]))/10.);
	  //if ( baseError < 1e-4) { baseError = 1e-4; }
	  double baseError = fPhred2Err[cQuals[i]-33];

	  if ( cBases[i] == a1 ) {
	    //probReadIBD[0] = probReadnonIBD = (1.-alleleFreq)*(1.-baseError) + alleleFreq * (baseError/3.);
	    probReadIBD[0] = 1.-baseError; 
	    probReadIBD[1] = 0.5-baseError/3.;
	    probReadIBD[2] = baseError/3.;
	    //++nRGRefCnts[nRGIndices[i]];
	    //++nSMRefCnt;
	  }
	  else if ( cBases[i] == a2 ) {
	    //probReadIBD[0] = probReadnonIBD = (alleleFreq)*(1.-baseError) + (1.-alleleFreq) * (baseError/3.);
	    probReadIBD[0] = baseError/3.;
	    probReadIBD[1] = 0.5-baseError/3.;
	    probReadIBD[2] = 1.-baseError;
	    //++nRGAltCnts[nRGIndices[i]];
	    //++nSMAltCnt;
	  }
	  else {
	    //probReadIBD[0] = probReadnonIBD
	    probReadIBD[0] = probReadIBD[1] = probReadIBD[2] = baseError/3.;
	    //++nRGEtcCnts[nRGIndices[i]];
	    //++nSMEtcCnt;
	  }
	  ++nRGMarkerBases[nRGIndices[i]];
	  ++nSMBases;

	  for(uint32_t j=0; j < 3; ++j) {
	    if ( bPrecise ) {
	      fRGGenoIBDLLKs[nRGIndices[i]*3+j] += log(probReadIBD[j]);
	      fSMGenoIBDLLKs[j] += log(probReadIBD[j]);
	    }
	    else {
	      fRGGenoIBDLLKs[nRGIndices[i]*3+j] *= probReadIBD[j];
	      fSMGenoIBDLLKs[j] *= probReadIBD[j];
	    }
	  }
	  //if ( bPrecise ) {
	  //  fRGnonIBDLLKs[nRGIndices[i]] += log(probReadnonIBD);
	  //  fSMnonIBDLLK += log(probReadnonIBD);
	  //}
	  //else {
	  //  fRGnonIBDLLKs[nRGIndices[i]] *= probReadnonIBD;
	  //  fSMnonIBDLLK *= probReadnonIBD;
	  //}

	  for(uint32_t j=0; j < nCompInds; ++j) {
	    if ( bedGenos[j] > 0 ) {
	      ++nRGIndGenoMarkerBases[nRGIndices[i]*nCompInds*3 + j*3 + bedGenos[j]-1];
	      ++nSMIndGenoMarkerBases[j*3 + bedGenos[j]-1];
	      if ( cBases[i] == a1 ) {
		++nRGIndGenoMarkerRefs[nRGIndices[i]*nCompInds*3 + j*3 + bedGenos[j]-1];
		++nSMIndGenoMarkerRefs[j*3 + bedGenos[j]-1];
	      }
	      else if ( cBases[i] == a2 ) {
		++nRGIndGenoMarkerAlts[nRGIndices[i]*nCompInds*3 + j*3 + bedGenos[j]-1];
		++nSMIndGenoMarkerAlts[j*3 + bedGenos[j]-1];
	      }
	    }
	  }
	}

	for(uint32_t i=0; i < nRGs; ++i) {
	  if ( nRGMarkerBases[i] > 0 ) {
	    ++(nRGMarkers[i]);
	    nRGBases[i] += nRGMarkerBases[i];
	    fRGSumPriorHets[i] += (2.*alleleFreq*(1.-alleleFreq));

	    if ( bPrecise ) {
	      // based on log-likeligood
	      fRGSumPosteriorHets[i] +=
		( (2.*alleleFreq*(1.-alleleFreq))/
		( (1.-alleleFreq)*(1.-alleleFreq)*exp(fRGGenoIBDLLKs[i*3 + 0] - fRGGenoIBDLLKs[i*3 + 1]) +
		  (2.*alleleFreq*(1.-alleleFreq)) +
		  alleleFreq*alleleFreq*exp(fRGGenoIBDLLKs[i*3 + 2] - fRGGenoIBDLLKs[i*3 + 1]) ) );
	    }
	    else {
	      fRGSumPosteriorHets[i] +=
		(2.*alleleFreq*(1.-alleleFreq))*fRGGenoIBDLLKs[i*3 + 1]/
		( (1.-alleleFreq)*(1.-alleleFreq)*fRGGenoIBDLLKs[i*3 + 0] +
		  (2.*alleleFreq*(1.-alleleFreq))*fRGGenoIBDLLKs[i*3 + 1] +
		  alleleFreq*alleleFreq*fRGGenoIBDLLKs[i*3 + 2] );
	    }
	  }
	}
	++nSMMarkers;
	fSMSumPriorHets += (2.*alleleFreq*(1.-alleleFreq));
	if ( bPrecise ) {
	  // based on log-likeligood
	  fSMSumPosteriorHets +=
	    ( (2.*alleleFreq*(1.-alleleFreq))/
	      ( (1.-alleleFreq)*(1.-alleleFreq)*exp(fSMGenoIBDLLKs[0] - fSMGenoIBDLLKs[1]) +
		(2.*alleleFreq*(1.-alleleFreq)) +
		alleleFreq*alleleFreq*exp(fSMGenoIBDLLKs[2] - fSMGenoIBDLLKs[1]) ) );
	}
	else {
	  fSMSumPosteriorHets +=
	    (2.*alleleFreq*(1.-alleleFreq))*fSMGenoIBDLLKs[1]/
	    ( (1.-alleleFreq)*(1.-alleleFreq)*fSMGenoIBDLLKs[0] +
	      (2.*alleleFreq*(1.-alleleFreq))*fSMGenoIBDLLKs[1] +
	      alleleFreq*alleleFreq*fSMGenoIBDLLKs[2] );
	}

	// assume uniform prior on each genotype and calculate the
	// expected number of heterozygous alleles
	// e.g. when DP=2
	// Pr(AB|Reads) = Pr(Reads|AB)*Pr(AB)/Pr(Reads)
	//  = Pr(Reads|AB)*Pr(AB)/(Pr(Reads|AB)Pr(AB)+Pr(Reads|AA)Pr(AA)+Pr(Reads|BB)+Pr(BB)
	//  = Pr(Reads|AB)/(Pr(Reads|AA)+Pr(Reads|AB)+Pr(Reads|BB))
	//  if DP=1, the value will be close to 0.5 for every marker

	//for(uint32_t j=0; j < nRGs; ++j) {
	  //uint32_t k = j*(maxDepth+1)+nRGRefCnts[j]+nRGAltCnts[j];
	  //if ( nRGRefCnts[j] > 0 ) {
	  //  if ( nRGAltCnts[j] > 0 ) {
	  //    ++nRGMarkerDepthRefAlt[k];
	  //  }
	  //  else {
	  //    ++nRGMarkerDepthRefOnly[k];
	  //  }
	  //  ++nRGMarkerDepth[k];
	  //  fRGSumEHetDepth[k] += (2.*alleleFreq*(1.-alleleFreq));
	  //}
	  //else if ( nRGAltCnts[j] > 0 ) {
	  //  ++nRGMarkerDepth[k];
	  //  ++nRGMarkerDepthAltOnly[k];
	  //  fRGSumEHetDepth[k] += (2.*alleleFreq*(1.-alleleFreq));
	  //}
	  //else { // skip
	  //}

	  //nSMRefCnt += nRGRefCnts[j];
	  //nSMAltCnt += nRGAltCnts[j];
	  //nSMEtcCnt += nRGEtcCnts[j];
	//}

	//if ( nSMRefCnt > 0 ) {
	//  if ( nSMAltCnt > 0 ) {
	//    ++nSMMarkerDepthRefAlt[nSMRefCnt+nSMAltCnt];
	//  }
	//  else {
	//    ++nSMMarkerDepthRefOnly[nSMRefCnt+nSMAltCnt];
	//  }
	//  ++nSMMarkerDepth[nSMRefCnt+nSMAltCnt];
	// fSMSumEHetDepth[nSMRefCnt+nSMAltCnt] += (2.*alleleFreq*(1.-alleleFreq));
	//}
	//else if ( nSMAltCnt > 0 ) {
	//  ++nSMMarkerDepth[nSMRefCnt+nSMAltCnt];
	//  ++nSMMarkerDepthAltOnly[nSMRefCnt+nSMAltCnt];
	//  fSMSumEHetDepth[nSMRefCnt+nSMAltCnt] += (2.*alleleFreq*(1.-alleleFreq));
	//}

	// calculating per-RG LLK
	for(uint32_t i=0; i < nRGs; ++i) {
	  double fRGNonIBDLLK; // likelihood of reads given non-IBDs 
	  if ( bPrecise ) {
	    uint32_t maxIdx = 0;
	    for(uint32_t j=1; j < 3; ++j) {
	      if ( fRGGenoIBDLLKs[i*3+j] > fRGGenoIBDLLKs[i*3+maxIdx] ) {
		maxIdx = j;
	      }
	    }
	    fRGNonIBDLLK = fRGGenoIBDLLKs[i*3+maxIdx] + log((1.-alleleFreq)*(1.-alleleFreq)*exp(fRGGenoIBDLLKs[i*3+0]-fRGGenoIBDLLKs[i*3+maxIdx]) + 2.*alleleFreq*(1.-alleleFreq)*exp(fRGGenoIBDLLKs[i*3+1]-fRGGenoIBDLLKs[i*3+maxIdx]) + alleleFreq*alleleFreq*exp(fRGGenoIBDLLKs[i*3+2]-fRGGenoIBDLLKs[i*3+maxIdx]));
	  }
	  else {
	    fRGNonIBDLLK = (1.-alleleFreq)*(1.-alleleFreq)*fRGGenoIBDLLKs[i*3+0] + 2.*alleleFreq*(1.-alleleFreq)*fRGGenoIBDLLKs[i*3+1] + alleleFreq*alleleFreq*fRGGenoIBDLLKs[i*3+2];
	  }

	  for(uint32_t j=0; j < nCompInds; ++j) {
	    if ( bPrecise ) {
	      double probReadsIBD = (bedGenos[j] > 0) ? ( (1.-1.5*genoError) * exp(fRGGenoIBDLLKs[ i*3 + bedGenos[j] - 1] - fRGNonIBDLLK)  + 0.5*genoError*(exp(fRGGenoIBDLLKs[i*3+0] - fRGNonIBDLLK) + exp(fRGGenoIBDLLKs[i*3+1] - fRGNonIBDLLK) + exp(fRGGenoIBDLLKs[i*3+2] - fRGNonIBDLLK) ) ) : 1.;
	      for(uint32_t k=0; k < numIBDs; ++k) {
		double llk = fRGNonIBDLLK + log( vfIBDs[k]*probReadsIBD + (1.-vfIBDs[k])*1. );
		
		fSumRGIndLLKs[i*nCompInds*numIBDs+j*numIBDs+k] += llk;
	      }
	    }
	    else {
	      double probReadsIBD = (bedGenos[j] > 0) ? ( (1.-1.5*genoError) * ( fRGGenoIBDLLKs[ i*3 + bedGenos[j] - 1 ] )  + 0.5*genoError*(fRGGenoIBDLLKs[i*3+0] + fRGGenoIBDLLKs[i*3+1] + fRGGenoIBDLLKs[i*3+2]) ) : fRGNonIBDLLK;
	      for(uint32_t k=0; k < numIBDs; ++k) {
		double llk = log( vfIBDs[k]*probReadsIBD + (1.-vfIBDs[k])*fRGNonIBDLLK );
		fSumRGIndLLKs[i*nCompInds*numIBDs+j*numIBDs+k] += llk;
	      }
	    }
	  }
	}

	// calcuating per-SM LLK
	double fSMNonIBDLLK;
	if ( bPrecise ) {
	  uint32_t maxIdx = 0;
	  for(uint32_t j=1; j < 3; ++j) {
	    if ( fSMGenoIBDLLKs[j] > fSMGenoIBDLLKs[maxIdx] ) {
	      maxIdx = j;
	    }
	  }
	  fSMNonIBDLLK = fSMGenoIBDLLKs[maxIdx] + log((1.-alleleFreq)*(1.-alleleFreq)*exp(fSMGenoIBDLLKs[0]-fSMGenoIBDLLKs[maxIdx]) + 2.*alleleFreq*(1.-alleleFreq)*exp(fSMGenoIBDLLKs[1]-fSMGenoIBDLLKs[maxIdx]) + alleleFreq*alleleFreq*exp(fSMGenoIBDLLKs[2]-fSMGenoIBDLLKs[maxIdx]));
	}
	else {
	  fSMNonIBDLLK = (1.-alleleFreq)*(1.-alleleFreq)*fSMGenoIBDLLKs[0] + 2.*alleleFreq*(1.-alleleFreq)*fSMGenoIBDLLKs[1] + alleleFreq*alleleFreq*fSMGenoIBDLLKs[2];
	}

	for(uint32_t j=0; j < nCompInds; ++j) {
	  if ( bPrecise ) {
	    double probReadsIBD = (bedGenos[j] > 0) ? ( (1.-1.5*genoError) * exp(fSMGenoIBDLLKs[ bedGenos[j] - 1] - fSMNonIBDLLK)  + 0.5*genoError*(exp(fSMGenoIBDLLKs[0] - fSMNonIBDLLK) + exp(fSMGenoIBDLLKs[1] - fSMNonIBDLLK) + exp(fSMGenoIBDLLKs[2] - fSMNonIBDLLK) ) ) : 1.;
	    for(uint32_t k=0; k < numIBDs; ++k) {
	      double llk = fSMNonIBDLLK + log( vfIBDs[k]*probReadsIBD + (1.-vfIBDs[k])*1. );
	      fSumSMIndLLKs[j*numIBDs+k] += llk;
	    }
	  }
	  else {
	    double probReadsIBD = (bedGenos[j] > 0) ? ( (1.-1.5*genoError) * ( fSMGenoIBDLLKs[ bedGenos[j] - 1 ] )  + 0.5*genoError*(fSMGenoIBDLLKs[0] + fSMGenoIBDLLKs[1] + fSMGenoIBDLLKs[2]) ) : fSMNonIBDLLK;
	    for(uint32_t k=0; k < numIBDs; ++k) {
	      double llk = log( vfIBDs[k]*probReadsIBD + (1.-vfIBDs[k])*fSMNonIBDLLK );
	      fSumSMIndLLKs[j*numIBDs+k] += llk;
	    }
	  }
	}

	// Accounting for library complexity
	// alpha - % MIX - excessive heterozygosity
	// beta - % HOM - excessive homozygosity
	// R - observed reads
	// Z - complexity bottlneck status ( beta % of bottlneck ) - only one base sampled from AF
	// G - True underlying genotypes of the individuals - 2 * 2 alleles

	// Goo's routine
	// Given : 
	//   - alleleFrequency
	//   - a1 : reference allele
	//   - a2 : alternative allele
	//   - nRGIndices
	//   - cBases
	//   - cQuals

	//  (HOM,MIX) = (numHomFracs-1,0),...(-1,0),(0,0),(0,1),..,(0,numMixFracs-1)
	for(int32_t l = 0-static_cast<int32_t>(numHomFracs)+1; l < static_cast<int32_t>(numMixFracs); ++l) {
	  uint32_t h = (l < 0) ? static_cast<uint32_t>(0-l) : 0;
	  uint32_t j = (l > 0) ? static_cast<uint32_t>(l) : 0;
	  uint32_t lu = static_cast<uint32_t>( l + numHomFracs -1 );
	  
	  double alpha = vfMixFracs[j];
	  double beta = vfHomFracs[h];

	  double genoFreq[3] = {(1.-beta)*(1.-alleleFreq)*(1.-alleleFreq) + beta*(1.-alleleFreq), (1.-beta)*2.*alleleFreq*(1.-alleleFreq), (1.-beta)*alleleFreq*alleleFreq + beta*alleleFreq};

	  for(uint32_t i=0; i < 9; ++i) {
	    fSMMixMarkerLLKs[i] = (bPrecise ? 0. : 1.);
	  }
	  for(uint32_t i=0; i < 9*nRGs; ++i) {
	    fRGMixMarkerLLKs[i] = (bPrecise ? 0. : 1.);
	  }
	  //double valuesToMultiply[9] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
	  double lks[9];
	  for(uint32_t i=0; i < nRGIndices.size(); ++i) {
	    double baseError = fPhred2Err[cQuals[i]-33];
	    double baseMatch = 1.-baseError;
	    if ( cBases[i] == a1 ) {
	      lks[0*3+0] = baseMatch;
	      lks[0*3+1] = (0.5+0.5*alpha)*baseMatch + ((1.-alpha)/6.)*baseError;
	      lks[0*3+2] = alpha*baseMatch + ((1.-alpha)/3.)*baseError;
	      lks[1*3+0] = (1.-0.5*alpha)*baseMatch + (alpha/6.)*baseError;
	      lks[1*3+1] = 0.5*baseMatch + 1./6.*baseError;
	      lks[1*3+2] = (0.5*alpha)*baseMatch + (1./3.-alpha/6.)*baseError;
	      lks[2*3+0] = (1.-alpha)*baseMatch + (alpha/3.)*baseError;
	      lks[2*3+1] = (0.5*(1.-alpha))*baseMatch + (1./6.+alpha/6.)*baseError;
	      lks[2*3+2] = baseError/3.;
	    }
	    else if ( cBases[i] == a2 ) {
	      lks[2*3+2] = baseMatch;
	      lks[2*3+1] = (0.5+0.5*alpha)*baseMatch + ((1.-alpha)/6.)*baseError;
	      lks[2*3+0] = alpha*baseMatch + ((1.-alpha)/3.)*baseError;
	      lks[1*3+2] = (1.-0.5*alpha)*baseMatch + (alpha/6.)*baseError;
	      lks[1*3+1] = 0.5*baseMatch + 1./6.*baseError;
	      lks[1*3+0] = (0.5*alpha)*baseMatch + (1./3.-alpha/6.)*baseError;
	      lks[0*3+2] = (1.-alpha)*baseMatch + (alpha/3.)*baseError;
	      lks[0*3+1] = (0.5*(1.-alpha))*baseMatch + (1./6.+alpha/6.)*baseError;
	      lks[0*3+0] = baseError/3.;
	    }
	    else {
	      lks[0*3+0] = baseError/3.;
	      lks[0*3+1] = baseError/3.;
	      lks[0*3+2] = baseError/3.;
	      lks[1*3+0] = baseError/3.;
	      lks[1*3+1] = baseError/3.;
	      lks[1*3+2] = baseError/3.;
	      lks[2*3+0] = baseError/3.;
	      lks[2*3+1] = baseError/3.;
	      lks[2*3+2] = baseError/3.;
	    }
	    
	    for(uint32_t k=0; k < 9; ++k) {
	      if ( bPrecise ) {
		//valuesToMultiply[k] *= lks[k];
		fSMMixMarkerLLKs[k] += log(lks[k]);
		fRGMixMarkerLLKs[nRGIndices[i] * 9 + k] += log(lks[k]);
	      }
	      else {
		fSMMixMarkerLLKs[k] *= lks[k];
		fRGMixMarkerLLKs[nRGIndices[i] * 9 + k] *= lks[k];
	      }
	    }
	  }
	  
	  double perMarkerProb = 0;
	  if ( bPrecise ) {
	    uint32_t maxIdx = 0;
	    for(uint32_t k = 1; k < 9; ++k) {
	      if ( fSMMixMarkerLLKs[maxIdx] < fSMMixMarkerLLKs[k] ) {
		maxIdx = k;
	      }
	    }
	    for(uint32_t k1 = 0; k1 < 3; ++k1) {
	      for(uint32_t k2 = 0; k2 < 3; ++k2) {
		perMarkerProb += (exp(fSMMixMarkerLLKs[k1*3+k2]-fSMMixMarkerLLKs[maxIdx])*genoFreq[k1]*genoFreq[k2]);
	      }
	    }
	    fSMHomMixLLKs[lu] += (fSMMixMarkerLLKs[maxIdx]+log(perMarkerProb));
	  }
	  else {
	    for(uint32_t k1 = 0; k1 < 3; ++k1) {
	      for(uint32_t k2 = 0; k2 < 3; ++k2) {
		perMarkerProb += (fSMMixMarkerLLKs[k1*3+k2]*genoFreq[k1]*genoFreq[k2]);
	      }
	    }
	    fSMHomMixLLKs[lu] += log(perMarkerProb);
	  }
	  
	  for(uint32_t k=0; k < nRGs; ++k) {
	    perMarkerProb = 0;
	    if ( bPrecise ) {
	      uint32_t maxIdx = 0;
	      for(uint32_t k1 = 1; k1 < 9; ++k1) {
		if ( fRGMixMarkerLLKs[k*9+maxIdx] < fRGMixMarkerLLKs[k*9+k1] ) {
		  maxIdx = k1;
		}
	      }
	      
	      for(uint32_t k1 = 0; k1 < 3; ++k1) {
		for(uint32_t k2 = 0; k2 < 3; ++k2) {
		  //perMarkerProb += (fRGMixMarkerLKs[k*9+k1*3+k2]*genoFreq[k1]*genoFreq[k2]);
		  perMarkerProb += (exp(fRGMixMarkerLLKs[k*9+k1*3+k2]-fRGMixMarkerLLKs[k*9+maxIdx])*genoFreq[k1]*genoFreq[k2]);
		}
	      }
	      //fRGMixLLKs[k*numMixFracs + j] += log(perMarkerProb);
	      fRGHomMixLLKs[k*(numHomFracs+numMixFracs-1) + lu] += (fRGMixMarkerLLKs[k*9+maxIdx]+log(perMarkerProb));
	    }
	    else {
	      for(uint32_t k1 = 0; k1 < 3; ++k1) {
		for(uint32_t k2 = 0; k2 < 3; ++k2) {
		  perMarkerProb += (fRGMixMarkerLLKs[k*9+k1*3+k2]*genoFreq[k1]*genoFreq[k2]);
		}
	      }
	      fRGHomMixLLKs[k*(numMixFracs+numHomFracs-1) + lu] += log(perMarkerProb);
	    }
	  }
	}
      }
      else {
      Logger::gLogger->writeLog("Skipping marker %s at %s:%d due to excessive depth %d",pBedMatrix->currentMarker()->sMarkerID.c_str(), bUCSC ? pBedMatrix->currentMarker()->getUCSCChromosomeName().c_str() : pBedMatrix->currentMarker()->getNCBIChromosomeName().c_str(), pBedMatrix->currentMarker()->nBasePosition, nRGIndices.size());
      }
    }
    else {
      //abort();
    }
  }
  //Logger::gLogger->writeLog("Finished Reading %u records",numSectionRecords);
  

  Logger::gLogger->writeLog("Writing output file %s",sOutFile.c_str());

  double* bestRGAllIBDs = new double[nRGs]();
  uint32_t* bestRGAllSMIndices = new uint32_t[nRGs]();
  double* bestRGAllLLKs = new double[nRGs]();

  double* bestRGSelfIBDs = new double[nRGs]();
  double* bestRGSelfLLKs = new double[nRGs]();

  FILE* fpBest = NULL;
  double bestSMAllIBD = 0;
  uint32_t bestSMAllSMIndex = 0;
  double bestSMAllLLK = 0-DBL_MAX;
  double bestSMSelfLLK = 0-DBL_MAX;
  double bestSMSelfIBD = 0;
  
  if ( nCompInds > 0 ) {
    for(uint32_t i=0; i < nRGs; ++i) {
      bestRGAllLLKs[i] = 0-DBL_MAX;
      bestRGSelfLLKs[i] = 0-DBL_MAX;
      for(uint32_t j=0; j < nCompInds; ++j) {
	for(uint32_t k=0; k < numIBDs; ++k) {
	  if ( fSumRGIndLLKs[i*nCompInds*numIBDs+j*numIBDs+k] > bestRGAllLLKs[i] ) {
	    bestRGAllLLKs[i] = fSumRGIndLLKs[i*nCompInds*numIBDs+j*numIBDs+k];
	    bestRGAllSMIndices[i] = j;
	    bestRGAllIBDs[i] = vfIBDs[k];
	  }
	}
      }
      
      int32_t matchIndex = -1;
      for(uint32_t k=0; k < nCompInds; ++k) {
	if ( vsSMIDs[i].compare(pBedMatrix->vpIndividuals[k]->sIndID) == 0 ) {
	  matchIndex = k;
	  break;
	}
      }
      if ( matchIndex < 0 ) {
	bestRGSelfLLKs[i] = 0-DBL_MAX;
	bestRGSelfIBDs[i] = 0-DBL_MAX;
      }
      else {
	for(uint32_t k=0; k < numIBDs; ++k) {
	  if ( fSumRGIndLLKs[i*nCompInds*numIBDs+matchIndex*numIBDs+k] > bestRGSelfLLKs[i] ) {
	    bestRGSelfLLKs[i] = fSumRGIndLLKs[i*nCompInds*numIBDs+matchIndex*numIBDs+k];
	    bestRGSelfIBDs[i] = vfIBDs[k];
	  }
	}
      }
    }

    {
      for(uint32_t j=0; j < nCompInds; ++j) {
	for(uint32_t k=0; k < numIBDs; ++k) {
	  if ( fSumSMIndLLKs[j*numIBDs+k] > bestSMAllLLK ) {
	    bestSMAllLLK = fSumSMIndLLKs[j*numIBDs+k];
	    bestSMAllSMIndex = j;
	    bestSMAllIBD = vfIBDs[k];
	  }
	}
      }
      
      int32_t matchIndex = -1;
      for(uint32_t k=0; k < nCompInds; ++k) {
	if ( vsSMIDs[0].compare(pBedMatrix->vpIndividuals[k]->sIndID) == 0 ) {
	  matchIndex = k;
	  break;
	}
      }
      if ( matchIndex < 0 ) {
	bestSMSelfLLK = 0-DBL_MAX;
	bestSMSelfIBD = 0-DBL_MAX;
      }
      else {
	for(uint32_t k=0; k < numIBDs; ++k) {
	  if ( fSumSMIndLLKs[matchIndex*numIBDs+k] > bestSMSelfLLK ) {
	    bestSMSelfLLK = fSumSMIndLLKs[matchIndex*numIBDs+k];
	    bestSMSelfIBD = vfIBDs[k];
	  }
	}
      }
    }

    fpBest = fopen((sOutFile+".bestRG").c_str(),"w");
    if ( fpBest == NULL ) {
      Logger::gLogger->error("Cannot open output file %s",sOutFile.c_str());        
    }
    fprintf(fpBest,"SEQ_SM\tRG\tBEST_SM\tBESTIBD\tBESTIBDLLK\tBESTIBDLLK-\t#GENOS\t#BASES\t%%GENREF\t%%GENHET\t%%GENALT\tDPREF\tRDPHET\tRDPALT\tREF-A1%%\tREF-A2%%\tHET-A1%%\tHET-A2%%\tALT-A1%%\tALT-A2%%\t#DP\t%%HETAF\t%%HETSEQ\tEXHET\t%%MIX\t%%HOM\tBESTHOMMIXLLK\tBESTHOMMIXLLK-\n");
  }

  FILE* fpSelf = fopen((sOutFile+".selfRG").c_str(),"w");
  if ( fpSelf == NULL ) {
    Logger::gLogger->error("Cannot open output file %s",sOutFile.c_str());        
  }
  fprintf(fpSelf,"SEQ_SM\tRG\tSELF_SM\tSELFIBD\tSELFIBDLLK\tSELFIBDLLK-\t#GENOS\t#BASES\t%%GENREF\t%%GENHET\t%%GENALT\tDPREF\tRDPHET\tRDPALT\tREF-A1%%\tREF-A2%%\tHET-A1%%\tHET-A2%%\tALT-A1%%\tALT-A2%%\t#DP\t%%HETAF\t%%HETSEQ\tEXHET\t%%MIX\t%%HOM\tBESTHOMMIXLLK\tBESTHOMMIXLLK-\n");

  Logger::gLogger->writeLog("-----------------------------------------------------------------------");
  Logger::gLogger->writeLog("RG\t%HOM\t%MIX\tLLK");
  Logger::gLogger->writeLog("-----------------------------------------------------------------------");
  for(uint32_t i=0; i < nRGs; ++i) {
    for(int32_t l = 0-static_cast<int32_t>(numHomFracs)+1; 
	l < static_cast<int32_t>(numMixFracs); ++l) {
      uint32_t h = (l < 0) ? static_cast<uint32_t>(0-l) : 0;
      uint32_t k = (l > 0) ? static_cast<uint32_t>(l) : 0;
      uint32_t lu = static_cast<uint32_t>( l + numHomFracs -1 );
      Logger::gLogger->writeLog("%s\t%.3lf\t%.3lf\t%.3le",vsRGIDs[i].c_str(),vfHomFracs[h],vfMixFracs[k],fRGHomMixLLKs[i*(numMixFracs+numHomFracs-1) + lu]);
    }
  }

  for(uint32_t i=0; i < nRGs; ++i) {
    // check if the SMIDs exist in the individuals
    int32_t matchIndex = -1;
    uint32_t j = bestRGAllSMIndices[i];
    for(uint32_t k=0; k < nCompInds; ++k) {
      if ( vsSMIDs[i].compare(pBedMatrix->vpIndividuals[k]->sIndID) == 0 ) {
	matchIndex = k;
	break;
      }
    }

    double bestHetLLK = 0-DBL_MAX;
    double pureHetLLK, bestMixFrac = 0, bestHomFrac = 0;

    for(int32_t l = 0-static_cast<int32_t>(numHomFracs)+1; 
	l < static_cast<int32_t>(numMixFracs); ++l) {
      uint32_t h = (l < 0) ? static_cast<uint32_t>(0-l) : 0;
      uint32_t k = (l > 0) ? static_cast<uint32_t>(l) : 0;
      uint32_t lu = static_cast<uint32_t>( l + numHomFracs -1 );

      if ( bestHetLLK < fRGHomMixLLKs[i*(numMixFracs+numHomFracs-1) + lu] ) {
	bestHetLLK = fRGHomMixLLKs[i*(numMixFracs+numHomFracs-1) + lu];
	bestMixFrac = vfMixFracs[k];
	bestHomFrac = vfHomFracs[h];
      }
    }
    pureHetLLK = fRGHomMixLLKs[i*(numMixFracs+numHomFracs-1) + numHomFracs-1];

    if ( nCompInds > 0 ) {
      fprintf(fpBest,"%s",vsSMIDs[i].c_str());  // SEQ_SM
      fprintf(fpBest,"\t%s",vsRGIDs[i].c_str());  // RG
      fprintf(fpBest,"\t%s\t%.3lf\t%.3le",pBedMatrix->vpIndividuals[j]->sIndID.c_str(),bestRGAllIBDs[i],bestRGAllLLKs[i]); // BEST_SM, BESTIBD, BESTLLK
      if ( matchIndex < 0 ) {
	fprintf(fpBest,"\tN/A");
      }
      else {
	fprintf(fpBest,"\t%.3le",bestRGAllLLKs[i]-fSumRGIndLLKs[i*nCompInds*numIBDs+matchIndex*numIBDs+numIBDs-1]); // BESTLLK-
      }
    }
    
    fprintf(fpSelf,"%s",vsSMIDs[i].c_str());  // SEQ_SM
    fprintf(fpSelf,"\t%s",vsRGIDs[i].c_str());  // RG
    if ( matchIndex < 0 ) {
      fprintf(fpSelf,"\tN/A\tN/A\tN/A\tN/A");
    }
    else {
      fprintf(fpSelf,"\t%s\t%.3lf\t%.3le",vsSMIDs[i].c_str(),bestRGSelfIBDs[i],bestRGSelfLLKs[i]); // SELF_SM, SELFIBD, SELFLLK
      fprintf(fpSelf,"\t%.3le",bestRGSelfLLKs[i]-fSumRGIndLLKs[i*nCompInds*numIBDs+matchIndex*numIBDs+numIBDs-1]); // BESTLLK-
    }
      
    uint32_t nGenos, nBases;
    double refDepth, hetDepth, altDepth;
    if ( nCompInds > 0 ) {
      nGenos = nIndGenoMarkers[3*j+0]+nIndGenoMarkers[3*j+1]+nIndGenoMarkers[3*j+2];
      fprintf(fpBest,"\t%u",nGenos); // #GENOS
      nBases = nRGIndGenoMarkerBases[i*nCompInds*3+j*3+0]+nRGIndGenoMarkerBases[i*nCompInds*3+j*3+1]+nRGIndGenoMarkerBases[i*nCompInds*3+j*3+2];
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%u",nBases); // #BASES
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+0])/static_cast<double>(nGenos)); // %GENREF
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+1])/static_cast<double>(nGenos)); // %GENHET
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+2])/static_cast<double>(nGenos)); // %GENALT
      refDepth = static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+0])/static_cast<double>(nIndGenoMarkers[3*j+0]);
      hetDepth = static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+1])/static_cast<double>(nIndGenoMarkers[3*j+1]);
      altDepth = static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+2])/static_cast<double>(nIndGenoMarkers[3*j+2]);
      fprintf(fpBest,"\t%.4lf\t%.5lf\t%.5lf",refDepth,hetDepth/refDepth,altDepth/refDepth);
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nRGIndGenoMarkerRefs[i*nCompInds*3+j*3+0])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+0])); //REF_A1%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nRGIndGenoMarkerAlts[i*nCompInds*3+j*3+0])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+0])); //REF_A2%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nRGIndGenoMarkerRefs[i*nCompInds*3+j*3+1])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+1])); //HET_A1%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nRGIndGenoMarkerAlts[i*nCompInds*3+j*3+1])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+1])); //HET_A2%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nRGIndGenoMarkerRefs[i*nCompInds*3+j*3+2])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+2])); //ALT_A1%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nRGIndGenoMarkerAlts[i*nCompInds*3+j*3+2])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+2])); //ALT_A2%
    }
      
    if ( matchIndex < 0 ) {
      fprintf(fpSelf,"\t%u",nRGMarkers[i]); // #GENOS
      fprintf(fpSelf,"\t%u",nRGBases[i]); // #BASES
      fprintf(fpSelf,"\tN/A"); // %GENREF
      fprintf(fpSelf,"\tN/A"); // %GENHET
      fprintf(fpSelf,"\tN/A"); // %GENALT
      fprintf(fpSelf,"\tN/A"); // DPREF
      fprintf(fpSelf,"\tN/A"); // RDPHET
      fprintf(fpSelf,"\tN/A"); // RDPALT
      fprintf(fpSelf,"\tN/A"); // REF-A1%
      fprintf(fpSelf,"\tN/A"); // REF-A2%
      fprintf(fpSelf,"\tN/A"); // HET-A1%
      fprintf(fpSelf,"\tN/A"); // HET-A2%
      fprintf(fpSelf,"\tN/A"); // ALT-A1%
      fprintf(fpSelf,"\tN/A"); // ALT-A2%
    }
    else {
      j = matchIndex;
      nGenos = nIndGenoMarkers[3*j+0]+nIndGenoMarkers[3*j+1]+nIndGenoMarkers[3*j+2];
      fprintf(fpSelf,"\t%u",nGenos); // #GENOS
      nBases = nRGIndGenoMarkerBases[i*nCompInds*3+j*3+0]+nRGIndGenoMarkerBases[i*nCompInds*3+j*3+1]+nRGIndGenoMarkerBases[i*nCompInds*3+j*3+2];
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%u",nBases); // #BASES
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+0])/static_cast<double>(nGenos)); // %GENREF
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+1])/static_cast<double>(nGenos)); // %GENHET
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+2])/static_cast<double>(nGenos)); // %GENALT

      refDepth = static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+0])/static_cast<double>(nIndGenoMarkers[3*j+0]);
      hetDepth = static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+1])/static_cast<double>(nIndGenoMarkers[3*j+1]);
      altDepth = static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+2])/static_cast<double>(nIndGenoMarkers[3*j+2]);

      fprintf(fpSelf,"\t%.4lf\t%.5lf\t%.5lf",refDepth,hetDepth/refDepth,altDepth/refDepth);
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nRGIndGenoMarkerRefs[i*nCompInds*3+j*3+0])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+0])); //REF_A1%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nRGIndGenoMarkerAlts[i*nCompInds*3+j*3+0])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+0])); //REF_A2%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nRGIndGenoMarkerRefs[i*nCompInds*3+j*3+1])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+1])); //HET_A1%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nRGIndGenoMarkerAlts[i*nCompInds*3+j*3+1])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+1])); //HET_A2%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nRGIndGenoMarkerRefs[i*nCompInds*3+j*3+2])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+2])); //ALT_A1%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nRGIndGenoMarkerAlts[i*nCompInds*3+j*3+2])/static_cast<double>(nRGIndGenoMarkerBases[i*nCompInds*3+j*3+2])); //ALT_A2%
    }
    
    //uint32_t nDP2 = 0;
    //double fEHetDP2AF = 0.;
    //double fEHetDP2SEQ = 0.;
    //for(uint32_t k=2; k <= maxDepth; ++k) {
    //  nDP2 += nRGMarkerDepth[i*(maxDepth+1)+k];
    //  fEHetDP2AF += fRGSumEHetDepth[i*(maxDepth+1)+k];
    //  fEHetDP2SEQ += (static_cast<double>(nRGMarkerDepthRefAlt[i*(maxDepth+1)+k])/(1.-pow(0.5,static_cast<double>(k-1))));
    //}
    //fEHetDP2AF /= static_cast<double>(nDP2);
    //fEHetDP2SEQ /= static_cast<double>(nDP2);

    double DP = static_cast<double>(nRGBases[i])/static_cast<double>(nMarkerCnt);
    double fHetAF = fRGSumPriorHets[i]/static_cast<double>(nRGMarkers[i]);
    double fHetSEQ = fRGSumPosteriorHets[i]/static_cast<double>(nRGMarkers[i]);
    double fHetRatio = fHetSEQ/fHetAF;

    if ( nCompInds > 0 ) {
      //fprintf(fpBest,"\t%u\t%.5lf\t%.5lf\t%.5lf",nDP2,fEHetDP2AF,fEHetDP2SEQ,fEHetDP2SEQ/fEHetDP2AF);
      fprintf(fpBest,"\t%.3lf\t%.5lf\t%.5lf\t%.5lf",DP,fHetAF,fHetSEQ,fHetRatio);
      fprintf(fpBest,"\t%.3lf\t%.3lf\t%.3le\t%.3le",bestMixFrac,bestHomFrac,bestHetLLK,bestHetLLK-pureHetLLK);
      fprintf(fpBest,"\n");
    }
    
    fprintf(fpSelf,"\t%.3lf\t%.5lf\t%.5lf\t%.5lf",DP,fHetAF,fHetSEQ,fHetRatio);
    fprintf(fpSelf,"\t%.3lf\t%.3lf\t%.3le\t%.3le",bestMixFrac,bestHomFrac,bestHetLLK,bestHetLLK-pureHetLLK);
    fprintf(fpSelf,"\n");
  }
  if ( nCompInds > 0 ) fclose(fpBest);
  fclose(fpSelf);

  if ( sameSMFlag ) {
    if ( nCompInds > 0 ) {
      fpBest = fopen((sOutFile+".bestSM").c_str(),"w");
      if ( fpBest == NULL ) {
	Logger::gLogger->error("Cannot open output file %s",sOutFile.c_str());    
      }
      fprintf(fpBest,"SEQ_SM\tRG\tBEST_SM\tBESTIBD\tBESTIBDLLK\tBESTIBDLLK-\t#GENOS\t#BASES\t%%GENREF\t%%GENHET\t%%GENALT\tDPREF\tRDPHET\tRDPALT\tREF-A1%%\tREF-A2%%\tHET-A1%%\tHET-A2%%\tALT-A1%%\tALT-A2%%\t#DP\t%%HETAF\t%%HETSEQ\tEXHET\t%%MIX\t%%HOM\tBESTHOMMIXLLK\tBESTHOMMIXLLK-\n");
    }
    
    fpSelf = fopen((sOutFile+".selfSM").c_str(),"w");
    if ( fpSelf == NULL ) {
      Logger::gLogger->error("Cannot open output file %s",sOutFile.c_str());    
    }
    fprintf(fpSelf,"SEQ_SM\tRG\tSELF_SM\tSELFIBD\tSELFIBDLLK\tSELFIBDLLK-\t#GENOS\t#BASES\t%%GENREF\t%%GENHET\t%%GENALT\tDPREF\tRDPHET\tRDPALT\tREF-A1%%\tREF-A2%%\tHET-A1%%\tHET-A2%%\tALT-A1%%\tALT-A2%%\t#DP\t%%HETAF\t%%HETSEQ\tEXHET\t%%MIX\t%%HOM\tBESTHOMMIXLLK\tBESTHOMMIXLLK-\n");

    Logger::gLogger->writeLog("-----------------------------------------------------------------------");
    Logger::gLogger->writeLog("SM\t%HOM\t%MIX\tLLK");
    Logger::gLogger->writeLog("-----------------------------------------------------------------------");

    for(int32_t l = 0-static_cast<int32_t>(numHomFracs)+1; 
	l < static_cast<int32_t>(numMixFracs); ++l) {
      uint32_t h = (l < 0) ? static_cast<uint32_t>(0-l) : 0;
      uint32_t k = (l > 0) ? static_cast<uint32_t>(l) : 0;
      uint32_t lu = static_cast<uint32_t>( l + numHomFracs -1 );
      Logger::gLogger->writeLog("%s\t%.3lf\t%.3lf\t%.3le",sBamSMID.c_str(),vfHomFracs[h],vfMixFracs[k],fSMHomMixLLKs[lu]);
    }

    // check if the SMIDs exist in the individuals
    int32_t matchIndex = -1;
    uint32_t j = bestSMAllSMIndex;
    for(uint32_t k=0; k < nCompInds; ++k) {
      if ( vsSMIDs[0].compare(pBedMatrix->vpIndividuals[k]->sIndID) == 0 ) {
	matchIndex = k;
	break;
      }
    }

    if ( nCompInds > 0 ) {
      fprintf(fpBest,"%s",vsSMIDs[0].c_str());  // SEQ_SM
      fprintf(fpBest,"\tN/A");  // RG = N/A
      fprintf(fpBest,"\t%s\t%.3lf\t%.3le",pBedMatrix->vpIndividuals[j]->sIndID.c_str(),bestSMAllIBD,bestSMAllLLK); // BEST_SM, BESTIBD, BESTLLK
      if ( matchIndex < 0 ) {
	fprintf(fpBest,"\tN/A");
      }
      else {
	fprintf(fpBest,"\t%.3le",bestSMAllLLK-fSumSMIndLLKs[matchIndex*numIBDs+numIBDs-1]); // BESTLLK-
      }
    }
    
    fprintf(fpSelf,"%s",vsSMIDs[0].c_str());  // SEQ_SM
    fprintf(fpSelf,"\tN/A");  // RG
    if ( matchIndex < 0 ) {
      fprintf(fpSelf,"\tN/A\tN/A\tN/A\tN/A");
    }
    else {
      fprintf(fpSelf,"\t%s\t%.3lf\t%.3le",vsSMIDs[0].c_str(),bestSMSelfIBD,bestSMSelfLLK); // SELF_SM, SELFIBD, SELFLLK
      fprintf(fpSelf,"\t%.3le",bestSMSelfLLK-fSumSMIndLLKs[matchIndex*numIBDs+numIBDs-1]); // BESTLLK-
    }
      
    uint32_t nGenos, nBases;
    double refDepth, hetDepth, altDepth;
    if ( nCompInds > 0 ) {
      nGenos = nIndGenoMarkers[3*j+0]+nIndGenoMarkers[3*j+1]+nIndGenoMarkers[3*j+2];
      fprintf(fpBest,"\t%u",nGenos); // #GENOS
      nBases = nSMIndGenoMarkerBases[j*3+0]+nSMIndGenoMarkerBases[j*3+1]+nSMIndGenoMarkerBases[j*3+2];
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%u",nBases); // #BASES
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+0])/static_cast<double>(nGenos)); // %GENREF
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+1])/static_cast<double>(nGenos)); // %GENHET
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+2])/static_cast<double>(nGenos)); // %GENALT
      refDepth = static_cast<double>(nSMIndGenoMarkerBases[j*3+0])/static_cast<double>(nIndGenoMarkers[3*j+0]);
      hetDepth = static_cast<double>(nSMIndGenoMarkerBases[j*3+1])/static_cast<double>(nIndGenoMarkers[3*j+1]);
      altDepth = static_cast<double>(nSMIndGenoMarkerBases[j*3+2])/static_cast<double>(nIndGenoMarkers[3*j+2]);
      fprintf(fpBest,"\t%.4lf\t%.5lf\t%.5lf",refDepth,hetDepth/refDepth,altDepth/refDepth);
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nSMIndGenoMarkerRefs[j*3+0])/static_cast<double>(nSMIndGenoMarkerBases[j*3+0])); //REF_A1%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nSMIndGenoMarkerAlts[j*3+0])/static_cast<double>(nSMIndGenoMarkerBases[j*3+0])); //REF_A2%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nSMIndGenoMarkerRefs[j*3+1])/static_cast<double>(nSMIndGenoMarkerBases[j*3+1])); //HET_A1%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nSMIndGenoMarkerAlts[j*3+1])/static_cast<double>(nSMIndGenoMarkerBases[j*3+1])); //HET_A2%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nSMIndGenoMarkerRefs[j*3+2])/static_cast<double>(nSMIndGenoMarkerBases[j*3+2])); //ALT_A1%
      fprintf(fpBest,"\t");
      fprintf(fpBest,"%.5lf",static_cast<double>(nSMIndGenoMarkerAlts[j*3+2])/static_cast<double>(nSMIndGenoMarkerBases[j*3+2])); //ALT_A2%
    }

    if ( matchIndex < 0 ) {
      fprintf(fpSelf,"\t%u",nSMMarkers); // #GENOS
      fprintf(fpSelf,"\t%u",nSMBases); // #BASES
      fprintf(fpSelf,"\tN/A"); // %GENREF
      fprintf(fpSelf,"\tN/A"); // %GENHET
      fprintf(fpSelf,"\tN/A"); // %GENALT
      fprintf(fpSelf,"\tN/A"); // DPREF
      fprintf(fpSelf,"\tN/A"); // RDPHET
      fprintf(fpSelf,"\tN/A"); // RDPALT
      fprintf(fpSelf,"\tN/A"); // REF-A1%
      fprintf(fpSelf,"\tN/A"); // REF-A2%
      fprintf(fpSelf,"\tN/A"); // HET-A1%
      fprintf(fpSelf,"\tN/A"); // HET-A2%
      fprintf(fpSelf,"\tN/A"); // ALT-A1%
      fprintf(fpSelf,"\tN/A"); // ALT-A2%
    }
    else {
      j = matchIndex;
      nGenos = nIndGenoMarkers[3*j+0]+nIndGenoMarkers[3*j+1]+nIndGenoMarkers[3*j+2];
      fprintf(fpSelf,"\t%u",nGenos); // #GENOS
      nBases = nSMIndGenoMarkerBases[j*3+0]+nSMIndGenoMarkerBases[j*3+1]+nSMIndGenoMarkerBases[j*3+2];
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%u",nBases); // #BASES
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+0])/static_cast<double>(nGenos)); // %GENREF
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+1])/static_cast<double>(nGenos)); // %GENHET
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nIndGenoMarkers[3*j+2])/static_cast<double>(nGenos)); // %GENALT
      refDepth = static_cast<double>(nSMIndGenoMarkerBases[j*3+0])/static_cast<double>(nIndGenoMarkers[3*j+0]);
      hetDepth = static_cast<double>(nSMIndGenoMarkerBases[j*3+1])/static_cast<double>(nIndGenoMarkers[3*j+1]);
      altDepth = static_cast<double>(nSMIndGenoMarkerBases[j*3+2])/static_cast<double>(nIndGenoMarkers[3*j+2]);
      fprintf(fpSelf,"\t%.4lf\t%.5lf\t%.5lf",refDepth,hetDepth/refDepth,altDepth/refDepth);
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nSMIndGenoMarkerRefs[j*3+0])/static_cast<double>(nSMIndGenoMarkerBases[j*3+0])); //REF_A1%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nSMIndGenoMarkerAlts[j*3+0])/static_cast<double>(nSMIndGenoMarkerBases[j*3+0])); //REF_A2%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nSMIndGenoMarkerRefs[j*3+1])/static_cast<double>(nSMIndGenoMarkerBases[j*3+1])); //HET_A1%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nSMIndGenoMarkerAlts[j*3+1])/static_cast<double>(nSMIndGenoMarkerBases[j*3+1])); //HET_A2%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nSMIndGenoMarkerRefs[j*3+2])/static_cast<double>(nSMIndGenoMarkerBases[j*3+2])); //ALT_A1%
      fprintf(fpSelf,"\t");
      fprintf(fpSelf,"%.5lf",static_cast<double>(nSMIndGenoMarkerAlts[j*3+2])/static_cast<double>(nSMIndGenoMarkerBases[j*3+2])); //ALT_A2%
    }
    
    //uint32_t nDP2 = 0;
    //double fEHetDP2AF = 0.;
    //double fEHetDP2SEQ = 0.;
    //for(uint32_t k=2; k <= maxDepth; ++k) {
    //  nDP2 += nSMMarkerDepth[k];
    //  fEHetDP2AF += fSMSumEHetDepth[k];
    //  fEHetDP2SEQ += (static_cast<double>(nSMMarkerDepthRefAlt[k])/(1.-pow(0.5,static_cast<double>(k-1))));
    //}
    //fEHetDP2AF /= static_cast<double>(nDP2);
    //fEHetDP2SEQ /= static_cast<double>(nDP2);

    double DP = static_cast<double>(nSMBases)/static_cast<double>(nMarkerCnt);
    double fHetAF = fSMSumPriorHets/static_cast<double>(nSMMarkers);
    double fHetSEQ = fSMSumPosteriorHets/static_cast<double>(nSMMarkers);
    double fHetRatio = fHetSEQ/fHetAF;

    //if ( nCompInds > 0 ) fprintf(fpBest,"\t%u\t%.5lf\t%.5lf\t%.5lf",nDP2,fEHetDP2AF,fEHetDP2SEQ,fEHetDP2SEQ/fEHetDP2AF);
    //fprintf(fpSelf,"\t%u\t%.5lf\t%.5lf\t%.5lf",nDP2,fEHetDP2AF,fEHetDP2SEQ,fEHetDP2SEQ/fEHetDP2AF);
    if ( nCompInds > 0 ) fprintf(fpBest,"\t%.3lf\t%.5lf\t%.5lf\t%.5lf",DP,fHetAF,fHetSEQ,fHetRatio);
    fprintf(fpSelf,"\t%.3lf\t%.5lf\t%.5lf\t%.5lf",DP,fHetAF,fHetSEQ,fHetRatio);

    double bestHetLLK = 0-DBL_MAX;
    double pureHetLLK, bestMixFrac = 0, bestHomFrac = 0;
    for(int32_t l = 0-static_cast<int32_t>(numHomFracs)+1; 
	l < static_cast<int32_t>(numMixFracs); ++l) {
      uint32_t h = (l < 0) ? static_cast<uint32_t>(0-l) : 0;
      uint32_t k = (l > 0) ? static_cast<uint32_t>(l) : 0;
      uint32_t lu = static_cast<uint32_t>( l + numHomFracs -1 );

      if ( bestHetLLK < fSMHomMixLLKs[lu] ) {
	bestHetLLK = fSMHomMixLLKs[lu];
	bestMixFrac = vfMixFracs[k];
	bestHomFrac = vfHomFracs[h];
      }
    }
    pureHetLLK = fSMHomMixLLKs[numHomFracs-1];
    
    if ( nCompInds > 0 ) fprintf(fpBest,"\t%.3lf\t%.3lf\t%.3le\t%.3le",bestMixFrac,bestHomFrac,bestHetLLK,bestHetLLK-pureHetLLK);
    fprintf(fpSelf,"\t%.3lf\t%.3lf\t%.3le\t%.3le",bestMixFrac,bestHomFrac,bestHetLLK,bestHetLLK-pureHetLLK);

    if ( nCompInds > 0 ) fprintf(fpBest,"\n");
    fprintf(fpSelf,"\n");
  }
  if ( nCompInds > 0 ) fclose(fpBest);
  fclose(fpSelf);

  Logger::gLogger->writeLog("Finished writing output files %s.{bestRG,selfRG,bestSM,selfSM}",sOutFile.c_str());

  /*
  delete [] nIndGenoMarkers;
  delete [] nRGIndGenoMarkerBases;
  delete [] nRGIndGenoMarkerRefs;
  delete [] nRGIndGenoMarkerAlts;
  delete [] fSumRGIndLLKs;
  delete [] fRGGenoIBDLKs;
  delete [] fRGnonIBDLKs;
  delete [] bestIBDs;
  delete [] bestSMIndices;
  delete [] bestLLKs;
  delete [] nRGMarkerDepth;
  delete [] nRGMarkerDepthRefOnly;
  delete [] nRGMarkerDepthRefAlt;
  delete [] nRGMarkerDepthAltOnly;
  delete [] nTotalMarkerDepth;
  delete [] nTotalMarkerDepthRefOnly;
  delete [] nTotalMarkerDepthRefAlt;
  delete [] nTotalMarkerDepthAltOnly;
  delete [] nRGRefCnts;
  delete [] nRGAltCnts;
  delete [] nRGEtcCnts;*/

  return 0;
}
