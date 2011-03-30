#include <math.h>
//#include <boost/thread/mutex.hpp>

#include "VcfFile.h"
#include "FilterStat.h"
#include "Logger.h"

bool FilterStat::loadAnchorVcf(const char* file) {
  if ( !sAnchorVcf.IsEmpty() ) {
    Logger::gLogger->error("Anchor VCF is specified twice");
  }
  sAnchorVcf = file;
  VcfFile vcf;
  vcf.setSiteOnly(true);
  vcf.openForRead(file,1);

  VcfMarker* pMarker;
  for( int i=0; vcf.iterateMarker(); ++i ) {
    pMarker = vcf.getLastMarker();

    if ( i == 0 ) {
      sChrom = pMarker->sChrom;
    }
    else {
      if ( sChrom.Compare(pMarker->sChrom) != 0 ) {
	Logger::gLogger->error("Only a single-chromosome VCF is currently supported");
      }
    }

    vPos.push_back(pMarker->nPos);
    if ( pMarker->asAlts.Length() > 1 ) {
      vAl1.push_back(pMarker->asAlts[0].c_str());
      vAl2.push_back(pMarker->asAlts[1].c_str());
    }
    else {
      vAl1.push_back(pMarker->sRef.c_str());
      vAl2.push_back(pMarker->asAlts[0].c_str());
    }

    for(int j=0; j < pMarker->asInfoKeys.Length(); ++j) {
      if ( pMarker->asInfoKeys[j].Compare("AF") == 0 ) {
	vAFs.push_back(atof(pMarker->asInfoValues[j]));
	break;
      }
    }

    for(int j=0; j < FILTER_STAT_COUNTS; ++j) {
      vCounts.push_back(0);
    }
  }
  return true;
}

bool FilterStat::appendStatVcf(const char* file) {
  VcfFile vcf;
  vcf.setSiteOnly(false);
  vcf.setParseValues(true);
  vcf.setParseGenotypes(false);
  vcf.setParseDosages(false);
  vcf.openForRead(file,1);
  
  VcfMarker* pMarker;
  StringArray tok;
  for( int i=0, j=0; vcf.iterateMarker(); ++i, ++j ) {
    pMarker = vcf.getLastMarker();
    if ( sChrom.Compare(pMarker->sChrom) != 0 ) {
      Logger::gLogger->error("Chromosome name does not match - %s vs %s",sChrom.c_str(),pMarker->sChrom.c_str());
    }

    while ( vPos[j] < pMarker->nPos ) { ++j; }

    if ( vPos[j] > pMarker->nPos ) {
      Logger::gLogger->error("Position %s:%d is not observed in the anchor VCF",sChrom.c_str(),pMarker->nPos);
    }

    std::vector<int> vAlleles;
    std::vector<int> vStrands;

    //fprintf(stderr,"%s:%d\t%s\n",pMarker->sChrom.c_str(),pMarker->nPos,pMarker->asFormatKeys[0].c_str());

    for(int k=0; k < pMarker->asFormatKeys.Length(); ++k) {
      if ( pMarker->asFormatKeys[k].Compare("BASE") == 0 ) {
	tok.ReplaceColumns(pMarker->asSampleValues[k],',');
	for(int l=0; l < tok.Length(); ++l) {
	  if ( tok[l].Compare(vAl1[j].c_str()) == 0 ) {
	    vAlleles.push_back(0);
	  }
	  else if ( tok[l].Compare(vAl2[j].c_str()) == 0 ) {
	    vAlleles.push_back(1);
	  }
	  else {
	    vAlleles.push_back(2);
	  }
	}
      }
      else if ( pMarker->asFormatKeys[k].Compare("STRAND") == 0 ) {
	tok.ReplaceColumns(pMarker->asSampleValues[k],',');
	for(int l=0; l < tok.Length(); ++l) {
	  if ( tok[l].Compare("F") == 0 ) {
	    vStrands.push_back(0);
	  }
	  else {
	    vStrands.push_back(1);
	  }
	}
      }
    }

    //fprintf(stderr,"%s:%d\t%d",pMarker->sChrom.c_str(),pMarker->nPos,(int)vAlleles.size());
    //for(int k=0; k < (int) vAlleles.size(); ++k) {
    //  fprintf(stderr,"\t%d",vAlleles[k]*2+vStrands[k]);
    //}
    //fprintf(stderr,"\n");

    if ( vAlleles.size() != vStrands.size() ) {
      Logger::gLogger->error("Alleles and Strands do not match in size at %s:%d, in %s",pMarker->sChrom.c_str(), pMarker->nPos, file);
    }

    // updates the counts - needs synchronization
    {
      //boost::mutex::scoped_lock lock(mutex);
      for(int k=0; k < (int) vAlleles.size(); ++k) {
	++(vCounts[FILTER_STAT_COUNTS*j + vAlleles[k]*2 + vStrands[k]]);
      }
    }
  }
  return true;
}

bool FilterStat::writeMergedVcf(const char* outFile) {
  IFILE oFile = ifopen(outFile,"wb");
  if ( oFile == NULL ) {
    Logger::gLogger->error("Cannot open output file %s",outFile);
  }

  VcfFile vcf;
  vcf.setSiteOnly(false);
  vcf.setParseValues(true);
  vcf.openForRead(sAnchorVcf.c_str(),1);  

  vcf.printVCFHeader(oFile);

  VcfMarker* pMarker;
  String STC, STR;
  for( int i=0; vcf.iterateMarker(); ++i ) {
     pMarker = vcf.getLastMarker();

     int c[FILTER_STAT_COUNTS];
     for(int j=0; j < FILTER_STAT_COUNTS; ++j) {
       c[j] = vCounts[FILTER_STAT_COUNTS*i+j];
     }
     STC.printf("%d,%d,%d,%d,%d,%d",c[0],c[1],c[2],c[3],c[4],c[5]);

     if ( ( c[0]+c[1] > 4 ) && ( c[1]+c[3] > 4 ) && ( c[0]+c[2] > 4 ) && ( c[1]+c[3] > 4 ) ) { 
       STR.printf("%.2lf",((c[0]+.5)*(c[3]+.5)-(c[1]+.5)*(c[2]+.5))/sqrt((c[0]+c[1]+1.)*(c[2]+c[3]+1.)*(c[0]+c[2]+1.)*(c[1]+c[3]+1.)));
     }
     else {
       STR = "0";
     }
     pMarker->asInfoKeys.Add("STC");
     pMarker->asInfoKeys.Add("STR");
     pMarker->asInfoValues.Add(STC);
     pMarker->asInfoValues.Add(STR);

     pMarker->printVCFMarker(oFile,false);
  }
  ifclose(oFile);
  return true;
}
