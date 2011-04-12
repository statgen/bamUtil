////////////////////////////////////////////////////////////////////// 
// thunder/Main.cpp 
// (c) 2000-2008 Goncalo Abecasis
// 
// This file is distributed as part of the MaCH source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile MaCH.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Saturday April 12, 2008
// 

//static int NUM_NON_GLF_FIELDS = 9;

#include <vector>
#include <map>
 
#include "ShotgunHaplotyper.h"
#include "ShotgunManners.h"
#include "OutputHandlers.h"
#include "DosageCalculator.h"
#include "MergeHaplotypes.h"
#include "HaplotypeLoader.h"
#include "Parameters.h"
#include "InputFile.h"
#include "Error.h"
#include "VcfFile.h"

float * thetas = NULL;
int     nthetas = 0;

float * error_rates = NULL;
int     nerror_rates = 0;

// print output files directly in VCF format
// inVcf contains skeleton of VCF information to copy from
// consensus contains the haplotype information to replace GT field
// dosage contains dosage information to be added as DS field
// thetas contains recombination rate information between markers
// error-rates contains per-marker error rates
// rsqs contains rsq_hat estimates ??
void OutputVCFConsensus(const String& inVcf, Pedigree & ped, ConsensusBuilder & consensus, DosageCalculator &doses, const String& filename, float* thetas, float* error_rates)
{
  consensus.Merge(); // calculate consensus sequence

  if (consensus.stored)
    printf("Merged sampled haplotypes to generate consensus\n"
	   "   Consensus estimated to have %.1f errors in missing genotypes and %.1f flips in haplotypes\n\n",
	   consensus.errors, consensus.flips);

  // read and write VCF inputs
  try {

    fprintf(stderr,"Outputing VCF file %s\n",filename.c_str());

    VcfFile* pVcf = new VcfFile;
    IFILE outVCF = ifopen(filename.c_str(),"wb");
    if ( outVCF == NULL ) {
      error("Cannot open output file %s for writing",filename.c_str());
      exit(-1);
    }

    pVcf->bSiteOnly = false;
    pVcf->bParseGenotypes = false;
    pVcf->bParseDosages = false;
    pVcf->bParseValues = true;
    pVcf->openForRead(inVcf.c_str());
    pVcf->printVCFHeader(outVCF); // print header file

    char** haplotypes = consensus.consensus;

    // check the sanity of data
    if ( pVcf->getSampleCount() == 0 ) {
      throw VcfFileException("No individual genotype information exist in the input VCF file %s",filename.c_str());
    }

    int nSamples = pVcf->getSampleCount();

    // build map of personID -> sampleIndex
    std::map<std::string,int> pedMap;
    for (int i=0; i < ped.count; ++i) {
      pedMap[ped[i].pid.c_str()] = i;
      //fprintf(stderr,"Adding (%s,%d)\n", ped[i].pid.c_str(), i);
    }

    std::vector<int> vcf2ped;
    for (int i=0; i < nSamples; ++i) {
      std::map<std::string,int>::iterator found = pedMap.find(pVcf->vpVcfInds[i]->sIndID.c_str());
      if ( found == pedMap.end() ) {
	error("Cannot find individual ID %s",pVcf->vpVcfInds[i]->sIndID.c_str());
	exit(-1);
      }
      else {
	//fprintf(stderr,"Found (%s,%d)\n", pVcf->vpVcfInds[i]->sIndID.c_str(), found->second);
	vcf2ped.push_back(found->second);
      }
    }

    // read VCF lines
    VcfMarker* pMarker = new VcfMarker;
    
    char sDose[255];
    double freq, maf, avgPost, rsq;
    
    for(int m = 0; pVcf->iterateMarker(); ++m ) {
      //fprintf(stderr,"m=%d\n",m);

      pMarker = pVcf->getLastMarker();

      doses.CalculateMarkerInfo(m, freq, maf, avgPost, rsq);

      //fprintf(stderr,"foo1\n");

      sprintf(sDose,"%.4lf",1.-freq);
      pMarker->asInfoKeys.Add("LDAF");
      pMarker->asInfoValues.Add(sDose);
      sprintf(sDose,"%.4lf",avgPost);
      pMarker->asInfoKeys.Add("AVGPOST");
      pMarker->asInfoValues.Add(sDose);
      sprintf(sDose,"%.4lf",rsq);
      pMarker->asInfoKeys.Add("RSQ");
      pMarker->asInfoValues.Add(sDose);
      sprintf(sDose,"%.4lf",error_rates[m]);
      pMarker->asInfoKeys.Add("ERATE");
      pMarker->asInfoValues.Add(sDose);

      //fprintf(stderr,"foo2\n");
      
      int GTidx = pMarker->asFormatKeys.Find("GT");
      if ( GTidx < -1 ) {
	throw VcfFileException("Cannot recognize GT key in FORMAT field");
      }
      pMarker->asFormatKeys.InsertAt(GTidx+1, "DS");
      int DSidx = GTidx + 1;
      
      int nFormats = pMarker->asFormatKeys.Length();

      //fprintf(stderr,"foo3\n");
      
      for(int i=0; i < nSamples; ++i) {
	int pi = vcf2ped[i];
	// modify GT values;
	//fprintf(stderr,"i=%d, pi=%d, GTidx = %d, nFormats = %d, asSampleValues.Length() = %d, haplotypes = %x\n",i,pi,GTidx,nFormats,pMarker->asSampleValues.Length(), haplotypes);
	if ( pMarker->asAlts.Length() == 1 ) {
	  pMarker->asSampleValues[nFormats*i + GTidx].printf("%d|%d",haplotypes[pi*2][m],haplotypes[pi*2+1][m]);
	}
	else {
	  pMarker->asSampleValues[nFormats*i + GTidx].printf("%d|%d",haplotypes[pi*2][m]+1,haplotypes[pi*2+1][m]+1);
	}
	// add DS values
	sprintf(sDose,"%.3lf",2-doses.GetDosage(pi,m));
	pMarker->asSampleValues.InsertAt(nFormats*i + DSidx, sDose);
      }
      //fprintf(stderr,"foo4\n");
      pMarker->printVCFMarker(outVCF,false); // print marker to output file
    }
    delete pVcf;
    //delete pMarker;
    ifclose(outVCF);
  }
  catch ( VcfFileException e ) {
    error ( e.what() );
  }
}

void UpdateVector(float * current, float * & vector, int & n, int length)
   {
   if (n++ == 0)
      {
      vector = new float[length];

      for (int i = 0; i < length; i++)
         vector[i] = current[i];
      }
   else
      for (int i = 0; i < length; i++)
         vector[i] += current[i];
   }

void UpdateErrorRates(Errors * current, float * & vector, int & n, int length)
   {
   if (n++ == 0)
      {
      vector = new float[length];

      for (int i = 0; i < length; i++)
         vector[i] = current[i].rate;
      }
   else
      for (int i = 0; i < length; i++)
         vector[i] += current[i].rate;
   }

int MemoryAllocationFailure()
   {
   printf("FATAL ERROR - Memory allocation failed\n");
   return -1;
   }

void LoadShotgunSamples(Pedigree &ped, const String & filename) {
  printf("starting LoadShotgunSamples\n\n");
   try {
     VcfFile* pVcf = new VcfFile;
     pVcf->bSiteOnly = true;
     pVcf->bParseGenotypes = false;
     pVcf->bParseDosages = false;
     pVcf->bParseValues = false;
     pVcf->openForRead(filename.c_str());

     // check the sanity of data
     if ( pVcf->getSampleCount() == 0 ) {
       throw VcfFileException("No individual genotype information exist in the input VCF file %s",filename.c_str());
     }

     for(int i = 0; i < pVcf->getSampleCount(); ++i) {
       ped.AddPerson(pVcf->vpVcfInds[i]->sIndID, pVcf->vpVcfInds[i]->sIndID, "0", "0", 1);
     }

     delete pVcf;
   }
   catch ( VcfFileException e ) {
     error ( e.what() );
   }
   
   ped.Sort();
   printf("Loaded %d individuals from shotgun sequence data\n\n", ped.count);
}

void LoadPolymorphicSites(const String& filename) {
   try {
     VcfFile* pVcf = new VcfFile;
     pVcf->bSiteOnly = true;
     pVcf->bParseGenotypes = false;
     pVcf->bParseDosages = false;
     pVcf->bParseValues = false;
     pVcf->openForRead(filename.c_str());

     VcfMarker* pMarker = new VcfMarker;

     StringArray altalleles;
     String markerName;

     while( pVcf->iterateMarker() ) {
       int markers = Pedigree::markerCount;
       pMarker = pVcf->getLastMarker();

       markerName.printf("%s:%d",pMarker->sChrom.c_str(),pMarker->nPos);
       int marker = Pedigree::GetMarkerID(markerName);
       int al1, al2;

       //printf("Re-opening VCF file\n");

       if ( pMarker->asAlts.Length() == 2 ) {
	 al1 = Pedigree::LoadAllele(marker, pMarker->asAlts[0]);
	 al2 = Pedigree::LoadAllele(marker, pMarker->asAlts[1]);
       }
       else {
	 al1 = Pedigree::LoadAllele(marker, pMarker->sRef);
	 al2 = Pedigree::LoadAllele(marker, pMarker->asAlts[0]);
       }

       //printf("Re-opening VCF file\n");

       if ( markers != marker ) {
	 error("Each polymorphic site should only occur once, but site %s is duplicated\n", markerName.c_str());
       }


       if (al1 != 1 || al2 != 2) {
	 error("Allele labels '%s' and '%s' for polymorphic site '%s' are not valid\n", (const char *) altalleles[0], (const char *) altalleles[1], markerName.c_str());   
       }
     }
     delete pVcf;
     //delete pMarker;
   }
   catch ( VcfFileException e ) {
     error(e.what());
   }
}

void LoadShotgunResults(Pedigree &ped, char** genotypes, char* refalleles, double* freq1s, const String & filename) {
  printf("starting LoadShotgunResults\n\n");

   try {
     VcfFile* pVcf = new VcfFile;
     pVcf->bSiteOnly = false;
     pVcf->bParseGenotypes = false;
     pVcf->bParseDosages = false;
     pVcf->bParseValues = true;
     pVcf->openForRead(filename.c_str());

     // check the sanity of data
     if ( pVcf->getSampleCount() == 0 ) {
       throw VcfFileException("No individual genotype information exist in the input VCF file %s",filename.c_str());
     }

     int nSamples = pVcf->getSampleCount();
     int * personIndices = new int [ped.count];
     StringIntHash originalPeople; // key: famid+subID, value: original order (0 based);
     int person = 0;          
     for (int i = 0; i < nSamples; i++) {
       originalPeople.Add(pVcf->vpVcfInds[i]->sIndID+"."+pVcf->vpVcfInds[i]->sIndID, person);
       person++;
     }
     
     for (int i = 0; i < ped.count; i++) {
       personIndices[originalPeople.Integer(ped[i].famid+"."+ped[i].pid)] = i;
     }
     
     int markerindex = 0;
     
     printf("starting LoadPolymorphicSites\n\n");
     VcfMarker* pMarker = new VcfMarker;

     while( pVcf->iterateMarker() ) {
       refalleles[markerindex] = 1;
     
       pMarker = pVcf->getLastMarker();

       int AFidx = pMarker->asInfoKeys.Find("AF");
       if ( AFidx == -1 ) {
	 throw VcfFileException("Cannot recognize AF key in FORMAT field");
       }
       if ( pMarker->asAlts.Length() == 1 ) {
	 freq1s[markerindex] = (1.-pMarker->asInfoValues[AFidx].AsDouble());
       }
       else {
	 freq1s[markerindex] = pMarker->asInfoValues[AFidx].AsDouble();
       }

      StringArray phred;
      int PLidx = pMarker->asFormatKeys.Find("PL");
      int formatLength = pMarker->asFormatKeys.Length();
      if ( PLidx < 0 ) {
	PLidx = pMarker->asFormatKeys.Find("GL");
      }
      if ( PLidx < 0 ) {
	fprintf(stderr,"Missing PL key in FORMAT field\n");
      }
      int genoindex = markerindex*3;

      for (int i = 0; i < ped.count; i++)
	{
	  phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i*formatLength], ",");
	  int phred11 = phred[0].AsInteger();
	  int phred12 = phred[1].AsInteger();
	  int phred22 = phred[2].AsInteger();
	  //printf("phred scores are %d, %d, %d\n", phred11, phred12, phred22);
	  
	  if (phred11 > 255) phred11 = 255;
          if (phred12 > 255) phred12 = 255;
          if (phred22 > 255) phred22 = 255;

	  genotypes[personIndices[i]][genoindex] = phred11;
	  genotypes[personIndices[i]][genoindex+1] = phred12;
	  genotypes[personIndices[i]][genoindex+2] = phred22;
	}
       ++markerindex;
     }
     
     delete pVcf;
     //delete pMarker;
   }
   catch ( VcfFileException e ) {
     error ( e.what() );
   }
}

int main(int argc, char ** argv)
   {
   String shotgunfile, mapfile, outfile("mach1.out");
   String crossFile, errorFile;

   double errorRate = 0.01;
   int seed = 123456, warmup = 0, states = 0;
   int burnin = 0, rounds = 0, polling = 0, samples = 0;
   bool weighted = false, compact = false;
   bool mle = false, mledetails = false, uncompressed = false;
   bool inputPhased = false;
   bool phaseByRef = false;
   bool randomPhase = true;
   bool noWeight = true;
   bool likeWeight = false;
   bool matchWeight = false;

   SetupCrashHandlers();
   SetCrashExplanation("reading command line options");

   printf("Thunder_Glf 1.0.9 -- Markov Chain Haplotyping for Shotgun Sequence Data\n"
          "(c) 2005-2007 Goncalo Abecasis, Yun Li, with thanks to Paul Scheet\n\n");

   ParameterList pl;

BEGIN_LONG_PARAMETERS(longParameters)
   LONG_PARAMETER_GROUP("Shotgun Sequences")
      LONG_STRINGPARAMETER("shotgun", &shotgunfile)
   LONG_PARAMETER_GROUP("Optional Files")
      LONG_STRINGPARAMETER("crossoverMap", &crossFile)
      LONG_STRINGPARAMETER("errorMap", &errorFile)
      LONG_STRINGPARAMETER("physicalMap", &mapfile)
   LONG_PARAMETER_GROUP("Markov Sampler")
      LONG_INTPARAMETER("seed", &seed)
      LONG_INTPARAMETER("burnin", &burnin)
      LONG_INTPARAMETER("rounds", &rounds)
   LONG_PARAMETER_GROUP("Haplotyper")
      LONG_INTPARAMETER("states", &states)
      LONG_DOUBLEPARAMETER("errorRate", &errorRate)
      LONG_PARAMETER("weighted", &weighted)
      LONG_PARAMETER("compact", &compact)
  LONG_PARAMETER_GROUP("Phasing")
      EXCLUSIVE_PARAMETER("randomPhase", &randomPhase)
      EXCLUSIVE_PARAMETER("inputPhased", &inputPhased)
      EXCLUSIVE_PARAMETER("refPhased",  &phaseByRef)
  LONG_PARAMETER_GROUP("Weighting")
      EXCLUSIVE_PARAMETER("noWeight", &noWeight)
      EXCLUSIVE_PARAMETER("likeWeight", &likeWeight)
      EXCLUSIVE_PARAMETER("matchWeight",  &matchWeight)
   LONG_PARAMETER_GROUP("Imputation")
      LONG_PARAMETER("geno", &OutputManager::outputGenotypes)
      LONG_PARAMETER("quality", &OutputManager::outputQuality)
      LONG_PARAMETER("dosage", &OutputManager::outputDosage)
      LONG_PARAMETER("probs", &OutputManager::outputProbabilities)
      LONG_PARAMETER("mle", &mle)
   LONG_PARAMETER_GROUP("Output Files")
      LONG_STRINGPARAMETER("prefix", &outfile)
      LONG_PARAMETER("phase", &OutputManager::outputHaplotypes)
      LONG_PARAMETER("uncompressed", &OutputManager::uncompressed)
      LONG_PARAMETER("mldetails", &mledetails)
   LONG_PARAMETER_GROUP("Interim Output")
      LONG_INTPARAMETER("sampleInterval", &samples)
      LONG_INTPARAMETER("interimInterval", &polling)
END_LONG_PARAMETERS();

   pl.Add(new LongParameters("Available Options", longParameters));

   pl.Add(new HiddenString('m', "Map File", mapfile));
   pl.Add(new HiddenString('o', "Output File", outfile));
   pl.Add(new HiddenInteger('r', "Haplotyping Rounds", rounds));
   pl.Add(new HiddenDouble('e', "Error Rate", errorRate));

   pl.Read(argc, argv);
   pl.Status();

   if ( OutputManager::outputDosage  == false ) { // hmkang 
     error("--dosage flag must be set in this implementation");
   }

   // Setup random seed ...
   globalRandom.Reset(seed);

   SetCrashExplanation("loading information on polymorphic sites");

   // Setup and load a list of polymorphic sites, each with two allele labels ...
   Pedigree ped;

   SetCrashExplanation("loading shotgun data - first pass");

   LoadShotgunSamples(ped, shotgunfile);
   
   LoadPolymorphicSites(shotgunfile);
   
   SetCrashExplanation("loading map information for polymorphic sites");

   printf("Loaded information on %d polymorphic sites\n\n", Pedigree::markerCount);

   Pedigree::LoadMarkerMap(mapfile);

   // Check if physical map is available
   bool positionsAvailable = true;

   for (int i = 0; i < ped.markerCount; i++)
      if (Pedigree::GetMarkerInfo(i)->chromosome < 0)
         {
         positionsAvailable = false;
         break;
         }

   if (positionsAvailable)
      {
      printf("    Physical map will be used to improve crossover rate estimates.\n");

      for (int i = 1; i < ped.markerCount; i++)
         if (ped.GetMarkerInfo(i)->position <= ped.GetMarkerInfo(i - 1)->position ||
             ped.GetMarkerInfo(i)->chromosome != ped.GetMarkerInfo(i - 1)->chromosome)
            {
            printf("    FATAL ERROR -- Problems with physical map ...\n\n"
                   "    Before continuing, check the following:\n"
                   "    * All markers are on the same chromosome\n"
                   "    * All marker positions are unique\n"
                   "    * Markers in pedigree and haplotype files are ordered by physical position\n\n");
            return -1;
            }
      }

   printf("\n");
   
   printf("Processing input files and allocating memory for haplotyping\n");

   SetCrashExplanation("allocating memory for haplotype engine and consensus builder");

   ShotgunHaplotyper engine;

   engine.economyMode = compact;

   engine.EstimateMemoryInfo(ped.count, ped.markerCount, states, compact, false);
   engine.AllocateMemory(ped.count, states, ped.markerCount);

   // Copy genotypes into haplotyping engine
   if (engine.readyForUse)
      LoadShotgunResults(ped, engine.genotypes, engine.refalleles, engine.freq1s, shotgunfile);

   // printf("Done loading shotgun file\n\n");
   
   // Copy phased haplotypes into haplotyping engine
   engine.phased = 0;

   if (engine.readyForUse == false || engine.ForceMemoryAllocation() == false)
      return MemoryAllocationFailure();

   if (positionsAvailable && engine.AllocateDistances())
      {
      for (int i = 1; i < ped.markerCount; i++)
         engine.distances[i - 1] = ped.GetMarkerInfo(i)->position -
                                   ped.GetMarkerInfo(i-1)->position;
      }

   engine.ShowMemoryInfo();

   if (mle)
      {
      engine.ShowMLEMemoryInfo();
      if (!engine.AllocateMLEMemory())
         return MemoryAllocationFailure();
      }

   ConsensusBuilder::EstimateMemoryInfo(rounds - burnin, ped.count * 2, ped.markerCount);
   ConsensusBuilder consensus(rounds - burnin, ped.count * 2, ped.markerCount);

   if (consensus.readyForUse == false)
      return MemoryAllocationFailure();

   DosageCalculator::storeDistribution = OutputManager::outputDosage ||
                                         OutputManager::outputQuality ||
                                         OutputManager::outputGenotypes;

   DosageCalculator::EstimateMemoryInfo(rounds - burnin, ped.count, ped.markerCount);
   DosageCalculator doses(rounds - burnin, ped.count, ped.markerCount);

   if (doses.readyForUse == false)
      return MemoryAllocationFailure();

   /*
   if (weighted) {
     engine.weightByMismatch = true;
     printf("Using weighting approach\n\n");
     }*/
   if ( likeWeight ) {
     engine.weightByLikelihood = true;
   }
   else if ( matchWeight ) {
     engine.weightByLongestMatch = true;
   }
   //   engine.CalculateWeights();

   printf("Memory allocated successfully\n\n");

   SetCrashExplanation("loading error rate and cross over maps");

   engine.SetErrorRate(errorRate);

   bool newline = engine.LoadCrossoverRates(crossFile);
   newline |= engine.LoadErrorRates(errorFile);
   if (newline) printf("\n");

   SetCrashExplanation("searching for initial haplotype set");

   if ( inputPhased ) {
     printf("Loading phased information from the input VCF file\n\n");
     engine.LoadHaplotypesFromVCF(shotgunfile);
   }
   else if ( phaseByRef ) {
     printf("Assigning haplotypes based on reference genome\n\n");
     engine.PhaseByReferenceSetup();
   }
   else {
     printf("Assigning random set of haplotypes\n\n");
     engine.RandomSetup();
   }
   printf("Found initial haplotype set\n\n");
//OutputManager::WriteHaplotypes(outfile, ped, engine.haplotypes);
//return 0;

   SetCrashExplanation("revving up haplotyping engine");

   SetCrashExplanation("interating through markov chain haplotyping procedure");

   for (int i = 0; i < rounds; i++)
      {
      engine.LoopThroughChromosomes();

      engine.UpdateThetas();
      errorRate = engine.UpdateErrorRate();

      printf("Markov Chain iteration %d [%d mosaic crossovers]\n",
             i + 1, engine.TotalCrossovers() );

      if (i < burnin)
         continue;

      if (OutputManager::outputHaplotypes)
         consensus.Store(engine.haplotypes);

      if (doses.storeDosage || doses.storeDistribution)
         doses.Update(engine.haplotypes);

      UpdateVector(engine.thetas, thetas, nthetas, engine.markers - 1);
      UpdateErrorRates(engine.error_models, error_rates, nerror_rates, engine.markers);

      if (polling > 0 && ((i - burnin) % polling) == 0) {
	OutputVCFConsensus(shotgunfile, ped, consensus, doses, outfile + ".prelim" + (i+1) + ".vcf.gz", thetas, error_rates);
      	//OutputManager::OutputConsensus(ped, consensus, doses, outfile + ".prelim" + (i + 1));
      }

      if (samples > 0 && ((i - burnin) % samples) == 0)
         OutputManager::WriteHaplotypes(outfile + ".sample" + (i + 1) + ".gz",  ped, engine.haplotypes);

      }

   if (rounds) printf("\n");

   SetCrashExplanation("estimating maximum likelihood solution, conditional on current state");

   if (mle)
      {
      // Use best available error and crossover rates for MLE
      if (nerror_rates)
         for (int i = 0; i < engine.markers; i++)
            engine.SetErrorRate(i, error_rates[i] / nerror_rates);

      if (nthetas)
         for (int i = 0; i < engine.markers - 1; i++)
            engine.thetas[i] = thetas[i] / nthetas;

      engine.OutputMLEs(ped, outfile, mledetails);
      }

//   ParseHaplotypes(engine.haplotypes, engine.individuals * 2 - 2, engine.markers, 32);

   SetCrashExplanation("outputing solution");

   // If we did multiple rounds of haplotyping, then generate consensus
   if (rounds > 1)
     OutputVCFConsensus(shotgunfile, ped, consensus, doses, outfile + ".vcf.gz", thetas, error_rates);
     //OutputManager::OutputConsensus(ped, consensus, doses, outfile);
   else
      OutputManager::WriteHaplotypes(outfile, ped, engine.haplotypes);

   printf("Estimated mismatch rate in Markov model is: %.5f\n", errorRate);
   }




 
