////////////////////////////////////////////////////////////////////// 
// mach1/Main.cpp 
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
 
#include "AssociationAnalysis.h"
#include "CostCalculator.h"
#include "OutputHandlers.h"
#include "MergeHaplotypes.h"
#include "HaplotypeLoader.h"
#include "Parameters.h"
#include "Manners.h"
#include <ctime>

float * thetas = NULL;
int     nthetas = 0;

float * error_rates = NULL;
int     nerror_rates = 0;

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

void EvaluateHaplotypes(DosageCalculator & doses, Pedigree & ped, char ** genotypes)
   {
   int matches = 0;
   int mismatches = 0;
   int partialmatches = 0;

   for (int i = 0; i < ped.count; i++)
      for (int j = 0; j < ped.markerCount; j++)
         if (genotypes[i][j] == 0 && ped[i].markers[j].isKnown())
            {
            int actual = ped[i].markers[j].SequenceCoded() - 1;
            int imputed = doses.GetBestGenotype(i, j);

            if (actual == imputed)
               matches++;
            else if (actual == 1 || imputed == 1)
               partialmatches++;
            else
               mismatches++;
            }

   double total = matches + partialmatches + mismatches + 1e-30;

   printf("Comparing %.0f masked genotypes with estimates ...\n", total);
   printf("   Estimated per genotype error rate is  %.4f\n",
          (mismatches + partialmatches) / total);
   printf("   Estimated per allele error rate is %.4f\n\n",
          (mismatches + partialmatches * 0.5) / total);
   }

int MemoryAllocationFailure()
   {
   printf("FATAL ERROR - Memory allocation failed\n");
   return -1;
   }

int main(int argc, char ** argv)
   {
   clock_t startt, endt;
   startt = clock();
   String datfile, pedfile, hapfile, mapfile, hapsnps, outfile("mach1.out");
   String crossFile, errorFile;

   double errorRate = 0.001, mask = 0.0;
   int seed = 123456, warmup = 0, states = 0;
   int burnin = 0, rounds = 0, polling = 0, samples = 0;
   bool weighted = false, greedy = false, compact = false;
   bool association = false, quickNPL = false, mle = false, mledetails = false;

   SetupCrashHandlers();
   SetCrashExplanation("reading command line options");

#ifndef VERSION
   printf("Mach 1.0 -- Markov Chain Haplotyping\n"
#else
   printf("Mach " VERSION " -- Markov Chain Haplotyping\n"
#endif
          "(c) 2005-2007 Goncalo Abecasis, with thanks to Yun Li, Paul Scheet\n\n");

   ParameterList pl;

BEGIN_LONG_PARAMETERS(longParameters)
   LONG_PARAMETER_GROUP("Input Files")
      LONG_STRINGPARAMETER("datfile", &datfile)
      LONG_STRINGPARAMETER("pedfile", &pedfile)
      LONG_DOUBLEPARAMETER("mask", &mask)
   LONG_PARAMETER_GROUP("Optional Files")
      LONG_STRINGPARAMETER("crossoverMap", &crossFile)
      LONG_STRINGPARAMETER("errorMap", &errorFile)
      LONG_STRINGPARAMETER("physicalMap", &mapfile)
   LONG_PARAMETER_GROUP("Phased Data")
      LONG_STRINGPARAMETER("snps", &hapsnps)
      LONG_STRINGPARAMETER("haps", &hapfile)
      LONG_PARAMETER("hapmapFormat", &HaplotypeLoader::hapmapFormat)
      LONG_PARAMETER("autoFlip", &HaplotypeLoader::autoFlip)
      LONG_PARAMETER("greedy", &greedy)
   LONG_PARAMETER_GROUP("Markov Sampler")
      LONG_INTPARAMETER("seed", &seed)
      LONG_INTPARAMETER("burnin", &burnin)
      LONG_INTPARAMETER("rounds", &rounds)
   LONG_PARAMETER_GROUP("Mapping Options")
      LONG_PARAMETER("npl", &quickNPL)
      LONG_PARAMETER("association", &association)
   LONG_PARAMETER_GROUP("Haplotyper")
      LONG_INTPARAMETER("states", &states)
      LONG_DOUBLEPARAMETER("errorRate", &errorRate)
      LONG_PARAMETER("weighted", &weighted)
      LONG_PARAMETER("compact", &compact)
   LONG_PARAMETER_GROUP("Imputation")
      LONG_PARAMETER("geno", &OutputManager::outputGenotypes)
      LONG_PARAMETER("quality", &OutputManager::outputQuality)
      LONG_PARAMETER("dosage", &OutputManager::outputDosage)
      LONG_PARAMETER("probs", &OutputManager::outputProbabilities)
      LONG_PARAMETER("mle", &mle)
   LONG_PARAMETER_GROUP("Output Files")
      LONG_STRINGPARAMETER("prefix", &outfile)
      LONG_PARAMETER("phase", &OutputManager::outputHaplotypes)
      LONG_PARAMETER("mldetails", &mledetails)
   LONG_PARAMETER_GROUP("Interim Output")
      LONG_INTPARAMETER("sampleInterval", &samples)
      LONG_INTPARAMETER("interimInterval", &polling)
   BEGIN_LEGACY_PARAMETERS()
      LONG_STRINGPARAMETER("legend", &hapsnps)
END_LONG_PARAMETERS();

   pl.Add(new LongParameters("Available Options", longParameters));

   pl.Add(new HiddenString('d', "Data File", datfile));
   pl.Add(new HiddenString('p', "Pedigree File", pedfile));
   pl.Add(new HiddenString('m', "Map File", mapfile));
   pl.Add(new HiddenString('h', "Haplotype File", hapfile));
   pl.Add(new HiddenString('s', "Haplotype SNP File", hapsnps));
   pl.Add(new HiddenString('o', "Output File", outfile));
   pl.Add(new HiddenInteger('r', "Haplotyping Rounds", rounds));
   pl.Add(new HiddenDouble('e', "Error Rate", errorRate));

//   pl.Add(new IntParameter('w', "Warm-Up Sample", warmup));
//   pl.Add(new HiddenIntParameter('s', "Random Seed", seed));

   pl.Read(argc, argv);
   pl.Status();

   // Setup random seed ...
   globalRandom.Reset(seed);

   SetCrashExplanation("loading phased chromosomes");

   // Load phased haplotypes, if available
   HaplotypeLoader haps;

   haps.LoadMarkerList(hapsnps);
   haps.LoadHaplotypes(hapfile);

   if (haps.count || haps.markerCount)
      printf("\n");

   if (greedy)
      {
      if (haps.count == 0)
         {
         printf("   GREEDY SOLUTION NOT AVAILABLE. Although you requested a greedy\n"
                "   solution, no phased haplotypes were provided as input.\n\n");
         greedy = false;
         }
      else
         {
         printf("   GREEDY SOLUTION. Phased haplotypes will be used to resolve ambiguous\n"
                "   individuals and generate a greedy solution.\n\n");
         //states = haps.count;
         // default --greedy (using all reference haplotypes
         // 1.0.16.b allows --greedy together with --states SS (SS <= haps.count)
         if (states == 0 || states > haps.count) {states = haps.count;}
         }
      }

   SetCrashExplanation("loading pedigree file");

   // Setup and load pedigree ...
   Pedigree ped;

   ped.Prepare(datfile);
   ped.Load(pedfile);
   ped.LoadMarkerMap(mapfile);

   printf("Loaded pedigree with:\n"
          "    %d individuals to be haplotyped at %d markers\n",
          ped.count, ped.markerCount);

   if (mask > 0.0)
      printf("    %.0f%% of genotypes will be masked prior to haplotyping\n", mask * 100.0);

   // Check if physical map is available
   bool positionsAvailable = true;

   for (int i = 0; i < ped.markerCount; i++)
      if (ped.GetMarkerInfo(i)->chromosome < 0)
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

   // Check that haplotypes and pedigree are consistent
   haps.ConsistencyCheck(ped);

   printf("Formating genotypes and allocating memory for haplotyping\n");
   ped.ShowMemoryInfo();
   haps.ShowMemoryInfo();

   SetCrashExplanation("allocating memory for haplotype engine and consensus builder");

   Haplotyper engine;

   engine.economyMode = compact;

   engine.EstimateMemoryInfo(ped.count + haps.count / 2, ped.markerCount, states, compact, haps.count != 0);
   engine.AllocateMemory(ped.count + haps.count / 2, states, ped.markerCount);

   // Copy genotypes into haplotyping engine
   if (engine.readyForUse)
      for (int i = 0; i < ped.count; i++)
         for (int j = 0; j < ped.markerCount; j++)
            {
            if (mask == 0.0 || globalRandom.Next() > mask)
               engine.genotypes[i][j] = ped[i].markers[j].SequenceCoded();
            else
               engine.genotypes[i][j] = 0;
            }

   // Verify that no more than two alleles were present for each marker
   if (engine.readyForUse)
      {
      StringHash badMarkers;

      for (int i = 0; i < ped.count; i++)
         for (int j = 0; j < ped.markerCount; j++)
            if (engine.genotypes[i][j] > 3)
               badMarkers.Add(ped.markerNames[j]);

      if (badMarkers.Entries() > 0)
         {
         printf("\n\nFATAL ERROR:\n"
                "This version of MaCH is designed for bi-allelic markers\n"
                "However, the following marker(s) have >2 alleles:\n   ");

         int togo = badMarkers.Entries();
         for (int i = 0, new_line = 3, lines = 0; i < badMarkers.Capacity(); i++)
            if (badMarkers.SlotInUse(i))
               {
               if (new_line + badMarkers[i].Length() > 78)
                  printf("\n   "), new_line = 3, lines++;

               if (lines > 10 && togo > 5) break;

               printf("%s ", (const char *) badMarkers[i]);
               new_line += badMarkers[i].Length();
               togo--;
               }

         if (togo) printf("\n%d additional markers not listed", togo);

         printf("\n\nPlease remove or recode markers with more than 2 alleles\n\n");
         return(-1);
         }
      }

   // Copy phased haplotypes into haplotyping engine
   engine.phased = haps.count / 2;

   if (engine.readyForUse)
      for (int i = 0; i < (haps.count & ~1); i++)
         for (int j = 0; j < ped.markerCount; j++)
            engine.haplotypes[ped.count * 2 + i][j] = haps.haplotypes[i][j] - 1;

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

   if (quickNPL && Pedigree::affectionCount)
      {
      engine.EstimateDiseaseMemoryInfo(ped.count, ped.markerCount, ped.affectionCount);

      if (engine.AllocateDiseaseStatus(Pedigree::affectionCount))
         {
         for (int i = 0; i < ped.count; i++)
            for (int j = 0; j < ped.affectionCount; j++)
               engine.diseaseStatus[i][j] = ped[i].affections[j];
         }
      else
         return MemoryAllocationFailure();
      }
   else
      quickNPL = false;

   ConsensusBuilder::EstimateMemoryInfo(rounds - burnin, ped.count * 2, ped.markerCount);
   ConsensusBuilder consensus(rounds - burnin, ped.count * 2, ped.markerCount);

   if (consensus.readyForUse == false)
      return MemoryAllocationFailure();

   DosageCalculator::storeDistribution = OutputManager::outputDosage ||
                                         OutputManager::outputQuality ||
                                         OutputManager::outputGenotypes ||
					 OutputManager::outputProbabilities ||
                                         mask > 0.0;

   DosageCalculator::EstimateMemoryInfo(rounds - burnin, ped.count, ped.markerCount);
   DosageCalculator doses(rounds - burnin, ped.count, ped.markerCount);

   if (doses.readyForUse == false)
      return MemoryAllocationFailure();

   if (weighted)
      engine.CalculateWeights();

   printf("Memory allocated successfully\n\n");

   SetCrashExplanation("loading error rate and cross over maps");

   engine.SetErrorRate(errorRate);

   bool newline = engine.LoadCrossoverRates(crossFile);
   newline |= engine.LoadErrorRates(errorFile);
   if (newline) printf("\n");

   SetCrashExplanation("searching for initial haplotype set");

   engine.greedy = greedy;
   engine.RandomSetup();
   printf("Found initial haplotype set\n\n");

   SetCrashExplanation("revving up haplotyping engine");

   if (warmup)
      {
      engine.WarmUp(warmup, 5);

      printf("Warmed up haplotyping engine ...\n\n");
      }

//   ParseHaplotypes(engine.haplotypes, engine.individuals * 2 - 2, engine.markers, 32);

// The cost calculator uses heurestics to try and find faster haplotyping
// strategies -- however, these are not yet implemented!
// CostCalculator blueSky;

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

      // blueSky.OptimizeCost(engine.haplotypes, engine.individuals * 2, engine.markers);

      if (polling > 0 && ((i - burnin) % polling) == 0)
      {
         OutputManager::OutputConsensus(ped, consensus, doses, outfile + ".prelim" + (i + 1));

   FILE * file = fopen(outfile + ".prelim" + (i + 1) + ".rec", "wt");        

   if (file == NULL)
      printf("Error opening output file [%s.prelim.%d.rec]\n", (const char *) outfile, i+1);
   else
      {
      fprintf(file, "Interval AvgRate LastRate\n");
      for (int j = 0; j < engine.markers - 1; j++)
            fprintf(file, "%d-%d %.4f %.4f\n", j + 1, j + 2,
                          nthetas ? thetas[j] / nthetas : engine.thetas[j],
                          engine.thetas[j]);
      fclose(file);

      printf("Wrote out file [%s.prelim.%d.rec] with mosaic crossover rates ...\n", (const char *) outfile, i+1);
      }    

   file = fopen(outfile +  ".prelim" + (i + 1) + ".erate", "wt");

   if (file == NULL)
      printf("Error opening output file [%s.prelim.%d.erate]\n", (const char *) outfile, i+1);
   else
      {           
      fprintf(file, "Marker AvgRate LastRate\n");    
      for (int j = 0; j < engine.markers; j++)     
         fprintf(file, "%s %.4f %.4f\n", (const char *) ped.markerNames[j],
                       nerror_rates ? error_rates[j] / nerror_rates : engine.GetErrorRate(j),
                       engine.GetErrorRate(j));
      fclose(file);

      printf("Wrote out file [%s.prelim.%d.erate] with per marker error rates ...\n\n",
             (const char *) outfile, i+1);
      }


      }

      if (samples > 0 && ((i - burnin) % samples) == 0)
         OutputManager::WriteHaplotypes(outfile + ".sample" + (i + 1), ped, engine.haplotypes);

      UpdateVector(engine.thetas, thetas, nthetas, engine.markers - 1);
      UpdateErrorRates(engine.error_models, error_rates, nerror_rates, engine.markers);

      engine.updateDiseaseScores = quickNPL;
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
      OutputManager::OutputConsensus(ped, consensus, doses, outfile);
   else
      OutputManager::WriteHaplotypes(outfile, ped, engine.haplotypes);

   if (doses.storeDosage || doses.storeDistribution)
      doses.OutputMarkerInfo(outfile + ".info");

   FILE * file = fopen(outfile + ".rec", "wt");

   if (file == NULL)
      printf("Error opening output file [%s.rec]\n", (const char *) outfile);
   else
      {
      fprintf(file, "Interval AvgRate LastRate\n");
      for (int i = 0; i < engine.markers - 1; i++)
            fprintf(file, "%d-%d %.4f %.4f\n", i + 1, i + 2,
                          nthetas ? thetas[i] / nthetas : engine.thetas[i],
                          engine.thetas[i]);
      fclose(file);

      if (thetas != NULL) delete [] thetas;

      printf("Wrote out file [%s.rec] with mosaic crossover rates ...\n", (const char *) outfile);
      }

   file = fopen(outfile + ".erate", "wt");

   if (file == NULL)
      printf("Error opening output file [%s.erate]\n", (const char *) outfile);
   else
      {
      fprintf(file, "Marker AvgRate LastRate\n");
      for (int i = 0; i < engine.markers; i++)
         fprintf(file, "%s %.4f %.4f\n", (const char *) ped.markerNames[i],
                       nerror_rates ? error_rates[i] / nerror_rates : engine.GetErrorRate(i),
                       engine.GetErrorRate(i));
      fclose(file);

      if (error_rates != NULL) delete [] error_rates;

      printf("Wrote out file [%s.erate] with per marker error rates ...\n\n",
             (const char *) outfile);
      }

   if (quickNPL && rounds)
      AssociationAnalysis::ScoreNPL(outfile, ped, engine, rounds - burnin);

   if (association && rounds)
      AssociationAnalysis::ScoreMarkers(outfile, ped, doses);

   if (mask > 0.0 && rounds)
      EvaluateHaplotypes(doses, ped, engine.genotypes);

   printf("Estimated mismatch rate in Markov model is: %.5f\n\n", errorRate);

   endt = clock();
   int lapsetime = (int) ((double)(endt - startt) / CLOCKS_PER_SEC);
   printf("Analysis took %d seconds\n\n", lapsetime);

   }


 
