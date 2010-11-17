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
 
#include "ShotgunHaplotyper.h"
#include "ShotgunManners.h"
#include "OutputHandlers.h"
#include "MergeHaplotypes.h"
#include "HaplotypeLoader.h"
#include "Parameters.h"
#include "InputFile.h"
#include "Error.h"
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

int MemoryAllocationFailure()
   {
   printf("FATAL ERROR - Memory allocation failed\n");
   return -1;
   }

void LoadPolymorphicSites(const String & filename)
   {
   StringArray input, tokens;

   input.Read(filename);

   for (int i = 0; i < input.Length(); i++)
      {
      tokens.ReplaceTokens(input[i]);

      if (tokens.Length() < 3) continue;

      int markers = Pedigree::markerCount;
      int marker  = Pedigree::GetMarkerID(tokens[0]);
      int al1 = Pedigree::LoadAllele(marker, tokens[1]);
      int al2 = Pedigree::LoadAllele(marker, tokens[2]);

      if (markers != marker)
         error("Each polymorphic site should only occur once, but site %s is duplicated\n",
               (const char *) tokens[0]);

      if (al1 != 1 || al2 != 2)
         error("Allele labels '%s' and '%s' for polymorphic site '%s' are not valid\n",
               (const char *) tokens[1], (const char *) tokens[2], (const char *) tokens[0]);
      }

   if (Pedigree::markerCount == 0)
      error("No information on polymorphic sites available,\n"
            "please check you provided correct filename.\n");
   }

void LoadShotgunSamples(Pedigree & ped, const String & filename)
   {
   bool fail;

   String      input;
   StringArray tokens;

   IFILE f = ifopen(filename, "rb");

   if (f == NULL)
      error("Failed to open file with read count data,\n"
            "please check you provided correct filename.\n");

   while (!ifeof(f))
      {
      input.ReadLine(f);
      tokens.ReplaceTokens(input);

      if (tokens.Length() < 5) continue;

      ped.AddPerson(tokens[0], tokens[1], tokens[2], tokens[3], ped.TranslateSexCode(tokens[5], fail), true);
      }

   ifclose(f);

   ped.Sort();
   }

void LoadShotgunResults(Pedigree & ped, char ** genotypes, const String & filename)
   {
   bool fail;

   String      input;
   StringArray tokens;

   IFILE f = ifopen(filename, "rt");

   if (f == NULL)
      error("Failed to open file with read count data,\n"
            "please check you provided correct filename.\n");

   while (!ifeof(f))
      {
      input.ReadLine(f);
      tokens.ReplaceTokens(input);

      if (tokens.Length() == 0) continue;

      if (tokens.Length() != ped.markerCount * 2 + 5)
         error("Incorrect number of columns for line beggining:\n\n"
                "  %.70s\n\n"
                "Expecting ids for family, individual, parents and sex, followed\n"
                "by %d columns with read count summaries\n\n",
                (const char *) input, ped.markerCount * 2);

      int person = ped.FindPerson(tokens[0], tokens[1])->serial;

      for (int i = 0; i < ped.markerCount; i++)
         {
         int allele1 = tokens[5 + i * 2].AsInteger();
         int allele2 = tokens[5 + i * 2 + 1].AsInteger();

         if (allele1 >= 16) allele1 = 15;
         if (allele2 >= 16) allele2 = 15;

         genotypes[person][i] = allele2 * 16 + allele1;
         }
      }
   ifclose(f);
   }



int main(int argc, char ** argv)
   {
   String polymorphicSites, readCounts, mapfile, outfile("mach1.out");
   String crossFile, errorFile;

   double sequencingError = 0.005, errorRate = 0.01;
   int seed = 123456, warmup = 0, states = 0;
   int burnin = 0, rounds = 0, polling = 0, samples = 0;
   bool weighted = false, compact = false;
   bool mle = false, mledetails = false;

   SetupCrashHandlers();
   SetCrashExplanation("reading command line options");

   printf("Thunder 1.0.16.a -- Markov Chain Haplotyping for Shotgun Sequence Data\n"
          "(c) 2005-2007 Goncalo Abecasis, with thanks to Yun Li, Paul Scheet\n\n");

   ParameterList pl;
   clock_t startt, endt;
   startt = clock();

BEGIN_LONG_PARAMETERS(longParameters)
   LONG_PARAMETER_GROUP("Shotgun Sequences")
      LONG_STRINGPARAMETER("polymorphicSites", &polymorphicSites)
      LONG_STRINGPARAMETER("readCounts", &readCounts)
      LONG_DOUBLEPARAMETER("seqError", &sequencingError)
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
   LONG_PARAMETER_GROUP("Imputation")
      LONG_PARAMETER("geno", &OutputManager::outputGenotypes)
      LONG_PARAMETER("quality", &OutputManager::outputQuality)
      LONG_PARAMETER("dosage", &OutputManager::outputDosage)
      LONG_PARAMETER("mle", &mle)
   LONG_PARAMETER_GROUP("Output Files")
      LONG_STRINGPARAMETER("prefix", &outfile)
      LONG_PARAMETER("phase", &OutputManager::outputHaplotypes)
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

   // Setup random seed ...
   globalRandom.Reset(seed);

   SetCrashExplanation("loading information on polymorphic sites");

   // Setup and load a list of polymorphic sites, each with two allele labels ...
   Pedigree ped;

   SetCrashExplanation("loading shotgun data - first pass");

   LoadShotgunSamples(ped, readCounts);

   LoadPolymorphicSites(polymorphicSites);

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

   engine.EstimateMemoryInfo(ped.count, ped.markerCount, states, compact);
   engine.AllocateMemory(ped.count, states, ped.markerCount);

   // Copy genotypes into haplotyping engine
   if (engine.readyForUse)
      LoadShotgunResults(ped, engine.genotypes, readCounts);

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

   if (weighted)
      engine.CalculateWeights();

   printf("Memory allocated successfully\n\n");

   SetCrashExplanation("loading error rate and cross over maps");

   engine.SetErrorRate(errorRate);
   engine.SetShotgunError(sequencingError);

   bool newline = engine.LoadCrossoverRates(crossFile);
   newline |= engine.LoadErrorRates(errorFile);
   if (newline) printf("\n");

   SetCrashExplanation("searching for initial haplotype set");

   engine.RandomSetup();
   printf("Found initial haplotype set\n\n");

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

      if (polling > 0 && ((i - burnin) % polling) == 0)
         OutputManager::OutputConsensus(ped, consensus, doses, outfile + ".prelim" + (i + 1));

      if (samples > 0 && ((i - burnin) % samples) == 0)
         OutputManager::WriteHaplotypes(outfile + ".sample" + (i + 1), ped, engine.haplotypes);

      UpdateVector(engine.thetas, thetas, nthetas, engine.markers - 1);
      UpdateErrorRates(engine.error_models, error_rates, nerror_rates, engine.markers);
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

   printf("Estimated mismatch rate in Markov model is: %.5f\n", errorRate);
   endt = clock();
   int lapsetime = (int) ((double)(endt - startt) / CLOCKS_PER_SEC);
   printf("Analysis took %d seconds\n\n", lapsetime);

   }




 
