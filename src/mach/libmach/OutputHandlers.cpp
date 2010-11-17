////////////////////////////////////////////////////////////////////// 
// mach1/OutputHandlers.cpp 
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
 
#include "OutputHandlers.h"

bool OutputManager::outputHaplotypes = false;
bool OutputManager::outputGenotypes = false;
bool OutputManager::outputDosage = false;
bool OutputManager::outputProbabilities = false;
bool OutputManager::outputQuality = false;

void OutputManager::WriteHaplotypes(const char * outfile, Pedigree & ped, char ** haplotypes)
   {
   if (!OutputManager::outputHaplotypes)
      return;

   IFILE file = ifopen(outfile, "wt");

   if (file == NULL)
      {
      printf("Error opening output file [%s]\n", (const char *) outfile);
      exit(1);
      }

   for (int i = 0; i < ped.count * 2; i++)
      {
      ifprintf(file, "%s->%s HAPLO%d ",
                    (const char *) ped[i / 2].famid,
                    (const char *) ped[i / 2].pid, i % 2 + 1);

      for (int j = 0; j < ped.markerCount; j++)
         if (ped[i/2].markers[j].isKnown())
            ifprintf(file, "%s", (const char *) ped.GetMarkerInfo(j)->GetAlleleLabel(haplotypes[i][j] + 1).ToUpper());
         else
            ifprintf(file, "%s", (const char *) ped.GetMarkerInfo(j)->GetAlleleLabel(haplotypes[i][j] + 1).ToLower());

      ifprintf(file, "\n");
      }

   ifclose(file);

   printf("Wrote out file [%s] with phased chromosomes ...\n", (const char *) outfile);
   }

void OutputManager::OutputConsensus(Pedigree & ped, ConsensusBuilder & consensus,
                                    DosageCalculator & doses, String filename)
   {
   consensus.Merge();

   if (consensus.stored)
      printf("Merged sampled haplotypes to generate consensus\n"
             "   Consensus estimated to have %.1f errors in missing genotypes and %.1f flips in haplotypes\n\n",
             consensus.errors, consensus.flips);

   if (outputHaplotypes) WriteHaplotypes(filename + ".gz", ped, consensus.consensus);
   if (outputGenotypes)  WriteGenotypes(filename + ".geno.gz", ped, doses);
   if (outputQuality)    WriteQuality(filename + ".qc.gz", ped, doses);
   if (outputDosage)     WriteDosages(filename + ".dose.gz", ped, doses);
   if (outputProbabilities) WriteProbabilities(filename + ".prob.gz", ped, doses);
   }

void OutputManager::WriteGenotypes(const char * outfile, Pedigree & ped, DosageCalculator & doses)
   {
   IFILE file = ifopen(outfile, "wt");

   if (file == NULL)
      {
      printf("Error opening output file [%s]\n", (const char *) outfile);
      exit(1);
      }

   for (int i = 0; i < ped.count; i++)
      {
      ifprintf(file, "%s->%s GENO ", (const char *) ped[i].famid,
                                    (const char *) ped[i].pid);

      for (int j = 0; j < ped.markerCount; j++)
         {
         int best = doses.GetBestGenotype(i, j);

         MarkerInfo * info = ped.GetMarkerInfo(j);

         ifprintf(file, " %s/%s", (const char *) info->GetAlleleLabel((best + 2) / 2),
                                 (const char *) info->GetAlleleLabel((best + 3) / 2));
         }

      ifprintf(file, "\n");
      }

   ifclose(file);

   printf("Wrote out file [%s] with imputed genotypes ...\n", (const char *) outfile);
   }

void OutputManager::WriteDosages(const char * outfile, Pedigree & ped, DosageCalculator & doses)
   {
   IFILE file = ifopen(outfile, "wt");

   if (file == NULL)
      {
      printf("Error opening output file [%s]\n", (const char *) outfile);
      exit(1);
      }

   for (int i = 0; i < ped.count; i++)
      {
      ifprintf(file, "%s->%s ALLELE_DOSE ", (const char *) ped[i].famid,
                                           (const char *) ped[i].pid);

      for (int j = 0; j < ped.markerCount; j++)
         ifprintf(file, "%.3f ", doses.GetDosage(i, j));

      ifprintf(file, "\n");
      }

   ifclose(file);

   printf("Wrote out file [%s] with dosage information...\n", (const char *) outfile);
   }

void OutputManager::WriteQuality(const char * outfile, Pedigree & ped, DosageCalculator & doses)
   {
   IFILE file = ifopen(outfile, "wt");

   if (file == NULL)
      {
      printf("Error opening output file [%s]\n", (const char *) outfile);
      exit(1);
      }

   for (int i = 0; i < ped.count; i++)
      {
      ifprintf(file, "%s->%s GENO_QUALITY ", (const char *) ped[i].famid,
                                           (const char *) ped[i].pid);

      for (int j = 0; j < ped.markerCount; j++)
         ifprintf(file, "%.3f ", doses.GetQuality(i, j));

      ifprintf(file, "\n");
      }

   ifclose(file);

   printf("Wrote out file [%s] with quality scores...\n", (const char *) outfile);
   }

void OutputManager::WriteProbabilities(const char * outfile, Pedigree & ped, DosageCalculator & doses)
   {
   IFILE file = ifopen(outfile, "wt");

   if (file == NULL)
      {
      printf("Error opening output file [%s]\n", (const char *) outfile);
      exit(1);
      }

   for (int i = 0; i < ped.count; i++)
      {
      ifprintf(file, "%s->%s GENO_PROBS", (const char *) ped[i].famid,
                                          (const char *) ped[i].pid);

      unsigned int n0, n1, n2;

      for (int j = 0; j < ped.markerCount; j++)
         {
         doses.GetCounts(i, j, n0, n1, n2);

         double d = n0 + n1 + n2 + 1e-30;

         ifprintf(file, " %.3f %.3f", n0 / d, n1 / d);
         }

      ifprintf(file, "\n");
      }

   ifclose(file);

   printf("Wrote out file [%s] with genotype probabilities...\n", (const char *) outfile);
   }

 
