////////////////////////////////////////////////////////////////////// 
// mach1/MergeHaplotypes.cpp 
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
 
#include "MergeHaplotypes.h"
#include "MemoryAllocators.h"
#include "OutputHandlers.h"
#include "MemoryInfo.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef CHAR_BIT
#define    SEVEN     (CHAR_BIT - 1)
#define    EIGHT     (CHAR_BIT)
#else
#define    SEVEN      7
#define    EIGHT      8
#endif

ConsensusBuilder::ConsensusBuilder(int N, int H, int M)
   {
   if (OutputManager::outputHaplotypes == false)
      {
      stored = samples = haplotypes = markers = 0;

      sampledHaplotypes = NULL;
      consensus = NULL;

      readyForUse = true;
      
      return;
      }

   stored = 0;

   samples = N;
   haplotypes = H;
   markers = M;

   sampledHaplotypes = AllocateCharCube(samples, haplotypes, (markers + SEVEN) / EIGHT );
   consensus = AllocateCharMatrix(haplotypes, markers);

   readyForUse = consensus != NULL && sampledHaplotypes != NULL;
   }

ConsensusBuilder::~ConsensusBuilder()
   {


   if (sampledHaplotypes != NULL)
      FreeCharCube(sampledHaplotypes, samples, haplotypes);

   if (consensus != NULL)
      FreeCharMatrix(consensus, haplotypes);
   }

void ConsensusBuilder::Store(char ** newHaplotypes)
   {
   if (sampledHaplotypes == 0) return;

   for (int i = 0; i < haplotypes; i++)
      {
      int byte = 0, mask = 1;

      for (int j = 0; j < markers; j++)
         {
         if (newHaplotypes[i][j])
            sampledHaplotypes[stored][i][byte] |= mask;
         else
            sampledHaplotypes[stored][i][byte] &= ~mask;

         mask = (mask == (1 << SEVEN)) ? (byte++, 1) : mask * 2;
         }
      }

   stored++;
   }

void ConsensusBuilder::Merge()
   {
   // Don't try to build consensus with no data
   if (stored == 0) return;

   // Initialize haplotype quality scores
   flips = errors = 0;

   // The phase for each pair of haplotypes indicates their ordering
   // in relation to the consensus
   char * phase = new char [stored];

   // Loop through each haplotype in the set
   for (int h = 0; h < haplotypes; h += 2)
      {
      // Select phase based on the first heterozygous position for each haplotype
      for (int i = 0; i < stored; i++)
         {
         phase[i] = 0;

         for (int j = 0, byte = 0, mask = 1; j < markers; j++)
            {
            if ((sampledHaplotypes[i][h][byte] ^ sampledHaplotypes[i][h + 1][byte]) & mask)
               {
               phase[i] = (sampledHaplotypes[i][h][byte] & mask) != 0;
               break;
               }

            mask = (mask == (1 << SEVEN)) ? (byte++, 1) : mask * 2;
            }
         }

      // Build consensus one position at a time ...
      for (int i = 0, byte = 0, mask = 1; i < markers; i++)
         {
         int counts[4] = {0, 0, 0, 0};

         // Count the number of occurences for each genotype
         for (int j = 0; j < stored; j++)
            {
            int allele1 = (sampledHaplotypes[j][h + phase[j]][byte] & mask) != 0;
            int allele2 = (sampledHaplotypes[j][h + (phase[j] ^ 1)][byte] & mask) != 0;

            counts[allele1 * 2 + allele2]++;
            }

         // Record the expect number of copies of the common allele
         // dosage[h / 2][i] = (short)((counts[1] + counts[2] + 2 * counts[3]) * scale + 0.5);

         // Select the most likely genotype
         int best = 0;

         for (int j = 1; j < 4; j++)
            if (counts[j] > counts[best])
               best = j;

         // Assign it to the consensus
         consensus[h][i] = best / 2;
         consensus[h + 1][i] = best % 2;

         // Count number of samples with an alternative solutions
         int alternative_solution_weight = (best == 0 || best == 3) ?
             (stored - counts[best]) : (counts[0] + counts[3]);

         // Update estimated haplotype quality scores
         errors += alternative_solution_weight;

         // If a heterozygous genotype was selected, update the phase for other informative
         // haplotypes
         if (best != 0 && best != 3 && counts[best ^ 3] > 0)
            {
            // Update estimated flip scores
            flips += counts[best ^ 3];

            best = best ^ 3;

            for (int j = 0; j < stored; j++)
               {
               int allele1 = (sampledHaplotypes[j][h + phase[j]][byte] & mask) != 0;
               int allele2 = (sampledHaplotypes[j][h + (phase[j] ^ 1)][byte] & mask) != 0;

               if ((allele1 * 2 + allele2) == best)
                  phase[j] = phase[j] ^ 1;
               }
            }

         mask = (mask == (1 << SEVEN)) ? (byte++, 1) : mask * 2;
         }
      }

   delete [] phase;

   flips /= stored;
   errors /= stored;
   }

void ConsensusBuilder::ShowMemoryInfo()
   {
   EstimateMemoryInfo(samples, haplotypes, markers);
   }

void ConsensusBuilder::EstimateMemoryInfo(int Samples, int Haplotypes, int Markers)
   {
   if (OutputManager::outputHaplotypes == false)
      return;

   double bytes = sizeof(char) * (double) Samples * Haplotypes * (Markers + SEVEN) / EIGHT +
                  sizeof(char) * (double) Haplotypes * Markers;

   printf("   %40s %s\n", "Consensus builder ...", (const char *) MemoryInfo(bytes));
   }
 
