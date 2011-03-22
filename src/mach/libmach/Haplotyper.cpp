////////////////////////////////////////////////////////////////////// 
// mach1/Haplotyper.cpp 
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
 
#include "Haplotyper.h"
#include "StringArray.h"
#include "Pedigree.h"
#include "MemoryAllocators.h"
#include "MemoryInfo.h"
#include "Error.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef   __DOUBLE_HAPLOTYPING__
#define  float double
#endif

#ifndef min
#define min(a,b)     ((a)<(b)?(a):(b))
#endif

// Coding for individual genotypes used here is:
//
//    0 -- missing!
//    1 -- homozygous for allele 1
//    2 -- heterozygous
//    3 -- homozygous for allele 2
//
//    16+  -- ordered genotype
//    1    -- first allele is 1
//    2    -- first allele is 2
//    4    -- second allele is 1
//    8    -- second allele is 2
//
//    32+  -- allele choice depends on another individual
//

Haplotyper::Haplotyper()
   {
   phased = 0;
   individuals = 0;
   markers = 0;
   states = 0;

   genotypes = NULL;
   haplotypes = NULL;
   weights = NULL;

   thetas = NULL;
   distances = NULL;

   readyForUse = false;
   greedy = false;
   economyMode = false;

   memoryBlock = NULL;
   smallMemoryBlock = NULL;

   stack = NULL;
   stackPtr = 0;

   updateDiseaseScores = false;
   diseaseCount = 0;
   diseaseStatus = NULL;
   diseaseScores = NULL;

   posterior = NULL;
   rightMatrices = NULL;

   orderedGenotypes = false;
   orderedGenotypeFlags = NULL;

   skipSanityCheck = false;

   // Entry [x][y] in the penetrance matrix gives the probability of
   // observing genotype y when the true genotype is x. By default,
   // we set up the matrix to assume no errors.
   // SetErrorRate(0.0);
   }

Haplotyper::~Haplotyper()
   {
   if (individuals != 0 && markers != 0)
      {
      // In a efficient version, we should not allocate memory for phased haplotypes
      // for (int i = 0; i < individuals - phased; i++)
      //   delete [] genotypes[i];

      for (int i = 0; i < individuals; i++)
         {
         delete [] genotypes[i];
         delete [] haplotypes[i * 2];
         delete [] haplotypes[i * 2 + 1];
         }

      delete [] genotypes;
      delete [] haplotypes;
      delete [] marginals;

      delete [] leftMatrices;
      delete [] leftProbabilities;

      delete [] thetas;
      delete [] crossovers;

      delete [] error_models;
      delete [] penetrances;

      for (int i = 0; i < markers; i++)
         if (memoryBlock[i] != NULL)
            delete [] memoryBlock[i];

      delete [] stack;
      delete [] memoryBlock;
      delete [] smallMemoryBlock;

      if (distances != NULL)
         delete [] distances;

      if (orderedGenotypeFlags != NULL)
         delete [] orderedGenotypeFlags;
      }

   if (weights != NULL)
      delete [] weights;

   if (diseaseCount != 0)
      {
      FreeCharMatrix(diseaseStatus, individuals);
      FreeFloatMatrix(diseaseScores, individuals);
      FreeFloatMatrix(nplScores, 3);
      }

   if (rightMatrices != NULL)
      FreeFloatMatrix(rightMatrices, 2);

   if (posterior != NULL)
      {
      FreeFloatMatrix(posterior, 3);
      FreeFloatMatrix(mlinfo, 4);
      }
   }

void Haplotyper::AllocateWeights()
   {
   // Free original vector, if required
   if (weights != NULL) delete [] weights;

   // Allocate vector for weights
   weights = new float [individuals];

   // Trap out of memory conditions
   if (weights == NULL)
      error("Out of memory allocating weights for each individual\n");
   }

void Haplotyper::FreeWeights()
   {
   delete [] weights;
   weights = NULL;
   }

void Haplotyper::CalculateWeights()
   {
   AllocateWeights();

   // Calculate weights ...
   float sum = 0.0;
   for (int i = 0; i < individuals - phased; i++)
      {
      weights[i] = 0.0;

      for (int j = 0; j < markers; j++)
         if (genotypes[i][j] != GENOTYPE_MISSING)
            weights[i]++;

     if (weights[i] == 0.0)
      weights[i] = 1e-30;

     sum += weights[i];
     }

   // Phase known haplotypes get the maximum weight
   for (int i = individuals - phased; i < individuals; i++)
     sum += weights[i] = markers;

   // Give up if there are no genotyped individuals
   if (sum == 0.0)
     FreeWeights();
   }

void Haplotyper::RandomSetup(Random * rand)
   {
   if (rand == NULL)
      rand = &globalRandom;

   for (int j = 0; j < markers; j++)
      {
      int alleles = 0, mac = 0;

      for (int i = 0; i < individuals - phased; i++)
         if (genotypes[i][j] != GENOTYPE_MISSING)
            if ((genotypes[i][j] & GENOTYPE_ORDERED) == 0)
               alleles += 2, mac += genotypes[i][j] - 1;
            else
               {
               alleles += ((genotypes[i][j] & FIRST_ALLELE) != 0) +
                          ((genotypes[i][j] & SECOND_ALLELE) != 0);
               mac += ((genotypes[i][j] & FIRST_ALLELE_TWO) != 0) +
                      ((genotypes[i][j] & SECOND_ALLELE_TWO) != 0);
               }

      for (int i = individuals - phased; i < individuals; i++)
         {
         mac += haplotypes[i * 2][j] == 1;
         mac += haplotypes[i * 2 + 1][j] == 1;
         }
      alleles += phased * 2;

      if (alleles == 0)
         {
         for (int i = 0; i < individuals; i++)
            haplotypes[i * 2][j] = haplotypes[i * 2 + 1][j] = 0;
         continue;
         }

      float freq = mac / (float) alleles;

      for (int i = 0; i < individuals - phased; i++)
         switch (genotypes[i][j] & ~GENOTYPE_LINKED)
            {
            case 0:
              haplotypes[i * 2][j] = rand->Next() < freq;
              haplotypes[i * 2 + 1][j] = rand->Next() < freq;
              break;
            case 1:
              haplotypes[i * 2][j] = 0;
              haplotypes[i * 2 + 1][j] = 0;
              break;
            case 2:
              {
              int bit = rand->Binary();

              haplotypes[i * 2][j] = bit;
              haplotypes[i * 2 + 1][j] = bit ^ 1;
              }
              break;
            case 3:
              haplotypes[i * 2][j] = 1;
              haplotypes[i * 2 + 1][j] = 1;
              break;
            default:
               // Ordered genotype ...
               haplotypes[i * 2][j] = genotypes[i][j] & FIRST_ALLELE ?
                  (genotypes[i][j] & FIRST_ALLELE_TWO) > 0 : rand->Next() > freq;
               haplotypes[i * 2 + 1][j] = genotypes[i][j] & SECOND_ALLELE ?
                  (genotypes[i][j] & SECOND_ALLELE_TWO) > 0 : rand->Next() > freq;
            }
      }
   }

bool Haplotyper::AllocateMemory(int persons, int maxStates, int m)
   {
   individuals = persons;
   states = maxStates > 1 && maxStates < individuals * 2 ? maxStates & ~1: individuals * 2 - 2;
   markers = m;

   genotypes = AllocateCharMatrix(individuals, markers);
   haplotypes = AllocateCharMatrix(individuals * 2, markers);

   marginals = new float [states * 2];

   leftMatrices = new float * [markers];
   leftProbabilities = new float [markers];

   memoryBlock = new float * [markers];
   smallMemoryBlock = new float * [markers];
   smallFree = 0;

   stack = new int [markers];
   stackPtr = -1;

   thetas = new float [markers - 1];
   crossovers = new int [markers - 1];

   error_models = new Errors [markers];
   penetrances = new float [markers * 9];

   if (genotypes == NULL || haplotypes == NULL || marginals == NULL ||
       leftMatrices == NULL || leftProbabilities == NULL || thetas == NULL ||
       crossovers == NULL || error_models == NULL || penetrances == NULL)
       return readyForUse = false;

   for (int i = 0; i < markers; i++)
      memoryBlock[i] = smallMemoryBlock[i] = NULL;

   for (int i = 0; i < markers - 1; i++)
      thetas[i] = 0.01;

   gridSize = economyMode ? (int) sqrt((double)markers) : markers;

   orderedGenotypeFlags = new int [individuals];

   for (int i = 0; i < individuals; i++)
      orderedGenotypeFlags = 0;

   return readyForUse = true;
   }

bool Haplotyper::AllocateDiseaseStatus(int nDiseases)
   {
   diseaseCount = nDiseases;

   if (diseaseCount == 0)
      return true;

   readyForUse = AllocateRightMatrices();

   if (!readyForUse)
      return false;

   diseaseStatus = AllocateCharMatrix(individuals, nDiseases);
   diseaseScores = AllocateFloatMatrix(individuals, nDiseases * markers);
   nplScores = AllocateFloatMatrix(3, nDiseases);

   if (diseaseScores == NULL || diseaseStatus == NULL || nplScores == NULL)
      {
      diseaseCount = 0;
      if (diseaseStatus != NULL) FreeCharMatrix(diseaseStatus, individuals);
      if (diseaseScores != NULL) FreeFloatMatrix(diseaseScores, individuals);
      if (nplScores != NULL) FreeFloatMatrix(nplScores, 2);
      return false;
      }

   for (int i = 0; i < individuals; i++)
      for (int j = 0;  j < nDiseases * markers; j++)
         diseaseScores[i][j] = 0.0;

   for (int i = 0; i < individuals; i++)
      for (int j = 0; j < nDiseases; j++)
         diseaseStatus[i][j] = 0;

   return true;
   }

bool Haplotyper::AllocateMLEMemory()
   {
   readyForUse = AllocateRightMatrices();

   if (!readyForUse)
      return false;

   posterior = AllocateFloatMatrix(3, markers);
   mlinfo = AllocateFloatMatrix(4, markers);

   readyForUse = posterior != NULL && mlinfo != NULL;

   if (!readyForUse)
      {
      if (rightMatrices != NULL) FreeFloatMatrix(rightMatrices, 2);
      if (posterior != NULL) FreeFloatMatrix(posterior, 3);
      if (mlinfo != NULL) FreeFloatMatrix(mlinfo, 4);
      }

   return readyForUse;
   }

bool Haplotyper::AllocateRightMatrices()
   {
   int matrixSize = orderedGenotypes ? states * states : states * (states + 1) / 2;

   rightMatrices = AllocateFloatMatrix(2, matrixSize);

   return rightMatrices != NULL;
   }

bool Haplotyper::AllocateDistances()
   {
   distances = new float [markers - 1];

   return distances != NULL;
   }

bool Haplotyper::ForceMemoryAllocation()
   {
   // Cycle through individuals, with the exact same steps as the actual
   // haplotyper and request memory ... by requesting all memory upfront,
   // we force crashes to happen early.
   for (int i = 0; i < individuals - phased; i++)
      {
      ResetMemoryPool();
      GetMemoryBlock(0);

      if (leftMatrices[0] == NULL)
         return false;

      int skipped = 0;
      for (int j = 1; j < markers; j++)
         if (genotypes[i][j] != GENOTYPE_MISSING || j == markers - 1)
            {
            GetMemoryBlock(j);

            if (leftMatrices[j] == NULL)
               return false;
            }
         else
            skipped++;

      if (skipped == 0) break;
      }

   if (!phased)
      return true;

   ResetMemoryPool();
   for (int j = 0; j < markers; j++)
      {
      GetSmallMemoryBlock(j);

      if (leftMatrices[j] == NULL)
         return false;
      }

   return true;
   }

void Haplotyper::Transpose(float * source, float * dest, float theta)
   {
   if (theta == 0.0)
      {
      for (int i = 0; i < states; i++)
         for (int j = 0; j <= i; j++, dest++, source++)
            *dest = *source;

      return;
      }

   float sum = 0.0;
   float * probability = source;
   float * output = dest;

   for (int i = 0; i < states; i++)
      marginals[i] = 0.0;

   for (int i = 0; i < states; i++)
      {
      for (int j = 0; j < i; j++)
         {
         sum += *probability;
         marginals[i] += *probability;
         marginals[j] += *probability;
         probability++;
         }

      sum += *probability;
      marginals[i] += (*probability) * 2.0;
      probability++;
      }

   probability = source;

   float no_change = (1.0 - theta) * (1.0 - theta);
   float one_change = (1.0 - theta) * theta / states;
   float two_changes = sum * theta * theta / (states * states);

   // Automatically rescale likelihoods when they get too small
   if (sum < 1e-15)
      {
      no_change *= 1e30;
      one_change *= 1e30;
      two_changes *= 1e30;
      }

   // This final loop actually transposes the probabilities for each state
   if (weights == NULL)
      for (int i = 0; i < states; i++)
         {
         for (int j = 0; j < i; j++)
            {
            *output = *probability * no_change +
                       marginals[i] * one_change +
                       marginals[j] * one_change +
                       2 * two_changes;

            probability++;
            output++;
            }

         *output = *probability * no_change +
                    marginals[i] * one_change +
                    two_changes;

         probability++;
         output++;
         }
   else
      for (int i = 0; i < states; i++)
         {
         for (int j = 0; j < i; j++)
            {
            *output = *probability * no_change +
                       marginals[i] * one_change * weights[j / 2] +
                       marginals[j] * one_change * weights[i / 2]+
                       2 * two_changes * weights[i / 2] * weights[j / 2];

            probability++;
            output++;
            }

         *output = *probability * no_change +
                    marginals[i] * one_change * weights[i / 2] +
                    two_changes * weights[i / 2] * weights[i / 2];

         probability++;
         output++;
         }
   }

void Haplotyper::TransposeOrdered(float * source, float * dest, float theta)
   {
   if (theta == 0.0)
      {
      for (int i = 0; i < states; i++)
         for (int j = 0; j < states; j++, dest++, source++)
            *dest = *source;

      return;
      }

   float sum = 0.0;
   float * probability = source;
   float * output = dest;

   for (int i = 0; i < 2 * states; i++)
      marginals[i] = 0.0;

   for (int i = 0; i < states; i++)
      for (int j = 0; j < states; j++)
         {
         sum += *probability;
         marginals[i] += *probability;
         marginals[states + j] += *probability;
         probability++;
         }

   probability = source;

   float no_change = (1.0 - theta) * (1.0 - theta);
   float one_change = (1.0 - theta) * theta / states;
   float two_changes = sum * theta * theta / (states * states);

   // Automatically rescale likelihoods when they get too small
   if (sum < 1e-15)
      {
      no_change *= 1e30;
      one_change *= 1e30;
      two_changes *= 1e30;
      }

   // This final loop actually transposes the probabilities for each state
   if (weights == NULL)
      for (int i = 0; i < states; i++)
         for (int j = 0; j < i; j++, probability++, output++)
            *output = *probability * no_change +
                       marginals[i] * one_change +
                       marginals[states + j] * one_change +
                       2 * two_changes;
   else
      for (int i = 0; i < states; i++)
         for (int j = 0; j < states; j++, probability++, output++)
            *output = *probability * no_change +
                       marginals[i] * one_change * weights[j / 2] +
                       marginals[states + j] * one_change * weights[i / 2]+
                       2 * two_changes * weights[i / 2] * weights[j / 2];
   }

void Haplotyper::TransposeHaplotype(float * source, float * dest, float theta)
   {
   if (theta == 0.0)
      {
      for (int i = 0; i < states; i++)
         dest[i] = source[i];

      return;
      }

   float sum = 0.0;
   for (int i = 0; i < states; i++)
      sum += source[i];

   float no_change = 1.0 - theta;
   float one_change = sum * theta / states;

   // Automatically rescale likelihoods when they get too small
   if (sum < 1e-15)
      {
      no_change *= 1e30;
      one_change *= 1e30;
      }

   // This final loop actually transposes the probabilities for each state
   if (weights == NULL)
      for (int i = 0; i < states; i++)
         dest[i] = source[i] * no_change + one_change;
   else
      for (int i = 0; i < states; i++)
         dest[i] = source[i] * no_change + one_change * weights[i / 2];
   }

void Haplotyper::ConditionOnData(float * matrix, int marker, char genotype)
   {
   // We treat missing genotypes as uninformative about the mosaic's
   // underlying state. If we were to allow for deletions and the like,
   // that may no longer be true.
   if (genotype == GENOTYPE_MISSING)
      return;

   for (int i = 0; i < states; i++)
      {
      double factors[2];

      factors[0] = Penetrance(marker, haplotypes[i][marker], genotype - 1);
      factors[1] = Penetrance(marker, haplotypes[i][marker] + 1, genotype - 1);

      for (int j = 0; j <= i; j++, matrix++)
         *matrix *= factors[haplotypes[j][marker]];
      }
   }

void Haplotyper::ConditionOnOrderedData(float * matrix, int marker, char genotype)
   {
   if (genotype == GENOTYPE_MISSING)
      return;

   if ((genotype & GENOTYPE_ORDERED) == 0)
      {
      for (int i = 0; i < states; i++)
         {
         double factors[2];

         factors[0] = Penetrance(marker, haplotypes[i][marker], genotype);
         factors[1] = Penetrance(marker, haplotypes[i][marker] + 1, genotype);

         for (int j = 0; j < states; j++, matrix++)
            *matrix *= factors[haplotypes[j][marker]];
         }

      return;
      }

   double erate = GetErrorRate(marker);
   double complement = 1.0 - erate;

   double first_factor[2], second_factor[2];

   first_factor[0] = genotype & FIRST_ALLELE ?
      (genotype & FIRST_ALLELE_ONE ? complement : erate) : 1.0;
   first_factor[1] = genotype & FIRST_ALLELE ?
      (genotype & FIRST_ALLELE_ONE ? erate : complement) : 1.0;

   second_factor[0] = genotype & SECOND_ALLELE ?
      (genotype & SECOND_ALLELE_ONE ? complement : erate) : 1.0;
   second_factor[1] = genotype & SECOND_ALLELE ?
      (genotype & SECOND_ALLELE_ONE ? erate : complement) : 1.0;

   for (int i = 0; i < states; i++)
      {
      double factors[2];

      factors[0] = first_factor[haplotypes[i][marker]] * second_factor[0];
      factors[1] = first_factor[haplotypes[i][marker]] * second_factor[1];

      for (int j = 0; j < states; j++, matrix++)
         *matrix *= factors[haplotypes[j][marker]];
      }
   }

void Haplotyper::ConditionHaplotypeOnData(float * matrix, int marker, char allele)
   {
   double factors[2];

   factors[0] = GetErrorRate(marker);
   factors[1] = 1.0 - factors[0];

   for (int i = 0; i < states; i++)
      matrix[i] *= haplotypes[i][marker] == allele ? factors[1] : factors[0];
   }

bool Haplotyper::SanityCheck()
   {
   bool okay = true;

#ifdef _DEBUG

   // Current implementation ignores ordered genotype data
   if (skipSanityCheck)
      return true;

   for (int i = 0; i < markers; i++)
      if (haplotypes[states][i] + haplotypes[states + 1][i] + 1 !=
          genotypes[states / 2][i] && genotypes[states / 2][i] != GENOTYPE_MISSING)
          printf("Mismatch at marker %d\n", i + 1, okay = false);
#endif

   return okay;
   }

void Haplotyper::SetupPrior(float * matrix)
   {
   float prior = 1.0 / (states * states);

   if (weights == NULL)
      {
      for (int i = 0; i < states; i++)
         {
         for (int j = 0; j < i; j++)
            {
            *matrix = 2.0 * prior;
            matrix++;
            }
         *matrix = prior;
         matrix++;
         }
      }
   else
      for (int i = 0; i < states; i++)
         {
         for (int j = 0; j < i; j++)
            {
            *matrix = 2.0 * prior * weights[i / 2] * weights[j / 2];
            matrix++;
            }
         *matrix = prior * weights[i / 2] * weights[i / 2];
         matrix++;
         }
   }

void Haplotyper::SetupOrderedPrior(float * matrix)
   {
   float prior = 1.0 / (states * states);

   if (weights == NULL)
      for (int i = 0; i < states; i++)
         for (int j = 0; j < states; j++, matrix++)
            *matrix = prior;
   else
      for (int i = 0; i < states; i++)
         for (int j = 0; j < states; j++, matrix++)
            *matrix = prior * weights[i / 2] * weights[j / 2];
   }

void Haplotyper::ScoreLeftConditional()
   {
   ResetMemoryPool();
   GetMemoryBlock(0);

   SetupPrior(leftMatrices[0]);
   ConditionOnData(leftMatrices[0], 0, genotypes[states / 2][0]);

   double theta = 0.0;
   float *from = leftMatrices[0];
   for (int i = 1; i < markers; i++)
      {
      // Cumulative recombination fraction allows us to skip uninformative positions
      theta = theta + thetas[i - 1] - theta * thetas[i - 1];

      // Skip over uninformative positions to save time
      if (genotypes[states / 2][i] != GENOTYPE_MISSING || i == markers - 1)
         {
         GetMemoryBlock(i);

         Transpose(from, leftMatrices[i], theta);
         ConditionOnData(leftMatrices[i], i, genotypes[states / 2][i]);

         theta = 0;
         from = leftMatrices[i];
         }
      }

   MarkMemoryPool();
   }

void Haplotyper::ImputeGenotypes()
   {
   RewindMemoryPool();

   // Process the last position
   RetrieveMemoryBlock(markers - 1);
   ImputeGenotypes(leftMatrices[markers - 1], markers - 1);

   SetupPrior(rightMatrices[0]);
   ConditionOnData(rightMatrices[0], 0, genotypes[states / 2][0]);

   float *temp;
   float *from = rightMatrices[0];
   float *to = rightMatrices[1];

   for (int i = markers - 2; i >= 0; i--)
      {
      // Move things along
      Transpose(from, to, thetas[i]);

      // Find nearest informative marker
      double theta = 0.0;
      int left = i;

      while (left > 0 && genotypes[states / 2][left] == GENOTYPE_MISSING)
         {
         // Cumulative recombination fraction to nearest marker
         theta = theta + thetas[left - 1] - theta * thetas[left - 1];
         left--;
         }

      RetrieveMemoryBlock(left);
      float * leftMatrix = leftMatrices[left];

      if (left != i)
         {
         Transpose(leftMatrix, from, theta);
         leftMatrix = from;
         }

      ImputeGenotypes(leftMatrix, to, i);
      ConditionOnData(to, i, genotypes[states / 2][i]);

      temp = from;
      from = to;
      to = temp;
      }
   }

void Haplotyper::ImputeGenotypes(float * matrix, int marker)
   {
   posterior[0][marker] = posterior[1][marker] = posterior[2][marker] = 0.0;

   for (int i = 0; i < states; i++)
      if (haplotypes[i][marker])
         for (int j = 0; j <= i; j++, matrix++)
            posterior[haplotypes[j][marker] + 1][marker] += *matrix;
      else
         for (int j = 0; j <= i; j++, matrix++)
            posterior[haplotypes[j][marker]][marker] += *matrix;

   NormalizePosterior(marker);
   }

void Haplotyper::ImputeGenotypes(float * matrix1, float * matrix2, int marker)
   {
   posterior[0][marker] = posterior[1][marker] = posterior[2][marker] = 0.0;

   if (weights == NULL)
     for (int i = 0; i < states; i++, matrix1++, matrix2++)
       {
       if (haplotypes[i][marker])
         for (int j = 0; j < i; j++, matrix1++, matrix2++)
            posterior[haplotypes[j][marker] + 1][marker] += *matrix1 * *matrix2 * 0.5;
       else
         for (int j = 0; j < i; j++, matrix1++, matrix2++)
            posterior[haplotypes[j][marker]][marker] += *matrix1 * *matrix2 * 0.5;

       posterior[haplotypes[i][marker] * 2][marker] += *matrix1 * *matrix2;
       }
   else
     for (int i = 0; i < states; i++, matrix1++, matrix2++)
         {
         if (haplotypes[i][marker])
            for (int j = 0; j < i; j++, matrix1++, matrix2++)
               posterior[haplotypes[j][marker] + 1][marker] +=
              *matrix1 * *matrix2 * 0.5 / (weights[j / 2] * weights[i / 2] + 1e-30);
       else
         for (int j = 0; j < i; j++, matrix1++, matrix2++)
            posterior[haplotypes[j][marker]][marker] +=
              *matrix1 * *matrix2 * 0.5 / (weights[j / 2] * weights[i / 2]  + 1e-30);

       posterior[haplotypes[i][marker] * 2][marker] +=
            *matrix1 * *matrix2 / (weights[i / 2] * weights[i / 2] + 1e-30);
         }


   NormalizePosterior(marker);
   }

void Haplotyper::ScoreLeftConditionalForOrderedGenotypes()
   {
   ResetMemoryPool();
   GetMemoryBlock(0);

   SetupOrderedPrior(leftMatrices[0]);
   ConditionOnOrderedData(leftMatrices[0], 0, genotypes[states / 2][0]);

   double theta = 0.0;
   float *from = leftMatrices[0];
   for (int i = 1; i < markers; i++)
      {
      // Cumulative recombination fraction allows us to skip uninformative positions
      theta = theta + thetas[i - 1] - theta * thetas[i - 1];

      // Skip over uninformative positions to save time
      if (genotypes[states / 2][i] != GENOTYPE_MISSING || i == markers - 1)
         {
         GetMemoryBlock(i);

         TransposeOrdered(from, leftMatrices[i], theta);
         ConditionOnOrderedData(leftMatrices[i], i, genotypes[states / 2][i]);

         theta = 0;
         from = leftMatrices[i];
         }
      }

   MarkMemoryPool();
   }

void Haplotyper::ImputeGenotypesFromOrderedData()
   {
   RewindMemoryPool();

   // Process the last position
   RetrieveMemoryBlock(markers - 1);
   ImputeGenotypesFromOrderedData(leftMatrices[markers - 1], markers - 1);

   SetupOrderedPrior(rightMatrices[0]);
   ConditionOnOrderedData(rightMatrices[0], 0, genotypes[states / 2][0]);

   float *temp;
   float *from = rightMatrices[0];
   float *to = rightMatrices[1];

   for (int i = markers - 2; i >= 0; i--)
      {
      // Move things along
      TransposeOrdered(from, to, thetas[i]);

      // Find nearest informative marker
      double theta = 0.0;
      int left = i;

      while (left > 0 && genotypes[states / 2][left] == GENOTYPE_MISSING)
         {
         // Cumulative recombination fraction to nearest marker
         theta = theta + thetas[left - 1] - theta * thetas[left - 1];
         left--;
         }

      RetrieveMemoryBlock(left);
      float * leftMatrix = leftMatrices[left];

      if (left != i)
         {
         TransposeOrdered(leftMatrix, from, theta);
         leftMatrix = from;
         }

      ImputeGenotypesFromOrderedData(leftMatrix, to, i);
      ConditionOnOrderedData(to, i, genotypes[states / 2][i]);

      temp = from;
      from = to;
      to = temp;
      }
   }

void Haplotyper::ImputeGenotypesFromOrderedData(float * matrix, int marker)
   {
   posterior[0][marker] = posterior[1][marker] = posterior[2][marker] = 0.0;

   for (int i = 0; i < states; i++)
      if (haplotypes[i][marker])
         for (int j = 0; j < states; j++, matrix++)
            posterior[haplotypes[j][marker] + 1][marker] += *matrix;
      else
         for (int j = 0; j < states; j++, matrix++)
            posterior[haplotypes[j][marker]][marker] += *matrix;

   NormalizePosterior(marker);
   }

void Haplotyper::ImputeGenotypesFromOrderedData(float * matrix1, float * matrix2, int marker)
   {
   posterior[0][marker] = posterior[1][marker] = posterior[2][marker] = 0.0;

   if (weights == NULL)
      for (int i = 0; i < states; i++)
         if (haplotypes[i][marker])
            for (int j = 0; j < states; j++, matrix1++, matrix2++)
               posterior[haplotypes[j][marker] + 1][marker] += *matrix1 * *matrix2;
         else
            for (int j = 0; j < states; j++, matrix1++, matrix2++)
               posterior[haplotypes[j][marker]][marker] += *matrix1 * *matrix2;
   else
      for (int i = 0; i < states; i++)
         if (haplotypes[i][marker])
            for (int j = 0; j < states; j++, matrix1++, matrix2++)
               posterior[haplotypes[j][marker] + 1][marker] +=
              *matrix1 * *matrix2 / (weights[j / 2] * weights[i / 2] + 1e-30);
       else
         for (int j = 0; j < states; j++, matrix1++, matrix2++)
            posterior[haplotypes[j][marker]][marker] +=
              *matrix1 * *matrix2 / (weights[j / 2] * weights[i / 2] + 1e-30);

   NormalizePosterior(marker);
   }

void Haplotyper::ScoreLeftConditionalForHaplotype()
   {
   ResetMemoryPool();
   GetSmallMemoryBlock(0);

   float * matrix = leftMatrices[0];

   float prior = 1.0 / states;

   if (weights == NULL)
      for (int i = 0; i < states; i++)
         matrix[i] = prior;
   else
      for (int i = 0; i < states; i++)
         matrix[i] = prior * weights[i / 2];

   ConditionHaplotypeOnData(leftMatrices[0], 0, haplotypes[states][0]);

   for (int i = 1; i < markers; i++)
      {
      GetSmallMemoryBlock(i);

      TransposeHaplotype(leftMatrices[i - 1], leftMatrices[i], thetas[i - 1]);
      ConditionHaplotypeOnData(leftMatrices[i], i, haplotypes[states][i]);
      }

   MarkMemoryPool();
   }

void Haplotyper::SampleChromosomes(Random * rand)
   {
   // Print(markers - 1);
   RewindMemoryPool();
   RetrieveMemoryBlock(markers - 1);

   float * probability = leftMatrices[markers - 1];
   float sum = 0.0;

   // Calculate sum over all states
   for (int i = 0; i < states; i++)
      for (int j = 0; j <= i; j++)
         {
         sum += *probability;
         probability++;
         }

   // Sample number and select state
   float choice = rand->Uniform(0, sum);

   sum = 0.0;

   int first = 0, second = 0;
   for (probability = leftMatrices[markers - 1]; first < states; first++)
      {
      for (second = 0; second <= first; second++)
         {
         sum += *probability;
         probability++;

         if (sum >= choice) break;
         }

      if (second <= first) break;
      }

   // printf("Cumulative probability: %g\n", sum);
   // printf("           Random draw: %g\n", choice);
   // printf("        Selected state: %g\n", *(probability - 1));

   for (int j = markers - 2; j >= 0; j--)
      {
      // printf("Sum: %f, Chose (%d,%d)\n", sum, first, second);

      ImputeAlleles(j + 1, first, second, rand);

      // Starting marker for this iteration
      int   j0 = j;

      // Cumulative recombination fraction, skipping over uninformative
      // positions
      float theta = thetas[j];
      while (genotypes[states / 2][j] == GENOTYPE_MISSING && j > 0)
         {
         --j;
         theta = theta + thetas[j] - theta * thetas[j];
         }

      // When examining the previous location we consider three alternatives:
      // states that could be reached when both haplotypes recombine (11),
      // states that can be reached when the first (10) or second (01) haplotype recombines,
      // and the states that can be reached without recombination.

      float sum00 = 0.0, sum01 = 0.0, sum10 = 0.0, sum11 = 0.0;

      RetrieveMemoryBlock(j);
      probability = leftMatrices[j];

      for (int k = 0; k < states; k++)
         for (int l = 0; l <= k; l++, probability++)
            {
            sum11 += *probability;
            if (first == k || first == l) sum01 += *probability;
            if (second == k || second == l) sum10 += *probability;
            if (first == k && second == l || first == l && second == k) sum00 += *probability;
            }

      if (weights != NULL)
         {
         sum01 *= weights[second / 2];
         sum10 *= weights[first / 2];
         sum11 *= weights[second / 2] * weights[first / 2];
         }

      sum = sum11 * theta * theta / (states * states) +
            (sum10 + sum01) * theta * (1.0 - theta) / states +
            sum00 * (1.0 - theta) * (1.0 - theta);

      // Sample number and decide how many state changes occurred between the
      // two positions
      choice = rand->Uniform(0, sum);

      // The most likely outcome is that no changes occur ...
      choice -= sum00 * (1.0 - theta) * (1.0 - theta);
      if (choice <= 0.0)
         {
         // Record outcomes for intermediate, uninformative, positions
         FillPath(states, j, j0 + 1, first);
         FillPath(states + 1, j, j0 + 1, second);

         continue;
         }

      // But perhaps the first or second haplotype recombined
      probability = leftMatrices[j];

      choice -= sum10 * theta * (1.0 - theta) / states;
      if (choice <= 0.0)
         {
         // The first haplotype changed ...
         choice = choice * states / (theta * (1.0 - theta));

         // Record the original state
         int first0 = first;

         if (weights != NULL) choice /= weights[first / 2];

         for (first = 0; first < states; first++)
            {
            if (first >= second)
               choice += probability[first * (first + 1) / 2 + second];
            else
               choice += probability[second * (second + 1) / 2 + first];

            if (choice >= 0.0) break;
            }

         // Record outcomes for intermediate, uninformative, positions
         SamplePath(states, j, j0 + 1, first, first0, rand);
         FillPath(states + 1, j, j0 + 1, second);

         continue;
         }

      choice -= sum01 * theta * (1.0 - theta) / states;
      if (choice <= 0.0)
         {
         // The second haplotype changed ...
         choice = choice * states / (theta * (1.0 - theta));

         // Save the original state
         int second0 = second;

         if (weights != NULL) choice /= weights[second / 2];

         for (second = 0; second < states; second++)
            {
            if (first >= second)
               choice += probability[first * (first + 1) / 2 + second];
            else
               choice += probability[second * (second + 1) / 2 + first];

            if (choice >= 0.0) break;
            }

         // Record outcomes for intermediate, uninformative, positions
         FillPath(states, j, j0 + 1, first);
         SamplePath(states + 1, j, j0 + 1, second, second0, rand);

         continue;
         }

      // Try to select any other state
      choice *= states * states / (theta * theta);
      sum = 0.0;

      // Save the original states
      int first0 = first;
      int second0 = second;

      if (weights != NULL) choice /= weights[first / 2] * weights[second / 2];

      for (first = 0; first < states; first++)
         {
         for (second = 0; second <= first; second++, probability++)
            {
            sum += *probability;

            if (sum > choice) break;
            }

         if (second <= first) break;
         }

      if (rand->Binary())
         {
         int temp = first;
         first = second;
         second = temp;
         }

      // Record outcomes for intermediate, uninformative, positions
      SamplePath(states, j, j0 + 1, first, first0, rand);
      SamplePath(states + 1, j, j0 + 1, second, second0, rand);
      }

   ImputeAlleles(0, first, second, rand);
   }

void Haplotyper::SampleChromosomesFromOrderedData(Random * rand)
   {
   // Print(markers - 1);
   RewindMemoryPool();
   RetrieveMemoryBlock(markers - 1);

   float * probability = leftMatrices[markers - 1];
   float sum = 0.0;

   // Calculate sum over all states
   for (int i = 0; i < states; i++)
      for (int j = 0; j < states; j++, probability++)
         sum += *probability;

   // Sample number and select state
   float choice = rand->Uniform(0, sum);

   sum = 0.0;

   int first = 0, second = 0;
   for (probability = leftMatrices[markers - 1]; first < states; first++)
      {
      for (second = 0; second < states; second++, probability++)
         {
         sum += *probability;

         if (sum >= choice) break;
         }

      if (second < states) break;
      }

   for (int j = markers - 2; j >= 0; j--)
      {
      // TODO -- Impute alleles should take ordered data into account!
      ImputeAlleles(j + 1, first, second, rand);

      // Starting marker for this iteration
      int   j0 = j;

      // Cumulative recombination fraction, skipping over uninformative
      // positions
      float theta = thetas[j];
      while (genotypes[states / 2][j] == GENOTYPE_MISSING && j > 0)
         {
         --j;
         theta = theta + thetas[j] - theta * thetas[j];
         }

      // When examining the previous location we consider three alternatives:
      // states that could be reached when both haplotypes recombine (11),
      // states that can be reached when the first (10) or second (01) haplotype recombines,
      // and the states that can be reached without recombination.

      float sum00 = 0.0, sum01 = 0.0, sum10 = 0.0, sum11 = 0.0;

      RetrieveMemoryBlock(j);
      probability = leftMatrices[j];

      // CURRENT EDITING SPOT !!
      for (int k = 0; k < states; k++)
         for (int l = 0; l <= k; l++, probability++)
            {
            sum11 += *probability;
            if (first == k || first == l) sum01 += *probability;
            if (second == k || second == l) sum10 += *probability;
            if (first == k && second == l || first == l && second == k) sum00 += *probability;
            }

      if (weights != NULL)
         {
         sum01 *= weights[second / 2];
         sum10 *= weights[first / 2];
         sum11 *= weights[second / 2] * weights[first / 2];
         }

      sum = sum11 * theta * theta / (states * states) +
            (sum10 + sum01) * theta * (1.0 - theta) / states +
            sum00 * (1.0 - theta) * (1.0 - theta);

      // Sample number and decide how many state changes occurred between the
      // two positions
      choice = rand->Uniform(0, sum);

      // The most likely outcome is that no changes occur ...
      choice -= sum00 * (1.0 - theta) * (1.0 - theta);
      if (choice <= 0.0)
         {
         // Record outcomes for intermediate, uninformative, positions
         FillPath(states, j, j0 + 1, first);
         FillPath(states + 1, j, j0 + 1, second);

         continue;
         }

      // But perhaps the first or second haplotype recombined
      probability = leftMatrices[j];

      choice -= sum10 * theta * (1.0 - theta) / states;
      if (choice <= 0.0)
         {
         // The first haplotype changed ...
         choice = choice * states / (theta * (1.0 - theta));

         // Record the original state
         int first0 = first;

         if (weights != NULL) choice /= weights[first / 2];

         for (first = 0; first < states; first++)
            {
            if (first >= second)
               choice += probability[first * (first + 1) / 2 + second];
            else
               choice += probability[second * (second + 1) / 2 + first];

            if (choice >= 0.0) break;
            }

         // Record outcomes for intermediate, uninformative, positions
         SamplePath(states, j, j0 + 1, first, first0, rand);
         FillPath(states + 1, j, j0 + 1, second);

         continue;
         }

      choice -= sum01 * theta * (1.0 - theta) / states;
      if (choice <= 0.0)
         {
         // The second haplotype changed ...
         choice = choice * states / (theta * (1.0 - theta));

         // Save the original state
         int second0 = second;

         if (weights != NULL) choice /= weights[second / 2];

         for (second = 0; second < states; second++)
            {
            if (first >= second)
               choice += probability[first * (first + 1) / 2 + second];
            else
               choice += probability[second * (second + 1) / 2 + first];

            if (choice >= 0.0) break;
            }

         // Record outcomes for intermediate, uninformative, positions
         FillPath(states, j, j0 + 1, first);
         SamplePath(states + 1, j, j0 + 1, second, second0, rand);

         continue;
         }

      // Try to select any other state
      choice *= states * states / (theta * theta);
      sum = 0.0;

      // Save the original states
      int first0 = first;
      int second0 = second;

      if (weights != NULL) choice /= weights[first / 2] * weights[second / 2];

      for (first = 0; first < states; first++)
         {
         for (second = 0; second <= first; second++, probability++)
            {
            sum += *probability;

            if (sum > choice) break;
            }

         if (second <= first) break;
         }

      if (rand->Binary())
         {
         int temp = first;
         first = second;
         second = temp;
         }

      // Record outcomes for intermediate, uninformative, positions
      SamplePath(states, j, j0 + 1, first, first0, rand);
      SamplePath(states + 1, j, j0 + 1, second, second0, rand);
      }

   ImputeAlleles(0, first, second, rand);
   }

void Haplotyper::SampleHaplotypeSource(Random * rand)
   {
   RewindMemoryPool();

   float * probability = leftMatrices[markers - 1];
   float sum = 0.0;

   // Calculate sum over all states
   for (int i = 0; i < states; i++)
      sum += probability[i];

   // Sample number and select state
   float choice = rand->Uniform(0, sum);

   int haplotype;
   sum = 0.0;

   for (haplotype = 0; haplotype < states; haplotype++)
      {
      sum += probability[haplotype];

      if (sum >= choice) break;
      }

   for (int j = markers - 2; j >= 0; j--)
      {
      // Track whether imputed state matches observed allele
      if (haplotypes[haplotype][j+1] == haplotypes[states][j+1])
         error_models[j + 1].matches++;
      else
         error_models[j + 1].mismatches++;

      float theta = thetas[j];

      probability = leftMatrices[j];

      float nocross = probability[haplotype] * (1.0 - theta);

      float sum = 0.0;
      for (int k = 0; k < states; k++)
         sum += probability[k];

      float cross = sum * theta / states;

      if (weights != NULL)
         cross *= weights[haplotype / 2];

      // Sample number and decide how many state changes occurred between the
      // two positions
      choice = rand->Uniform(0, nocross + cross);

      // The most likely outcome is that no changes occur ...
      if (choice <= nocross)
         continue;

      crossovers[j]++;

      // If a crossover occured, we need to sample a state according to probability
      choice = rand->Uniform(0, sum);

      sum = 0.0;
      for (haplotype = 0; haplotype < states; haplotype++)
         {
         sum += probability[haplotype];

         if (sum >= choice) break;
         }
      }

      // Track whether imputed state matches observed allele
      if (haplotypes[haplotype][0] == haplotypes[states][0])
         error_models[0].matches++;
      else
         error_models[0].mismatches++;
   }

void Haplotyper::UpdateThetasWithDistances()
   {
   double scale = 1.0 / (individuals * 2);

   // First we estimate a base line rate to be applied to intervals with
   // 0 or 1 observed "crossovers"
   double base_crossovers = 1;
   double base_length = 0;

   for (int i = 0; i < markers - 1; i++)
      if (crossovers[i] <= 1)
         {
         base_crossovers += crossovers[i];
         base_length += distances[i];
         }

   double base_rate = base_crossovers * scale / (base_length ? base_length : 1);

   // Then we update the rate for each interval using either the number
   // of observed crossovers (if > 1) or the baseline rate
   for (int i = 0; i < markers - 1; i++)
      if (crossovers[i] > 1)
         thetas[i] = crossovers[i] * scale;
      else
         thetas[i] = base_rate * distances[i];
   }

void Haplotyper::UpdateThetas()
   {
   if (distances != NULL)
      {
      UpdateThetasWithDistances();
      return;
      }

   double scale = 1.0 / (individuals * 2);

   // First we estimate a base line rate to be applied to intervals with
   // 0 or 1 observed "crossovers"
   int base_count = 1, base_intervals = 0;

   for (int i = 0; i < markers - 1; i++)
      if (crossovers[i] <= 1)
         base_count += crossovers[i], base_intervals++;

   double base_rate = base_count * scale / (base_intervals ? base_intervals : 1);

   // Then we update the rate for each interval using either the number
   // of observed crossovers (if > 1) or the baseline rate
   for (int i = 0; i < markers - 1; i++)
      if (crossovers[i] > 1)
         thetas[i] = crossovers[i] * scale;
      else
         thetas[i] = base_rate;
   }

double Haplotyper::UpdateErrorRate()
   {
   // Group markers into those with low error rates, which are estimated as a
   // group, and those with high error rates, which are estimated individually
   Errors baseModel;

   for (int i = 0; i < markers; i++)
      if (error_models[i].mismatches <= 2)
         baseModel += error_models[i];
      else
         SetErrorRate(i, error_models[i].Update());

   baseModel.Update();

   for (int i = 0; i < markers; i++)
      if (error_models[i].mismatches <= 2)
         SetErrorRate(i, baseModel.rate);

   return baseModel.rate;
   }

double Haplotyper::GetErrorRate()
   {
   double average = 0.0;

   for (int i = 0; i < markers; i++)
      average += GetErrorRate(i);

   return average / (markers + 1e-30);
   }

int Haplotyper::TotalCrossovers()
   {
   int total = 0;

   for (int i = 0; i < markers - 1; i++)
      total += crossovers[i];

   return total;
   }

void Haplotyper::ResetCrossovers()
   {
   for (int i = 0; i < markers - 1; i++)
      crossovers[i] = 0;

   for (int i = 0; i < markers; i++)
      error_models[i].Reset();
   }

void Haplotyper::Print()
   {
   printf("Reference Haplotypes\n");

   for (int i = 0; i < states; i++)
      {
      printf("%3d ", i);

      for (int j = 0; j < markers; j++)
         printf("%d", haplotypes[i][j]);
      printf("\n");
      }

   printf("\nGenotypes to be phased\n%3s ", "");

   for (int j = 0; j < markers; j++)
      printf(genotypes[states / 2][j] == GENOTYPE_MISSING ? "." : "%d", genotypes[individuals - 1][j]);

   printf("\nSelected Haplotypes\n");

   for (int i = states; i < states + 2; i++)
      {
      printf("%3d ", i);

      for (int j = 0; j < markers; j++)
         printf("%d", haplotypes[i][j]);
      printf("\n");
      }

   printf("\n");

   if (individuals < 10)
      for (int i = 0; i < markers; i++)
         Print(i);

//   for (int i = 0; i < markers; i++)
//      {
//      printf("Left Conditional Matrix at marker %d", i);
//      Print(i);
//      printf("\n");
//      }
   }

void Haplotyper::Print(int marker)
   {
   float * pointer = leftMatrices[marker];

   for (int i = 0; i < individuals * 2 - 2; i++)
      for (int j = 0; j <= i; j++, pointer++)
         printf("state (%d,%d) = %g [geno = %d/%d] \n",
                 i, j, *pointer, haplotypes[i][marker], haplotypes[j][marker]);
   }

void Haplotyper::SetErrorRate(int marker, float rate)
   {
   // These are the penetrances for underlying homozygous genotypes
   Penetrance(marker, 0, 0) = Penetrance(marker, 2, 2) = square(1.0 - rate);
   Penetrance(marker, 0, 1) = Penetrance(marker, 2, 1) = 2 * (1.0  - rate) * rate;
   Penetrance(marker, 0, 2) = Penetrance(marker, 2, 0) = square(rate);

   // These are the penetrances for underlying heterozygous genotypes
   Penetrance(marker, 1, 0) = Penetrance(marker, 1, 2) = (1.0 - rate) * rate;
   Penetrance(marker, 1, 1) = square(1.0 - rate) + square(rate);

   // Save estimated error rate
   error_models[marker].rate = rate;
   }

void Haplotyper::SetErrorRate(float rate)
   {
   for (int i = 0; i < markers; i++)
      SetErrorRate(i, rate);
   }

void Haplotyper::UpdateDiseaseScores(int marker, int state)
   {
//   printf("Sampled state %d [status = %d, score = %.1f]\n",
//          state, diseaseStatus[state / 2][0], nplScores[diseaseStatus[state / 2][0]][0]);

   marker *= diseaseCount;

   for (int j = 0; j < diseaseCount; j++)
      diseaseScores[states / 2][marker + j] += nplScores[diseaseStatus[state / 2][j]][j];
   }

void Haplotyper::ImputeAlleles(int marker, int state1, int state2, Random * rand)
   {
   // if (updateDiseaseScores)
   //   {
   //   UpdateDiseaseScores(marker, state1);
   //   UpdateDiseaseScores(marker, state2);
   //   }

   int imputed1 = haplotypes[state1][marker];
   int imputed2 = haplotypes[state2][marker];

   int genotype = genotypes[states / 2][marker];

   if (genotype != GENOTYPE_HOMOZYGOUS_FOR_ONE &&
       genotype != GENOTYPE_HOMOZYGOUS_FOR_TWO)
      {
      haplotypes[states][marker] = imputed1;
      haplotypes[states + 1][marker] = imputed2;
      }

   if (genotype == GENOTYPE_MISSING) return;

   int differences = abs(genotype - imputed1 - imputed2 - 1);

   if (genotype == GENOTYPE_HETEROZYGOUS && differences == 0)
      error_models[marker].uncertain_pairs++;
   else
      {
      error_models[marker].matches += 2 - differences;
      error_models[marker].mismatches += differences;
      }

   if (genotype != GENOTYPE_HETEROZYGOUS) return;

   if (imputed1 == imputed2)
      if (rand->Binary())
         haplotypes[states][marker] = !imputed2;
      else
         haplotypes[states + 1][marker] = !imputed1;
   }

void Haplotyper::ImputeAllele(int haplotype, int marker, int state)
   {
   // if (updateDiseaseScores) UpdateDiseaseScores(marker, state);

   haplotypes[haplotype][marker] = haplotypes[state][marker];
   }

void Haplotyper::BuildConsensus(int samples)
   {
   char ** sample = AllocateCharMatrix(samples * 2, markers);
   char ** consensus = AllocateCharMatrix(individuals * 2, markers);

   ResetCrossovers();

   for (int i = 0, slot = individuals - 1; i < individuals - phased; i++)
      {
      SwapIndividuals(i, slot);

      // Initialize sampled haplotypes
      for (int j = 0; j < markers; j++)
         if (genotypes[slot][j] == GENOTYPE_HOMOZYGOUS_FOR_ONE ||
             genotypes[slot][j] == GENOTYPE_HOMOZYGOUS_FOR_TWO)
            for (int k = 0; k < samples; k++)
               sample[k * 2][j] = sample[k * 2 + 1][j] = genotypes[slot][j] / 2;

      ScoreLeftConditional();

      for (int j = 0; j < samples; j++)
         {
         Swap(haplotypes[slot * 2], sample[j * 2]);
         Swap(haplotypes[slot * 2 + 1], sample[j * 2 + 1]);

         SampleChromosomes(&globalRandom);

         Swap(haplotypes[slot * 2], sample[j * 2]);
         Swap(haplotypes[slot * 2 + 1], sample[j * 2 + 1]);
         }

      BuildConsensus(consensus + i * 2, sample, samples);

      SwapIndividuals(i, slot);
      }

   FreeCharMatrix(sample, samples * 2);
   FreeCharMatrix(haplotypes, individuals * 2);

   haplotypes = consensus;
   }

void Haplotyper::BuildConsensus(char ** consensus, char ** haplotypes, int count)
   {
   // The phase for each pair of haplotypes indicates their ordering
   // in relation to the consensus
   char * phase = new char [count];

   // Trap out of memory conditions
   if (phase == NULL)
      error("Out of memory allocating phase bit-array\n");

   // Select phase based on the first heterozygous position for each haplotype
   for (int i = 0; i < count; i++)
      {
      phase[i] = 0;

      for (int j = 0; j < markers; j++)
         if (haplotypes[i * 2][j] != haplotypes[i * 2 + 1][j])
            {
            phase[i] = haplotypes[i * 2][j] > haplotypes[i * 2 + 1][j];
            break;
            }
      }

   // Build consensus one position at a time ...
   for (int i = 0; i < markers; i++)
      {
      int counts[4] = {0, 0, 0, 0};

      // Count the number of occurences for each genotype
      for (int j = 0; j < count; j++)
         counts[haplotypes[j * 2 + phase[j]][i] * 2 + haplotypes[j * 2 + (phase[j] ^ 1)][i]]++;

      // Select the most likely genotype
      int best = 0;

      for (int j = 1; j < 4; j++)
         if (counts[j] > counts[best])
            best = j;

      // Assign it to the consensus
      consensus[0][i] = best / 2;
      consensus[1][i] = best % 2;

      // If a heterozygous genotype was selected, update the phase for other informative
      // haplotypes
      if (best == 0 || best == 3) continue;

      int complement = (best ^ 3);

      for (int j = 0; j < count; j++)
         if ((haplotypes[j * 2 + phase[j]][i] * 2 + haplotypes[j * 2 + (phase[j] ^ 1)][i]) == complement)
            phase[j] = phase[j] ^ 1;
      }

   delete [] phase;
   }

void Haplotyper::WarmUp(int seeds, int rounds)
   {
   if (seeds < 0 || seeds > individuals || rounds <= 0)
      return;

   int saved_individuals = individuals;

   individuals = seeds;

   for (int i = 0; i < rounds; i++)
      {
      LoopThroughChromosomes();
      UpdateThetas();
      UpdateErrorRate();
      }

   for (int i = seeds; i < saved_individuals; i++)
      {
      SwapIndividuals(i, individuals - 1);

      ScoreLeftConditional();
      SampleChromosomes(&globalRandom);

      SwapIndividuals(i, individuals - 1);
      }

   individuals = saved_individuals;
   }

void Haplotyper::SwapIndividuals(int a, int b)
   {
   // if (b < 0 || b >= individuals)
   //   printf("Bad Swap!");

   Swap(genotypes[a], genotypes[b]);
   Swap(haplotypes[a * 2], haplotypes[b * 2]);
   Swap(haplotypes[a * 2 + 1], haplotypes[b * 2 + 1]);

   if (diseaseCount)
      {
      Swap(diseaseStatus[a], diseaseStatus[b]);
      Swap(diseaseScores[a], diseaseScores[b]);
      }

   if (weights != NULL)
      {
      float temp = weights[a];
      weights[a] = weights[b];
      weights[b] = temp;
      }
   }

void Haplotyper::SwapHaplotypes(int a, int b)
   {
   Swap(haplotypes[a], haplotypes[b]);
   }

void Haplotyper::ScaleWeights()
   {
   float sum = 0.0;

   for (int i = 0; i < states / 2; i++)
      sum += weights[i];

   float scale = states / sum * 0.5;

   for (int i = 0; i < individuals; i++)
      weights[i] *= scale;
   }

void Haplotyper::SelectReferenceSet(int * array, int forWhom)
   {
   if (greedy)
      {
      // Sanity check
      // assert(states == phased * 2);

      if (states == phased * 2)
      // default greedy
         {
         // We exclude inferred haplotypes from the reference set
         for (int i = 0; i < individuals - phased; i++)
            array[i] = 0;

         // We include phased haplotypes as our reference set
         for (int i = individuals - phased; i < individuals - 1; i++)
            array[i] = 1;

         // For the last entry in the reference set, we may need to pick
         // a pair of inferred haplotypes
         if (forWhom < individuals - phased)
            array[forWhom] = 1;
         else
            array[globalRandom.NextInt() % (individuals - phased)] = 1;
         } // end of default greedy
      else
      // approximate greedy
         {
         for (int i = 0; i < individuals - phased; i++)
            array[i] = 0;

         if (forWhom >= individuals-phased)
            globalRandom.Choose( & array[individuals - phased], phased-1, states / 2);
         else
            {
            // total phased trials, states/2 successes, so could overwrite the original array[individuals - 1]
            globalRandom.Choose( & array[individuals - phased], phased, states / 2);
            if (array[individuals - 1] == 1) array[forWhom] = 1;
            array[individuals - 1] = 1;
            }
         } // end of approximate greedy
      }

   else if (weights != NULL)
      globalRandom.Choose(array, weights, individuals - 1, states / 2);
   else
      globalRandom.Choose(array, individuals - 1, states / 2);

   // Swap reference set into position
   for (int j = 0, out = 0; j < individuals; j++)
      if (array[j])
         SwapIndividuals(j, out++);
   }

void Haplotyper::LoopThroughChromosomes()
   {
   bool approximate = (states == individuals * 2 - 2) ? false : true;

   ResetCrossovers();

   int * array = NULL;

   if (approximate)
      {
      array = new int [individuals];

      if (array == NULL)
         error("Out of memory allocating array for sampling individuals\n");

      array[individuals - 1] = 1;
      }

   for (int i = individuals - 1; i >= 0; i--)
      {
      SwapIndividuals(i, individuals - 1);

      if (approximate)
         SelectReferenceSet(array, i);

      if (weights != NULL)
         ScaleWeights();

      if (updateDiseaseScores)
         ScoreNPL();

      if (i < individuals - phased)
         {
         ScoreLeftConditional();
         SampleChromosomes(&globalRandom);

         if (updateDiseaseScores && diseaseCount)
            IntegrateNPL();

#ifdef _DEBUG
         if (!SanityCheck())
            {
            printf("\nProblems above occurred haplotyping individual %d\n\n", i);
            Print();
            }
#endif
         }
      else
         {
         ScoreLeftConditionalForHaplotype();
         SampleHaplotypeSource(&globalRandom);
         SwapHaplotypes(states, states + 1);
         ScoreLeftConditionalForHaplotype();
         SampleHaplotypeSource(&globalRandom);
         SwapHaplotypes(states, states + 1);
         }

      if (approximate)
         for (int j = individuals - 1, out = states / 2; j >= 0; j--)
            if (array[j])
               SwapIndividuals(j, out--);

      SwapIndividuals(i, individuals - 1);
      }

   if (approximate)
      delete [] array;
   }

void Haplotyper::OutputMLEs(Pedigree & ped, const String & prefix, bool mldetails)
   {
   IFILE dose = ifopen(prefix + ".mldose.gz", "wt");
   IFILE geno = ifopen(prefix + ".mlgeno.gz", "wt");
   FILE * info = fopen(prefix + ".mlinfo", "wt");
   IFILE qc;
   IFILE prob;
   if (mldetails)
   {
     qc = ifopen(prefix + ".mlqc.gz", "wt");
     prob = ifopen(prefix + ".mlprob.gz", "wt");
   }

   printf("Estimating MLE for missing genotypes conditional on current state...\n");

   if (dose == NULL || info == NULL || geno == NULL)
      error("Failed to open output file for MLE estimates");

   if (mldetails && (qc == NULL | prob == NULL))
      error ("Failed to open output file for detailed MLE estimates");

   ResetMarkerInfo();

   bool approximate = (states == individuals * 2 - 2) ? false : true;

   ResetCrossovers();

   int * array = NULL;

   if (approximate)
      {
      array = new int [individuals];

      if (array == NULL)
         printf("Out of memory allocating array for sampling individuals\n");

      array[individuals - 1] = 1;
      }

   int matches = 0, partialmatches = 0, mismatches = 0;
   for (int i = 0; i < individuals - phased; i++)
      {
      SwapIndividuals(i, individuals - 1);

      if (approximate)
         SelectReferenceSet(array, i);

      if (weights != NULL)
         ScaleWeights();

      ScoreLeftConditional();
      ImputeGenotypes();

      if (approximate)
         for (int j = individuals - 1, out = states / 2; j >= 0; j--)
            if (array[j])
               SwapIndividuals(j, out--);

      SwapIndividuals(i, individuals - 1);

      ifprintf(dose, "%s->%s ML_DOSE", (const char *) ped[i].famid, (const char *) ped[i].pid);
      ifprintf(geno, "%s->%s ML_GENO", (const char *) ped[i].famid, (const char *) ped[i].pid);

      if (qc)
         ifprintf(qc, "%s->%s ML_QC", (const char *) ped[i].famid,  (const char *) ped[i].pid);

      if (prob)
         ifprintf(prob, "%s->%s ML_PROB", (const char *) ped[i].famid, (const char *) ped[i].pid);

      for (int marker = 0; marker < markers; marker++)
         {
         int best = 0;
         if (posterior[1][marker] > posterior[best][marker]) best = 1;
         if (posterior[2][marker] > posterior[best][marker]) best = 2;

         MarkerInfo * info = ped.GetMarkerInfo(marker);

         ifprintf(dose, " %.3f", posterior[0][marker] * 2.0 + posterior[1][marker]);
         ifprintf(geno, " %s/%s", (const char *) info->GetAlleleLabel((best + 2) / 2),
                                 (const char *) info->GetAlleleLabel((best + 3) / 2));

         if (qc)
            ifprintf(qc, " %.3f", posterior[best][marker]);

         if (prob)
            ifprintf(prob, " %.3f %.3f", posterior[0][marker], posterior[1][marker]);

         if (genotypes[i][marker] == 0 && ped[i].markers[marker].isKnown())
            {
            int observed = ped[i].markers[marker].SequenceCoded() - 1;

            if (observed == best)
               matches++;
            else if (observed == 1 || best == 1)
               partialmatches++;
            else mismatches++;
            }
         }

      UpdateMarkerInfo();

      ifprintf(dose, "\n");
      ifprintf(geno, "\n");

      if (qc) ifprintf(qc, "\n");
      if (prob) ifprintf(prob, "\n");
      }

   OutputMarkerInfo(info);

   fclose(info);
   ifclose(dose);
   ifclose(geno);
   if (qc) ifclose(qc);
   if (prob) ifclose(prob);

   printf("          File [%s.mlinfo] contains marker summary information ...\n"
          "          File [%s.mldose] contains MLE for dosage ...\n"
          "          File [%s.mlgeno] contains MLE for most likely genotype\n",
          (const char *) prefix, (const char *) prefix, (const char *) prefix);

   printf(!mldetails ? "\n" :
          "          File [%s.mlqc] contains MLE for quality score ...\n"
          "          File [%s.mlprobs] contains MLE probabilities for each genotype ...\n\n",
          (const char *) prefix, (const char *) prefix);

   if (matches + mismatches + partialmatches)
      {
      double total = matches + mismatches + partialmatches;

      printf("Comparing %.0f masked genotypes with MLE estimates ...\n", total);
      printf("   Estimated per genotype error rate is %.4f\n",
            (mismatches + partialmatches) / total);
      printf("   Estimated per allele error rate is %.4f\n\n",
            (mismatches + partialmatches * 0.5) / total);
      }

   if (approximate)
      delete [] array;
   }

void Haplotyper::ShowMemoryInfo()
   {
   int blocks = 0;

   for (int i = 0; i < markers; i++)
      if (memoryBlock[i])
         blocks++;

   if (states <= 0 || states > individuals * 2 - 2)
      states = 2 * individuals - 2;

   double bytes = sizeof(char) * (double) individuals * markers * 3 // Genotypes, Haplotypes
                + sizeof(float) * (double) states * 2               // Marginals
                + sizeof(float) * (double) blocks * states * (states + 1) / 2  // matrices
                + sizeof(float) * (double) markers * 11             // penetrances, probabilities, thetas
                + sizeof(int) * (double) markers                    // crossover counts
                + sizeof(Errors) * (double) markers;             // error model information

   printf("   %40s %s\n", "Haplotyping engine (actual) ...", (const char *) MemoryInfo(bytes));
   }

void Haplotyper::EstimateMemoryInfo(int Individuals, int Markers, int States, bool Compact, bool Phased)
   {
   if (States <= 0 || States > Individuals * 2 - 2)
      States = 2 * Individuals - 2;

   int positions = Compact ? 2 * (int) sqrt((double)Markers) + 1 : Markers;

   if (Phased)
      if (Markers / ((States + 1) / 2) + 1 > positions)
         positions = Markers / ((States + 1) / 2) + 1;

   double bytes = sizeof(char) * (double) Individuals * Markers * 3 // Genotypes, Haplotypes
                + sizeof(float) * (double) States * 2               // Marginals
                + sizeof(float) * (double) positions * States * (States + 1) / 2  // matrices
                + sizeof(float) * (double) Markers * 11             // penetrances, probabilities, thetas
                + sizeof(int) * (double) Markers                    // crossover counts
                + sizeof(Errors) * (double) Markers;             // error model information

   printf("   %40s %s\n", "Haplotyping engine (max) ...", (const char *) MemoryInfo(bytes));
   }

void Haplotyper::ShowMLEMemoryInfo()
   {
   EstimateMLEMemoryInfo(individuals, markers, states);
   }

void Haplotyper::EstimateMLEMemoryInfo(int Individuals, int Markers, int States)
   {
   if (States <= 0 || States > Individuals * 2 - 2)
      States = 2 * Individuals - 2;

   double bytes = sizeof(float) * (double) 2 * States * (States + 1) / 2 +
                  sizeof(float) * (double) 7 * Markers;

   printf("   %40s %s\n", "MLE Estimator ...", (const char *) MemoryInfo(bytes));
   }

void Haplotyper::EstimateDiseaseMemoryInfo(int Individuals, int Markers, int Diseases)
   {
   // TODO -- Disease memory info estimate should take into account the two right matrices

   double bytes = sizeof(short) * (double) Diseases * Markers * Individuals +
                  sizeof(char) * (double) Diseases * Individuals;

   printf("   %40s %s\n", "Non-parametric scores ...", (const char *) MemoryInfo(bytes));
   }

void Haplotyper::FillPath(int haplotype, int fromMarker, int toMarker, int state)
   {
   fromMarker++;

   while (fromMarker < toMarker)
      ImputeAllele(haplotype, fromMarker++, state);
   }

void Haplotyper::SamplePath(int haplo, int fromMarker, int toMarker, int fromState, int toState, Random * rand)
   {
   double theta = 0.0;

   // Calculate overall recombination fraction for the interval
   for (int i = fromMarker; i < toMarker; i++)
      theta = thetas[i] + theta - theta * thetas[i];

   // Impute a path between the two end markers, assuming no genotypes
   // are observed -- the only constraint is that we must start at
   // fromState and end at toState with at least one intervening recombinant
   while (fromMarker < toMarker - 1)
      {
      double r = rand->Uniform(0.0, theta);

      double theta1 = thetas[fromMarker];

      if (theta < 0.9)
         // Fast closed formula
         theta = (theta - theta1) / (1.0 - theta1);
      else
         {
         theta = 0.0;

         // More accurate, iterative formula
         for (int i = fromMarker + 1; i < toMarker; i++)
            theta = thetas[i] + theta - theta * thetas[i];
         }

      if (r > theta1)
         {
         // No recombinant in the in first interval
         ImputeAllele(haplo, ++fromMarker, fromState);
         continue;
         }

      crossovers[fromMarker]++;
      if (r < theta1 * (1.0 - theta))
         {
         // No recombinant in the second interval
         FillPath(haplo, fromMarker, toMarker, toState);
         return;
         }
      else if (weights != NULL)
         {
         // Recombinants in both intervals, so we must sample
         // an intervening state -- potentially taking weights
         // into account
         double sum = 0.0;

         for (int i = 0; i < states; i++)
            sum += weights[i];

         r = rand->Uniform(0, sum);

         sum = weights[0];
         fromState = 0;

         for (int i = 1; i < states, sum < r; i++)
            {
            sum += weights[i];
            fromState++;
            }

         ImputeAllele(haplo, ++fromMarker, fromState);
         }
      else
         ImputeAllele(haplo, ++fromMarker, fromState = (int) (rand->Next() * states));
      }

   // If we get here, record obligate recombinant between two consecutive markers
   crossovers[fromMarker]++;
   }

// Memory management functions
//

void Haplotyper::GetMemoryBlock(int marker)
   {
   if (!economyMode || marker == 0 || marker > stack[stackPtr] + gridSize)
      {
      stack[++stackPtr] = marker;
      leftMatrices[marker] = GetLargeBlock();

      ResetReuseablePool();
      }
   else
      leftMatrices[marker] = GetReuseableBlock();
   }

void Haplotyper::GetSmallMemoryBlock(int marker)
   {
   leftMatrices[marker] = GetSmallBlock();
   }

void Haplotyper::RetrieveMemoryBlock(int marker)
   {
   if (stack[stackPtr] <= marker)
      return;
   else
      {
      ResetReuseablePool();

      double theta = 0.0;
      float *from = leftMatrices[stack[--stackPtr]];

      for (int i = stack[stackPtr] + 1; i <= marker; i++)
         {
         // Cumulative recombination fraction allows us to skip uninformative positions
         theta = theta + thetas[i - 1] - theta * thetas[i - 1];

         // Skip over uninformative positions to save time
         if (genotypes[states / 2][i] != GENOTYPE_MISSING || i == markers - 1)
            {
            leftMatrices[i] = GetReuseableBlock();

            Transpose(from, leftMatrices[i], theta);
            ConditionOnData(leftMatrices[i], i, genotypes[states / 2][i]);

            theta = 0;
            from = leftMatrices[i];
            }
         }
      }
   }

float * Haplotyper::AllocateMemoryBlock()
   {
   int blockSize = orderedGenotypes ? states * states : states * (states + 1) / 2;

   float * block = new float [blockSize];

   for (int i = 0; i < blockSize / states; i++)
      if (smallFree < markers)
         smallMemoryBlock[smallFree++] = block + i * states;

   return block;
   }

float * Haplotyper::GetLargeBlock()
   {
   if (memoryBlock[nextAvailable] == NULL)
      memoryBlock[nextAvailable] = AllocateMemoryBlock();

   return memoryBlock[nextAvailable++];
   }

float * Haplotyper::GetSmallBlock()
   {
   if (smallMemoryBlock[nextSmallAvailable] == NULL)
      {
      while (memoryBlock[nextAvailable] != NULL)
         nextAvailable++;

      memoryBlock[nextAvailable++] = AllocateMemoryBlock();
      }

   return smallMemoryBlock[nextSmallAvailable++];
   }

float * Haplotyper::GetReuseableBlock()
   {
   if (memoryBlock[nextReuseable] == NULL)
      memoryBlock[nextReuseable] = AllocateMemoryBlock();

   return memoryBlock[nextReuseable--];
   }

void Haplotyper::ResetMemoryPool()
   {
   nextAvailable = nextSmallAvailable = 0;
   nextReuseable = markers - 1;
   stackPtr = -1;
   }

void Haplotyper::ResetReuseablePool()
   {
   nextReuseable = markers - 1;
   }

void Haplotyper::ScoreNPL()
   {
   if (diseaseCount == 0)
      return;

   for (int i = 0; i < 3; i++)
      for (int j = 0; j < diseaseCount; j++)
         nplScores[i][j] = 0.0;

   for (int i = 0;  i < states / 2; i++)
      for (int j = 0; j < diseaseCount; j++)
         if (diseaseStatus[i][j])
            nplScores[3 - diseaseStatus[i][j]][j] += weights == NULL ? 1 : weights[i];

   for (int j = 0; j < diseaseCount; j++)
      nplScores[1][j] = -nplScores[1][j];
   }

void Haplotyper::MarkMemoryPool()
   {
   savedStackPtr = stackPtr;
   }

void Haplotyper::RewindMemoryPool()
   {
   if (stackPtr != savedStackPtr)
      {
      stackPtr = savedStackPtr;

      if (stack[stackPtr] != markers - 1)
         stack[++stackPtr] = markers - 1;
      }
   }

bool Haplotyper::LoadCrossoverRates(const char * filename)
   {
   IFILE f = ifopen(filename, "rb");

   if (f == NULL)
   {
      printf("Warning: crossover rate map file [%s] not available\n", (const char *) filename);
      return false;               
   }


   LoadCrossoverRates(f);
   ifclose(f);

   return true;
   }

void Haplotyper::LoadCrossoverRates(IFILE file)
   {
   String      buffer;
   StringArray tokens;

   printf("Loading crossover rates ...\n");
   buffer.ReadLine(file);

   int interval = 0;
   while (!ifeof(file) && interval < markers - 1)
      {
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (tokens.Length() == 0) continue;

      if (tokens.Length() < 2)
         error("The following line could not be parsed:\n\n%s\n"
               "Each line should list an interval followed by the corresponding\n"
               "crossover rate. Intervals should be in map order.", (const char *) buffer);

      thetas[interval++] = tokens[1].AsDouble();
      }
   }

bool Haplotyper::LoadErrorRates(const char * filename)
   {
   IFILE f = ifopen(filename, "rb");

   if (f == NULL)
   {
      printf("Warning: error rate map file [%s] not available\n", (const char *) filename);
      return false;               
   }

   LoadErrorRates(f);
   ifclose(f);

   return true;
   }

void Haplotyper::LoadErrorRates(IFILE file)
   {
   String      buffer;
   StringArray tokens;

   printf("Loading mosaic error rates ...\n");
   buffer.ReadLine(file);

   int marker = 0;
   while (!ifeof(file) && marker < markers)
      {
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (tokens.Length() == 0) continue;

      if (tokens.Length() < 2)
         error("The following line could not be parsed:\n\n%s\n"
               "Each line should list a marker followed by the corresponding\n"
               "error rate. Markers should be in map order.", (const char *) buffer);

      SetErrorRate(marker++, tokens[1].AsDouble());
      }
   }

void Haplotyper::ResetMarkerInfo()
   {
   for (int i = 0; i < markers; i++)
      mlinfo[0][i] = mlinfo[1][i] = mlinfo[2][i] = mlinfo[3][i] = 0.0;
   }

#ifndef square
#define square(x)    ((x)*(x))
#endif

void Haplotyper::UpdateMarkerInfo()
   {
   for (int i = 0; i < markers; i++)
      {
      int best = 0;
      if (posterior[1][i] > posterior[best][i]) best = 1;
      if (posterior[2][i] > posterior[best][i]) best = 2;

      mlinfo[0][i] += posterior[0][i];
      mlinfo[1][i] += posterior[1][i];
      mlinfo[2][i] += square(posterior[0][i] + posterior[1][i] * 0.50);
      mlinfo[3][i] += posterior[best][i];
      }
   }

void Haplotyper::OutputMarkerInfo(FILE * output)
   {
   fprintf(output, "SNP\tAl1\tAl2\tFreq1\tMAF\tQuality\tRsq\n");

   for (int i = 0; i < markers; i++)
      {
      double p0 = mlinfo[0][i] / (individuals - phased + 1e-30);
      double p1 = mlinfo[1][i] / (individuals - phased + 1e-30);
      double sumsq = mlinfo[2][i] / (individuals - phased + 1e-30);
      double qc = mlinfo[3][i] / (individuals - phased + 1e-30);

      double freq = p0 + p1 * 0.50;
      double var1 = max(p0 + p1 * 0.25 - square(freq), 0);
      double var2 = max(sumsq - square(freq), 0);

      // To avoid problems due to rounding in calculation of sumsq - square(freq)
      if (var2 < 1e-7) var2 = 0.0;

      MarkerInfo * info = Pedigree::GetMarkerInfo(i);

      fprintf(output, "%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n",
              (const char *) info->name,
              (const char *) info->GetAlleleLabel(1),
              info->CountAlleles() > 1 ? (const char *) info->GetAlleleLabel(2) : "-",
              freq, freq > 0.50 ? 1.0 - freq : freq, qc, var2 / (var1 + 1e-30));
      }
   }

void Haplotyper::IntegrateNPL(float * matrix, int marker)
   {
   float * matrixStart = matrix;
   marker *= diseaseCount;

   for (int disease = 0; disease < diseaseCount; disease++)
      {
      double source_status[3] = { 0.0, 0.0, 0.0 };

      matrix = matrixStart;
      for (int i = 0; i < states; i++)
         {
         int i_status = diseaseStatus[i/2][disease];

         for (int j = 0; j <= i; j++, matrix++)
            {
            source_status[i_status] += *matrix;
            source_status[diseaseStatus[j/2][disease]] += *matrix;
            }
         }

      double sum = source_status[1] + source_status[2];

      if (sum > 0.0)
         diseaseScores[states/2][marker + disease] +=
            nplScores[1][disease] * source_status[1] / sum +
            nplScores[2][disease] * source_status[2] / sum;
      }
   }

void Haplotyper::IntegrateNPL(float * matrix1, float * matrix2, int marker)
   {
   float  * matrixStart1 = matrix1;
   float  * matrixStart2 = matrix2;
   marker *= diseaseCount;

   if (weights == NULL)
      for (int disease = 0; disease < diseaseCount; disease++)
         {
         double source_status[3] = { 0.0, 0.0, 0.0 };

         matrix1 = matrixStart1; matrix2 = matrixStart2;
         for (int i = 0; i < states; i++, matrix1++, matrix2++)
            {
            int i_status = diseaseStatus[i/2][disease];

            for (int j = 0; j < i; j++, matrix1++, matrix2++)
               {
               source_status[i_status] += *matrix1 * *matrix2;
               source_status[diseaseStatus[j/2][disease]] += *matrix1 * *matrix2;
               }

            source_status[i_status] += *matrix1 * *matrix2 * 4.0;
            }

         double sum = source_status[1] + source_status[2];

         if (sum > 0.0)
            diseaseScores[states/2][marker + disease] +=
               nplScores[1][disease] * source_status[1] / sum +
               nplScores[2][disease] * source_status[2] / sum;
         }
   else
      {
      for (int disease = 0; disease < diseaseCount; disease++)
         {
         double source_status[3] = { 0.0, 0.0, 0.0 };

         matrix1 = matrixStart1; matrix2 = matrixStart2;
         for (int i = 0; i < states; i++, matrix1++, matrix2++)
            {
            int i_status = diseaseStatus[i/2][disease];

            for (int j = 0; j <= i; j++, matrix1++, matrix2++)
               {
               double cell = *matrix1 * *matrix2 / (weights[j / 2] * weights[i / 2] + 1e-30);

               source_status[i_status] += cell;
               source_status[diseaseStatus[j/2][disease]] += cell;
               }

         source_status[i_status] += *matrix1 * *matrix2 * 4.0 / (weights[i / 2] * weights[i / 2] + 1e-30);
            }

         double sum = source_status[1] + source_status[2];

         if (sum > 0.0)
            diseaseScores[states/2][marker + disease] +=
               nplScores[1][disease] * source_status[1] / sum +
               nplScores[2][disease] * source_status[2] / sum;
         }
      }
   }

void Haplotyper::IntegrateNPL()
   {
   RewindMemoryPool();

   // Process the last position
   RetrieveMemoryBlock(markers - 1);
   IntegrateNPL(leftMatrices[markers - 1], markers - 1);

   SetupPrior(rightMatrices[0]);
   ConditionOnData(rightMatrices[0], 0, genotypes[states / 2][0]);

   float *temp;
   float *from = rightMatrices[0];
   float *to = rightMatrices[1];

   for (int i = markers - 2; i >= 0; i--)
      {
      // Move things along
      Transpose(from, to, thetas[i]);

      // Find nearest informative marker
      double theta = 0.0;
      int left = i;

      while (left > 0 && genotypes[states / 2][left] == GENOTYPE_MISSING)
         {
         // Cumulative recombination fraction to nearest marker
         theta = theta + thetas[left - 1] - theta * thetas[left - 1];
         left--;
         }

      RetrieveMemoryBlock(left);
      float * leftMatrix = leftMatrices[left];

      if (left != i)
         {
         Transpose(leftMatrix, from, theta);
         leftMatrix = from;
         }

      IntegrateNPL(leftMatrix, to, i);
      ConditionOnData(to, i, genotypes[states / 2][i]);

      temp = from;
      from = to;
      to = temp;
      }
   }

 
