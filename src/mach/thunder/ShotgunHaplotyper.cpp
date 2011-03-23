////////////////////////////////////////////////////////////////////// 
// thunder/ShotgunHaplotyper.cpp 
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
#include "MemoryAllocators.h"

#include <math.h>

ShotgunHaplotyper::ShotgunHaplotyper()
   {
   skipSanityCheck = true;

   shotgunErrorMatrix = AllocateFloatMatrix(3, 256);

   SetShotgunError(0.005);
   }

ShotgunHaplotyper::~ShotgunHaplotyper()
   {
   FreeFloatMatrix(shotgunErrorMatrix, 3);
   }

void ShotgunHaplotyper::CalculateWeights()
   {
   AllocateWeights();

   // Calculate weights ...
   float sum = 0.0;
   for (int i = 0; i < individuals - phased; i++)
      {
      weights[i] = 0.0;

      for (int j = 0; j < markers; j++)
         weights[i] += (genotypes[i][j] % 16) + (genotypes[i][j] / 16);

      sum += weights[i];
      }

   // Give up if there are no genotyped individuals
   if (sum == 0.0)
      FreeWeights();
   }

void ShotgunHaplotyper::RandomSetup(Random * rand)
   {
   if (rand == NULL)
      rand = &globalRandom;

   for (int j = 0; j < markers; j++)
      {
      int alleles = 0, mac = 0;

      for (int i = 0; i < individuals; i++)
         {
         int g = (unsigned char) genotypes[i][j];

         alleles += (g % 16);
         mac += (g / 16);
         }
      alleles += mac;

      if (alleles == 0)
         {
         for (int i = 0; i < individuals; i++)
            haplotypes[i * 2][j] = haplotypes[i * 2 + 1][j] = 0;
         continue;
         }

      double freq = mac / (double) alleles;

      double prior_11 = (1.0 - freq) * (1.0 - freq);
      double prior_12 = 2.0 * freq * (1.0 - freq);
      double prior_22 = freq * freq;

      for (int i = 0; i < individuals; i++)
         {
         int observed = (unsigned char) (genotypes[i][j]);

         double posterior_11 = prior_11 * shotgunErrorMatrix[0][observed];
         double posterior_12 = prior_12 * shotgunErrorMatrix[1][observed];
         double posterior_22 = prior_22 * shotgunErrorMatrix[2][observed];
         double sum = posterior_11 + posterior_12 + posterior_22;

         if (sum == 0)
            printf("Problem!\n");

         posterior_11 /= sum;
         posterior_12 /= sum;

         double r = rand->Next();

         if (r < posterior_11)
            {
            haplotypes[i * 2][j] = 0;
            haplotypes[i * 2 + 1][j] = 0;
            }
         else if (r < posterior_11 + posterior_12)
            {
            bool bit = rand->Binary();

            haplotypes[i * 2][j] = bit;
            haplotypes[i * 2 + 1][j] = bit ^ 1;
            }
         else
            {
            haplotypes[i * 2][j] = 1;
            haplotypes[i * 2 + 1][j] = 1;
            }
         }
      }
   }

void ShotgunHaplotyper::SetShotgunError(double rate)
   {
   // Store the background rate
   shotgunError = rate;

   // First calculate binomial coefficients
   int binomial[33][33];

   binomial[0][0] = 1;
   binomial[1][0] = binomial[1][1] = 1;

   for (int i = 2; i < 32; i++)
      {
      binomial[i][0] = binomial[i][i] = 1;

      for (int j = 1; (j < i) && (j < 16); j++)
         binomial[i][j] = binomial[i-1][j] + binomial[i-1][j-1];
      }

   // Next setup the error matrices for each possible true genotype
   for (int i = 0; i < 16; i++)
      for (int j = 0; j < 16; j++)
         if (rate == 0)
            {
            shotgunErrorMatrix[0][j*16 + i] = j == 0 ? 1.0 : 0.0;
            shotgunErrorMatrix[1][j*16 + i] = pow(0.5, i + j) * binomial[i+j][j];
            shotgunErrorMatrix[2][j*16 + i] = i == 0 ? 1.0 : 0.0;
            }
         else
            {
            shotgunErrorMatrix[0][j*16 + i] = pow(1.0 - rate, i) * pow(rate, j) * binomial[i+j][j];
            shotgunErrorMatrix[1][j*16 + i] = pow(0.5, i + j) * binomial[i+j][j];
            shotgunErrorMatrix[2][j*16 + i] = pow(rate, i) * pow(1.0 - rate, j) * binomial[i+j][j];
            }
   }

void ShotgunHaplotyper::ConditionOnData(float * matrix, int marker, char genotype)
   {
   // We treat missing genotypes as uninformative about the mosaic's
   // underlying state. If we were to allow for deletions and the like,
   // that may no longer be true.
   if (genotype == GENOTYPE_MISSING)
      return;

   int g = (unsigned char) genotype;

   double conditional_probs[3];

   for (int i = 0; i < 3; i++)
      conditional_probs[i] = Penetrance(marker, i, 0) *  shotgunErrorMatrix[0][g] +
                             Penetrance(marker, i, 1) *  shotgunErrorMatrix[1][g] +
                             Penetrance(marker, i, 2) *  shotgunErrorMatrix[2][g];

   for (int i = 0; i < states; i++)
      {
      double factors[2];

      factors[0] = conditional_probs[haplotypes[i][marker]];
      factors[1] = conditional_probs[haplotypes[i][marker] + 1];

      for (int j = 0; j <= i; j++, matrix++)
         *matrix *= factors[haplotypes[j][marker]];
      }
   }

void ShotgunHaplotyper::ImputeAlleles(int marker, int state1, int state2, Random * rand)
   {
   int copied1 = haplotypes[state1][marker];
   int copied2 = haplotypes[state2][marker];

   int genotype = (unsigned char) genotypes[states / 2][marker];

   double posterior_11 = Penetrance(marker, copied1 + copied2, 0) * shotgunErrorMatrix[0][genotype];
   double posterior_12 = Penetrance(marker, copied1 + copied2, 1) * shotgunErrorMatrix[1][genotype];
   double posterior_22 = Penetrance(marker, copied1 + copied2, 2) * shotgunErrorMatrix[2][genotype];
   double sum = posterior_11 + posterior_12 + posterior_22;

   posterior_11 /= sum;
   posterior_22 /= sum;

   double r = rand->Next();

   if (r < posterior_11)
      {
      haplotypes[states][marker] = 0;
      haplotypes[states + 1][marker] = 0;
      }
   else if (r < posterior_11 + posterior_22)
      {
      haplotypes[states][marker] = 1;
      haplotypes[states + 1][marker] = 1;
      }
   else if (copied1 != copied2)
      {
      double rate = GetErrorRate(marker);

      if (rand->Next() < rate * rate / ((rate * rate) + (1.0 - rate) * (1.0 - rate)))
         {
         copied1 = !copied1;
         copied2 = !copied2;
         }

      haplotypes[states][marker] = copied1;
      haplotypes[states + 1][marker] = copied2;
      }
   else
      {
      bool bit = rand->Binary();

      haplotypes[states][marker] = bit;
      haplotypes[states + 1][marker] = bit ^ 1;
      }

   int imputed1 = haplotypes[states][marker];
   int imputed2 = haplotypes[states + 1][marker];

   int differences = abs(copied1 - imputed1) + abs(copied2 - imputed2);

   error_models[marker].matches += 2 - differences;
   error_models[marker].mismatches += differences;
   }


 
