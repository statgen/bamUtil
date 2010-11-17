////////////////////////////////////////////////////////////////////// 
// mach1/DosageCalculator.cpp 
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
 
#include "DosageCalculator.h"
#include "MemoryAllocators.h"
#include "Pedigree.h"
#include "MemoryInfo.h"

#include <math.h>
#include <limits.h>

#ifndef UCHAR_MAX
#define UCHAR_MAX  255
#endif

#ifndef USHRT_MAX
#define USHRT_MAX  65535
#endif

#ifndef UINT_MAX
#define UINT_MAX   4294967295U
#endif

bool DosageCalculator::storeDosage = false;
bool DosageCalculator::storeDistribution = false;

DosageCalculator::DosageCalculator(int N, int G, int M)
   {
   storeDosage |= storeDistribution;

   stored = 0;
   samples = N;
   genotypes = G;
   markers = M;

   cTwo = cDosage = NULL;
   sDosage = sDosage = NULL;
   iDosage = iDosage = NULL;

   if (!storeDosage)
      {
      readyForUse = true;
      wordSize = 0;
      return;
      }

   wordSize = N <= (UCHAR_MAX / 2) ? 1 : N < (USHRT_MAX / 2) ? 2 : 4;

   switch (wordSize)
      {
      case 1 : cDosage = AllocateMatrix<unsigned char>(G, M, 0); break;
      case 2 : sDosage = AllocateMatrix<unsigned short>(G, M, 0); break;
      case 4 : iDosage = AllocateMatrix<unsigned int>(G, M, 0); break;
      }

   readyForUse = cDosage != NULL || sDosage != NULL || iDosage != NULL;

   if (!storeDistribution | !readyForUse) return;

   switch (wordSize)
      {
      case 1 : cTwo = AllocateMatrix<unsigned char>(G, M, 0); break;
      case 2 : sTwo = AllocateMatrix<unsigned short>(G, M, 0); break;
      case 4 : iTwo = AllocateMatrix<unsigned int>(G, M, 0); break;
      }

   readyForUse = cTwo != NULL || sTwo != NULL || iTwo != NULL;
   }

DosageCalculator::~DosageCalculator()
   {
   switch (wordSize)
      {
      case 1 :
         if (cDosage != NULL) FreeMatrix(cDosage, genotypes);
         if (cTwo != NULL) FreeMatrix(cTwo, genotypes);
         break;
      case 2 :
         if (sDosage != NULL) FreeMatrix(sDosage, genotypes);
         if (sTwo != NULL) FreeMatrix(sTwo, genotypes);
         break;
      case 4 :
         if (iDosage != NULL) FreeMatrix(iDosage, genotypes);
         if (iTwo != NULL) FreeMatrix(iTwo, genotypes);
         break;
      }
   }

double DosageCalculator::GetDosage(int individual, int genotype)
   {
   if (stored == 0) return 0.0;

   switch (wordSize)
      {
      case 1 : return (2 * stored - cDosage[individual][genotype]) / (double) stored;
      case 2 : return (2 * stored - sDosage[individual][genotype]) / (double) stored;
      case 4 : return (2 * stored - iDosage[individual][genotype]) / (double) stored;
      }

   return 0.0;
   }

void DosageCalculator::GetCounts(int individual, int genotype,
                       unsigned int & n0, unsigned int & n1, unsigned int & n2)
   {
   n0 = n1 = n2 = 0;

   if (stored == 0)
      return;

   switch (wordSize)
      {
      case 1 :
         n2 = cTwo[individual][genotype];
         n1 = cDosage[individual][genotype] - 2 * n2;
         break;
      case 2 :
         n2 = sTwo[individual][genotype];
         n1 = sDosage[individual][genotype] - 2 * n2;
         break;
      case 4 :
         n2 = iTwo[individual][genotype];
         n1 = iDosage[individual][genotype] - 2 * n2;
         break;
      }

   n0 = stored - n2 - n1;
   }

double DosageCalculator::GetQuality(int individual, int genotype)
   {
   if (stored == 0) return 0.0;

   unsigned int n2, n1, n0;

   GetCounts(individual, genotype, n0, n1, n2);

   if (n0 >= n2 && n0 >= n1)
      return n0 / (double) stored;

   if (n1 >= n2)
      return n1 / (double) stored;

   return n2 / (double) stored;
   }

int DosageCalculator::GetBestGenotype(int individual, int genotype)
   {
   if (stored == 0) return 1;

   unsigned int n2, n1, n0;

   GetCounts(individual, genotype, n0, n1, n2);

   if (n0 >= n2 && n0 >= n1)
      return 0;

   if (n1 >= n2)
      return 1;

   return 2;
   }

void DosageCalculator::EstimateMemoryInfo(int Samples, int Genotypes, int Markers)
   {
   if (storeDosage == false && storeDistribution == false)
      return;

   int bytesPerItem = Samples < (UCHAR_MAX / 2) ? sizeof(unsigned char) :
                      Samples < (USHRT_MAX / 2) ? sizeof(unsigned short) :
                                                  sizeof(unsigned int);

   double bytes = bytesPerItem * (double) Genotypes * Markers;

   if (storeDistribution) bytes *= 2;

   printf("   %40s %s\n", "Dosage Calculator ...", (const char *) MemoryInfo(bytes));
   }

void DosageCalculator::ShowMemoryInfo()
   {
   EstimateMemoryInfo(samples, genotypes, markers);
   }

void DosageCalculator::Update(char ** newHaplotypes)
   {
   if (storeDosage == false)
      return;

   if (wordSize == 1)
      for (int i = 0; i < genotypes * 2; i += 2)
         for (int j = 0; j < markers; j++)
            cDosage[i/2][j] += newHaplotypes[i][j] + newHaplotypes[i + 1][j];
   else if (wordSize == 2)
      for (int i = 0; i < genotypes * 2; i += 2)
         for (int j = 0; j < markers; j++)
            sDosage[i/2][j] += newHaplotypes[i][j] + newHaplotypes[i + 1][j];
   else if (wordSize == 4)
      for (int i = 0; i < genotypes * 2; i += 2)
         for (int j = 0; j < markers; j++)
            iDosage[i/2][j] += newHaplotypes[i][j] + newHaplotypes[i + 1][j];

   if (storeDistribution)
      {
      if (wordSize == 1)
         for (int i = 0; i < genotypes * 2; i += 2)
            for (int j = 0; j < markers; j++)
               cTwo[i/2][j] += (newHaplotypes[i][j] + newHaplotypes[i + 1][j]) == 2;
      else if (wordSize == 2)
         for (int i = 0; i < genotypes * 2; i += 2)
            for (int j = 0; j < markers; j++)
               sTwo[i/2][j] += (newHaplotypes[i][j] + newHaplotypes[i + 1][j]) == 2;
      else if (wordSize == 4)
         for (int i = 0; i < genotypes * 2; i += 2)
            for (int j = 0; j < markers; j++)
               iTwo[i/2][j] += (newHaplotypes[i][j] + newHaplotypes[i + 1][j]) == 2;
      }


   stored++;
   }

#ifndef square
#define square(x) ((x)*(x))
#endif

#ifndef max
#define max(a,b)     ((a)>(b)?(a):(b))
#endif

void DosageCalculator::OutputMarkerInfo(const char * filename)
   {
   FILE * f = fopen(filename, "wt");

   if (f == NULL)
      {
      printf("Failed to open output file [%s]\n", filename);
      return;
      }

   OutputMarkerInfo(f);

   fclose(f);

   printf("Wrote out file [%s] with marker information\n", filename);
   }

void DosageCalculator::OutputMarkerInfo(FILE * output)
   {
   if (stored == 0) return;

   if (!storeDistribution)
      {
      OutputBasicMarkerInfo(output);
      return;
      }

   fprintf(output, "SNP\tAl1\tAl2\tFreq1\tMAF\tQuality\tRsq\n");

   double scale_sg = 1.0 / (samples * genotypes + 1e-30);
   double scale_g  = 1.0 / (genotypes + 1e-30);
   double scale_ss = 1.0 / (samples * samples + 1e-30);

   for (int marker = 0; marker < markers; marker++)
      {
      double p0 = 0.0, p1 = 0.0;
      double qc = 0.0, sumsq = 0.0;

      for (int sample = 0; sample < genotypes; sample++)
         {
         unsigned int n0, n1, n2;

         GetCounts(sample, marker, n0, n1, n2);

         p0 += n0; p1 += n1;
         qc += (n0 > n1 && n0 > n2) ? n0 : (n1 > n2) ? n1 : n2;
         sumsq += square(n0 + n1 * 0.5) * scale_ss;
         }

      p0 *= scale_sg; p1 *= scale_sg;
      qc *= scale_sg; sumsq *= scale_g;

      double freq = p0 + p1 * 0.50;
      double var1 = max(p0 + p1 * 0.25 - square(freq), 0);
      double var2 = max(sumsq - square(freq), 0);

      MarkerInfo * info = Pedigree::GetMarkerInfo(marker);

      fprintf(output, "%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n",
                     (const char *) info->name,
                     (const char *) info->GetAlleleLabel(1),
                     info->CountAlleles() > 1 ? (const char *) info->GetAlleleLabel(2) : "-",
                     freq, freq > 0.50 ? 1.0 - freq : freq, qc, var2 / (var1 + 1e-30));
      }
   }

void DosageCalculator::OutputBasicMarkerInfo(FILE * output)
   {
   fprintf(output, "SNP\tAl1\tAl2\tFreq\tRsq_hat\n");

   double scale  = 1.0 / (genotypes + 1e-30);

   for (int marker = 0; marker < markers; marker++)
      {
      double sum = 0.0, sumsq = 0.0;

      for (int sample = 0; sample < genotypes; sample++)
         {
         double dose = GetDosage(sample, marker);

         sum += dose;
         sumsq += dose * dose;
         }

      sum *= scale;
      sumsq *= scale;

      double freq = sum * 0.50;
      double var1 = 2 * freq * (1.0 - freq);
      double var2 = max(sumsq - sum * sum, 0.0);

      MarkerInfo * info = Pedigree::GetMarkerInfo(marker);

      fprintf(output, "%s\t%s\t%s\t%.4f\t%.4f\n",
                     (const char *) info->name,
                     (const char *) info->GetAlleleLabel(1),
                     info->CountAlleles() > 1 ? (const char *) info->GetAlleleLabel(2) : "-",
                     freq > 0.50 ? 1.0 - freq : freq, var2 > var1 ? 1.0 : var2 / (var1 + 1e-30));
      }
   }


 
