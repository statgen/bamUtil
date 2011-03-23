////////////////////////////////////////////////////////////////////// 
// mach1/CostCalculator.cpp 
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
 
#include "CostCalculator.h"
#include "MemoryAllocators.h"
#include "MathStats.h"
#include "LongLongCounter.h"

#include <math.h>

CostCalculator::CostCalculator()
   {
   cost = NULL;
   path = NULL;
   }

CostCalculator::~CostCalculator()
   {
   if (cost != NULL) delete [] cost;
   if (path != NULL) delete [] path;
   }

void CostCalculator::OptimizeCost(char ** haplotypes, int count, int markers)
   {
   LongCounter uniqueHaplotypes;

   int bits  = sizeof(long long) * 8 - 1;
   int limit = count / 2;

   long long * individualHaplotypes = new long long [count];
   int ** uniqueCounts = AllocateIntMatrix(markers, bits);

   // First we construct a matrix with the number of unique haplotypes
   // along different portions of the current solution
   for (int i = 0; i < markers; i++)
      {
      uniqueHaplotypes.Clear();
      // Retrieve one marker haplotypes
      for (int j = 0; j < count; j++)
         {
         individualHaplotypes[j] = haplotypes[j][i];
         uniqueHaplotypes.IncrementCount(individualHaplotypes[j]);
         }

      // Count the number of unique haplotypes
      uniqueCounts[i][0] = uniqueHaplotypes.Entries();

      for (int j = 1; j < bits; j++)
         {
         if (uniqueHaplotypes.Entries() > limit || i + j >= markers)
            {
            uniqueCounts[i][j] = count;
            continue;
            }

         uniqueHaplotypes.Clear();

         for (int k = 0; k < count; k++)
            {
            individualHaplotypes[k] = individualHaplotypes[k] * 2 + haplotypes[k][i+j];
            uniqueHaplotypes.IncrementCount(individualHaplotypes[k]);
            }

         uniqueCounts[i][j] = uniqueHaplotypes.Entries();
         }
      }

   // Finally, we use dynamic programming to find the best cost path
   if (cost != NULL) delete [] cost;
   if (path != NULL) delete [] path;

   cost = new double [markers];
   path = new int [markers];

   cost[0] = BasicCost(count);
   path[0] = 0;

   for (int i = 1; i < markers; i++)
      {
      cost[i] = cost[i - 1] + BasicCost(count);
      path[i] = i;

      for (int j = 1; j < bits; j++)
         if (i - j >= 0)
            {
            double alternate_cost = cost[i - j]
                        + ReducedCost(uniqueCounts[i-j][j], count) * j
                        + TranslationCost(count);

            if (alternate_cost < cost[i])
               {
               cost[i] = alternate_cost;
               path[i] = i - j;
               }
            }
      }

#ifdef _DEBUG
   // Estimate savings
   double naiveCost = markers * BasicCost(count);

   printf("Reduced state space optimization would result in a speedup of about %.1f-fold\n",
           naiveCost / cost[markers-1]);

   int position = markers - 1;
   printf("  Optimal Path: ");
   while (position != 0)
      if (path[position] == position)
         printf(" Step(%d) ", position--);
      else
         {
         printf(" Condense(%d -> %d)", position, path[position]);
         position = path[position];
         }
   printf("\n");
#endif

   // We are all done, free memory
   FreeIntMatrix(uniqueCounts, markers);
   delete [] individualHaplotypes;
   }
 
