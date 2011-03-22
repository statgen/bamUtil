////////////////////////////////////////////////////////////////////// 
// mach1/CostCalculator.h 
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
 
#ifndef __COSTCALCULATOR_H__
#define __COSTCALCULATOR_H__

class CostCalculator
   {
   public:
      CostCalculator();
      ~CostCalculator();

      void OptimizeCost(char ** haplotypes, int count, int markers);

      double * cost;
      int *    path;

   protected:
      double BasicCost(int count)
         { return count * count; }

      double ReducedCost(int unique, int count)
         { return 4.0 * unique * unique + unique * count; }

      double TranslationCost(int count)
         { return 2.0 * count * count; }
   };

#endif

 
