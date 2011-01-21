////////////////////////////////////////////////////////////////////// 
// mach1/DosageCalculator.h 
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
 
#ifndef __DOSAGECALCULATOR_H__
#define __DOSAGECALCULATOR_H__

#include <stdio.h>

class DosageCalculator
   {
   public:
      bool readyForUse;

      static bool storeDosage;
      static bool storeDistribution;

      DosageCalculator(int samples, int genotypes, int markers);
      ~DosageCalculator();

      void Update(char ** newHaplotypes);

      void GetCounts(int individual, int genotype,
                     unsigned int & n0, unsigned int & n1, unsigned int & n2);

      double GetDosage(int individual, int marker);
      double GetQuality(int individual, int marker);
      int    GetBestGenotype(int individual, int marker);

      void   CalculateMarkerInfo(int marker, double& freq, double& maf, double& avgPost, double& rsq);
      void   OutputMarkerInfo(const char * filename);
      void   OutputMarkerInfo(FILE * output);
      void   OutputBasicMarkerInfo(FILE * output);

      static void EstimateMemoryInfo(int samples, int genotypes, int markers);
      void ShowMemoryInfo();

   private:
      int wordSize;

      unsigned char **   cDosage;
      unsigned short **  sDosage;
      unsigned int **    iDosage;

      unsigned char **   cTwo;
      unsigned short **  sTwo;
      unsigned int **    iTwo;

      int samples;
      int genotypes;
      int markers;

      int stored;
   };


#endif

 
