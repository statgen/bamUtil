////////////////////////////////////////////////////////////////////// 
// thunder/ShotgunHaplotyper.h 
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

class ShotgunHaplotyper : public Haplotyper
   {
   public:
      ShotgunHaplotyper();
      ~ShotgunHaplotyper();

      virtual void CalculateWeights();
      virtual void RandomSetup(Random * rand = NULL);
      virtual void ConditionOnData(float * matrix, int marker, char genotype);
      virtual void ImputeAlleles(int marker, int state1, int state2, Random * rand);

      void   SetShotgunError(double rate);
      double GetShotgunError() { return shotgunError; }

      static void EstimateMemoryInfo(int Individuals, int Markers, int States, bool Compact)
         {
         Haplotyper::EstimateMemoryInfo(Individuals, Markers, States, Compact, false);
         }

   protected:
      float ** shotgunErrorMatrix;
      float    shotgunError;
   };

 
