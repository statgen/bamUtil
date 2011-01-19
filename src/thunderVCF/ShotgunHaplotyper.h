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
      
      char   * refalleles;
      double * freq1s;

      virtual void CalculateWeights();
      virtual void RandomSetup(Random * rand = NULL);
      virtual void LoadHaplotypesFromVCF(String& fileName);
      virtual void ConditionOnData(float * matrix, int marker, 
                                   char phred11, char phred12, char phred22);
      virtual void ImputeAlleles(int marker, int state1, int state2, Random * rand);
      
      virtual bool AllocateMemory(int nIndividuals, int maxHaplos, int nMarkers);
      virtual void EstimateMemoryInfo(int Individuals, int Markers, int States, bool Compact, bool Phased);
      virtual void ScoreLeftConditional();
      
      virtual void SampleChromosomes(Random * rand);
      virtual bool ForceMemoryAllocation();
      
      
      virtual void RetrieveMemoryBlock(int marker);
      
      //void   SetShotgunError(double rate);
      //double GetShotgunError() { return shotgunError; }
      void   CalculatePhred2Prob();

   protected:
      //float ** shotgunErrorMatrix;
      //float    shotgunError;
      bool     phred2probAvailable;
      double * phred2prob;
   };

 
