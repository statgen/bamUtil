////////////////////////////////////////////////////////////////////// 
// mach1/MergeHaplotypes.h 
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
 
#ifndef __MERGE_HAPLOTYPES__
#define __MERGE_HAPLOTYPES__

class ConsensusBuilder
   {
   public:
      int stored;

      int samples;
      int haplotypes;
      int markers;

      char *** sampledHaplotypes;
      char **  consensus;

      bool readyForUse;

      ConsensusBuilder(int samples, int haplotypes, int markers);
      ~ConsensusBuilder();

      void Store(char ** newHaplotypes);
      void Merge();

      // Quality scores for estimated haplotypes
      float flips;
      float errors;

      // Report memory usage
      void ShowMemoryInfo();
      static void EstimateMemoryInfo(int Samples, int Haplotypes, int Markers);

   private:
   };


#endif


 
