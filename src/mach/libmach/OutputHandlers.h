////////////////////////////////////////////////////////////////////// 
// mach1/OutputHandlers.h 
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
 
#ifndef __OUTPUTHANDLERS_H__
#define __OUTPUTHANDLERS_H__

#include "MergeHaplotypes.h"
#include "DosageCalculator.h"
#include "Pedigree.h"
#include "InputFile.h"

class OutputManager
   {
   public:
      static void WriteHaplotypes(const char * outfile, Pedigree & ped, char ** haplotypes);
      static void WriteQuality(const char * outfile, Pedigree & ped, short ** matrix);
      static void WriteGenotypes(const char * outfile, Pedigree & ped, DosageCalculator & doses);
      static void WriteDosages(const char * outfile, Pedigree & ped, DosageCalculator & doses);
      static void WriteProbabilities(const char * outfile, Pedigree & ped, DosageCalculator & doses);
      static void WriteQuality(const char * outfile, Pedigree & ped, DosageCalculator & doses);

      static void OutputConsensus(Pedigree & ped, ConsensusBuilder & consensus, DosageCalculator & doses, String filename);

      static bool outputHaplotypes;
      static bool outputGenotypes;
      static bool outputDosage;
      static bool outputProbabilities;
      static bool outputQuality;
   };


#endif

 
