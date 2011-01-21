////////////////////////////////////////////////////////////////////// 
// mach1/AssociationAnalysis.h 
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
 
#ifndef __ASSOCIATIONANALYSIS_H__
#define __ASSOCIATIONANALYSIS_H__

#include "Pedigree.h"
#include "Haplotyper.h"
#include "DosageCalculator.h"

class AssociationAnalysis
   {
   public:
      static void ScoreNPL(String & prefix, Pedigree & ped, Haplotyper & engine, int rounds);
      static void ScoreMarkers(String & prefix, Pedigree & ped, DosageCalculator & doses);

   private:
      static void OutputPValue(FILE * output, double t, double df);
   };

#endif


 
