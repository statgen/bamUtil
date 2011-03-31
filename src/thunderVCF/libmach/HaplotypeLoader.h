////////////////////////////////////////////////////////////////////// 
// mach1/HaplotypeLoader.h 
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
 
#ifndef __HAPLOTYPE_LOADER_H__
#define __HAPLOTYPE_LOADER_H__

#include "Pedigree.h"

#include <stdio.h>

class HaplotypeLoader
   {
   public:
      HaplotypeLoader();
      ~HaplotypeLoader();

      static bool hapmapFormat;
      static bool autoFlip;
      static bool loadPositions;

      void LoadMarkerList(const char * filename);
      void LoadMarkerList(IFILE file);
      void LoadHapMapLegendFile(IFILE file);

      void LoadHaplotypes(const char * filename);
      void LoadHaplotypes(IFILE file);
      void LoadHapMapHaplotypes(IFILE file);

      void WriteMarkerList(const char * filename, int from = -1, int to = -1);
      void WriteMarkerList(FILE * file, int from = -1, int to = -1);

      void WriteHaplotypes(const char * filename, int from = -1, int to = -1);
      void WriteHaplotypes(FILE * file, int from = -1, int to = -1);

      void ConsistencyCheck(Pedigree & ped);
      bool RenameAlleles(Pedigree & ped, int marker);
      bool FixStrand(Pedigree & ped, int marker);

      int     markerCount;
      int     count;
      char ** haplotypes;

      void ShowMemoryInfo();

   private:
      double square(double x)
         { return x * x; }

      bool FlipAllele(Pedigree & ped, int marker, int allele);
      bool FlipAlleles(Pedigree & ped, int marker, int allele1, int allele2);

      const char * FlipAllele(String & allele1);
   };

#endif
 
