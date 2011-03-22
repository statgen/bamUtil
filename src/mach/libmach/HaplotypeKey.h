////////////////////////////////////////////////////////////////////// 
// mach1/HaplotypeKey.h 
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
 
#ifndef __HAPLOTYPE_KEY__
#define __HAPLOTYPE_KEY__

class HaplotypeHash
   {
   public:
      int *  codes;
      int *  counts;
      int *  extras;

      int    size;
      int    count;

      HaplotypeHash(int size = 257);
      ~HaplotypeHash();

      bool   IncrementCount(int code);
      bool   DecrementCount(int code);
      int    GetPosition(int code);

      bool   SlotIsEmpty(int slot)
               { return counts[slot] == 0; }

      void   Clear();

   private:
      void   InsertAtSlot(int slot, int code)
         {
         codes[slot] = code;
         counts[slot]++;
         count++;
         }
   };

class HaplotypeKey
   {
   public:
      int * codes;
      int * map;
      int   count;
      int   from, to;

      int   unique;
      int * haplotypes;
      int * counts;

      HaplotypeKey(int haplotypeCount, int maxUnique);
      ~HaplotypeKey();

      void  Initialize(char ** haplotypes, int start, int stop);

      void  Clear();

      // Extend the current encoding by adding one marker
      void  ExtendCodes(char ** haplotypes);

      // Or trim it by removing the last marker
      void  TrimCodes();

      // List all unique haplotypes and map each haplotype to an id
      void  BuildMap();
      void  HashCodes();

      // Transform haplotype into a numeric code
      int   EncodeHaplotype(char * haplotype);

      // Replace haplotype at a specific slot
      void  ReplaceHaplotype(int slot, char * haplotype);

   private:
      HaplotypeHash hash;

      int    maxUnique;
   };

void ParseHaplotypes(char ** haplotypes, int count, int markers, int threshold);

#endif


 
