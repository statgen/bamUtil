////////////////////////////////////////////////////////////////////// 
// mach1/HaplotypeKey.cpp 
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
 
#include "HaplotypeKey.h"
#include "Error.h"

#include <stdio.h>

HaplotypeHash::HaplotypeHash(int SIZE)
   {
   size = SIZE;
   count = 0;

   extras = new int [size];
   counts = new int [size];
   codes = new int [size];
   }

HaplotypeHash::~HaplotypeHash()
   {
   delete [] extras;
   delete [] counts;
   delete [] codes;
   }

void HaplotypeHash::Clear()
   {
   for (int i = 0; i < size; i++)
      counts[i] = 0;
   count = 0;
   }

int HaplotypeHash::GetPosition(int code)
   {
   int h = code % size;

   while (true)
      {
      // Slot is at the correct position
      if (codes[h] == code)
         return h;

      // Rehashing is required ...
      h++;

      if (h == size) h = 0;
      }
   }

bool HaplotypeHash::IncrementCount(int code)
   {
   int h = code % size;

   while (true)
      {
      // Slot is empty
      if (counts[h] == 0)
         {
         InsertAtSlot(h, code);
         return true;
         }

      // Slot is at the correct position
      if (codes[h] == code)
         {
         counts[h]++;
         return false;
         }

      // Rehashing is required ...
      h++;

      if (h == size) h = 0;
      }
   }

bool HaplotypeHash::DecrementCount(int code)
   {
   int h = code % size;

   while (true)
      {
      // Slot is at the correct position
      if (codes[h] == code)
         {
         counts[h]--;

         // Entry was deleted ...
         if (counts[h] == 0)
            {
            count--;

            // Rehash subsequent entries, as necessary
            for (int rehash = h + 1; rehash != h; rehash == rehash + 1 == size ? 0 : rehash + 1)
               {
               // Done rehashing when an empty slot is reached
               if (counts[rehash] == 0) break;

               // Find out new position for the code
               int h = codes[rehash] % size;

               // If the same as current position, nothing to do
               if (h == rehash) continue;

               // Otherwise, bubble it up to an appropriate slot
               while (counts[h] != 0 && h != rehash)
                  if (++h == size)
                     h = 0;

               // Move entry if required
               if (h != rehash)
                  {
                  counts[h] = counts[rehash];
                  codes[h] = codes[rehash];

                  counts[rehash] = 0;
                  }
               }

            return true;
            }

         return false;
         }

      // Rehashing is required ...
      h++;

      if (h == size) h = 0;
      }
   }

HaplotypeKey::HaplotypeKey(int haplotypeCount, int maximumUnique)
   {
   codes = new int [count = haplotypeCount];
   map = new int [count];

   haplotypes = new int [2 * (maxUnique = maximumUnique)];
   counts = new int [2 * maxUnique];
   }

HaplotypeKey::~HaplotypeKey()
   {
   delete [] codes;
   delete [] map;
   delete [] haplotypes;
   delete [] counts;
   }

void HaplotypeKey::Clear()
   {
   for (int i = 0; i < count; i++)
      codes[i] = 0;
   }

void HaplotypeKey::ExtendCodes(char ** haplotypes)
   {
   for (int i = 0; i < count; i++)
      codes[i] = codes[i] * 2 + haplotypes[i][to];
   to++;
   }

void HaplotypeKey::TrimCodes()
   {
   for (int i = 0; i < count; i++)
      codes[i] /= 2;
   to--;
   }

void HaplotypeKey::HashCodes()
   {
   hash.Clear();

   for (int i = 0; i < count; i++)
      hash.IncrementCount(codes[i]);
   }

void HaplotypeKey::BuildMap()
   {
   unique = 0;

   for (int i = 0; i < hash.size; i++)
      if (!hash.SlotIsEmpty(i))
         {
         hash.extras[i] = unique;
         haplotypes[unique] = hash.codes[i];
         counts[unique] = hash.counts[i];

         unique++;
         }

   for (int i = 0; i < count; i++)
      map[i] = hash.extras[hash.GetPosition(codes[i])];
   }

int HaplotypeKey::EncodeHaplotype(char * haplotype)
   {
   int code = 0;

   for (int i = from; i < to; i++)
      code = code * 2 + haplotype[i];

   return code;
   }

void HaplotypeKey::ReplaceHaplotype(int slot, char * haplotype)
   {
   int new_code = EncodeHaplotype(haplotype);

   if (new_code == codes[slot])
      return;

   if (hash.DecrementCount(codes[slot]))
      if (--unique > map[slot])
         {
         for (int i = 0; i < hash.size; i++)
            if (hash.extras[i] > map[slot])
               hash.extras[i]--;

         for (int i = 0; i < count; i++)
            if (map[i] > map[slot])
               map[i]--;
         }

   if (hash.IncrementCount(codes[slot] = new_code))
      {
      map[slot] = hash.extras[hash.GetPosition(new_code)] = unique;
      haplotypes[unique] = new_code;
      counts[unique] = new_code;

      unique++;

      if (unique >= maxUnique * 2)
         error("Haplotype hashing failed -- too many distinct haplotypes between"
               "markers %d and %d\n", from + 1, to + 1);
      }
   else
      map[slot] = hash.extras[hash.GetPosition(new_code)];
   }

void HaplotypeKey::Initialize(char ** haplotypes, int start, int stop)
   {
   from = to = start;

   int max = 1;

   // The first few markers are guaranteed to keep us under the complexity
   // limit
   Clear();
   while (max < maxUnique && to != stop)
      {
      ExtendCodes(haplotypes);
      max *= 2;
      }

   // After enough markers have been added, we must monitor further extensions
   HashCodes();
   while (hash.count <= maxUnique && from - to < 31 && to != stop)
      {
      ExtendCodes(haplotypes);
      HashCodes();
      }

   if (hash.count > maxUnique)
      {
      TrimCodes();
      HashCodes();
      }

   BuildMap();
   }

void ParseHaplotypes(char ** haplotypes, int count, int markers, int threshold)
   {
   HaplotypeKey key(count, threshold);

   key.to = -1;

   int blocks = 0;
   while (key.to < markers - 1)
      {
      key.Initialize(haplotypes, key.to + 1, markers);

      printf("Block from %5d - %5d :   %3d haplotypes\n", key.from, key.to, key.unique);

      blocks++;
      }

   printf("Region can be parsed into %d blocks\n\n", blocks);
   }
 
