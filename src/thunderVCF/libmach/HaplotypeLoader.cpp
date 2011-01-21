////////////////////////////////////////////////////////////////////// 
// mach1/HaplotypeLoader.cpp 
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
 
#include "HaplotypeLoader.h"
#include "MemoryAllocators.h"
#include "MemoryInfo.h"
#include "Error.h"

#include <ctype.h>

#ifndef ZEPS
#define ZEPS 1e-30
#endif

bool HaplotypeLoader::hapmapFormat = false;
bool HaplotypeLoader::autoFlip = false;
bool HaplotypeLoader::loadPositions = false;

HaplotypeLoader::HaplotypeLoader()
   {
   haplotypes = NULL;
   markerCount = 0;
   count = 0;
   }

HaplotypeLoader::~HaplotypeLoader()
   {
   if (haplotypes != NULL)
      FreeCharMatrix(haplotypes, count + 1);
   }

void HaplotypeLoader::LoadMarkerList(const char * filename)
   {
   IFILE f = ifopen(filename, "rb");

   if (f == NULL)
      return; // error("Marker list [%s] could not be opened\n", filename);

   LoadMarkerList(f);
   ifclose(f);
   }

void HaplotypeLoader::LoadHaplotypes(const char * filename)
   {
   IFILE f = ifopen(filename, "rb");

   if (f == NULL)
      {
      if (Pedigree::markerCount)
         error("File [%s] with phased haplotypes could not be opened\n", filename);

      return;
      }

   LoadHaplotypes(f);
   ifclose(f);
   }

void HaplotypeLoader::LoadMarkerList(IFILE file)
   {
   if (hapmapFormat)
      {
      LoadHapMapLegendFile(file);
      return;
      }

   String      buffer;
   StringArray tokens;

   printf("Loading list of markers in phased haplotype ...\n");

   while (!ifeof(file))
      {
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (tokens.Length() == 0) continue;

      if (tokens.Length() > 1)
         error("Each line should list exactly one marker name, but the\n"
               "following line appears to include extra information:\n\n"
               "%s\n%s",
               (const char *) buffer,
               tokens.Length() != 4 ? "" : "\n"
               "If you are using a HapMap-style legend file, remember to\n"
               "use the --hapmapFormat command line option.\n\n");

      int markerId = Pedigree::GetMarkerID(tokens[0]);

      if (markerCount++ != markerId)
         error("Marker %s is duplicated.\n\n"
               "Every marker should have a unique name.\n",
               (const char *) tokens[0]);
      }
   }

void HaplotypeLoader::LoadHapMapLegendFile(IFILE file)
   {
   String      buffer;
   StringArray tokens;

   printf("Loading HapMap-style legend file ...\n");

   buffer.ReadLine(file);
   while (!ifeof(file))
      {
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (tokens.Length() == 0) continue;

      if (tokens.Length() != 4)
         error("Each line should list the marker name, position and alleles,\n"
               "but the following line includes %d items (instead of 4):\n\n"
               "%s\n", tokens.Length(), (const char *) buffer);

      int markerId = Pedigree::GetMarkerID(tokens[0]);

      if (markerCount++ != markerId)
         error("Marker %s is duplicated.\n\n"
               "Every marker should have a unique name.\n",
               (const char *) tokens[0]);

      MarkerInfo * info = Pedigree::GetMarkerInfo(markerId);

      info->NewAllele(tokens[2]);
      info->NewAllele(tokens[3]);

      if (!loadPositions) continue;
    
      info->position = tokens[1].AsDouble();
      }
   }

void HaplotypeLoader::LoadHaplotypes(IFILE file)
   {
   if (hapmapFormat)
      {
      LoadHapMapHaplotypes(file);
      return;
      }

   // Don't load haplotypes unless we have a marker list
   if (markerCount == 0)
      {
      printf("  WARNING -- Since no marker list was provided, haplotype file will be ignored\n\n");
      return;
      }

   printf("Loading phased haplotypes ...\n");

   String      buffer;
   StringArray tokens;

   // In the first pass, we simply count the number of non-blank lines
   while (!ifeof(file))
      {
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (!tokens.Length()) continue;

      count++;
      }

   // Check if we got some valid input
   if (count == 0 || markerCount == 0)
      return;

   // Then, we allocate memory for storing the phased haplotypes
   haplotypes = AllocateCharMatrix(count + 1, Pedigree::markerCount);

   // And finally, we load the data in a second pass
   ifrewind(file);

   int line = 0, index = 0;
   while (!ifeof(file))
      {
      line++;
      buffer.ReadLine(file);
      tokens.ReplaceTokens(buffer);

      if (tokens.Length() == 0) continue;

      int hapstart = tokens.Length() - 1;
      int offset = markerCount;

      while ((offset -= tokens[hapstart].Length()) > 0 && hapstart > 0)
         hapstart--;

      if (offset != 0)
         error("The haplotype file format was not recognized\n"
               "(Problem occured reading haplotype #%d in line #%d)\n\n"
               "Check that the number of markers matches the SNPs list\n",
               ++line, ++index);

      for (int i = 0; i < markerCount; i++)
         {
         MarkerInfo * info = Pedigree::GetMarkerInfo(i);

         if (offset == tokens[hapstart].Length())
            {
            offset = 0;
            hapstart++;
            }

         int al = info->GetAlleleNumber(tokens[hapstart][offset++]);

         if (al == -1)
            al = info->NewAllele(tokens[hapstart][offset-1]);

         if (al == 0)
            error("Missing data in haplotype %d at position %d (marker %s)\n",
                   ++index, ++i, (const char *) Pedigree::markerNames[i]);

         if (al > 2)
            error("More than 2 alleles at position %d (marker %s)\n",
                   ++i, (const char *) Pedigree::markerNames[i]);

         haplotypes[index][i] = al;
         }
      index++;
      }

   for (int i = 0; i < markerCount; i++)
      {
      MarkerInfo * info = Pedigree::GetMarkerInfo(i);

      // The last row in the haplotypes file stores the number of distinct
      // alleles observed in the input files
      haplotypes[count][i] = info->CountAlleles();
      }
   }

void HaplotypeLoader::LoadHapMapHaplotypes(IFILE file)
   {
   // Don't load haplotypes unless we have a marker list
   if (markerCount == 0)
      {
      printf("  WARNING -- Since no legend file was provided, haplotype file will be ignored\n\n");
      return;
      }

   printf("Loading HapMap-style phased haplotypes ...\n");

   String      buffer;
   StringArray tokens;

   // In the first pass, we simply count the number of non-blank lines
   while (!ifeof(file))
      {
      buffer.ReadLine(file);
      buffer.Trim();

      if (buffer.Length() != 2 * markerCount - 1)
         {
         if (buffer.Length())
            error("According to the legend file, there should be %d alleles per haplotype.\n"
                  "However, some lines have an unexpected character count\n", markerCount);

         continue;
         }

      count++;
      }

   // Check if we got some valid input
   if (count == 0 || markerCount == 0)
      return;

   // Then, we allocate memory for storing the phased haplotypes
   haplotypes = AllocateCharMatrix(count + 1, Pedigree::markerCount);

   // And finally, we load the data in a second pass
   ifrewind(file);

   int line = 0, index = 0;
   while (!ifeof(file))
      {
      line++;
      buffer.ReadLine(file);
      buffer.Trim();

      if (buffer.Length() != 2 * markerCount - 1) continue;

      bool badchar = false;
      for (int i = 0; i < markerCount; i++)
         {
         if (buffer[i * 2] == '0')
            haplotypes[index][i] = 1;
         else if (buffer[i * 2] == '1')
            haplotypes[index][i] = 2;
         else
            badchar = true;

         if (badchar || buffer[i * 2 + 1] != ' ' && buffer[i * 2 + 1] != 0)
            error("Haplotype file should include a series of '0's and '1's,\n"
                  "separated by spaces. However, an unexpected character was\n"
                  "encountered in line %d.\n", line);
          }
      index++;
      }

   for (int i = 0; i < markerCount; i++)
      haplotypes[count][i] = 2;
   }

void HaplotypeLoader::ConsistencyCheck(Pedigree & ped)
   {
   if (count == 0 || markerCount == 0)
      return;

   if (markerCount != Pedigree::markerCount)
      {
      printf("The following markers appear in the pedigree, but not in the phased haplotypes:");

      int skipped_markers = 0;
      for (int i = markerCount, line = 80, lines = 0; i < Pedigree::markerCount; i++)
         if (lines < 10)
            {
            if (line + Pedigree::markerNames[i].Length() + 1 > 79)
               if (lines == 9)
                  {
                  printf("\n");
                  lines++;
                  skipped_markers++;
                  continue;
                  }
               else
                  printf("\n   "), line = 3, lines++;

            printf("%s ", (const char *) Pedigree::markerNames[i]);
            line += Pedigree::markerNames[i].Length() + 1;
            }
         else
            skipped_markers++;

      if (skipped_markers)
         printf("These %d markers and %d other%s (%d total) will be ignored\n\n",
                Pedigree::markerCount - markerCount - skipped_markers,
                skipped_markers, skipped_markers == 1 ? "" : "s",
                Pedigree::markerCount - markerCount);
      else
         printf("These %d markers will be ignored\n", Pedigree::markerCount - markerCount);

      Pedigree::markerCount = markerCount;
      }

   bool warnings = false;
   bool errors = false;
   int numbersToBases = 0;

   for (int i = 0; i < markerCount; i++)
      {
      bool bad_marker = false;
      MarkerInfo * info = Pedigree::GetMarkerInfo(i);

      if (autoFlip && info->CountAlleles() > haplotypes[count][i])
          if (RenameAlleles(ped, i))
             numbersToBases++;

      // The last row in the haplotypes table stores the number of distinct alleles
      // observed in the input files
      if (info->CountAlleles() > haplotypes[count][i] && info->CountAlleles() > 2)
         {
         if (autoFlip && FixStrand(ped, i))
            printf("Fixed alleles for marker %s ... ", (const char *) info->name);
         else
            {
            printf("Mismatched alleles for marker %s ... ", (const char *) info->name);
            errors |= info->CountAlleles() > 2;
            bad_marker = info->CountAlleles() > 2;
            }

         printf("Phased Haps: [%s", (const char *) info->GetAlleleLabel(1));
         for (int j = 2; j < haplotypes[count][i] + 1; j++)
            printf(",%s", (const char *) info->GetAlleleLabel(j));
         printf("]  Pedigree: [%s", (const char *) info->GetAlleleLabel(haplotypes[count][i]+1));
         for (int j = haplotypes[count][i] + 2; j < info->CountAlleles() + 1; j++)
            printf(",%s", (const char *) info->GetAlleleLabel(j));
         printf("]\n");

         warnings = true;
         }

      if (bad_marker)
         {
         if (autoFlip)
            {
            printf("   Genotypes for marker %s will be discarded\n", (const char *) info->name);

            for (int j = 0; j < ped.count; j++)
               ped[j].markers[i][0] = ped[j].markers[i][1] = 0;
            }
         }

      // We do a final sanity check to see if the allele frequencies are
      // similar in the pedigree to be phased in the haplotype file
      int hapCount[2] = {0, 0};

      for (int j = 0; j < count; j++)
         hapCount[haplotypes[j][i] - 1]++;

      int pedCount[2] = {0, 0};

      for (int j = 0; j < ped.count; j++)
         if (ped[j].markers[i].isKnown())
            {
            pedCount[ped[j].markers[i][0] - 1]++;
            pedCount[ped[j].markers[i][1] - 1]++;
            }

      double pedTotal = pedCount[0] + pedCount[1];

      if (pedCount[0] + pedCount[1] <= 1) continue;

      double freq   = (pedCount[0] + hapCount[0]) / (count + pedTotal);
      double offset = square(hapCount[0] - freq * count);
      double chisq  = offset / (count * freq + ZEPS) +
                      offset / (count * (1.0 - freq) + ZEPS) +
                      offset / (pedTotal * (1.0 - freq) + ZEPS) +
                      offset / (pedTotal * freq + ZEPS);

      // Should only be exceed about 1/10,000 tries
      if (chisq > 15.13)
         {
         printf("Warning: Allele %s (at %s) has frequency %f in phased haplos, but %f in the sample\n",
                (const char *) info->GetAlleleLabel(1),
                (const char *) Pedigree::markerNames[i],
                hapCount[0] / (count + ZEPS), pedCount[0] / (pedTotal + ZEPS));
         warnings = true;
         }

      // TODO -- we should probably calculate Fst to compare the two
      // sets of haplotypes ...
      }

   if (numbersToBases)
      printf("Numeric labels converted to bases at %d markers...\n", numbersToBases);

   if (errors && !autoFlip)
      error("Please ensure that allele labels in pedigree are consistent with haplotype file\n");

   if (warnings)
      printf("\n");
   }

void HaplotypeLoader::ShowMemoryInfo()
   {
   if (count == 0) return;

   int bytes = markerCount * count * sizeof(char);

   printf("   %40s %s\n", "Phase known haplotypes ...", (const char *) MemoryInfo(bytes));
   }

// The function below converts numeric allele labels to base-pairs
// (assuming 1,2,3,4 match A,C,G,T)

bool HaplotypeLoader::RenameAlleles(Pedigree & ped, int marker)
    {
    int base_alleles = haplotypes[count][marker];
    int total_alleles = ped.CountAlleles(marker);

    MarkerInfo * info = ped.GetMarkerInfo(marker);

    // Only apply fix to two alleles at a time
    if (total_alleles - base_alleles > 2 || total_alleles > 4)
        return false;

    // Abort unless all old alleles are labeled as bases
    for (int i = 1; i <= base_alleles; i++)
      if (info->alleleLabels[i].Length() != 1 ||
          toupper(info->alleleLabels[i][0]) != 'A' &&
          toupper(info->alleleLabels[i][0]) != 'C' &&
          toupper(info->alleleLabels[i][0]) != 'G' &&
          toupper(info->alleleLabels[i][0]) != 'T')
          continue;

    // Abort unless all new alleles are named as numbers
    for (int i = base_alleles + 1; i <= total_alleles; i++)
        if (info->alleleLabels[i].Length() != 1 ||
            info->alleleLabels[i][0] < '1' ||
            info->alleleLabels[i][0] > '4')
            return false;

    StringArray newLabels;
    StringIntHash newNumbers;

    newLabels.Push("");
    for (int i = 1; i <= base_alleles; i++)
        {
        newLabels.Push(info->alleleLabels[i]);
        newNumbers.Add(info->alleleLabels[i], i);
        }

    int  nextAllele = base_alleles + 1;
    int  rename[5] = {0, 0, 0, 0, 0};
    const char * bases[4] = {"A", "C", "G", "T"};

    // Next, generate a new set of allele labels and an appropriate index
    for (int i = base_alleles + 1; i <= total_alleles; i++)
        {
        int base = info->alleleLabels[i][0] - '1';

        int newNumber = newNumbers.Integer(bases[base]);

        if (newNumber > 0)
           {
           rename[i] = newNumber;
           continue;
           }

        newLabels.Push(bases[base]);
        newNumbers.Add(newLabels.Last(), nextAllele);
        rename[i] = nextAllele++;
        }

   // Finally, apply the renaming filter to the rest of the pedigree
   for (int i = 0; i < ped.count; i++)
      {
      if (rename[ped[i].markers[marker][0]])
        ped[i].markers[marker][0] = rename[ped[i].markers[marker][0]];
      if (rename[ped[i].markers[marker][1]])
        ped[i].markers[marker][1] = rename[ped[i].markers[marker][1]];
      }

   info->alleleLabels = newLabels;
   info->alleleNumbers = newNumbers;

   return true;
   }

// The code below tries to automatically fix allele flips, but is not
// enabled by default -- I am not convinced this is a good idea??

bool HaplotypeLoader::FixStrand(Pedigree & ped, int marker)
   {
   int base_alleles = haplotypes[count][marker];

   for (int i = 0; i < ped.count; i++)
      if (ped[i].markers[marker].isKnown())
         if (ped[i].markers[marker].Lo() <= base_alleles ||
             ped[i].markers[marker].Hi() >  base_alleles + 2)
            return false;

   MarkerInfo * info = ped.GetMarkerInfo(marker);

   if (info->CountAlleles() == base_alleles + 1)
      return FlipAllele(ped, marker, base_alleles + 1);

   if (info->CountAlleles() == base_alleles + 2)
      return FlipAlleles(ped, marker, base_alleles + 1, base_alleles + 2);

   return false;
   }

bool HaplotypeLoader::FlipAllele(Pedigree & ped, int marker, int al1)
   {
   MarkerInfo * info = Pedigree::GetMarkerInfo(marker);
   String label1 = info->GetAlleleLabel(al1);

   if (label1.Length() != 1)
      return false;

   String flip1 = FlipAllele(label1);

   int flipped1 = info->GetAlleleNumber(flip1);

   if (flipped1 >= al1 || flipped1 < 1)
      return false;

   for (int i = 0; i < ped.count; i++)
      for (int j = 0; j < 2; j++)
         if (ped[i].markers[marker][j] == al1)
            ped[i].markers[marker][j] = flipped1;

   return true;
   }

bool HaplotypeLoader::FlipAlleles(Pedigree & ped, int marker, int al1, int al2)
   {
   MarkerInfo * info = Pedigree::GetMarkerInfo(marker);
   String label1 = info->GetAlleleLabel(al1);
   String label2 = info->GetAlleleLabel(al2);

   if (label1.Length() != 1 || label2.Length() != 1)
      return false;

   String flip1 = FlipAllele(label1);
   String flip2 = FlipAllele(label2);

   int flipped1 = info->GetAlleleNumber(flip1);
   int flipped2 = info->GetAlleleNumber(flip2);

   if (flipped1 == al2)
      return false;

   if (flipped1 > 2 || flipped2 > 2 || flipped1 < 1 || flipped2 < 1)
      return false;

   for (int i = 0; i < ped.count; i++)
      for (int j = 0; j < 2; j++)
      {
      if (ped[i].markers[marker][j] == al1)
         ped[i].markers[marker][j] = flipped1;
      else if (ped[i].markers[marker][j] == al2)
         ped[i].markers[marker][j] = flipped2;
      }

   return true;
   }

const char * HaplotypeLoader::FlipAllele(String & allele)
   {
   static const char * flips[4] = {"A", "C", "G", "T"};

   if (allele.Length() != 1)
      return "";

   switch (allele[0])
      {
      case 'A': case 'a':
         return flips[3];
      case 'C': case 'c':
         return flips[2];
      case 'G': case 'g':
         return flips[1];
      case 'T': case 't':
         return flips[0];
      default:
         return "";
      }
   }

void HaplotypeLoader::WriteMarkerList(const char * filename, int from, int to)
   {
   FILE * f = fopen(filename, "wb");

   if (f == NULL)
      return; // error("Marker list [%s] could not be opened\n", filename);

   WriteMarkerList(f, from, to);
   fclose(f);
   }

void HaplotypeLoader::WriteHaplotypes(const char * filename, int from, int to)
   {
   FILE * f = fopen(filename, "wb");

   if (f == NULL)
      return; // error("Marker list [%s] could not be opened\n", filename);

   WriteHaplotypes(f, from, to);
   fclose(f);
   }

void HaplotypeLoader::WriteMarkerList(FILE * output, int from, int to)
   {
   if (from == -1) from = 0;
   if (to == -1) to = count;

   if (hapmapFormat)
      {
      fprintf(output, "rs\tposition\t0\t1\n");
      for (int i = from; i < to; i++)
          {
          MarkerInfo * info = Pedigree::GetMarkerInfo(i);

          fprintf(output, "%s\t%d\t%s\t%s\n", 
                          (const char *) info->name, info->position > 0 ? (int) info->position : i+1,
                          (const char *) info->GetAlleleLabel(1),
                          (const char *) info->GetAlleleLabel(2));
          }
      }
   else
      {
      for (int i = from; i < to; i++)
          fprintf(output, "%s\n", (const char *) Pedigree::markerNames[i]);
      }
   }

void HaplotypeLoader::WriteHaplotypes(FILE * output, int from, int to)
   {
   if (from == -1) from = 0;
   if (to == -1) to = count;

   if (hapmapFormat)
      for (int i = 0; i < count; i++)
          {
          for (int j = from; j < to; j++)
               fprintf(output, "%d ", haplotypes[i][j] - 1); 
          fprintf(output, "\n");
          }
   else
      for (int i = 0; i < count; i++)
          {
          printf("HAP%d ", count);
          for (int j = from; j < to; j++)
               fprintf(output, "%s ", (const char *) Pedigree::GetMarkerInfo(i)->alleleLabels[haplotypes[i][j]]); 
          fprintf(output, "\n");
          }
   }
 
