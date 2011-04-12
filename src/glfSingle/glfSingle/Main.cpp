////////////////////////////////////////////////////////////////////// 
// glfSingle/Main.cpp 
// (c) 2000-2010 Goncalo Abecasis
// 
// This file is distributed as part of the Goncalo source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile Goncalo.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Thursday March 25, 2010
// 
 
#include "BaseQualityHelper.h"
#include "Parameters.h"
#include "glfHandler.h"
#include "IntArray.h"
#include "Error.h"

#include <math.h>
#include <time.h>

FILE * baseCalls = NULL;

void ReportHeader(const char * label)
   {
   if (baseCalls == NULL) return;

   fprintf(baseCalls, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", label);
   }

void ReportDate(FILE * output)
   {
   if (output == NULL) return;

   time_t systemTime;
   time(&systemTime);

   tm * digestedTime;
   digestedTime = gmtime(&systemTime);

   fprintf(output, "##filedate=%04d%02d%02d\n", digestedTime->tm_year + 1900,
                                                digestedTime->tm_mon + 1,
                                                digestedTime->tm_mday);
   }

void ReportSNP(glfHandler & glf, const char * filter, int allele1, int allele2, double posterior)
   {
   if (baseCalls == NULL) return;

   int refBase = glf.data.refBase;

   char alleles[] = { 0, 'a', 'c', 'g', 't' };

   // #Chrom   Pos   Id
   fprintf(baseCalls, "%s\t%d\t.\t", (const char *) glf.label, glf.position + 1);

   // Reference allele
   int nalleles = 1;
   fprintf(baseCalls, "%c\t", alleles[refBase]);

   // Other alleles
   if (allele1 != refBase)
      fprintf(baseCalls, "%c", alleles[allele1]), nalleles++;

   if (allele2 != refBase && allele2 != allele1)
      fprintf(baseCalls, "%s%c", nalleles > 1 ? "," : "", alleles[allele2]), nalleles++;

   if (nalleles == 1)
      fprintf(baseCalls, ".");

   fprintf(baseCalls, "\t");

   int quality = posterior > 0.9999999999 ? 100 : -10 * log10(1.0 - posterior);

   // Quality for this call
   fprintf(baseCalls, "%d\t", quality);

   // Filter for this call
   fprintf(baseCalls, "%s\t", filter == NULL || filter[0] == 0 ? "." : filter);

   // Information for this call
   fprintf(baseCalls, "%s\t", allele1 == allele2 ? "Hom=1" : "Het=1");

   // Format for this call
   fprintf(baseCalls, "%s\t", "GT:GD");

   // Genotype for this call
   fprintf(baseCalls, "%d/%d:%d\n",
            allele1 == refBase ? 0 : 1,
            (allele1 == allele2) ?
              (allele1 == refBase ? 0 : 1) :
              (allele2 == refBase ? 0 : nalleles - 1),
            glf.data.depth);
   }

int main(int argc, char ** argv)
   {
   String genofile, callfile, label;
   ParameterList pl;

   double posterior = 0.50;
   int    mapQuality = 60;
   int    minDepth = 0;
   int    maxDepth = 1000;
   bool   referenceCalls = false;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Map Quality Filter")
         LONG_INTPARAMETER("minMapQuality", &mapQuality)
      LONG_PARAMETER_GROUP("Total Depth Filter")
         LONG_INTPARAMETER("minDepth", &minDepth)
         LONG_INTPARAMETER("maxDepth", &maxDepth)
      LONG_PARAMETER_GROUP("Output")
         LONG_PARAMETER("reference", &referenceCalls)
   END_LONG_PARAMETERS();

   pl.Add(new StringParameter('g', "Genotype Likelihood File", genofile));
   pl.Add(new StringParameter('b', "Base Call File", callfile));
   pl.Add(new StringParameter('l', "Sample Label", label));
   pl.Add(new DoubleParameter('p', "Posterior Threshold", posterior));
   pl.Add(new LongParameters("Additional Options", longParameters));
   pl.Read(argc, argv);
   pl.Status();

   if (posterior < 0.50)
      error("Posterior threshold for genotype calls (-p option) must be > 0.50.");

   glfHandler glf;

   if (!glf.Open(genofile))
      error("Failed to open genotype likelihood file\n");

   baseCalls = fopen(callfile, "wt");

   if (baseCalls)
      {
      fprintf(baseCalls, "##fileformat=VCFv3.3\n");
      ReportDate(baseCalls);
      fprintf(baseCalls, "##source=glfSingle\n");
      fprintf(baseCalls, "##input=%s\n", (const char *) genofile);
      fprintf(baseCalls, "##minDepth=%d\n", minDepth);
      fprintf(baseCalls, "##maxDepth=%d\n", maxDepth);
      fprintf(baseCalls, "##minMapQuality=%d\n", mapQuality);
      fprintf(baseCalls, "##minPosterior=%.4f\n", posterior);
      fprintf(baseCalls, "##FILTER=mapQ,\"RMS Map Quality\"\n");
      fprintf(baseCalls, "##FILTER=depth,\"Read Depth\"\n");
      fprintf(baseCalls, "##FORMAT=GT,String,1,\"Genotype Call\"\n");
      fprintf(baseCalls, "##FORMAT=GQ,Integer,1,\"Genotype Call Quality\"\n");
      fprintf(baseCalls, "##FORMAT=DP,Integer,1,\"Read Depth\"\n");
      }

   ReportHeader(label.IsEmpty() ? genofile : label);

   IntArray coverage;
   IntArray target;

   while (glf.NextSection())
      {
      printf("Processing section %s with %d entries\n", (const char *) glf.label, glf.maxPosition);

      int actualBases = 0;
      int lowMapQuality = 0;
      int inadequateCoverage = 0;
      int homozygousReference = 0;
      int heterozygousPolymorphism = 0;
      int homozygousPolymorphism = 0;
      int baseCounts[5] = {0, 0, 0, 0, 0};

      coverage.Dimension(glf.maxPosition);
      coverage.Zero();
      target.Dimension(glf.maxPosition);
      target.Zero();

      String filter;

      while (glf.NextBaseEntry())
         {
         if (glf.position > glf.maxPosition)
            break;

         coverage[glf.position] = glf.data.depth;

         baseCounts[glf.data.refBase]++;

         if (glf.data.refBase == 0) continue;

         filter.Clear();

         actualBases++;

         // Check if we have good enough mapping quality
         if (glf.data.mapQuality < mapQuality)
            {
            filter.catprintf("%smapQ<%d", filter.Length() ? ";" : "", mapQuality);
            if (filter.Length() == 0) lowMapQuality++;
            }

         if (int(glf.data.depth) < minDepth)
            {
            filter.catprintf("%sdepth<%d", filter.Length() ? ";" : "", minDepth);
            if (filter.Length() == 0) inadequateCoverage++;
            }

         if (int(glf.data.depth) > maxDepth)
            {
            filter.catprintf("%sdepth>%d", filter.Length() ? ";" : "", maxDepth);
            if (filter.Length() == 0) inadequateCoverage++;
            }

         target[glf.position] = 1;

         // Create convenient aliases for each base
         unsigned char refBase = glf.data.refBase;
         unsigned char transition = (((glf.data.refBase - 1) ^ 2) + 1);
         unsigned char transvers1 = (((glf.data.refBase - 1) ^ 3) + 1);
         unsigned char transvers2 = (((glf.data.refBase - 1) ^ 1) + 1);

         int homRef = glf.GenotypeIndex(refBase, refBase);

         // Calculate partial likelihoods for a series of alternative genotypes
         double refTransition = 2./3. * glf.likelihoods[glf.GenotypeIndex(refBase, transition)];
         double refTransvers1 = 1./6. * glf.likelihoods[glf.GenotypeIndex(refBase, transvers1)];
         double refTransvers2 = 1./6. * glf.likelihoods[glf.GenotypeIndex(refBase, transvers2)];

         double homTransition = 2./3. * glf.likelihoods[glf.GenotypeIndex(transition, transition)];
         double homTransvers1 = 1./6. * glf.likelihoods[glf.GenotypeIndex(transvers1, transvers1)];
         double homTransvers2 = 1./6. * glf.likelihoods[glf.GenotypeIndex(transvers2, transvers2)];

         double transitiontv1 = 1./3. * glf.likelihoods[glf.GenotypeIndex(transition, transvers1)];
         double transitiontv2 = 1./3. * glf.likelihoods[glf.GenotypeIndex(transition, transvers2)];
         double transvers1tv2 = 1./3. * glf.likelihoods[glf.GenotypeIndex(transvers1, transvers2)];

         // Calculate likelihoods for three genotyping groupings
         double lRef = 0.9985 * glf.likelihoods[homRef];
         double lHet = 0.0010 * (refTransition + refTransvers1 + refTransvers2);
         double lHom = 0.0005 * (
            0.99975 * (homTransition + homTransvers1 + homTransvers2) +
            0.00025 * (transitiontv1 + transitiontv2 + transvers1tv2));

         double sum = lRef + lHet + lHom;

         double threshold = sum * posterior;

         if (lRef > threshold)
            {
            if (referenceCalls) ReportSNP(glf, filter, refBase, refBase, lRef / sum);
            if (filter.Length() == 0) homozygousReference++;
            }
         else if (lHet > threshold)
            {
            if (0.001 * refTransition > threshold)
               ReportSNP(glf, filter, refBase, transition, 0.001 * refTransition / sum); else
            if (0.001 * refTransvers1 > threshold)
               ReportSNP(glf, filter, refBase, transvers1, 0.001 * refTransvers1 / sum); else
            if (0.001 * refTransvers2 > threshold)
               ReportSNP(glf, filter, refBase, transvers2, 0.001 * refTransvers2 / sum); else
            continue;

            if (filter.Length() == 0) heterozygousPolymorphism++;
            }
         else if (lHom > threshold)
            {
            if (0.0005 * 0.99975 * homTransition > threshold)
               ReportSNP(glf, filter, transition, transition, 0.0005 * 0.99975 * homTransition / sum); else
            if (0.0005 * 0.99975 * homTransvers1 > threshold)
               ReportSNP(glf, filter, transvers1, transvers1, 0.0005 * 0.99975 * homTransvers1 / sum); else
            if (0.0005 * 0.99975 * homTransvers2 > threshold)
               ReportSNP(glf, filter, transvers2, transvers2, 0.0005 * 0.99975 * homTransvers2 / sum); else
            if (0.0005 * 0.00025 * transitiontv1 > threshold)
               ReportSNP(glf, filter, transition, transvers1, 0.0005 * 0.00025 * transitiontv1 / sum); else
            if (0.0005 * 0.00025 * transitiontv2 > threshold)
               ReportSNP(glf, filter, transition, transvers2, 0.0005 * 0.00025 * transitiontv2 / sum); else
            if (0.0005 * 0.00025 * transvers1tv2 > threshold)
               ReportSNP(glf, filter, transvers1, transvers2, 0.0005 * 0.00025 * transvers1tv2 / sum); else
            continue;
            if (filter.Length() == 0) homozygousPolymorphism++;
            }
         }

      printf("    Total Aligned Bases = %10.0f (%.1f raw depth)\n", coverage.dSum(), coverage.dSum() / (coverage.Length() + 1e-10));
      printf("Depth at Analysis Sites = %10.1f\n", coverage.dSumProduct(target) / (target.Sum() + 1e-10));

      printf("          Missing bases = %9d (%.3f%%)\n",
            baseCounts[0], (glf.maxPosition - actualBases) * 100. / glf.maxPosition);
      printf("        Reference bases = %9d (%.3f%%)\n",
            actualBases, actualBases * 100. / glf.maxPosition);

      printf("              A/T bases = %9d (%.3f%%, %d A, %d T)\n",
             baseCounts[1] + baseCounts[4],
            (baseCounts[1] + baseCounts[4]) * 100. / actualBases,
             baseCounts[1], baseCounts[4]);

      printf("              G/C bases = %9d (%.3f%%, %d G, %d C)\n",
             baseCounts[2] + baseCounts[3],
            (baseCounts[2] + baseCounts[3]) * 100. / actualBases,
             baseCounts[2], baseCounts[3]);

      printf("    Low Mapping Quality = %9d bases (%.3f%%)\n",
             lowMapQuality, lowMapQuality * 100. / actualBases);

      printf("       Inadequate Depth = %9d bases (%.3f%%)\n",
             inadequateCoverage, inadequateCoverage * 100. / actualBases);

      printf("  Reference Homozygotes = %9d bases (%.3f%%)\n",
             homozygousReference, homozygousReference * 100. / actualBases);

      printf("     Heterozygous Sites = %9d bases (%.3f%%)\n",
             heterozygousPolymorphism, heterozygousPolymorphism * 100. / actualBases);

      printf("            Homozygotes = %9d bases (%.3f%%)\n",
             homozygousPolymorphism, homozygousPolymorphism * 100. / actualBases);

      int noCalls = actualBases - homozygousReference - heterozygousPolymorphism - homozygousPolymorphism;
      printf("                No call = %9d bases (%.3f%%)\n",
            noCalls, noCalls * 100. / actualBases);
      }

   if (baseCalls != NULL)
      fclose(baseCalls);

   printf("RUN COMPLETED!\n");
   }


 
