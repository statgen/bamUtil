////////////////////////////////////////////////////////////////////// 
// glfMultiples/Main.cpp 
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
// Wednesday June 16, 2010
// 
 
#include "BaseQualityHelper.h"
#include "StringArray.h"
#include "StringAlias.h"
#include "Parameters.h"
#include "glfHandler.h"
#include "MathGold.h"
#include "Error.h"

#include <math.h>
#include <time.h>

class GenotypeLikelihood : public ScalarMinimizer
   {
   public:
      int n;
      int position;
      glfHandler * glf;

      void SetAlleles(int al1, int al2)
         {
         allele1 = al1;
         allele2 = al2;

         geno11 = glfHandler::GenotypeIndex(allele1, allele1);
         geno12 = glfHandler::GenotypeIndex(allele1, allele2);
         geno22 = glfHandler::GenotypeIndex(allele2, allele2);
         }

      virtual double f(double freq) { return -Evaluate(freq); }

      double Evaluate(double freq)
         {
         double prior11 = freq * freq;
         double prior12 = freq * (1.0 - freq) * 2.0;
         double prior22 = (1.0 - freq) * (1.0 - freq);

         double likelihood = 0.0;

         for (int i = 0; i < n; i++)
            likelihood += log(prior11 * glf[i].GetLikelihoods(position)[geno11] +
                              prior12 * glf[i].GetLikelihoods(position)[geno12] +
                              prior22 * glf[i].GetLikelihoods(position)[geno22] +
                              1e-30);

         return likelihood;
         }

      void GetPriors(double * priors, double freq)
         {
         for (int i = 0; i < 10; i++)
            priors[i] = 0.0;

         priors[geno11] = freq * freq;
         priors[geno12] = freq * (1.0  - freq) * 2.0;
         priors[geno22] = (1.0 - freq) * (1.0 - freq);
         }

      double OptimizeFrequency()
         {
         a = 0.00001; fa = f(a);
         b = 0.4; fb = f(b);
         c = 0.99999; fc = f(c);

         Brent(0.0001);

         return min;
         }

   protected:
      int allele1, allele2;
      int geno11, geno12, geno22;
   };

FILE * baseCalls = NULL;

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

void DumpDetails(glfHandler * glf, int n, int position, char refBase)
   {
   char alleles[] = { 0, 'a', 'c', 'g', 't' };

   printf("Dump for section %s, position %d [%c]\n",
          (const char *) glf[0].label, position, alleles[refBase]);

   printf("Depth");
   for (int i = 0; i < n; i++)
      printf("\t%d", glf[i].GetDepth(position));
   printf("\n");

   printf("MapQ");
   for (int i = 0; i < n; i++)
      printf("\t%d", glf[i].GetMapQuality(position));
   printf("\n");

   for (int i = 1, index = 0; i <= 4; i++)
      for (int j = i; j <= 4; j++, index++)
         {
         printf("%c/%c", alleles[i], alleles[j]);
         for (int k = 0; k < n; k++)
            printf("\t%d", glf[k].GetLogLikelihoods(position)[index]);
         printf("\n");
         }
   }

int GetBestGenotype(const double likelihoods[], const double priors[])
   {
   int best = 0;

   for (int i = 1; i < 10; i++)
      if (likelihoods[i] * priors[i] > likelihoods[best] * priors[best])
         best = i;

   return best;
   }

const char * GetBestGenotypeLabel(const double likelihoods[], const double priors[])
   {
   const char * genotypeLabel[10] = {"A/A", "A/C", "A/G", "A/T", "C/C", "C/G", "C/T", "G/G", "G/T", "T/T"};

   return genotypeLabel[GetBestGenotype(likelihoods, priors)];
   }

int GetBestRatio(const double likelihoods[], const double priors[])
   {
   double sum = 0.0;
   int best = 0;

   for (int i = 1; i < 10; i++)
      if (likelihoods[i] * priors[i] > likelihoods[best] * priors[best])
         best = i;

   for (int i = 0; i < 10; i++)
      sum += likelihoods[i] * priors[i];

   double error = 1.0 - likelihoods[best] * priors[best]/sum;

   if (error < 0.0000000001)
      return 100;

   return int (-log10(error) * 10 + 0.5);
   }

void ReportGenotypes(GenotypeLikelihood & lk, glfHandler * glf, int n, int position, int refAllele, int al1, int al2)
   {
   if (baseCalls == NULL)
      return;

   double priors[10];
   lk.GetPriors(priors, lk.min);

   int geno11 = glfHandler::GenotypeIndex(al1, al1);
   int geno12 = glfHandler::GenotypeIndex(al1, al2);
   int geno22 = glfHandler::GenotypeIndex(al2, al2);

   int label1 = al1 == refAllele ? 0 : 1;
   int label2 = al2 == refAllele ? 0 : al1 == al2 ? label1 : label1 + 1;

   String label11, label12, label22;
   label11.printf("%d/%d", label1, label1);
   label12.printf("%d/%d", label1, label2);
   label22.printf("%d/%d", label2, label2);

   // fprintf(baseCalls, "\t%.3f", lk.min);

   for (int i = 0; i < n; i++)
      {
      // Report on the best genotype for the current SNP model
      int quality = GetBestRatio(glf[i].GetLikelihoods(position), priors);
      int best = GetBestGenotype(glf[i].GetLikelihoods(position), priors);

      fprintf(baseCalls, "\t%s:%d:%d",
               (const char *)
               (best == geno11 ? label11 : best == geno12 ? label12 : label22),
               glf[i].GetDepth(position), quality);

      fprintf(baseCalls, ":%d,%d,%d", glf[i].GetLogLikelihoods(position)[geno11],
                                      glf[i].GetLogLikelihoods(position)[geno12],
                                      glf[i].GetLogLikelihoods(position)[geno22]);
      }

   fprintf(baseCalls, "\n");
   }

void ReportSNP(glfHandler * glf, int n, int position,
                                 int refBase, int allele1, int allele2,
                                 const char * filter,
                                 int totalCoverage, int averageMapQuality, double posterior)
   {
   if (baseCalls == NULL) return;

   char alleles[] = { 0, 'a', 'c', 'g', 't' };

   // #Chrom   Pos   Id
   fprintf(baseCalls, "%s\t%d\t.\t", (const char *) glf[0].label, position + 1);

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

   GenotypeLikelihood lk;

   // Find best frequency
   lk.glf = glf;
   lk.n = n;
   lk.position = position;
   lk.SetAlleles(allele1, allele2);
   lk.OptimizeFrequency();

   //double maf = lk.min > 0.50 ? 1.0 - lk.min : lk.min;
   double af = 1.0 - lk.min;

   // Information for this call
   fprintf(baseCalls, "depth=%d;mapQ=%d;AF=%.6lf\t", totalCoverage, averageMapQuality, af);

   // Format for this call
   fprintf(baseCalls, "%s", "GT:DP:GQ:GL");

   ReportGenotypes(lk, glf, n, position, refBase, allele1, allele2);
   }

double PolymorphismLikelihood
         (glfHandler * glf, int n, int position, int refAllele, int mutAllele)
   {
   GenotypeLikelihood lk;

   lk.glf = glf;
   lk.n = n;
   lk.position = position;
   lk.SetAlleles(refAllele, mutAllele);
   lk.OptimizeFrequency();

   return -lk.fmin;
   }

double SinkLikelihood
       (glfHandler * glf, int n, int position)
   {
   double lk = 0.0;
   double scale = -log(10.0) / 10;

   for  (int r = 1; r <= 4; r++)
      for (int m = r + 1; m <= 4; m++)
          {
          int geno = glfHandler::GenotypeIndex(r, m);

          double partial = log(1.0 / 6.0);
          for (int i = 0; i < n; i++)
             partial += glf[i].GetLogLikelihoods(position)[geno] * scale;

          if (lk == 0.0)
            {
            lk = partial;
            continue;
            }

          double max = partial > lk ? partial : lk;
          double min = partial > lk ? lk : partial;

          if (max - min > 50)
            lk = max;
          else
            lk = log(exp(partial - min) + exp(lk - min)) + min;
          }

   return lk;
   }

int main(int argc, char ** argv)
   {
   printf("glfMultiples -- SNP calls based on .glf or .glz files\n");
   printf("(c) 2008-2010 Goncalo Abecasis, Sebastian Zoellner, Yun Li\n\n");

   String positionfile;
   String callfile;
   String glfAliases;
   String glfPrefix;
   String glfSuffix;
   ParameterList pl;

   double posterior = 0.50;
   int    mapQuality = 60;
   int    minTotalDepth = 1;
   int    maxTotalDepth = 1000;
   bool   verbose = false;
   bool   mapQualityStrict = false;
   bool   hardFilter = false;
   bool   softFilter = true;
   bool   robustPrior = true;
   bool   uniformPrior = false;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Map Quality Filter")
         LONG_INTPARAMETER("minMapQuality", &mapQuality)
         LONG_PARAMETER("strict", &mapQualityStrict)
      LONG_PARAMETER_GROUP("Total Depth Filter")
         LONG_INTPARAMETER("minDepth", &minTotalDepth)
         LONG_INTPARAMETER("maxDepth", &maxTotalDepth)
      LONG_PARAMETER_GROUP("Position Filter")
         LONG_STRINGPARAMETER("positionFile", &positionfile)
      LONG_PARAMETER_GROUP("Filtering Options")
         EXCLUSIVE_PARAMETER("hardFilter", &hardFilter)
         EXCLUSIVE_PARAMETER("softFilter", &softFilter)
      LONG_PARAMETER_GROUP("Prior Options")
         EXCLUSIVE_PARAMETER("uniformPrior", &uniformPrior)
         EXCLUSIVE_PARAMETER("robustPrior", &robustPrior)
      LONG_PARAMETER_GROUP("Output")
         LONG_PARAMETER("verbose", &verbose)
      LONG_PARAMETER_GROUP("Sample Names")
         LONG_STRINGPARAMETER("glfAliases", &glfAliases)
      LONG_PARAMETER_GROUP("Prefixes and Suffixes")
         LONG_STRINGPARAMETER("glfPrefix",&glfPrefix)
         LONG_STRINGPARAMETER("glfSuffix",&glfSuffix)
   END_LONG_PARAMETERS();

   pl.Add(new StringParameter('b', "Base Call File", callfile));
   pl.Add(new DoubleParameter('p', "Posterior Threshold", posterior));
   pl.Add(new LongParameters("Additional Options", longParameters));
   int argstart = pl.ReadWithTrailer(argc, argv) + 1;
   pl.Status();

   if (posterior < 0.50)
      error("Posterior threshold for genotype calls (-p option) must be > 0.50.");

   time_t t;
   time(&t);

   printf("Analysis started on %s\n", ctime(&t));
   fflush(stdout);

   int n = argc - argstart;
   argv += argstart;

   if (n == 0)
      error("No glf files listed at the end of command line\n");

   // Prior for finding difference from the reference at a particular site
   double prior = 0.0;
   for (int i = 1; i <= 2 * n; i++)
      prior += 1.0 / i;
   prior *= 0.001;

   glfHandler * glf = new glfHandler[n];

   for (int i = 0; i < n; i++) {
     String glfName = glfPrefix + String(argv[i]) + glfSuffix;
     if (!glf[i].Open(glfName))
       error("Failed to open genotype likelihood file [%s]\n", glfName.c_str());
   }

   StringAlias aliases;
   aliases.ReadFromFile(glfAliases);

   printf("Calling genotypes for files ...\n");
   for (int i = 0; i < n; i++)
      printf("%s\n", argv[i]);
   printf("\n");

   baseCalls = fopen(callfile, "wt");

   if (baseCalls != NULL)
      {
      fprintf(baseCalls, "##fileformat=VCFv3.3\n");
      ReportDate(baseCalls);
      fprintf(baseCalls, "##source=glfMultiples\n");
      fprintf(baseCalls, "##minDepth=%d\n", minTotalDepth);
      fprintf(baseCalls, "##maxDepth=%d\n", maxTotalDepth);
      fprintf(baseCalls, "##minMapQuality=%d\n", mapQuality);
      fprintf(baseCalls, "##minPosterior=%.4f\n", posterior);
      fprintf(baseCalls, "##FILTER=mapQ,\"RMS Map Quality\"\n");
      fprintf(baseCalls, "##FILTER=depth,\"Read Depth\"\n");
      fprintf(baseCalls, "##FORMAT=GT,String,1,\"Genotype Call\"\n");
      fprintf(baseCalls, "##FORMAT=GQ,Integer,1,\"Genotype Call Quality\"\n");
      fprintf(baseCalls, "##FORMAT=DP,Integer,1,\"Read Depth\"\n");
      fprintf(baseCalls, "##FORMAT=GL,Integer,3,\"Genotype Likelihoods for 0/0,0/1,1/1 or 1/1,1/2,2/2\"\n");
      fprintf(baseCalls, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
      for (int i = 0; i < n; i++)
         fprintf(baseCalls, "\t%s", (const char *) aliases.GetAlias(argv[i]));
      fprintf(baseCalls, "\n");
      }

   StringArray positions;
   positions.Read(positionfile);

   int positionIndex = 0;
   int nextPosition = positions.Length() ? positions[positionIndex++].AsInteger() : 0;

   while (glf[0].NextSection())
      {
      for (int i = 1; i < n; i++)
         {
         glf[i].NextSection();

         if (glf[0].maxPosition != glf[i].maxPosition || glf[0].label != glf[i].label)
            {
            error("Genotype files '%s' and '%s' are not compatible ...\n"
                "    File '%s' has section %s with %d entries ...\n"
                "    File '%s' section %s with %d entries ...\n",
                argv[0], argv[i],
                argv[0], (const char *) glf[0].label, glf[0].maxPosition,
                argv[i], (const char *) glf[i].label, glf[i].maxPosition);
            }
         }

      printf("Processing section %s with %d entries\n", (const char *) glf[0].label, glf[0].maxPosition);

      char refBase = 0;
      int position = 0;
      int mapQualityFilter = 0;
      int depthFilter = 0;
      int homozygousReference = 0;
      int transitions = 0;
      int transversions = 0;
      int otherPolymorphisms = 0;
      int sinkFilter = 0;
      int baseCounts[5] = {0, 0, 0, 0, 0};

      String filter;

      while (true)
         {
         if (position > 0)
            {
            // Check whether we have reached the end of the current chromosome
            bool done = true;
            for (int i = 0; i < n; i++)
               if (glf[i].data.recordType != 0)
                  done = false;
            if (done) break;
            }

         // Advance to the next position where needed
         for (int i = 0; i < n; i++)
            if (glf[i].position == position)
               glf[i].NextBaseEntry();

         // Figure out the current analysis position
         refBase = glf[0].data.refBase;
         position = glf[0].position;
         for (int i = 0; i < n; i++)
            if (position > glf[i].position)
               {
               position = glf[i].position;
               refBase = glf[i].data.refBase;
               }

         // Avoid alignments that extend past the end of the chromosome
         if (position >= glf[0].maxPosition)
            break;

         baseCounts[refBase]++;

         // These lines can be uncommented for debugging purposes
         // for (int i = 0; i < n; i++)
         //   printf("GLF %d : position %d, refBase %d\n", i, position, refBase);
         // printf("Position: %d, refBase: %d\n", position, refBase);

         if (positions.Length())
            if (nextPosition != position)
               continue;
            else
               nextPosition = positions[positionIndex++].AsInteger();

         if (refBase == 0) continue;

         filter.Clear();

         int totalDepth = glf[0].GetDepth(position);
         for (int i = 1; i < n; i++)
            totalDepth += glf[i].GetDepth(position);

         if (totalDepth < minTotalDepth)
            {
            if (filter.Length() == 0) depthFilter++;
            if (hardFilter) continue;
            filter.catprintf("%sdepth<=%d", filter.Length() ? ";" : "", totalDepth);
            }

         if (totalDepth > maxTotalDepth)
            {
            if (filter.Length() == 0) depthFilter++;
            if (hardFilter) continue;
            filter.catprintf("%sdepth>=%d", filter.Length() ? ";" : "", totalDepth);
            }

         double averageMapQuality = glf[0].GetMapQuality(position);
         for (int i = 1; i < n; i++)
            averageMapQuality += glf[i].GetMapQuality(position);
         averageMapQuality /= n;

         bool passMapQualityFilter = false;

         // Check if we have at least one sample with good quality data
         for (int i = 0; i < n; i++)
            if (glf[i].GetMapQuality(position) >= mapQuality)
               passMapQualityFilter = true;

         if (!passMapQualityFilter)
            {
            if (filter.Length() == 0) mapQualityFilter++;
            if (hardFilter) continue;
            filter.catprintf("%smapQ<%d", filter.Length() ? ";" : "", mapQuality);
            }

         // Create convenient aliases for each base
         unsigned char transition = (((refBase - 1) ^ 2) + 1);
         unsigned char transvers1 = (((refBase - 1) ^ 3) + 1);
         unsigned char transvers2 = (((refBase - 1) ^ 1) + 1);

         int homRef = glf[0].GenotypeIndex(refBase, refBase);

         // Calculate likelihood assuming every is homozygous for the reference
         double lRef = log(1.0 - prior);
	 double pTs = uniformPrior ? 1./3. : 2./3.;
	 double pTv = uniformPrior ? 1./3. : 1./6.;
         for (int i = 0; i < n; i++)
            lRef += log(glf[i].GetLikelihoods(position)[homRef]);

         // Calculate likelihoods for the most likelily SNP configurations
         double refTransition = log(prior * pTs) + PolymorphismLikelihood(glf, n, position, refBase, transition);
         double refTransvers1 = log(prior * pTv) + PolymorphismLikelihood(glf, n, position, refBase, transvers1);
         double refTransvers2 = log(prior * pTv) + PolymorphismLikelihood(glf, n, position, refBase, transvers2);

         // Calculate likelihoods for less likely SNP configurations
         double transitiontv1 = log(prior * 0.001) + PolymorphismLikelihood(glf, n, position, transition, transvers1);
         double transitiontv2 = log(prior * 0.001) + PolymorphismLikelihood(glf, n, position, transition, transvers2);
         double transvers1tv2 = log(prior * 0.001) + PolymorphismLikelihood(glf, n, position, transvers1, transvers2);

         // Calculate the likelihood for unusual configurations where everyone is heterozygous ...
         double sink = n > 10 ? log(prior * 1e-8) + SinkLikelihood(glf, n, position) : -1e100;

         double lmax = max(
               max(max(lRef, refTransition),max(refTransvers1, refTransvers2)),
               max(max(transitiontv1, transitiontv2), max(transvers1tv2, sink)));

         double sum = exp(lRef - lmax) + exp(refTransition -lmax) +
                      exp(refTransvers1 - lmax) + exp(refTransvers2 - lmax) +
                      exp(transitiontv1 - lmax) + exp(transitiontv2 - lmax) +
                      exp(transvers1tv2 - lmax) + exp(sink - lmax);

         if (sum == 0.0) continue;

         if (exp(lRef - lmax)/sum > 1.0 - prior)
            {
            if (filter.Length() == 0) homozygousReference++;

            if (positions.Length())
               ReportSNP(glf, n, position, refBase, refBase, refBase, filter, totalDepth, averageMapQuality, lRef / sum);

            continue;
            }

         double quality = 1.0 - exp(lRef - lmax) / sum;

         if (verbose)
            {
            DumpDetails(glf, n, position, refBase);

            printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
                 lRef, refTransition, refTransvers1, refTransvers2,
                 transitiontv1, transitiontv2, transvers1tv2);
            }

         if (exp(refTransition - lmax)/sum > posterior)
            {
            ReportSNP(glf, n, position, refBase, refBase, transition,
                      filter, totalDepth, averageMapQuality, quality /* refTransition/sum */);
            if (filter.Length() == 0) transitions++;
            }
         else if (exp(refTransvers1 - lmax)/sum > posterior)
            {
            ReportSNP(glf, n, position, refBase, refBase, transvers1,
                      filter, totalDepth, averageMapQuality, quality /* refTransvers1/sum */);
            if (filter.Length() == 0) transversions++;
            }
         else if (exp(refTransvers2 - lmax)/sum > posterior)
            {
            ReportSNP(glf, n, position, refBase, refBase, transvers2,
                      filter, totalDepth, averageMapQuality, quality /* refTransvers2/sum */);
            if (filter.Length() == 0) transversions++;
            }
         else if (exp(transitiontv1 - lmax)/sum > posterior)
            {
            ReportSNP(glf, n, position, refBase, transition, transvers1,
                      filter, totalDepth, averageMapQuality, quality /* transitiontv1/sum */);
            if (filter.Length() == 0) otherPolymorphisms++;
            }
         else if (exp(transitiontv2 - lmax)/sum > posterior)
            {
            ReportSNP(glf, n, position, refBase, transition, transvers2,
                      filter, totalDepth, averageMapQuality, quality /* transitiontv2/sum */);
            if (filter.Length() == 0) otherPolymorphisms++;
            }
         else if (exp(transvers1tv2 - lmax)/sum > posterior)
            {
            ReportSNP(glf, n, position, refBase, transvers1, transvers2,
                      filter, totalDepth, averageMapQuality, quality /* transvers1tv2/sum */);
            if (filter.Length() == 0) otherPolymorphisms++;
            }
         else if (exp(sink - lmax)/sum > posterior)
            sinkFilter++;
         }

      int actualBases = glf[0].maxPosition - baseCounts[0];

      printf("          Missing bases = %9d (%.3f%%)\n",
            baseCounts[0], baseCounts[0] * 100. / glf[0].maxPosition);
      printf("        Reference bases = %9d (%.3f%%)\n",
            glf[0].maxPosition - baseCounts[0], (glf[0].maxPosition - baseCounts[0]) * 100. / glf[0].maxPosition);

      printf("              A/T bases = %9d (%.3f%%, %d A, %d T)\n",
             baseCounts[1] + baseCounts[4],
            (baseCounts[1] + baseCounts[4]) * 100. / actualBases,
             baseCounts[1], baseCounts[4]);

      printf("              G/C bases = %9d (%.3f%%, %d G, %d C)\n",
             baseCounts[2] + baseCounts[3],
            (baseCounts[2] + baseCounts[3]) * 100. / actualBases,
             baseCounts[2], baseCounts[3]);

      printf("           Depth Filter = %9d bases (%.3f%%)\n",
             depthFilter, depthFilter * 100. / actualBases);

      printf("     Map Quality Filter = %9d bases (%.3f%%)\n",
             mapQualityFilter, mapQualityFilter * 100. / actualBases);

      printf("        Non-Polymorphic = %9d bases (%.3f%%)\n",
             homozygousReference, homozygousReference * 100. / actualBases);

      printf("            Transitions = %9d bases (%.3f%%)\n",
             transitions, transitions * 100. / actualBases);

      printf("          Transversions = %9d bases (%.3f%%)\n",
             transversions, transversions * 100. / actualBases);

      printf("    Other Polymorphisms = %9d bases (%.3f%%)\n",
             otherPolymorphisms, otherPolymorphisms * 100. / actualBases);

      if (n > 10)
          printf("          Homology Sink = %9d bases (%.3f%%)\n",
                 sinkFilter, sinkFilter * 100. / actualBases);

      int noCalls = actualBases - homozygousReference - transitions - transversions - otherPolymorphisms - sinkFilter;
      printf("                No call = %9d bases (%.3f%%)\n",
            noCalls, noCalls * 100. / actualBases);

      fflush(stdout);
      }

   if (baseCalls != NULL)
      fclose(baseCalls);
   }


 
