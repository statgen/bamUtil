////////////////////////////////////////////////////////////////////// 
// mach1/AssociationAnalysis.cpp 
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
 
#include "AssociationAnalysis.h"
#include "MathStats.h"

#include <math.h>

void AssociationAnalysis::ScoreNPL(String & prefix, Pedigree & ped, Haplotyper & engine, int rounds)
   {
   FILE * file = fopen(prefix + ".npl", "wt");

   if (file == NULL)
      printf("Error opening output file [%s.npl]\n", (const char *) prefix);
   else
      {
      fprintf(file, "Marker");
      for (int j = 0; j < ped.affectionCount; j++)
         fprintf(file,
                 "\tstat(U)\ts.d.\tstat(A)\ts.d.\tt(%s)\tdf\tp-val"
                 "\tt(%s=1)\tdf\tp-val\tt(%s=2)\tdf\tp-val",
                 (const char *) ped.affectionNames[j],
                 (const char *) ped.affectionNames[j],
                 (const char *) ped.affectionNames[j]);
      fprintf(file, "\n");

      for (int i = 0; i < ped.markerCount; i++)
         {
         fprintf(file, "%s", (const char *) ped.markerNames[i]);

         for (int j = 0; j < ped.affectionCount; j++)
            {
            double sum[2] = {0.0, 0.0}, sumsq[2] = {0.0, 0.0};
            int    counts[2] = {0, 0};

            double scale = 1.0 / (rounds * engine.states);

            for (int k = 0; k < ped.count; k++)
               if (ped[k].affections[j] != 0)
                  {
                  double score = engine.diseaseScores[k][i * ped.affectionCount + j] * scale;

                  sum[ped[k].affections[j] - 1] += score;
                  sumsq[ped[k].affections[j] - 1] += score * score;
                  counts[ped[k].affections[j] - 1] ++;

#if _DEBUG
                  printf("Person %s->%s: %.3f\n",
                        (const char *) ped[k].famid, (const char *) ped[k].pid, score);
#endif
                  }

            if (counts[0] <= 1 || counts[1] <= 1)
               {
               fprintf(file, "\t-\t-\t-\t-\t-\t-\t-");
               continue;
               }

            sum[0] /= counts[0]; sumsq[0] /= counts[0];
            sum[1] /= counts[1]; sumsq[1] /= counts[1];

            sumsq[0] = (sumsq[0] - (sum[0] * sum[0])) * counts[0] / (counts[0] - 1);
            sumsq[1] = (sumsq[1] - (sum[1] * sum[1])) * counts[1] / (counts[1] - 1);

            double s0n = sumsq[0] / counts[0];
            double s1n = sumsq[1] / counts[1];

            double t = (sum[1] - sum[0]) / sqrt(s0n + s1n + 1e-10);
            double df = (s0n + s1n) * (s0n + s1n) /
                        (s0n * s0n / (counts[0] - 1) + s1n * s1n / (counts[1] - 1) + 1e-30);

            fprintf(file, "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t",
                          sum[0], sqrt(sumsq[0]), sum[1], sqrt(sumsq[1]), t, df);

            OutputPValue(file, t, df);

            double tunaff  = sum[0] / sqrt(s0n + 1e-10);
            double dfunaff = counts[0] - 1;

            fprintf(file, "\t%.3f\t%.0f\t", tunaff, dfunaff);

            OutputPValue(file, -tunaff, dfunaff);

            double taff  = sum[1] / sqrt(s1n + 1e-10);
            double dfaff = counts[1] - 1;

            fprintf(file, "\t%.3f\t%.0f\t", taff, dfaff);

            OutputPValue(file, taff, dfaff);
            }

         fprintf(file, "\n");
         }

      fclose(file);
      }

   printf("Wrote out file [%s.npl] with NPL statistics ...\n\n", (const char *) prefix);
   }

void AssociationAnalysis::ScoreMarkers(String & prefix, Pedigree & ped, DosageCalculator & doses)
   {
   FILE * file = fopen(prefix + ".assoc", "wt");

   if (file == NULL)
      printf("Error opening output file [%s.assoc]\n", (const char *) prefix);
   else
      {
      fprintf(file, "Marker");
      for (int j = 0; j < ped.affectionCount; j++)
         fprintf(file, "\tdose(U)\ts.d.\tdose(A)\ts.d.\tt(%s)\tdf\tp-val", (const char *) ped.affectionNames[j]);
      fprintf(file, "\n");

      for (int i = 0; i < ped.markerCount; i++)
         {
         fprintf(file, "%s", (const char *) ped.markerNames[i]);

         for (int j = 0; j < ped.affectionCount; j++)
            {
            double sum[2] = {0.0, 0.0}, sumsq[2] = {0.0, 0.0};
            int    counts[2] = {0, 0};

            for (int k = 0; k < ped.count; k++)
               if (ped[k].affections[j] != 0)
                  {
                  double score = doses.GetDosage(k, i);

                  sum[ped[k].affections[j] - 1] += score;
                  sumsq[ped[k].affections[j] - 1] += score * score;
                  counts[ped[k].affections[j] - 1] ++;
                  }

            if (counts[0] <= 1 || counts[1] <= 1)
               {
               fprintf(file, "\t-\t-\t-\t-\t-\t-\t-");
               continue;
               }

            sum[0] /= counts[0]; sumsq[0] /= counts[0];
            sum[1] /= counts[1]; sumsq[1] /= counts[1];

            sumsq[0] = (sumsq[0] - (sum[0] * sum[0])) * counts[0] / (counts[0] - 1);
            sumsq[1] = (sumsq[1] - (sum[1] * sum[1])) * counts[1] / (counts[1] - 1);

            double s0n = sumsq[0] / counts[0];
            double s1n = sumsq[1] / counts[1];

            double t = (sum[1] - sum[0]) / sqrt(s0n + s1n + 1e-10);
            double df = (s0n + s1n) * (s0n + s1n) /
                        (s0n * s0n / (counts[0] - 1) + s1n * s1n / (counts[1] - 1));

            fprintf(file, "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t",
                          sum[0], sqrt(sumsq[0]), sum[1], sqrt(sumsq[1]), t,
                          df);

            if (df < 1.0)
               fprintf(file, "-");
            else
               {
               double pvalue = tdist(t, df);

               fprintf(file, "%#.2g", pvalue);
               }
            }

         fprintf(file, "\n");
         }

      fclose(file);
      }

   printf("Wrote out file [%s.assoc] with marker-by-marker t statistics ...\n\n",
          (const char *) prefix);
   }

void AssociationAnalysis::OutputPValue(FILE * file, double t, double df)
   {
   if (df < 1.0)
      fprintf(file, "-");
   else
      {
      double pvalue = t >= 0.0 ? tdist(t, df) * 0.5 : 1.0 - tdist(t, df) * 0.5;

      fprintf(file, "%#.2g", pvalue);
      }
   }
 
