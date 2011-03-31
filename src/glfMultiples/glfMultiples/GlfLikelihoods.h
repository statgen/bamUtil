#ifndef __GLFLIKELIHOODS_H__
#define __GLFLIKELIHOODS_H__

#include "glfHandler.h"
#include "MathGold.h"
#include "Pedigree.h"

#include <math.h>

// Chromosome types
#define CT_AUTOSOME     0
#define CT_CHRX         1
#define CT_CHRY         2
#define CT_MITO         3

class GenotypeLikelihood : public ScalarMinimizer
   {
   public:
      int n;
      int position;
      glfHandler * glf;

      void SetAlleles(int al1, int al2)
         {
         chromosomeType = CT_AUTOSOME;

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

         double prior1 = freq;
         double prior2 = 1.0 - freq;

         double likelihood = 0.0;

         switch (chromosomeType)
            {
            case CT_MITO :
               prior11 = prior1;
               prior12 = 0.0;
               prior22 = prior2;
            case CT_AUTOSOME :
               for (int i = 0; i < n; i++)
                  likelihood += log(prior11 * glf[i].GetLikelihoods(position)[geno11] +
                                    prior12 * glf[i].GetLikelihoods(position)[geno12] +
                                    prior22 * glf[i].GetLikelihoods(position)[geno22] +
                                    1e-30);
                  break;
            case CT_CHRY :
               for (int i = 0; i < n; i++)
                  if (sex[i] == SEX_MALE)
                     likelihood += log(prior1 * glf[i].GetLikelihoods(position)[geno11] +
                                       prior2 * glf[i].GetLikelihoods(position)[geno22] +
                                       1e-30);
                  break;
            case CT_CHRX :
               for (int i = 0; i < n; i++)
                  if (sex[i] == SEX_MALE)
                     likelihood += log(prior1 * glf[i].GetLikelihoods(position)[geno11] +
                                       prior2 * glf[i].GetLikelihoods(position)[geno22] +
                                       1e-30);
                  else
                     likelihood += log(prior11 * glf[i].GetLikelihoods(position)[geno11] +
                                       prior12 * glf[i].GetLikelihoods(position)[geno12] +
                                       prior22 * glf[i].GetLikelihoods(position)[geno22] +
                                       1e-30);
            }


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
      int chromosomeType;
      int allele1, allele2;
      int geno11, geno12, geno22;

      char * sex;
   };

class FilterLikelihood : public ScalarMinimizer
   {
   public:
      int n;
      int position;
      glfHandler * glf;

      void SetReferenceAllele(int ref)
         {
         chromosomeType = CT_AUTOSOME;

         reference = ref;

         for (int a1 = 1; a1 <= 4; a1++)
            for (int a2 = a1; a2 <= 4; a2++)
               {
               int index = glfHandler::GenotypeIndex(a1, a2);

               if (a1 == a2)
                  group[index] = a1 == ref ? 0 : 2;
               else
                  group[index] = (a1 == ref || a2 == ref) ? 1 : 2;
               }
         }

      virtual double f(double freq) { return -Evaluate(freq); }

      double Evaluate(double freq)
         {
         double prior11 = freq * freq;
         double prior12 = freq * (1.0 - freq) * 2.0;
         double prior22 = (1.0 - freq) * (1.0 - freq);

         double prior1 = freq;
         double prior2 = 1.0 - freq;

         double likelihood = 0.0;

         double lk[3];
         for (int i = 0; i < n; i++)
            {
            for (int j = 0; j < 3; j++)
               lk[j] = 0.0;

            for (int j = 0; j < 10; j++)
               lk[group[j]] = max(lk[group[j]], glf[i].GetLikelihoods(position)[j]);

            switch (chromosomeType)
               {
               case CT_MITO :
                  prior11 = prior1;
                  prior12 = 0.0;
                  prior22 = prior2;
               case CT_AUTOSOME :
                  likelihood += log(prior11 * lk[0] + prior12 * lk[1] +
                                    prior22 * lk[2] + 1e-30);
                  break;
               case CT_CHRY :
                  if (sex[i] == SEX_MALE)
                     likelihood += log(prior1 * lk[0] + prior2 * lk[2] + 1e-30);
                  break;
               case CT_CHRX :
                  if (sex[i] == SEX_MALE)
                     likelihood += log(prior1 * lk[0] + prior2 * lk[2] + 1e-30);
                  else
                     likelihood += log(prior11 * lk[0] + prior12 * lk[1] +
                                       prior22 * lk[2] + 1e-30);
               }
            }

         return likelihood;
         }

      double OptimizeFrequency()
         {
         a = 0.00001; fa = f(a);
         b = 0.4; fb = f(b);
         c = 1.0 - 0.5 / n; fc = f(c);

         Brent(0.0001);

         return min;
         }

   protected:
      int chromosomeType;
      int reference;
      int group[10];

      char * sex;
   };


#endif

