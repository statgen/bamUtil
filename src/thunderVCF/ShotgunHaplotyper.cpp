////////////////////////////////////////////////////////////////////// 
// thunder/ShotgunHaplotyper.cpp 
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
 
#include "ShotgunHaplotyper.h"
#include "MemoryAllocators.h"
#include "MemoryInfo.h"
#include "Pedigree.h"
#include "VcfFile.h"
#include "Error.h"
#include <math.h>
#include <limits.h>

ShotgunHaplotyper::ShotgunHaplotyper()
   {
   phred2probAvailable = false;
   
   skipSanityCheck = true;

   //shotgunErrorMatrix = AllocateFloatMatrix(3, 256);

   //SetShotgunError(0.005);
   
   refalleles = NULL;
   freq1s = NULL;
   //weightByMismatch = false;
   weightByLikelihood = false;
   weightByLongestMatch = false;
   }

ShotgunHaplotyper::~ShotgunHaplotyper()
   {
   //FreeFloatMatrix(shotgunErrorMatrix, 3);
   
   if (refalleles != NULL)
      delete [] refalleles;
   
   if (freq1s != NULL)
      delete [] freq1s;
   
   if (phred2probAvailable)
      delete [] phred2prob;
   }

void ShotgunHaplotyper::CalculatePhred2Prob()
{
  if (phred2probAvailable)
     return;
     
  phred2prob = new double [500];
  for (int i = 0; i < 500; i++)
  {
    phred2prob[i] = pow(10.0, (-0.1 * i));
  }
  
  phred2probAvailable = true;
  
}

void ShotgunHaplotyper::WeightByLikelihood()
{
  // Calculate \logPr(Data|h1,h2), weight by Pr(Data|)
  AllocateWeights();
  int minPhred = INT_MAX;
  int* phreds = new int[individuals-1];
  for(int i=0; i < individuals-1; ++i) {
    int phred = 0;
    for(int l=0; l < markers; ++l) {
      char* h1 = haplotypes[2*i];
      char* h2 = haplotypes[2*i+1];
      phred += genotypes[i][ 3*l + h1[l] + h2[l] ];
    }
    phreds[i] = phred;
    if ( minPhred > phred ) minPhred = phred;
  }
  for(int i=0; i < individuals-1; ++i) {
    int phredDiff = phreds[i]-minPhred;
    weights[i] = 1./phredDiff;
    //weights[i] = (phredDiff > 60) ? 0.000001 : phred2prob[phredDiff];
  }
  delete [] phreds;  
}

void ShotgunHaplotyper::ChooseByLongestMatch(int* array) {
  int* sumDiff = new int[markers]();
  int* longestMatches = new int[individuals-1];
  int nDiff2, nMaxDiff2, nPrev, nLongestMatch;
  int min = INT_MAX, max = 0, sum = 0, sumsq = 0, minIdx = -1;

  for(int i=0; i < individuals-1; ++i) {
    // find longest match
    array[i] = 0; // initial do not deselect it
    nLongestMatch = 0;
    for(int j=0; j < 2; ++j) {
      char* ha  = haplotypes[2*(individuals-1)+j];
      for(int k=0; k < 2; ++k) {
	char* hb  = haplotypes[ 2*i + k ];
	for(int l=0; l < markers; ++l) {
	  sumDiff[l] = ( l == 0 ? (ha[l] ^ hb[l]) : sumDiff[l-1] + ha[l] ^ hb[l]);
	}
	nPrev = 0;
	nDiff2 = 0;
	nMaxDiff2 = 0;
	for(int l=1; l < markers; ++l) {
	  if (sumDiff[l] != sumDiff[l-1] ) {
	    if ( l - nPrev > nMaxDiff2 ) {
	      nMaxDiff2 = l-nPrev;
	    }
	    nPrev = l;
	  }
	}
	if ( nLongestMatch < nMaxDiff2 ) 
	  nLongestMatch = nMaxDiff2;
      }
    }
    if ( min > nLongestMatch ) min = nLongestMatch;
    if ( max < nLongestMatch ) max = nLongestMatch;
    sum += nLongestMatch;
    sumsq += (nLongestMatch*nLongestMatch);

    longestMatches[i] = nLongestMatch;
    if ( i < states / 2 ) {
      array[i] = 1;
    }
    else { // using n^2 algorithm for selection -- could be better, but should be fine
      int minVal = INT_MAX;
      int minIdx = -1;
      for(int j=0; j <= i; ++j) {
	if ( array[j] == 1 ) {
	  if ( minVal > longestMatches[j] ) {
	    minVal = longestMatches[j];
	    minIdx = j;
	  }
	}
      }
      if ( minIdx < 0 )
	error("Cannot find minimum longest match");

      if ( minIdx != i ) {
	array[minIdx] = 0;
	array[i] = 1;
      }
    }
  }
  //fprintf(stderr,"%d\t%d\t%.2lf\t%.2lf\n",min,max,(double)sum/(individuals-1),sqrt((double)sumsq/(individuals-1.)-(double)sum*sum/(individuals-1.)/(individuals-1.)));
  //for(int i=0; i < individuals-1; ++i) {
  //  fprintf(stderr,"%d\t%d\t%d\n",i,array[i],longestMatches[i]);
  //}
  delete [] sumDiff;  
  delete [] longestMatches;
}

void ShotgunHaplotyper::WeightByLongestMatch() 
{
  AllocateWeights();
  int* sumDiff = new int[markers]();
  int nDiff2, nMaxDiff2, nPrev, nLongestMatch;
  int min = INT_MAX, max = 0, sum = 0, sumsq = 0;

  for(int i=0; i < individuals-1; ++i) {
    // find longest match
    nLongestMatch = 0;
    for(int j=0; j < 2; ++j) {
      char* ha  = haplotypes[2*(individuals-1)+j];
      for(int k=0; k < 2; ++k) {
	char* hb  = haplotypes[ 2*i + k ];
	for(int l=0; l < markers; ++l) {
	  sumDiff[l] = ( l == 0 ? (ha[l] ^ hb[l]) : sumDiff[l-1] + ha[l] ^ hb[l]);
	}
	nPrev = 0;
	nDiff2 = 0;
	nMaxDiff2 = 0;
	for(int l=1; l < markers; ++l) {
	  if (sumDiff[l] != sumDiff[l-1] ) {
	    if ( l - nPrev > nMaxDiff2 ) {
	      nMaxDiff2 = l-nPrev;
	    }
	    nPrev = l;
	  }
	}
	if ( nLongestMatch < nMaxDiff2 ) 
	  nLongestMatch = nMaxDiff2;
      }
    }
    if ( min > nLongestMatch ) min = nLongestMatch;
    if ( max < nLongestMatch ) max = nLongestMatch;
    sum += nLongestMatch;
    sumsq += (nLongestMatch*nLongestMatch);
    weights[i] = (float)nLongestMatch;
  }
  fprintf(stderr,"%d\t%d\t%.2lf\t%.2lf\n",min,max,(double)sum/(individuals-1),sqrt((double)sumsq/(individuals-1.)-(double)sum*sum/(individuals-1.)/(individuals-1.)));
  delete [] sumDiff;
}

void ShotgunHaplotyper::WeightByMismatch()
{
  // assume that the individual to 
  AllocateWeights();

  float sum = 0.0;
  for (int i = 0; i < individuals-1; ++i) {
    // Calculate the mismatch between the haplotypes
    int minMismatches = INT_MAX;
    int mismatches = 1;
    for(int j=0; j < 2; ++j) {
      char* ha  = haplotypes[2*(individuals-1)+j];
      for(int k=0; k < 2; ++k) {
	char* hb  = haplotypes[ 2*i + k ];
	for(int l=0; l < markers; ++l) {
	  mismatches += (ha[l] ^ hb[l]);
	}
      }
    }
    if ( mismatches < minMismatches ) 
      minMismatches = mismatches;
    weights[i] = 1.0/minMismatches;
    sum += weights[i];
  }
  weights[individuals-1] = 0.0;
  //fprintf(stderr,"Weighting by mismatch count, sum = %f\n",sum);
}

void ShotgunHaplotyper::LoadHaplotypesFromVCF(String& fileName)
{
  //  if (rand == NULL)
  //   rand = &globalRandom;

  printf("starting LoadHaplotypesFromVCF\n");
  bool warningsPrinted = false;
  try {
    VcfFile* pVcf = new VcfFile;
    pVcf->bSiteOnly = false;
    pVcf->bParseGenotypes = true;
    pVcf->bParseDosages = false;
    pVcf->bParseValues = false;

    VcfMarker* pMarker = new VcfMarker;
    //CalculatePhred2Prob();

    pVcf->openForRead(fileName.c_str());

    for(int j=0; pVcf->iterateMarker(); ++j) {
      //fprintf(stderr,"j=%d\n",j);
      pMarker = pVcf->getLastMarker();
      for(int i=0; i < individuals; ++i) {
	//fprintf(stderr,"i=%d\n",j);
	unsigned short g = pMarker->vnSampleGenotypes[i];
	char g1, g2;

	// genotype is missing
	if ( g == 0xffff ) {
	  fprintf(stderr,"ERROR: Observed Missing genotypes");
	  abort();
	}
	else {
	  // genotype is unphased
	  if ( (g & 0x8000) == 0 ) {
	    if ( !warningsPrinted ) {
	      fprintf(stderr,"ERROR: Observed unphased genotypes %x",g);
	      abort();
	    }
	  }
	  g1 = (((g & 0x7f00) >> 8) & 0xff);
	  g2 = (g & 0x7f);
	}

	if ( pMarker->asAlts.Length() > 1 ) {
	  if ( g1 == 0 || g2 == 0 ) {
	    fprintf(stderr,"ERROR: TriAllelic Site, but '0' genotype is observed");
	    abort();
	  }
	  --g1;
	  --g2;
	}
	haplotypes[ i * 2 ][ j ] = g1; 
	haplotypes[ i * 2 + 1][ j ] = g2; 
      }
    }
    delete pVcf;
    //delete pMarker;
  }
  catch ( VcfFileException e ) {
    error( e.what() );
  }
}

void ShotgunHaplotyper::RandomSetup(Random * rand)
   {
   if (rand == NULL)
      rand = &globalRandom;

   CalculatePhred2Prob();
   
   for (int j = 0; j < markers; j++)
      {
      double mac = 0;
      int markerindex = 3*j;
      
      double hyperprior11 = freq1s[j] * freq1s[j];
      double hyperprior12 = 2.0 * freq1s[j] * (1.0 - freq1s[j]);
      double hyperprior22 = (1.0 - freq1s[j]) * (1.0 - freq1s[j]);
      
      for (int i = 0; i < individuals; i++)
         {
         double post11 = hyperprior11 * phred2prob[genotypes[i][markerindex]];
         double post12 = hyperprior12 * phred2prob[genotypes[i][markerindex+1]];
         double post22 = hyperprior22 * phred2prob[genotypes[i][markerindex+2]];
         double sumpost = post11 + post12 + post22;
         post11 /= sumpost;
         post12 /= sumpost;
         post22 /= sumpost;

         // estimated counts of AL2
         mac += post12+ 2*post22;
         }

      //here, each person contributes two alleles
      double freq = 0.5 * mac / (double) individuals;
      
      double prior_11 = (1.0 - freq) * (1.0 - freq);
      double prior_12 = 2.0 * freq * (1.0 - freq);
      double prior_22 = freq * freq;

      for (int i = 0; i < individuals; i++)
         {
         int observed = (unsigned char) (genotypes[i][j]);

         double posterior_11 = prior_11 * phred2prob[genotypes[i][markerindex]];
         double posterior_12 = prior_12 * phred2prob[genotypes[i][markerindex+1]];
         double posterior_22 = prior_22 * phred2prob[genotypes[i][markerindex+2]];
         double sum = posterior_11 + posterior_12 + posterior_22;

         if (sum == 0)
            printf("Problem!\n");

         posterior_11 /= sum;
         posterior_12 /= sum;

         double r = rand->Next();
         if (r < posterior_11)
            {
            haplotypes[i * 2][j] = 0;
            haplotypes[i * 2 + 1][j] = 0;
            }
         else if (r < posterior_11 + posterior_12)
            {
            bool bit = rand->Binary();

            haplotypes[i * 2][j] = bit;
            haplotypes[i * 2 + 1][j] = bit ^ 1;
            }
         else
            {
            haplotypes[i * 2][j] = 1;
            haplotypes[i * 2 + 1][j] = 1;
            }

         }
      }
   }


void ShotgunHaplotyper::PhaseByReferenceSetup(Random * rand)
   {
   if (rand == NULL)
      rand = &globalRandom;

   CalculatePhred2Prob();
   
   for (int j = 0; j < markers; j++)
      {
      double mac = 0;
      int markerindex = 3*j;
      
      double hyperprior11 = freq1s[j] * freq1s[j];
      double hyperprior12 = 2.0 * freq1s[j] * (1.0 - freq1s[j]);
      double hyperprior22 = (1.0 - freq1s[j]) * (1.0 - freq1s[j]);
      
      for (int i = 0; i < individuals; i++)
         {
         double post11 = hyperprior11 * phred2prob[genotypes[i][markerindex]];
         double post12 = hyperprior12 * phred2prob[genotypes[i][markerindex+1]];
         double post22 = hyperprior22 * phred2prob[genotypes[i][markerindex+2]];
         double sumpost = post11 + post12 + post22;
         post11 /= sumpost;
         post12 /= sumpost;
         post22 /= sumpost;

         // estimated counts of AL2
         mac += post12+ 2*post22;
         }

      //here, each person contributes two alleles
      double freq = 0.5 * mac / (double) individuals;
      
      double prior_11 = (1.0 - freq) * (1.0 - freq);
      double prior_12 = 2.0 * freq * (1.0 - freq);
      double prior_22 = freq * freq;

      for (int i = 0; i < individuals; i++)
         {
         int observed = (unsigned char) (genotypes[i][j]);

         double posterior_11 = prior_11 * phred2prob[genotypes[i][markerindex]];
         double posterior_12 = prior_12 * phred2prob[genotypes[i][markerindex+1]];
         double posterior_22 = prior_22 * phred2prob[genotypes[i][markerindex+2]];
         double sum = posterior_11 + posterior_12 + posterior_22;

         if (sum == 0)
            printf("Problem!\n");

         posterior_11 /= sum;
         posterior_12 /= sum;

         double r = rand->Next();
         if (r < posterior_11)
            {
            haplotypes[i * 2][j] = 0;
            haplotypes[i * 2 + 1][j] = 0;
            }
         else if (r < posterior_11 + posterior_12)
            {
	      //bool bit = rand->Binary();

	      haplotypes[i * 2][j] = 0;
	      haplotypes[i * 2 + 1][j] = 1;
            }
         else
            {
            haplotypes[i * 2][j] = 1;
            haplotypes[i * 2 + 1][j] = 1;
            }

         }
      }
   }

/*
void ShotgunHaplotyper::SetShotgunError(double rate)
   {
   // Store the background rate
   shotgunError = rate;

   // First calculate binomial coefficients
   int binomial[33][33];

   binomial[0][0] = 1;
   binomial[1][0] = binomial[1][1] = 1;

   for (int i = 2; i < 32; i++)
      {
      binomial[i][0] = binomial[i][i] = 1;

      for (int j = 1; (j < i) && (j < 16); j++)
         binomial[i][j] = binomial[i-1][j] + binomial[i-1][j-1];
      }

   // Next setup the error matrices for each possible true genotype
   for (int i = 0; i < 16; i++)
      for (int j = 0; j < 16; j++)
         if (rate == 0)
            {
            shotgunErrorMatrix[0][j*16 + i] = j == 0 ? 1.0 : 0.0;
            shotgunErrorMatrix[1][j*16 + i] = pow(0.5, i + j) * binomial[i+j][j];
            shotgunErrorMatrix[2][j*16 + i] = i == 0 ? 1.0 : 0.0;
            }
         else
            {
            shotgunErrorMatrix[0][j*16 + i] = pow(1.0 - rate, i) * pow(rate, j) * binomial[i+j][j];
            shotgunErrorMatrix[1][j*16 + i] = pow(0.5, i + j) * binomial[i+j][j];
            shotgunErrorMatrix[2][j*16 + i] = pow(rate, i) * pow(1.0 - rate, j) * binomial[i+j][j];
            }
   }
*/


void ShotgunHaplotyper::ConditionOnData(float * matrix, int marker, 
                                        char phred11, char phred12, char phred22)
   {
   // We treat missing genotypes as uninformative about the mosaic's
   // underlying state. If we were to allow for deletions and the like,
   // that may no longer be true.
   
   //if (genotype == GENOTYPE_MISSING)
      //return;

   double conditional_probs[3];
   int ph11 = (unsigned char) phred11;
   int ph12 = (unsigned char) phred12;
   int ph22 = (unsigned char) phred22;
   
   CalculatePhred2Prob();

   for (int i = 0; i < 3; i++)
      conditional_probs[i] = Penetrance(marker, i, 0) *  phred2prob[ph11] +
                             Penetrance(marker, i, 1) *  phred2prob[ph12] +
                             Penetrance(marker, i, 2) *  phred2prob[ph22];

   for (int i = 0; i < states; i++)
      {
      double factors[2];

      factors[0] = conditional_probs[haplotypes[i][marker]];
      factors[1] = conditional_probs[haplotypes[i][marker] + 1];

      for (int j = 0; j <= i; j++, matrix++)
         *matrix *= factors[haplotypes[j][marker]];
      }
   }

void ShotgunHaplotyper::ImputeAlleles(int marker, int state1, int state2, Random * rand)
   {
   int copied1 = haplotypes[state1][marker];
   int copied2 = haplotypes[state2][marker];

   int markerindex = marker * 3;
   int ph11 = (unsigned char) genotypes[states / 2][markerindex];
   int ph12 = (unsigned char) genotypes[states / 2][markerindex+1];
   int ph22 = (unsigned char) genotypes[states / 2][markerindex+2];

   CalculatePhred2Prob();

   double posterior_11 = Penetrance(marker, copied1 + copied2, 0) * phred2prob[ph11];
   double posterior_12 = Penetrance(marker, copied1 + copied2, 1) * phred2prob[ph12];
   double posterior_22 = Penetrance(marker, copied1 + copied2, 2) * phred2prob[ph22];
   double sum = posterior_11 + posterior_12 + posterior_22;

   posterior_11 /= sum;
   posterior_22 /= sum;

   double r = rand->Next();

   if (r < posterior_11)
      {
      haplotypes[states][marker] = 0;
      haplotypes[states + 1][marker] = 0;
      }
   else if (r < posterior_11 + posterior_22)
      {
      haplotypes[states][marker] = 1;
      haplotypes[states + 1][marker] = 1;
      }
   else if (copied1 != copied2)
      {
      double rate = GetErrorRate(marker);

      if (rand->Next() < rate * rate / ((rate * rate) + (1.0 - rate) * (1.0 - rate)))
         {
         copied1 = !copied1;
         copied2 = !copied2;
         }

      haplotypes[states][marker] = copied1;
      haplotypes[states + 1][marker] = copied2;
      }
   else
      {
      bool bit = rand->Binary();

      haplotypes[states][marker] = bit;
      haplotypes[states + 1][marker] = bit ^ 1;
      }

   int imputed1 = haplotypes[states][marker];
   int imputed2 = haplotypes[states + 1][marker];

   //int differences = abs(copied1 - imputed1) + abs(copied2 - imputed2);
   int differences = abs(copied1 + copied2 - imputed1 - imputed2);

   error_models[marker].matches += 2 - differences;
   error_models[marker].mismatches += differences;
   }


 
void ShotgunHaplotyper::EstimateMemoryInfo(int Individuals, int Markers, int States, bool Compact, bool Phased)
   {
   if (States <= 0 || States > Individuals * 2 - 2)
      States = 2 * Individuals - 2;

   int positions = Compact ? 2 * (int) sqrt((double)Markers) + 1 : Markers;

   if (Phased)
      if (Markers / ((States + 1) / 2) + 1 > positions)
         positions = Markers / ((States + 1) / 2) + 1;

   double bytes = sizeof(char) * (double) Individuals * Markers * 5 // Genotypes, Haplotypes
                + sizeof(float) * (double) States * 2               // Marginals
                + sizeof(float) * (double) positions * States * (States + 1) / 2  // matrices
                + sizeof(float) * (double) Markers * 11             // penetrances, probabilities, thetas
                + sizeof(int) * (double) Markers                    // crossover counts
                + sizeof(Errors) * (double) Markers             // error model information
                + sizeof(char) * (double) Markers;              // reference allele

   printf("   %40s %s\n", "Haplotyping engine (max) ...", (const char *) MemoryInfo(bytes));
   }

bool ShotgunHaplotyper::AllocateMemory(int persons, int maxStates, int m)
   {
   individuals = persons;
   states = maxStates > 1 && maxStates < individuals * 2 ? maxStates & ~1: individuals * 2 - 2;
   markers = m;

   genotypes = AllocateCharMatrix(individuals, markers*3);
   haplotypes = AllocateCharMatrix(individuals * 2, markers);
   
   refalleles = new char [markers];
   freq1s = new double[markers];

   marginals = new float [states * 2];

   leftMatrices = new float * [markers];
   leftProbabilities = new float [markers];

   memoryBlock = new float * [markers];
   smallMemoryBlock = new float * [markers];
   smallFree = 0;

   stack = new int [markers];
   stackPtr = -1;

   thetas = new float [markers - 1];
   crossovers = new int [markers - 1];

   error_models = new Errors [markers];
   penetrances = new float [markers * 9];

   if (genotypes == NULL || haplotypes == NULL || marginals == NULL ||
       leftMatrices == NULL || leftProbabilities == NULL || thetas == NULL ||
       crossovers == NULL || error_models == NULL || penetrances == NULL)
       return readyForUse = false;

   for (int i = 0; i < markers; i++)
      memoryBlock[i] = smallMemoryBlock[i] = NULL;

   for (int i = 0; i < markers - 1; i++)
      thetas[i] = 0.01;

   gridSize = economyMode ? (int) sqrt((double)markers) : markers;

   orderedGenotypeFlags = new int [individuals];

   for (int i = 0; i < individuals; i++)
      orderedGenotypeFlags = 0;

   return readyForUse = true;
   }



void ShotgunHaplotyper::ScoreLeftConditional()
   {
   ResetMemoryPool();
   GetMemoryBlock(0);

   SetupPrior(leftMatrices[0]);
   ConditionOnData(leftMatrices[0], 0, genotypes[states / 2][0], genotypes[states / 2][1], genotypes[states / 2][2]);

   double theta = 0.0;
   float *from = leftMatrices[0];
   for (int i = 1; i < markers; i++)
      {
      int markerindex = i*3;
      
      // Cumulative recombination fraction allows us to skip uninformative positions
      theta = theta + thetas[i - 1] - theta * thetas[i - 1];

      // Skip over uninformative positions to save time
      // maybe check for the difference between min(phred11, phred12, phred22) and the next smallest
      
      //if (genotypes[states / 2][i] != GENOTYPE_MISSING || i == markers - 1)
         {
         GetMemoryBlock(i);

         Transpose(from, leftMatrices[i], theta);
         ConditionOnData(leftMatrices[i], i, genotypes[states / 2][markerindex],
                         genotypes[states / 2][markerindex+1], genotypes[states / 2][markerindex+2]);

         theta = 0;
         from = leftMatrices[i];
         }
      }

   MarkMemoryPool();
   }



void ShotgunHaplotyper::RetrieveMemoryBlock(int marker)
   {
   if (stack[stackPtr] <= marker)
      return;
   else
      {
      ResetReuseablePool();

      double theta = 0.0;
      float *from = leftMatrices[stack[--stackPtr]];

      for (int i = stack[stackPtr] + 1; i <= marker; i++)
         {
         int markerindex = i * 3;
         // Cumulative recombination fraction allows us to skip uninformative positions
         theta = theta + thetas[i - 1] - theta * thetas[i - 1];

         // Skip over uninformative positions to save time
         // if (genotypes[states / 2][i] != GENOTYPE_MISSING || i == markers - 1)
            {
            leftMatrices[i] = GetReuseableBlock();

            Transpose(from, leftMatrices[i], theta);
            ConditionOnData(leftMatrices[i], i, genotypes[states / 2][markerindex],
                            genotypes[states / 2][markerindex+1], genotypes[states / 2][markerindex+2]);

            theta = 0;
            from = leftMatrices[i];
            }
         }
      }
   }

void ShotgunHaplotyper::SampleChromosomes(Random * rand)
   {
   // Print(markers - 1);
   RewindMemoryPool();
   RetrieveMemoryBlock(markers - 1);

   float * probability = leftMatrices[markers - 1];
   float sum = 0.0;

   // Calculate sum over all states
   for (int i = 0; i < states; i++)
      for (int j = 0; j <= i; j++)
         {
         sum += *probability;
         probability++;
         }

   // Sample number and select state
   float choice = rand->Uniform(0, sum);

   sum = 0.0;

   int first = 0, second = 0;
   for (probability = leftMatrices[markers - 1]; first < states; first++)
      {
      for (second = 0; second <= first; second++)
         {
         sum += *probability;
         probability++;

         if (sum >= choice) break;
         }

      if (second <= first) break;
      }

   // printf("Cumulative probability: %g\n", sum);
   // printf("           Random draw: %g\n", choice);
   // printf("        Selected state: %g\n", *(probability - 1));

   for (int j = markers - 2; j >= 0; j--)
      {
      // printf("Sum: %f, Chose (%d,%d)\n", sum, first, second);

      ImputeAlleles(j + 1, first, second, rand);

      // Starting marker for this iteration
      int   j0 = j;

      // Cumulative recombination fraction, skipping over uninformative
      // positions
      float theta = thetas[j];
      //while (genotypes[states / 2][j] == GENOTYPE_MISSING && j > 0)
         //{
         //--j;
         //theta = theta + thetas[j] - theta * thetas[j];
         //}

      // When examining the previous location we consider three alternatives:
      // states that could be reached when both haplotypes recombine (11),
      // states that can be reached when the first (10) or second (01) haplotype recombines,
      // and the states that can be reached without recombination.

      float sum00 = 0.0, sum01 = 0.0, sum10 = 0.0, sum11 = 0.0;

      RetrieveMemoryBlock(j);
      probability = leftMatrices[j];

      for (int k = 0; k < states; k++)
         for (int l = 0; l <= k; l++, probability++)
            {
            sum11 += *probability;
            if (first == k || first == l) sum01 += *probability;
            if (second == k || second == l) sum10 += *probability;
            if (first == k && second == l || first == l && second == k) sum00 += *probability;
            }

      if (weights != NULL)
         {
         sum01 *= weights[second / 2];
         sum10 *= weights[first / 2];
         sum11 *= weights[second / 2] * weights[first / 2];
         }

      sum = sum11 * theta * theta / (states * states) +
            (sum10 + sum01) * theta * (1.0 - theta) / states +
            sum00 * (1.0 - theta) * (1.0 - theta);

      // Sample number and decide how many state changes occurred between the
      // two positions
      choice = rand->Uniform(0, sum);

      // The most likely outcome is that no changes occur ...
      choice -= sum00 * (1.0 - theta) * (1.0 - theta);
      if (choice <= 0.0)
         {
         // Record outcomes for intermediate, uninformative, positions
         FillPath(states, j, j0 + 1, first);
         FillPath(states + 1, j, j0 + 1, second);

         continue;
         }

      // But perhaps the first or second haplotype recombined
      probability = leftMatrices[j];

      choice -= sum10 * theta * (1.0 - theta) / states;
      if (choice <= 0.0)
         {
         // The first haplotype changed ...
         choice = choice * states / (theta * (1.0 - theta));

         // Record the original state
         int first0 = first;

         if (weights != NULL) choice /= weights[first / 2];

         for (first = 0; first < states; first++)
            {
            if (first >= second)
               choice += probability[first * (first + 1) / 2 + second];
            else
               choice += probability[second * (second + 1) / 2 + first];

            if (choice >= 0.0) break;
            }

         // Record outcomes for intermediate, uninformative, positions
         SamplePath(states, j, j0 + 1, first, first0, rand);
         FillPath(states + 1, j, j0 + 1, second);

         continue;
         }

      choice -= sum01 * theta * (1.0 - theta) / states;
      if (choice <= 0.0)
         {
         // The second haplotype changed ...
         choice = choice * states / (theta * (1.0 - theta));

         // Save the original state
         int second0 = second;

         if (weights != NULL) choice /= weights[second / 2];

         for (second = 0; second < states; second++)
            {
            if (first >= second)
               choice += probability[first * (first + 1) / 2 + second];
            else
               choice += probability[second * (second + 1) / 2 + first];

            if (choice >= 0.0) break;
            }

         // Record outcomes for intermediate, uninformative, positions
         FillPath(states, j, j0 + 1, first);
         SamplePath(states + 1, j, j0 + 1, second, second0, rand);

         continue;
         }

      // Try to select any other state
      choice *= states * states / (theta * theta);
      sum = 0.0;

      // Save the original states
      int first0 = first;
      int second0 = second;

      if (weights != NULL) choice /= weights[first / 2] * weights[second / 2];

      for (first = 0; first < states; first++)
         {
         for (second = 0; second <= first; second++, probability++)
            {
            sum += *probability;

            if (sum > choice) break;
            }

         if (second <= first) break;
         }

      if (rand->Binary())
         {
         int temp = first;
         first = second;
         second = temp;
         }

      // Record outcomes for intermediate, uninformative, positions
      SamplePath(states, j, j0 + 1, first, first0, rand);
      SamplePath(states + 1, j, j0 + 1, second, second0, rand);
      }

   ImputeAlleles(0, first, second, rand);
   }

bool ShotgunHaplotyper::ForceMemoryAllocation()
   {
   // Cycle through individuals, with the exact same steps as the actual
   // haplotyper and request memory ... by requesting all memory upfront,
   // we force crashes to happen early.
   for (int i = 0; i < individuals - phased; i++)
      {
      ResetMemoryPool();
      GetMemoryBlock(0);

      if (leftMatrices[0] == NULL)
         return false;

      int skipped = 0;
      for (int j = 1; j < markers; j++)
         //if (genotypes[i][j] != GENOTYPE_MISSING || j == markers - 1)
            {
            GetMemoryBlock(j);

            if (leftMatrices[j] == NULL)
               return false;
            }
         //else
            //skipped++;

      if (skipped == 0) break;
      }

   if (!phased)
      return true;

   ResetMemoryPool();
   for (int j = 0; j < markers; j++)
      {
      GetSmallMemoryBlock(j);

      if (leftMatrices[j] == NULL)
         return false;
      }

   return true;
   }


void ShotgunHaplotyper::SelectReferenceSet(int * array, int forWhom) {
  //fprintf(stderr,"ShotgunHaplotyper::SelectReferenceSet() called\n");
  /*
  if ( weightByMismatch ) {
    WeightByMismatch();
    globalRandom.Choose(array, weights, individuals - 1, states / 2);
    }*/
  if ( weightByLikelihood ) {
    //fprintf(stderr,"Weighting by likelihood for %d\n",forWhom);
    WeightByLikelihood();
    globalRandom.Choose(array, weights, individuals - 1, states / 2);
  }
  else if ( weightByLongestMatch ) {
    ChooseByLongestMatch(array);
    //fprintf(stderr,"Weighting by longest match for %d\n",forWhom);
    //WeightByLongestMatch();
    //globalRandom.Choose(array, weights, individuals - 1, states / 2);
  }
  else {
    if (greedy)
      {
	// Sanity check
	assert(states == phased * 2);
	
	// We exclude inferred haplotypes from the reference set
	for (int i = 0; i < individuals - phased; i++)
	  array[i] = 0;
	
	// We include phased haplotypes as our reference set
	for (int i = individuals - phased; i < individuals - 1; i++)
	  array[i] = 1;
	
	// For the last entry in the reference set, we may need to pick
	// a pair of inferred haplotypes
	if (forWhom < individuals - phased)
	  array[forWhom] = 1;
	else
	  array[globalRandom.NextInt() % (individuals - phased)] = 1;
      }
    else if (weights != NULL)
      globalRandom.Choose(array, weights, individuals - 1, states / 2);
    else
      globalRandom.Choose(array, individuals - 1, states / 2);
  }
  
  // Swap reference set into position
  for (int j = 0, out = 0; j < individuals; j++)
    if (array[j])
      SwapIndividuals(j, out++);
}
