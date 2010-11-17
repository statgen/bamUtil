////////////////////////////////////////////////////////////////////// 
// mach1/Haplotyper.h 
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
 
#ifndef __HAPLOTYPER_H__
#define __HAPLOTYPER_H__

#include "Random.h"
#include "StringBasics.h"
#include "ErrorRate.h"
#include "InputFile.h"

#ifdef   __DOUBLE_HAPLOTYPING__
#define  float double
#endif

#include <stdlib.h>
#include <stdio.h>

class Pedigree;
class String;

class Haplotyper
   {
   public:
      // These are the basic variables that store information
      // about the underlying state

      int  phased, individuals, states, markers;

      char ** genotypes;
      char ** haplotypes;
      float * thetas;
      float * distances;

      // Store information about estimated error rates
      // and about observed mismatches between mosaic and
      // actual genotypes
      Errors * error_models;

      // These additional variables store optional information
      // about disease status for non-parametric mapping

      int      diseaseCount;
      char **  diseaseStatus;
      float ** diseaseScores;
      float ** nplScores;

      // Determines whether NPL scores are updated simultaneously with
      // haplotype sources
      bool updateDiseaseScores;

      // These flags report the current status and determine the
      // analytical strategy

      bool readyForUse;
      bool greedy;
      bool economyMode;

      // These variables store information about the presence of partially
      // ordered genotypes, both globally and at the individual level

      bool  orderedGenotypes;
      int * orderedGenotypeFlags;

      Haplotyper();
      ~Haplotyper();

      // Set the error rate. This is an omnibus rate which combines
      // the effects of mutation, gene-conversion and genotyping error.
      // With current genotyping technologies a value of 0.001 or higher
      // is reasonable
      void SetErrorRate(float rate);
      void SetErrorRate(int marker, float rate);

      double GetErrorRate(int marker)
         { return error_models[marker].rate; }
      double GetErrorRate();

      // Initialization function, allocates memory
      // necessary for tracking haplotypes and
      // genotypes
      bool AllocateMemory(int nIndividuals, int maxHaplos, int nMarkers);
      bool AllocateDiseaseStatus(int nDiseases);
      bool AllocateDistances();
      bool AllocateMLEMemory();
      bool AllocateRightMatrices();

      // By default, some memory allocation is delayed and carried out
      // on an as needed basis. This routine forces these allocations
      // to happen, based on current genotype data
      bool ForceMemoryAllocation();

      // Functions for initializing haplotype list
      // based on observed haplotypes
      virtual void RandomSetup(Random * rand = NULL);

      // Basic markov chain functionality (for genotypes)
      void Transpose(float * source, float * dest, float theta);
      void Transpose(float * source, float * dest, float * priors, float theta);

      virtual void ConditionOnData(float * matrix, int marker, char genotype);
      void ConditionOnData(int * stateKey, float * matrix, char genotype);

      void ScoreLeftConditional();
      void SampleChromosomes(Random * rand);

      void SetupPrior(float * matrix);

      // Routines for estimating MLEs for NPL scores
      void IntegrateNPL();
      void IntegrateNPL(float * matrix, int marker);
      void IntegrateNPL(float * left, float * right, int marker);

      // Routines for estimating MLEs for missing genotypes
      void ImputeGenotypes();
      void ImputeGenotypes(float * matrix, int marker);
      void ImputeGenotypes(float * left, float * right, int marker);

      // Analogous estimation routines that deal with partially ordered genotypes
      void ImputeGenotypesFromOrderedData();
      void ImputeGenotypesFromOrderedData(float * matrix, int marker);
      void ImputeGenotypesFromOrderedData(float * left, float * right, int marker);

      void OutputMLEs(Pedigree & ped, const String & prefix, bool detailed);

      //  Basic markov chain functionality (for phase known haplotypes)
      void TransposeHaplotype(float * source, float * dest, float theta);
      void ConditionHaplotypeOnData(float * matrix, int marker, char allele);
      void ScoreLeftConditionalForHaplotype();
      void SampleHaplotypeSource(Random * rand);

      // Markov chain functionality for dealing partially phased data
      void TransposeOrdered(float * source, float * dest, float theta);

      void ConditionOnOrderedData(float * matrix, int marker, char genotype);

      void ScoreLeftConditionalForOrderedGenotypes();
      void SampleChromosomesFromOrderedData(Random * rand);

      void SetupOrderedPrior(float * matrix);

      // Higher level markov chain functionality
      void WarmUp(int seeds, int rounds);
      void LoopThroughChromosomes();

      // Build a set of consensus haplotypes
      void BuildConsensus(int samples);

      // These functions update parameters based on the last iteration
      void   UpdateThetas();
      void   UpdateThetasWithDistances();
      double UpdateErrorRate();

      // These functions allow for different weights to be placed on each
      // possible haplotype. Currently these weights are simply based on
      // the number of positions where genotype data is available.
      virtual void CalculateWeights();
      void   AllocateWeights();
      void   FreeWeights();
      void   ScaleWeights();

      int  TotalCrossovers();

      // Report memory used by haplotyping engine
      void ShowMemoryInfo();
      void ShowMLEMemoryInfo();

      static void EstimateMemoryInfo(int Individuals, int Markers, int States, bool Compact, bool Phased);
      static void EstimateDiseaseMemoryInfo(int Individuals, int Markers, int Diseases);
      static void EstimateMLEMemoryInfo(int Individuals, int Markers, int Diseases);

      // Retrieve parameters from file
      bool LoadCrossoverRates(const char * filename);
      void LoadCrossoverRates(IFILE file);
      bool LoadErrorRates(const char * filename);
      void LoadErrorRates(IFILE file);

   protected:
      float *  marginals;
      float ** leftMatrices;
      float *  leftProbabilities;

      // These are used to impute missing genotypes in a forward-backward algorithm
      float ** rightMatrices;
      float ** posterior;
      float ** mlinfo;

      float *  weights;

      int   *  crossovers;

      // Impute two alleles at a particular marker, given sampled state
      virtual void  ImputeAlleles(int marker, int state1, int state2, Random * rand);
      void  ImputeAllele(int haplotype, int marker, int state);

      // Merge a set of sampled haplotypes into a consensus pair
      void  BuildConsensus(char ** consensus, char ** haplotypes, int count);

      float & Penetrance(int marker, int true_genotype, int observed_genotype)
         {
         return penetrances[marker * 9 + true_genotype * 3 + observed_genotype];
         }

      void FillPath(int haplotype, int fromMarker, int toMarker, int state);
      void SamplePath(int haplotype, int fromMarker, int toMarker, int fromState, int toState, Random * rand);

      // Routines for producing summary information about MLE estimates of genotypes
      void ResetMarkerInfo();
      void UpdateMarkerInfo();
      void OutputMarkerInfo(FILE * output);

      void NormalizePosterior(int marker)
         {
         double sum = posterior[0][marker] + posterior[1][marker] + posterior[2][marker];

         if (sum > 0.0)
            {
            posterior[0][marker] /= sum;
            posterior[1][marker] /= sum;
            posterior[2][marker] /= sum;
            }
         }

      bool skipSanityCheck;

   private:
      void Swap(char * & array1, char * & array2)
         { char * temp = array1; array1 = array2; array2 = temp; }

      void Swap(float * & array1, float * & array2)
         { float * temp = array1; array1 = array2; array2 = temp; }

      void SelectReferenceSet(int * choices, int forWhom);
      void SwapIndividuals(int a, int b);
      void SwapHaplotypes(int a, int b);

      void ScoreNPL();
      void UpdateDiseaseScores(int marker, int state);

      bool SanityCheck();

      void Print(int marker);
      void Print();

      float max(float a, float b)  { return a > b ? a : b; };
      float max(float a, float b, float c) { return max(max(a, b), c); }
      float max(float a, float b, float c, float d) { return max(max(a, b, c), d); }

      float square(float a) { return a * a; }

      float * penetrances;

      void  ResetCrossovers();

      // A series of memory management functions lets us delay allocation
      // of big blocks of memory until they are needed. Even more importantly
      // it allows us to reuse memory as needed.
      float ** memoryBlock;
      float ** smallMemoryBlock;
      int   *  stack, stackPtr, savedStackPtr;
      int      smallFree, nextAvailable, nextSmallAvailable, nextReuseable;
      int      gridSize;

      // This is the low level allocator
      float * AllocateMemoryBlock();

      // This retrieves a large block, used for modeling unphased genotypes
      float * GetLargeBlock();
      float * GetReuseableBlock();

      // This retrieves a smaller block, used for modeling a single haplotype
      float * GetSmallBlock();

      // This resets the memory pool
      void ResetMemoryPool();
      void ResetReuseablePool();

      // These commands allow us to run multiple passes through the data
      void MarkMemoryPool();
      void RewindMemoryPool();

      // These are the high level interfaces ...
      void GetMemoryBlock(int marker);
      void GetSmallMemoryBlock(int marker);
      void RetrieveMemoryBlock(int marker);
   };

#ifdef   __DOUBLE_HAPLOTYPING__
#undef   float
#endif

#define    GENOTYPE_MISSING                0
#define    GENOTYPE_HOMOZYGOUS_FOR_ONE     1
#define    GENOTYPE_HETEROZYGOUS           2
#define    GENOTYPE_HOMOZYGOUS_FOR_TWO     3

#define    GENOTYPE_ORDERED               16
#define    FIRST_ALLELE_ONE               1
#define    FIRST_ALLELE_TWO               2
#define    FIRST_ALLELE                   (1 | 2)
#define    SECOND_ALLELE_ONE              4
#define    SECOND_ALLELE_TWO              8
#define    SECOND_ALLELE                  (4 | 8)

#define    GENOTYPE_LINKED                32

#endif

 
