#ifndef __VCF_ENTRY_H
#define __VCF_ENTRY_H

#include "StringArray.h"
#include "Error.h"

class VcfEntry {
 public:
  static double phred2Err[256];
  String chrom;
  int pos;
  int ref;
  int al1;
  int al2;
  int qual;
  String filter;
  String info;
  String format;
  int num;
  StringArray genotypes;

 VcfEntry(const char* c, int p, int r, int a1, int a2, int q, const char* f, int n) :
  chrom(c), pos(p), ref(r), al1(a1), al2(a2), qual(q), filter(f), num(n) {
    //genotypes.Grow(n);
  }

  static bool initPhred2Error() {
    for(int i=0; i < 256; ++i) 
      phred2Err[i] = pow(.1, i/10.);
    return true;
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

  void fillGenotypes(GenotypeLikelihood& lk, glfHandler* glf, int totalCoverage, int averageMapQuality) {
    //fprintf(stderr,"starting fillGenotypes()\n");

    double priors[10];
    lk.GetPriors(priors, lk.min);

    int geno11 = glfHandler::GenotypeIndex(al1, al1);
    int geno12 = glfHandler::GenotypeIndex(al1, al2);
    int geno22 = glfHandler::GenotypeIndex(al2, al2);

    int label1 = al1 == ref ? 0 : 1;
    int label2 = al2 == ref ? 0 : al1 == al2 ? label1 : label1 + 1;

    String label11, label12, label22;
    label11.printf("%d/%d", label1, label1);
    label12.printf("%d/%d", label1, label2);
    label22.printf("%d/%d", label2, label2);

    int ACs[3] = {0,0,0};
    double ABnum = 0.5;
    double ABden = 1.0;
    double af = 1.0 - lk.min;
    double GPs[3] = {(1.-af)*(1.-af), 2.*af*(1.-af), af*af};

    for (int i = 0; i < num; i++)
    {
      //fprintf(stderr,"i = %d, num = %d, ACs={%d,%d,%d}\n",i,num,ACs[0],ACs[1],ACs[2]);

      // Report on the best genotype for the current SNP model
      int quality = GetBestRatio(glf[i].GetLikelihoods(pos), priors);
      int best = GetBestGenotype(glf[i].GetLikelihoods(pos), priors);
      int depth = glf[i].GetDepth(pos);
      int PLs[3] = {glf[i].GetLogLikelihoods(pos)[geno11],glf[i].GetLogLikelihoods(pos)[geno12],glf[i].GetLogLikelihoods(pos)[geno22]};

      int nrefDen = PLs[2]+PLs[0]-2*PLs[1]+6*depth;
      if ( nrefDen < 4 ) nrefDen = 4;
      if ( nrefDen < abs(PLs[0]-PLs[2]) )
	nrefDen = abs(PLs[0]-PLs[2]);

      double nref = 0.5 * depth * (1.0 + (double)(PLs[2]-PLs[0])/(double)nrefDen);
      double pHet = phred2Err[PLs[1]]*GPs[1]/(GPs[0]*phred2Err[PLs[0]] + GPs[1]*phred2Err[PLs[1]] + GPs[2]*phred2Err[PLs[2]]);
      ABnum += (pHet * nref);
      ABden += (pHet * depth);

      //fprintf(stderr,"foo %d\n",genotypes.Length());
      genotypes.Add("");
      //genotypes[i].Clear();
      if ( best == geno11 ) {
	genotypes[i] += label11;
	if ( depth > 0 ) ++ACs[0];
      }
      else if ( best == geno12 ) {
	genotypes[i] += label12;
	if ( depth > 0 ) ++ACs[1];
      }
      else {
	genotypes[i] += label22;
	if ( depth > 0 ) ++ACs[2];
      }
      //fprintf(stderr,"goo\n");

      genotypes[i].catprintf(":%d:%d",depth,quality);


      if (label1 == 0 && label2 == 0)
	continue;
      if (label2 < 2)
	genotypes[i].catprintf(":%d,%d,%d",PLs[0],PLs[1],PLs[2]);
      else {
	int genoRR, genoR1, genoR2;
        genoRR = glfHandler::GenotypeIndex(ref, ref);
        genoR1 = glfHandler::GenotypeIndex(ref, al1);
        genoR2 = glfHandler::GenotypeIndex(ref, al2);
        genotypes[i].catprintf(":%d,%d,%d,%d,%d,%d",
                    glf[i].GetLogLikelihoods(pos)[genoRR],
                    glf[i].GetLogLikelihoods(pos)[genoR1],
                    PLs[0],
                    glf[i].GetLogLikelihoods(pos)[genoR2],
                    PLs[1],
		    PLs[2]);
      }
    }
    double AB = ABnum/ABden;

    // fill the INFO field with
    // DP,MQ,NS,AN,AC,AF
    info.Clear();
    if ( al1 == ref ) {
      info.catprintf("DP=%d;MQ=%d;NS=%d;AN=%d;AC=%d;AF=%.6lf;AB=%.4lf",totalCoverage,averageMapQuality,ACs[0]+ACs[1]+ACs[2],2*(ACs[0]+ACs[1]+ACs[2]),ACs[1]+2*ACs[2],af,AB);
    }
    else {
      info.catprintf("DP=%d;MQ=%d;NS=%d;AN=%d;AC=%d,%d;AF=%.6lf,%.6lf;AB=%.4lf",totalCoverage,averageMapQuality,ACs[0]+ACs[1]+ACs[2],2*(ACs[0]+ACs[1]+ACs[2]),2*ACs[0]+ACs[1],ACs[1]+2*ACs[2],1-af,af,AB);
    }
  }

  void printSNP(FILE* f) {
    char alleles[] = { 0, 'A', 'C', 'G', 'T' };
    // CHROM POS ID
    fprintf(f, "%s\t%d\t.\t", chrom.c_str(), pos+1);
    // REF
    fprintf(f, "%c\t", alleles[ref]);
    // ALT
    int nalleles = 1;
    if (al1 != ref)
        fprintf(f, "%c", alleles[al1]), nalleles++;
    if (al2 != ref && al2 != al1)
        fprintf(f, "%s%c", nalleles > 1 ? "," : "", alleles[al2]), nalleles++;
    if (nalleles == 1)
      fprintf(f, ".");  // monomorphic SNP

    // QUAL
    fprintf(f, "\t%d\t", qual);
    // FILTER
    fprintf(f, "%s\t", filter.Length() == 0 ? "PASS" : filter.c_str());
    // INFO
    fprintf(f, "%s\t", info.Length() == 0 ? "." : info.c_str());
    // FORMAT
    fprintf(f, "GT:DP:GQ");
    if ((al2 != ref) || (al1 != ref))
    {
        fprintf(f, ":PL%s", al1 == ref ? "" : "3");
    }
    // GENOTYPES
    if ( num != genotypes.Length() )
      error("Internal error - num and genotypes.Length() does not match");

    for(int i=0; i < num; ++i) {
      fprintf(f,"\t%s",genotypes[i].c_str());
    }
    fprintf(f,"\n");
  }
};

#endif
