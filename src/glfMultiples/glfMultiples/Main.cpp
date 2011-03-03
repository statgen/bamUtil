/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "BaseQualityHelper.h"
#include "GlfLikelihoods.h"
#include "StringArray.h"
#include "StringAlias.h"
#include "Parameters.h"
#include "Pedigree.h"
#include "Error.h"
#include "VcfEntry.h"
//#include "BgzfFileType.h"

#include <math.h>
#include <time.h>
#include <limits.h>

FILE * baseCalls = NULL;

double VcfEntry::phred2Err[256];
bool dummy = VcfEntry::initPhred2Error();

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
    char alleles[] = { 0, 'A', 'C', 'G', 'T' };

    int firstGlf = 0;
    while (glf[firstGlf].isStub)
        firstGlf++;

    printf("Dump for section %s, position %d [%c]\n",
           (const char *) glf[firstGlf].label, position, alleles[(int)refBase]);

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

void ReportSNP(glfHandler * glf, int n, int position,
               int refBase, int allele1, int allele2,
               const char * filter,
               int totalCoverage, int averageMapQuality, double posterior)
{
    if (baseCalls == NULL) return;

    if (allele2 == refBase)
    {
        int swap = allele1;
        allele1 = allele2;
        allele2 = swap;
    }

    //char alleles[] = { 0, 'A', 'C', 'G', 'T' };

    // Find the location of the first non-stub glf
    int firstGlf = 0;
    while (glf[firstGlf].isStub)
        firstGlf++;

    int quality = posterior > 0.9999999999 ? 100 : -10 * log10(1.0 - posterior);

    VcfEntry entry(glf[firstGlf].label.c_str(), position, refBase, allele1, allele2, quality, filter, n);

    GenotypeLikelihood lk;
    // Find best frequency
    lk.glf = glf;
    lk.n = n;
    lk.position = position;
    lk.SetAlleles(allele1, allele2);
    lk.OptimizeFrequency();

    entry.fillGenotypes(lk, glf, totalCoverage, averageMapQuality);
    entry.printSNP(baseCalls);
}

double FilteringLikelihood
(glfHandler * glf, int n, int position, int refAllele)
{
    FilterLikelihood lk;

    lk.glf = glf;
    lk.n = n;
    lk.position = position;
    lk.SetReferenceAllele(refAllele);
    lk.OptimizeFrequency();

    return -lk.fmin;
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
    printf("(c) 2008-2011 Goncalo Abecasis, Sebastian Zoellner, Yun Li\n\n");

    String pedfile;
    String positionfile;
    String callfile;
    String glfAliases;
    String glfPrefix;
    String glfSuffix;
    ParameterList pl;

    double posterior = 0.50;
    int    mapQuality = 0;
    int    minTotalDepth = 1;
    int    maxTotalDepth = INT_MAX;
    bool   verbose = false;
    bool   mapQualityStrict = false;
    bool   hardFilter = false;
    bool   smartFilter = false;
    bool   softFilter = true;
    bool   robustPrior = true;
    bool   uniformPrior = false;
    String xLabel("X"), yLabel("Y"), mitoLabel("MT");
    int    xStart = 2699520, xStop = 154931044;

    BEGIN_LONG_PARAMETERS(longParameters)
        LONG_PARAMETER_GROUP("Pedigree File")
        LONG_STRINGPARAMETER("ped", &pedfile)
        LONG_PARAMETER_GROUP("Map Quality Filter")
        LONG_INTPARAMETER("minMapQuality", &mapQuality)
        LONG_PARAMETER("strict", &mapQualityStrict)
        LONG_PARAMETER_GROUP("Total Depth Filter")
        LONG_INTPARAMETER("minDepth", &minTotalDepth)
        LONG_INTPARAMETER("maxDepth", &maxTotalDepth)
        LONG_PARAMETER_GROUP("Position Filter")
        LONG_STRINGPARAMETER("positionFile", &positionfile)
        LONG_PARAMETER_GROUP("Chromosome Labels")
        LONG_STRINGPARAMETER("xChr", &xLabel)
        LONG_STRINGPARAMETER("yChr", &yLabel)
        LONG_STRINGPARAMETER("mito", &mitoLabel)
        LONG_INTPARAMETER("xStart", &xStart)
        LONG_INTPARAMETER("xStop", &xStop)
        LONG_PARAMETER_GROUP("Filtering Options")
        EXCLUSIVE_PARAMETER("hardFilter", &hardFilter)
        EXCLUSIVE_PARAMETER("smartFilter", &smartFilter)
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

    Pedigree ped;
    if (!pedfile.IsEmpty())
    {
        ped.pd.AddStringColumn("glfFile");
        ped.Load(pedfile);

        n = ped.count;
    }
    else
        if (n == 0)
            error("No pedigree file present and no glf files listed at the end of command line\n");

    // Prior for finding difference from the reference at a particular site
    //BgzfFileType::setRequireEofBlock(false);

    double prior = 0.0;
    for (int i = 1; i <= 2 * n; i++)
        prior += 1.0 / i;
    prior *= 0.001;

    glfHandler * glf = new glfHandler[n];

    bool firstGlf = n;
    if (ped.count)
    {
        bool warn = false;

        for (int i = n - 1; i > 0; i++)
        {
          
            if (!glf[i].Open(ped[i].strings[0]))
            {
                printf("Failed to open genotype likelihood file [%s] for individual %s:%s\n",
                       (const char *) ped[i].strings[0],
                       (const char *) ped[i].famid,
                       (const char *) ped[i].pid);

                glf[i].OpenStub();
                firstGlf = i;
            }

            if (warn)
                printf("\n");

            if (firstGlf == n)
                error("No genotype likelihood files could be opened");
        }
    }
    else
    {
        for (int i = firstGlf = 0; i < n; i++)
        {
            String glfName = glfPrefix + String(argv[i]) + glfSuffix;
            if (!glf[i].Open(glfName))
                error("Failed to open genotype likelihood file [%s]\n", glfName.c_str());
        }
    }

    StringAlias aliases;
    aliases.ReadFromFile(glfAliases);

    printf("Calling genotypes for files ...\n");
    for (int i = 0; i < n; i++)
        printf("%s\n", ped.count ? (const char *) ped[i].strings[0] : argv[i]);
    printf("\n");

    baseCalls = fopen(callfile, "wt");

    if (baseCalls != NULL)
    {
        fprintf(baseCalls, "##fileformat=VCFv4.0\n");
        ReportDate(baseCalls);
        fprintf(baseCalls, "##source=glfMultiples\n");
        fprintf(baseCalls, "##minDepth=%d\n", minTotalDepth);
        fprintf(baseCalls, "##maxDepth=%d\n", maxTotalDepth);
        fprintf(baseCalls, "##minMapQuality=%d\n", mapQuality);
        fprintf(baseCalls, "##minPosterior=%.4f\n", posterior);
        fprintf(baseCalls, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
        fprintf(baseCalls, "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root Mean Squared Mapping Quality\">\n");
        fprintf(baseCalls, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with coverage\">\n");
        fprintf(baseCalls, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles (with coverage)\">\n");
        fprintf(baseCalls, "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Alternative allele count (with coverage)\">\n");
        fprintf(baseCalls, "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Alternate allele frequency\">\n");
        fprintf(baseCalls, "##INFO=<ID=AB,Number=1,Type=Float,Description=\"Estimated allele balance between the alleles\">\n");
	if ( mapQuality > 0 ) {
	  fprintf(baseCalls, "##FILTER=<ID=mq%d,Description=\"Mapping Quality less than %d\">\n",mapQuality,mapQuality);
	}
	if ( minTotalDepth > 1 ) {
	  fprintf(baseCalls, "##FILTER=<ID=dp%d,Description=\"Total Read Depth less than %d\">\n",minTotalDepth,minTotalDepth);
	}
	if ( minTotalDepth < INT_MAX ) {
	  fprintf(baseCalls, "##FILTER=<ID=DP%d,Description=\"Total Read Depth greater than %d\">\n",maxTotalDepth,maxTotalDepth);
	}
        fprintf(baseCalls, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Most Likely Genotype\">\n");
        fprintf(baseCalls, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Call Quality\">\n");
        fprintf(baseCalls, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
        fprintf(baseCalls, "##FORMAT=<ID=GL,Number=3,Type=Integer,Description=\"Genotype Likelihoods for Genotypes 0/0,0/1,1/1\">\n");
        fprintf(baseCalls, "##FORMAT=<ID=GL3,Number=6,Type=Integer,Description=\"Genotype Likelihoods for Genotypes 0/0,0/1,1/1,0/2,1/2,2/2\">\n");
        fprintf(baseCalls, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (int i = 0; i < n; i++)
            fprintf(baseCalls, "\t%s", ped.count ?
                    (const char *) (ped[i].famid + ":" + ped[i].pid) :
                    (const char *) aliases.GetAlias(argv[i]));

        fprintf(baseCalls, "\n");
    }

    StringArray buffer, tokens;
    StringHash  positions;

    buffer.Read(positionfile);

    for (int i = 0; i < buffer.Length(); i++)
    {
        tokens.ReplaceTokens(buffer[i], " \t:");

        if (tokens.Length() != 2) continue;

        positions.Add(tokens[0] + ":" + (int(tokens[1].AsInteger() - 1)));
    }

    int chromosomeType = 0;

    while (glf[firstGlf].NextSection())
    {
        for (int i = firstGlf + 1; i < n; i++)
        {
            if (glf[i].isStub) continue;

            glf[i].NextSection();

            if (glf[firstGlf].maxPosition != glf[i].maxPosition || glf[firstGlf].label != glf[i].label)
            {
                error("Genotype files '%s' and '%s' are not compatible ...\n"
                      "    File '%s' has section %s with %d entries ...\n"
                      "    File '%s' section %s with %d entries ...\n",
                      ped.count ? (const char *) ped[firstGlf].strings[0] : argv[firstGlf],
                      ped.count ? (const char *) ped[i].strings[0] : argv[i],
                      ped.count ? (const char *) ped[firstGlf].strings[0] : argv[firstGlf],
                      (const char *) glf[firstGlf].label, glf[firstGlf].maxPosition,
                      ped.count ? (const char *) ped[i].strings[0] : argv[i],
                      (const char *) glf[i].label, glf[i].maxPosition);
            }
        }

        chromosomeType = CT_AUTOSOME;

        if (ped.count)
        {
            if (glf[firstGlf].label == xLabel) chromosomeType = CT_CHRX;
            if (glf[firstGlf].label == yLabel) chromosomeType = CT_CHRY;
            if (glf[firstGlf].label == mitoLabel) chromosomeType = CT_MITO;
        }

        printf("Processing section %s with %d entries\n",
               (const char *) glf[firstGlf].label, glf[firstGlf].maxPosition);

        char refBase = 0;
        int position = 0;
        int mapQualityFilter = 0;
        int depthFilter = 0;
        int homozygousReference = 0;
        int transitions = 0;
        int transversions = 0;
        int otherPolymorphisms = 0;
        int sinkFilter = 0;
        int smartFilterHits = 0;
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
            for (int i = 1; i < n; i++)
                if (position > glf[i].position)
                {
                    position = glf[i].position;
                    refBase = glf[i].data.refBase;
                }

            // Avoid alignments that extend past the end of the chromosome
            if (position >= glf[firstGlf].maxPosition)
                break;

            baseCounts[(int)refBase]++;

            // These lines can be uncommented for debugging purposes
            // for (int i = 0; i < n; i++)
            //   printf("GLF %d : position %d, refBase %d\n", i, position, refBase);
            // printf("Position: %d, refBase: %d\n", position, refBase);

            if (positions.Entries())
            {
                filter = glf[firstGlf].label + ":" + position;

                if (positions.Find(filter) < 0) continue;
            }

            if (refBase == 0) continue;

            // Corrected calculation of root-mean-square Map Quality score 
            // and check if we have at least one sample with good quality data

            int     currentDepth = 0, totalDepth = 0, numCovered = 0;
            double  currentQuality = 0.0, averageMapQuality = 0.0;
            bool    passMapQualityFilter = false;

            for (int i = 0; i < n; i++)
            {
                currentDepth = glf[i].GetDepth(position);
                if (currentDepth != 0)
                {
                    totalDepth += currentDepth;
                    numCovered++;     // not currently used -- will be "NS"
                    currentQuality = glf[i].GetMapQuality(position);
                    averageMapQuality += 
                        currentDepth * currentQuality * currentQuality;
                    if (currentQuality >= mapQuality)
                        passMapQualityFilter = true;
                }
            }
            averageMapQuality = sqrt(averageMapQuality / totalDepth);


            filter.Clear();

            if (!passMapQualityFilter)
            {
                if (filter.Length() == 0) mapQualityFilter++;
                if (hardFilter) continue;
                filter.catprintf("%smq%d", filter.Length() ? ";" : "", mapQuality);
            }

            if (totalDepth < minTotalDepth)
            {
                if (filter.Length() == 0) depthFilter++;
                if (hardFilter) continue;
                filter.catprintf("%sdp%d", filter.Length() ? ";" : "", minTotalDepth);
            }

            if (totalDepth > maxTotalDepth)
            {
                if (filter.Length() == 0) depthFilter++;
                if (hardFilter) continue;
                filter.catprintf("%sDP%d", filter.Length() ? ";" : "", maxTotalDepth);
            }

            // Create convenient aliases for each base
            unsigned char transition = (((refBase - 1) ^ 2) + 1);
            unsigned char transvers1 = (((refBase - 1) ^ 3) + 1);
            unsigned char transvers2 = (((refBase - 1) ^ 1) + 1);

            int homRef = glf[0].GenotypeIndex(refBase, refBase);

            // Calculate likelihood assuming every is homozygous for the reference
            double lRef = log(1.0 - prior);
            for (int i = 0; i < n; i++)
                lRef += log(glf[i].GetLikelihoods(position)[homRef]);

            // Figure out the correct type of analysis
            //int cType = chromosomeType != CT_CHRX ? chromosomeType :
            //    position >= xStart && position <= xStop ? CT_CHRX : CT_AUTOSOME;

            // Calculate maximum likelihood for a variant
            if (smartFilter)
            {
                double anyVariant = log(prior) + FilteringLikelihood(glf, n, position, refBase);

                if (exp(lRef - anyVariant) > (1.0 - posterior)/posterior)
                {
                    smartFilterHits++;
                    continue;
                }
            }

	    //fprintf(stderr,"position = %d\n",position);

            double pTs = uniformPrior ? 1./3. : 2./3.;
            double pTv = uniformPrior ? 1./3. : 1./6.;

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

                if (positions.Entries())
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

        int actualBases = glf[firstGlf].maxPosition - baseCounts[0];

        printf("          Missing bases = %9d (%.3f%%)\n",
               baseCounts[0], baseCounts[0] * 100. / glf[firstGlf].maxPosition);
        printf("        Reference bases = %9d (%.3f%%)\n",
               glf[firstGlf].maxPosition - baseCounts[0], (glf[firstGlf].maxPosition - baseCounts[0]) * 100. / glf[firstGlf].maxPosition);

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

        if (smartFilter)
            printf("           Smart Filter = %9d bases (%.3f%%)\n",
                   smartFilterHits, smartFilterHits * 100. / actualBases);

        int noCalls = actualBases - homozygousReference - transitions - transversions - otherPolymorphisms - sinkFilter;
        printf("                No call = %9d bases (%.3f%%)\n",
               noCalls, noCalls * 100. / actualBases);

        fflush(stdout);
    }

    if (baseCalls != NULL)
        fclose(baseCalls);

    time(&t);
    printf("\nAnalysis completed on %s\n", ctime(&t));
    fflush(stdout);
}
