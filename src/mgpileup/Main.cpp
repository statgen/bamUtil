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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>

#include "SamFile.h"
#include "PileupWithGenomeReference.h"
#include "PileupElementSummary.h"
#include "PileupElementBaseQual.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "InputFile.h"

int main(int argc, char ** argv)
{
    try 
    {
        std::string desc = "Example:\n\
./mgpileup -r /data/local/ref/karma.ref/human.g1k.v37 -b ../data/HG00160.chrom20.ILLUMINA.bwa.GBR.low_coverage.20100517.bam -v HG00160.chrom20.ILLUMINA.bwa.GBR.low_coverage.20100517.vcf -i -i ../data/LDL_b37.modified.genos.EUR.vcf -d \n\
Pileup takes in a BAM file and a genome reference file to return the following \n\
Alignment statistics are written to the VCF file specified, if the file name ends with .gz, the file output will be in gzip format.\n\
1. CHROM   : chromosome.\n\
2. POS     : position on chromosome.\n\
3. ID      : id.\n\
4. REF     : Reference base in the reference genome. A,C,T,G,N.\n\
5. ALT     : alt.\n\
6. QUAL    : quality score. \n\
7. FILTER  : filter.\n\
8. INFO    : info.\n\
9. FORMAT  : headers of custom data.\n\
   a. N         : number of contigs mapped at this locus.\n\
   b. BASE      : bases of respective contigs, '-' for deletions.\n\
   c. MAPQ      : quality score of the mapping for the contig, this is also shown for deletions.\n\
   d. BASEQ     : phred quality score, -1 for deletions.\n\
   e. STRAND    : F - forward, R - reverse.\n\
   f. CYCLE     : sequencing cycle, -1 for deletions.\n\
   g. GL        : genotype likelihood scores - AA,AC,AG,AT,CC,CG,CT,GG,GT,TT.\n\
10. <vcf output file name> : contains data described in FORMAT.\n";
   			
        std::string version = "1.0";
        TCLAP::CmdLine cmd(desc, ' ', version);
        TCLAP::ValueArg<std::string> argInputBAMFileName("b", "bam", "BAM file", true, "", "string");
        TCLAP::ValueArg<std::string> argRefSeqFileName("r", "reference", "Reference Sequence file", true, "", "string");
        TCLAP::ValueArg<std::string> argInputVCFFileName("i", "inputvcf", "VCF file listing the loci of interest (can be gzipped), bam index file is automatically assumed to be in the same location as the bam file.", false, "", "string");
        TCLAP::ValueArg<std::string> argPrevVCFFileName("p", "previousvcf", "VCF file from a prevoius pileup (can be gzipped)", false, "", "string");
        TCLAP::ValueArg<uint32_t> argMaxStoredLines("l", "maxStoredLines", "Maximum number of previous pileup lines to store before cutting off a bam region to analyze", false, 16000, "integer");
        TCLAP::ValueArg<uint32_t> argMaxRegionGap("g", "maxRegionGap", "Maximum gap between requested positions prior to breaking off a region", false, 16000, "integer");
        TCLAP::ValueArg<std::string> argOutputVCFFileName("v", "ouputvcf", "VCF file - if the extension is .gz, the written file will be a gzip file, (default is STDOUT)", false, "-", "string");
        TCLAP::SwitchArg argAddDelAsBase("d", "adddelasbase", "Adds deletions as base", cmd, false);
        TCLAP::SwitchArg argSummarize("s", "summarize", "Just print the summary of the allele counts for each position", cmd, false);

        cmd.add(argInputBAMFileName);
        cmd.add(argRefSeqFileName);
        cmd.add(argInputVCFFileName);
        cmd.add(argPrevVCFFileName);
        cmd.add(argMaxStoredLines);
        cmd.add(argMaxRegionGap);
        cmd.add(argOutputVCFFileName);
        cmd.parse(argc, argv);

        std::cerr << "Running mgpileup version " << version << std::endl; 
        std::string inputBAMFileName = argInputBAMFileName.getValue();
        std::cerr << "bam file                : " << inputBAMFileName << std::endl; 
        std::string refSeqFileName = argRefSeqFileName.getValue();	
        std::cerr << "reference sequence file : " << argRefSeqFileName.getValue() << std::endl; 
        std::string inputVCFFileName = argInputVCFFileName.getValue();
        uint32_t maxStoredLines = argMaxStoredLines.getValue();
        uint32_t maxRegionGap = argMaxRegionGap.getValue();
        std::string prevVCFFileName = argPrevVCFFileName.getValue();

        bool summarize = argSummarize.getValue();

        if(summarize)
        {
                 Pileup<PileupElementSummary> pileup(refSeqFileName);

                 pileup.processFile(inputBAMFileName);
        }
        else
        {
            if (inputVCFFileName != "")
            {
                std::cerr << "input VCF file          : " << inputVCFFileName << std::endl; 
            }
            if (prevVCFFileName != "")
            {
                std::cerr << "previous pileup VCF     : " << prevVCFFileName << std::endl; 
            }
            std::string outputVCFFileName = argOutputVCFFileName.getValue();
            if(outputVCFFileName[0] == '-')
            {
                std::cerr << "output VCF file         : STDOUT" << std::endl; 
            }
            else
            {
                std::cerr << "output VCF file         : " << outputVCFFileName << std::endl; 
            }

            std::cerr << "add deletions as bases  : " << (argAddDelAsBase.getValue()? "yes" : "no") << std::endl; 
            
            PileupWithGenomeReference<PileupElementBaseQual> pileup(1024, refSeqFileName, argAddDelAsBase.getValue());
            
            //process file with index    	
            if (inputVCFFileName != "")
            {
                pileup.processFile(inputBAMFileName, inputVCFFileName, prevVCFFileName, outputVCFFileName, maxStoredLines, maxRegionGap);
            }
            else
            {
                pileup.processFile(inputBAMFileName, outputVCFFileName);
            }
        }
    }
    catch (TCLAP::ArgException &e) 
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        abort();
    }

    return(0);
}
