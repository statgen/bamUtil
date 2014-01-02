/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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

//////////////////////////////////////////////////////////////////////////
#include "GapInfo.h"
#include "SamFile.h"
#include "BgzfFileType.h"
#include "SamFlag.h"

GapInfo::GapInfo()
    : BamExecutable()
{
}


void GapInfo::gapInfoDescription()
{
    std::cerr << " gapInfo - Print information on the gap between read pairs in a SAM/BAM File." << std::endl;
}


void GapInfo::description()
{
    gapInfoDescription();
}


void GapInfo::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam gapInfo --in <inputFile> --out <outputFile> [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in          : the SAM/BAM file to print read pair gap info for" << std::endl;
    std::cerr << "\t\t--out         : the output file to be written" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--refFile     : reference file, used to skip gaps that include reference base 'N' (for runs without --detailed)";
    std::cerr << "\t\t--detailed    : Print  the details for each read pair" << std::endl;
    std::cerr << "\tOptional Parameters for the Detailed Option:" << std::endl;
    std::cerr << "\t\t--checkFirst  : Check the first in pair flag and print \"NotFirst\" if it isn't first" << std::endl;
    std::cerr << "\t\t--checkStrand : Check the strand flag and print \"Reverse\" if it is reverse complimented" << std::endl;
    std::cerr << "\t\t--noeof       : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params      : Print the parameter settings to stderr" << std::endl;
}


int GapInfo::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    String refFile = "";
    bool detailed = false;
    bool checkFirst = false;
    bool checkStrand = false;
    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_PARAMETER("detailed", &detailed)
        LONG_PARAMETER_GROUP("Optional Detailed Parameters")
        LONG_PARAMETER("checkFirst", &checkFirst)
        LONG_PARAMETER("checkStrand", &checkStrand)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);

    // If no eof block is required for a bgzf file, set the bgzf file type to 
    // not look for it.
    if(noeof)
    {
        // Set that the eof block is not required.
        BgzfFileType::setRequireEofBlock(false);
    }

    // Check to see if the in file was specified, if not, report an error.
    if(inFile == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--in is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    // Check to see if the out file was specified, if not, report an error.
    if(outFile == "")
    {
        usage();
        inputParameters.Status();
        // Out file was not specified but it is mandatory.
        std::cerr << "--out is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(params)
    {
        inputParameters.Status();
    }

    return(processFile(inFile.c_str(), outFile.c_str(),
                       refFile, detailed, 
                       checkFirst, checkStrand));
}


int GapInfo::processFile(const char* inputFileName, const char* outputFileName,
                         const char* refFile, bool detailed,
                         bool checkFirst, bool checkStrand)
{
    // Open the file for reading.
    SamFile samIn;
    samIn.OpenForRead(inputFileName);

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    SamRecord samRecord;

    GenomeSequence* refPtr = NULL;
    if(strcmp(refFile, "") != 0)
    {
        refPtr = new GenomeSequence(refFile);
    }

    IFILE outFile = ifopen(outputFileName, "w");

    // Map for summary.
    std::map<int, int> gapInfoMap;


    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        uint16_t samFlags = samRecord.getFlag();

        if((!SamFlag::isMapped(samFlags)) || 
           (!SamFlag::isMateMapped(samFlags)) ||
           (!SamFlag::isPaired(samFlags)) ||
           (samFlags & SamFlag::SECONDARY_ALIGNMENT) || 
           (SamFlag::isDuplicate(samFlags)) ||
           (SamFlag::isQCFailure(samFlags)))
        {
            // unmapped, mate unmapped, not paired,
            // not the primary alignment,
            // duplicate, fails vendor quality check 
            continue;
        }

        // No gap info if the chromosome names are different or
        // are unknown.
        int32_t refID = samRecord.getReferenceID();
        if((refID != samRecord.getMateReferenceID()) || (refID == -1))
        {
            continue;
        }

        int32_t readStart = samRecord.get0BasedPosition();
        int32_t mateStart = samRecord.get0BasedMatePosition();

        // If the mate starts first, then the pair was processed by
        // the mate.
        if(mateStart < readStart)
        {
            continue;
        }
        if((mateStart == readStart) && (SamFlag::isReverse(samFlags)))
        {
            // read and mate start at the same position, so 
            // only process the forward strand.
            continue;
        }

        // Process this read pair.
        int32_t readEnd = samRecord.get0BasedAlignmentEnd();
        
        int32_t gapSize = mateStart - readEnd - 1;

        if(detailed)
        {
            // Output the gap info.
            ifprintf(outFile, "%s\t%d\t%d", 
                     samRecord.getReferenceName(), readEnd+1, gapSize);
            
            // Check if it is not the first or if it is not the forward strand.
            if(checkFirst && !SamFlag::isFirstFragment(samFlags))
            {
                ifprintf(outFile, "\tNotFirst");
            }
            if(checkStrand && SamFlag::isReverse(samFlags))
            {
                ifprintf(outFile, "\tReverse");
            }
            ifprintf(outFile, "\n");
        }
        else
        {
            // Summary.
            // Skip reads that are not the forward strand.
            if(SamFlag::isReverse(samFlags))
            {
                // continue
                continue;
            }

            // Forward.
            // Check the reference for 'N's.
            if(refPtr != NULL)
            {
                genomeIndex_t chromStartIndex = 
                    refPtr->getGenomePosition(samRecord.getReferenceName());
                if(chromStartIndex == INVALID_GENOME_INDEX)
                {
                    // Invalid position, so continue to the next one.
                    continue;
                }
                bool skipRead = false;
                for(int i = readEnd + 1; i < mateStart; i++)
                {
                    if((*refPtr)[i] == 'N')
                    {
                        // 'N' in the reference, so continue to the next read.
                        skipRead = true;
                        break;
                    }
                }
                if(skipRead)
                {
                    continue;
                }
            }
            
            // Update the gapInfo.
            gapInfoMap[gapSize]++;
        }
    }

    if(!detailed)
    {
        // Output the summary.
        ifprintf(outFile, "GapSize\tNumPairs\n");
        for(std::map<int,int>::iterator iter = gapInfoMap.begin(); 
            iter != gapInfoMap.end(); iter++)
        {
            ifprintf(outFile, "%d\t%d\n", (*iter).first, (*iter).second);
        }
    }
    

    SamStatus::Status returnStatus = samIn.GetStatus();
    if(returnStatus == SamStatus::NO_MORE_RECS)
    {
        return(SamStatus::SUCCESS);
    }
    return(returnStatus);
}
