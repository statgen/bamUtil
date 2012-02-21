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
    std::cerr << "\t\t--noeof       : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params      : Print the parameter settings to stderr" << std::endl;
}


int GapInfo::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    inputParameters.Read(argc-1, &(argv[1]));

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


    return(processFile(inFile.c_str(), outFile.c_str()));
}


int GapInfo::processFile(const char* inputFileName, const char* outputFileName)
{
    // Open the file for reading.
    SamFile samIn;
    samIn.OpenForRead(inputFileName);

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    SamRecord samRecord;

    IFILE outFile = ifopen(outputFileName, "w");

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
        if((mateStart == readStart) && (!SamFlag::isFirstFragment(samFlags)))
        {
            // read and mate start at the same position, so 
            // only process the first fragment.
            continue;
        }

        // Process this read pair.
        int32_t readEnd = samRecord.get0BasedAlignmentEnd();
        
        int32_t gapSize = mateStart - readEnd - 1;

        // Output the gap info.
        ifprintf(outFile, "%d\t%d", readEnd+1, gapSize);
        
        // Check if it is not the first or if it is not the forward strand.
        if(!SamFlag::isFirstFragment(samFlags))
        {
            ifprintf(outFile, "\tNotFirst");
        }
        if(SamFlag::isReverse(samFlags))
        {
            ifprintf(outFile, "\tReverse");
        }
        ifprintf(outFile, "\n");
    }

    SamStatus::Status returnStatus = samIn.GetStatus();
    if(returnStatus == SamStatus::NO_MORE_RECS)
    {
        return(SamStatus::SUCCESS);
    }
    return(returnStatus);
}
