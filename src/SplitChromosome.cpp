/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan
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
// This file contains the processing for the executable option "splitChromosome"
// which splits a sorted/indexed BAM file into one file per chromosome.
#include "SplitChromosome.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"

void SplitChromosome::splitChromosomeDescription()
{
    std::cerr << " splitChromosome - Split BAM by Chromosome" << std::endl;
}


void SplitChromosome::description()
{
    splitChromosomeDescription();
}


void SplitChromosome::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam splitChromosome --in <inputFilename>  --out <outputFileBaseName> [--noeof] [--bamout|--samout] [--params]"<< std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in       : the BAM file to be split" << std::endl;
    std::cerr << "\t\t--out      : the base filename for the SAM/BAM files to write into.  Does not include the extension.\n";
    std::cerr << "                 CHROM.bam or CHROM.sam will be appended to the basename where CHROM is the chromosome name.\n";
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--noeof  : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--bamout : write the output files in BAM format (default)." << std::endl;
    std::cerr << "\t\t--samout : write the output files in SAM format." << std::endl;
    std::cerr << "\t\t--params : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}




int SplitChromosome::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFileBase = "";
    bool noeof = false;
    bool bamOut = false;
    bool samOut = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFileBase)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PARAMETER_GROUP("Output Type")
           EXCLUSIVE_PARAMETER("bamout", &bamOut)
           EXCLUSIVE_PARAMETER("samout", &samOut)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);

    if(!samOut && !bamOut)
    {
        bamOut = true;
    }

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
    if(outFileBase == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--out is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(params)
    {
        inputParameters.Status();
    }

    SamFile samIn;
    // Open the file for reading.   
    samIn.OpenForRead(inFile);

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    samIn.setSortedValidation(SamFile::COORDINATE);

    // Read the sam records.
    SamRecord samRecord;

    // Track the status.
    SamStatus::Status status = SamStatus::SUCCESS;

    // Do not open the file until there is a record to write into it.
    SamFileWriter outFile;

    String outputName;
    
    int32_t prevRefID = -2;

    int numSectionRecords = 0;

    // Keep reading records until they aren't more.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Check if this is the first record in the section
        if(samRecord.getReferenceID() != prevRefID)
        {
            if(numSectionRecords != 0)
            {
                std::cerr << "Reference Name: "
                          << samHeader.getReferenceLabel(prevRefID)
                          << " has " << numSectionRecords << " records\n";
            }
            outputName = outFileBase;

            // Open a new file (will close the previous one)
            const char* refName = samRecord.getReferenceName();
            if(strcmp(refName, "*") == 0)
            {
                outputName += "UnknownChrom";
            }
            else
            {
                outputName += refName;
            }
            // Append the extension.
            if(bamOut)
            {
                outputName += ".bam";
            }
            else
            {
                outputName += ".sam";
            }
            outFile.OpenForWrite(outputName.c_str());
            outFile.WriteHeader(samHeader);
            numSectionRecords = 0;
            prevRefID = samRecord.getReferenceID();
        }

        numSectionRecords++;            

        // Successfully read a record from the file, so write it.
        outFile.WriteRecord(samHeader, samRecord);
    }

    // Output the last chromosome's info.
    if(numSectionRecords != 0)
    {
        std::cerr << "Reference Name: "
                  << samHeader.getReferenceLabel(prevRefID)
                  << " has " << numSectionRecords << " records\n";
    }

    std::cerr << "Number of records = " << samIn.GetCurrentRecordCount() 
              << std::endl;

    fprintf(stderr, "Returning: %d (%s)\n", status, SamStatus::getStatusString(status));
    return(status);
}

