/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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
// This file contains the processing for the executable option "stats"
// which generates some statistics for SAM/BAM files.
#include "Stats.h"
#include "SamFile.h"
#include "BgzfFileType.h"

void Stats::statsDescription()
{
    std::cerr << " stats - Stats a SAM/BAM File" << std::endl;
}


void Stats::description()
{
    statsDescription();
}


void Stats::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam stats --in <inputFile> [--disableStatistics] [--noeof] [--params] [--qual] [--maxNumQuals <maxNum>] [--unmappedQual] [--bamIndex <bamIndexFile>]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in : the SAM/BAM file to be statsd" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--disableStatistics : Turn off statistic generation" << std::endl;
    std::cerr << "\t\t--noeof             : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params            : Print the parameter settings" << std::endl;
    std::cerr << "\tQuality Statistics Parameters" << std::endl;
    std::cerr << "\t\t--qual              : Generate a count for each quality" << std::endl;
    std::cerr << "\t\t--maxNumQuals       : Maximum number of qualities to include in the counts" << std::endl;
    std::cerr << "\t\t                      Defaults to 100,000.  Use -1 to indicate all qualities." << std::endl;
    std::cerr << "\t\t--unmappedQual      : Only count qualities in unmapped reads" << std::endl;
    std::cerr << "\t\t                      Automatically sets qual & disableStatistics" << std::endl;
    std::cerr << "\t\t--bamIndex          : the path/name of the bam index file" << std::endl;
    std::cerr << "\t\t                      (if not specified, uses the --in value + \".bai\")" << std::endl;
    std::cerr << std::endl;
}


int Stats::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String indexFile = "";
    bool disableStatistics = false;
    bool noeof = false;
    bool params = false;
    bool qual = false;
    int maxNumQuals = 100000;
    bool unmappedQual = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_PARAMETER("disableStatistics", &disableStatistics)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PARAMETER_GROUP("Quality Histogram")
        LONG_PARAMETER("qual", &qual)
        LONG_INTPARAMETER("maxNumQuals", &maxNumQuals)
        LONG_PARAMETER("unmappedQual", &unmappedQual)
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    inputParameters.Read(argc-1, &(argv[1]));

    if(unmappedQual)
    {
        // If unmappedQual was set, set qual to true.
        qual = true;
        disableStatistics = true;
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
        std::cerr << "--in is a mandatory argument for stats, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    // IndexFile is only used if unmappedQual is specified
    if(unmappedQual && (indexFile == ""))
    {
        // In file was not specified, so set it to the in file
        // + ".bai"
        indexFile = inFile + ".bai";
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Open the file for reading.
    SamFile samIn;
    if(!samIn.OpenForRead(inFile))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    // Set that statistics should be generated.
    samIn.GenerateStatistics(!disableStatistics);

    // Read the sam header.
    SamFileHeader samHeader;
    if(!samIn.ReadHeader(samHeader))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    // Since the bamIndex file is not requried, only open
    // if necessary and catch if there is a failure.
    // Open the bam index file for reading.
    if(unmappedQual)
    {
        try
        {
            samIn.ReadBamIndex(indexFile);
            samIn.SetReadSection(-1);
            // Reset unmappedQual since nothing special needs
            // to be done from here on out.
            unmappedQual = false;
        }
        catch(std::exception& e)
        {
        }
    }

    // Read the sam records.
    SamRecord samRecord(ErrorHandler::RETURN);

    int numQual = 0;

    // Quality histogram.
    const int MAX_QUAL = 126;
    const int START_QUAL = 33;
    int hist[MAX_QUAL+1];
    for(int i = 0; i <= MAX_QUAL; i++)
    {
        hist[i] = 0;
    }

    // Keep reading records from the file until SamFile::ReadRecord
    // indicates to stop (returns false).
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // See if the quality histogram should be genereated.
        if(qual && ((maxNumQuals < 0) || (numQual < maxNumQuals)))
        {
            // Check if we need to check for unmapped.
            if(unmappedQual && (samRecord.getReferenceID() != -1))
            {
                continue;
            }
            // Get the quality.
            const char* qual = samRecord.getQuality();
            int index = 0;
            // Check for no quality ('*').
            if((qual[0] == '*') && (qual[1] == 0))
            {
                // This record does not have a quality string.
                continue;
            }

            while(((maxNumQuals < 0) || (numQual < maxNumQuals)) && (qual[index] != 0))
            {
                ++numQual;
                // Check for valid quality.
                if((qual[index] < START_QUAL) || (qual[index] > MAX_QUAL))
                {
                    std::cerr << "Invalid Quality found: " << qual[index] 
                              << ".  Must be between 0 and " << MAX_QUAL
                              << ".\n";
                    ++index;
                    continue;
                }
                
                // Increment the count for this quality.
                ++(hist[(int)(qual[index])]);
                ++index;
            }
        }
        if(disableStatistics && qual && (numQual >= maxNumQuals) && (maxNumQuals > 0))
        {
            // No more stats, so quit.
            break;
        }
    }

    std::cerr << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl << std::endl;

    if(!disableStatistics)
    {
        samIn.PrintStatistics();
        std::cerr << std::endl;
    }

    // Print the histogram.
    if(qual)
    {
        std::cerr << "Quality\tNum\n";
        for(int i = START_QUAL; i <= MAX_QUAL; i++)
        {
            std::cerr << i << "\t" << hist[i] << std::endl;
        }
        std::cerr << std::endl;
    }

    SamStatus::Status status = samIn.GetStatus();
    if(status == SamStatus::NO_MORE_RECS)
    {
        // A status of NO_MORE_RECS means that all reads were successful.
        status = SamStatus::SUCCESS;
    }

    fprintf(stderr, "Returning: %d (%s)\n", status, SamStatus::getStatusString(status));
    return(status);
}

