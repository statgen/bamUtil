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

//////////////////////////////////////////////////////////////////////////
#include "Filter.h"

#include "Parameters.h"
#include "BgzfFileType.h"
#include "GenomeSequence.h"
#include "SamFile.h"
#include "SamFilter.h"

void Filter::filterDescription()
{
    std::cerr << " filter - Filter reads by clipping ends with too high of a mismatch percentage and by marking reads unmapped if the quality of mismatches is too high" << std::endl;
}

void Filter::description()
{
    filterDescription();
}


void Filter::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam filter --in <inputFilename>  --refFile <referenceFilename>  --out <outputFilename> [--noeof] [--qualityThreshold <qualThresh>] [--defaultQualityInt <defaultQual>] [--mismatchThreshold <mismatchThresh>] [--params]"<< std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in       : the SAM/BAM file to be read" << std::endl;
    std::cerr << "\t\t--refFile  : the reference file" << std::endl;
    std::cerr << "\t\t--out      : the SAM/BAM file to write to" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--noeof             : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--qualityThreshold  : maximum sum of the mismatch qualities before marking\n"
              << "\t\t                      a read unmapped. (Defaults to 60)" << std::endl;
    std::cerr << "\t\t--defaultQualityInt : quality value to use for mismatches that do not have a quality\n" 
              << "\t\t                      (Defaults to 20)" << std::endl;
    std::cerr << "\t\t--mismatchThreshold : decimal value indicating the maximum ratio of mismatches to\n"
              << "\t\t                      matches and mismatches allowed before clipping from the ends\n"
              << "\t\t                      (Defaults to .10)" << std::endl;
    std::cerr << "\t\t--params            : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int Filter::execute(int argc, char **argv)
{
    String inFile = "";
    String refFile = "";
    String outFile = "-";
    bool noeof = false;
    bool params = false;

    uint32_t qualityThreshold = 60;
    uint32_t defaultQualityInt = 20;
    double mismatchThreshold = .10;

    // Read in the parameters.    
    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_PARAMETER("noeof", &noeof)
        LONG_INTPARAMETER("qualityThreshold", &qualityThreshold)
        LONG_INTPARAMETER("defaultQualityInt", &defaultQualityInt)
        LONG_DOUBLEPARAMETER("mismatchThreshold", &mismatchThreshold)
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
    if(outFile == "")
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

    // Open the reference.
    GenomeSequence reference(refFile);

    // Open the bam file.
    SamFile samIn;
    // Open the file for reading.   
    samIn.OpenForRead(inFile);

    // Open the output file.
    SamFile samOut;
    samOut.OpenForWrite(outFile);

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    // Write the header to the output file.
    samOut.WriteHeader(samHeader);

    // Read the sam records.
    SamRecord samRecord;

    SamStatus::Status returnStatus = SamStatus::SUCCESS;
    
    int clippedCount = 0;
    int mismatchThresholdFilterCount = 0;
    int qualityThresholdFilterCount = 0;

    // Keep reading records until they aren't anymore.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Store the original cigar and position for later reference.
        //        std::string origCigarString = samRecord.getCigar();
        //        int32_t origPosition = samRecord.get0BasedPosition();

        SamFilter::FilterStatus filterStatus = 
            SamFilter::clipOnMismatchThreshold(samRecord, reference,
                                               mismatchThreshold);

        // Check to see if the read was modified.
        if(filterStatus == SamFilter::FILTERED)
        {
            // The read was filtered, update the counter.
            ++ mismatchThresholdFilterCount;
        }
        else
        {
            // Check to see if the read was clipped.
            if(filterStatus == SamFilter::CLIPPED)
            {
                // The read was clipped, nothing to do other than
                // update the counter.
                ++clippedCount;
            }

            // Now filter on mismatch quality.
            filterStatus = 
                SamFilter::filterOnMismatchQuality(samRecord, reference,
                                                   qualityThreshold, 
                                                   defaultQualityInt);
            if(filterStatus == SamFilter::FILTERED)
            {
                // Filtered.
                ++qualityThresholdFilterCount;
            }
        }

        // Write the record.
        samOut.WriteRecord(samHeader, samRecord);

//         // If the read was clipped or filtered,
//         // write the read and the associated reference to a file.
//         if(filterStatus != SamFilter::NONE)
//         {
//             uint32_t chrStart = 
//                 reference.getGenomePosition(samRecord.getReferenceName());
//             std::string referenceString;
//             reference.getString(referenceString, 
//                                 chrStart + samRecord.get0BasedUnclippedStart(),
//                                 samRecord.getReadLength());
//             std:: cerr << samRecord.get0BasedPosition()
//                        << "\t"
//                        << samRecord.getCigar()
//                        << "\t"
//                        << samRecord.getSequence()
//                        << "\t"
//                        << samRecord.getQuality()
//                        << "\n"
//                        << origPosition
//                        << "\t"
//                        << origCigarString
//                        << "\t"
//                        << referenceString
//                        << "\t"
//                        << samRecord.getQuality()
//                        << "\n\n";
//         }
    }
    
    std::cerr << "Number of Reads Clipped by Filtering: "
              << clippedCount << "\n";
    std::cerr << "Number of Reads Filtered Due to MismatchThreshold: "
              << mismatchThresholdFilterCount << "\n";
    std::cerr << "Number of Reads Filtered Due to QualityThreshold: "
              << qualityThresholdFilterCount << "\n";

    return(returnStatus);
}
