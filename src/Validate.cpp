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
// This file contains the processing for the executable option "validate"
// which reads and validates SAM/BAM file and can generate some statistics
// from it.
#include "Validate.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"
#include "SamValidation.h"

void Validate::validateDescription()
{
    std::cerr << " validate - Validate a SAM/BAM File" << std::endl;
}


void Validate::description()
{
    validateDescription();
}


void Validate::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam validate --in <inputFile> [--noeof] [--so_flag|--so_coord|--so_query] [--maxErrors <numErrors>] [--verbose] [--printableErrors <numReportedErrors>] [--disableStatistics] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in : the SAM/BAM file to be validated" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--noeof             : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--refFile           : the reference file" << std::endl;
    std::cerr << "\t\t--so_flag           : validate the file is sorted based on the header's @HD SO flag." << std::endl;
    std::cerr << "\t\t--so_coord          : validate the file is sorted based on the coordinate." << std::endl;
    std::cerr << "\t\t--so_query          : validate the file is sorted based on the query name." << std::endl;
    std::cerr << "\t\t--maxErrors         : Number of records with errors/invalids to allow before quiting." << std::endl;
    std::cerr << "\t\t                      -1 (default) indicates to not quit until the entire file is validated." << std::endl;
    std::cerr << "\t\t                      0 indicates not to read/validate anything." << std::endl;
    std::cerr << "\t\t--verbose           : Print specific error details rather than just a summary" << std::endl;
    std::cerr << "\t\t--printableErrors   : Maximum number of records with errors to print the details of\n"
              << "\t\t                      before suppressing them when in verbose (defaults to 100)" 
              << std::endl;
    std::cerr << "\t\t--disableStatistics : Turn off statistic generation" << std::endl;
    std::cerr << "\t\t--params            : Print the parameter settings" << std::endl;
    //    std::cerr << "\t\t--quiet             : Suppress the display of errors and summary statistics" << std::endl;
    std::cerr << std::endl;
}


int Validate::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String refFile = "";
    int maxErrors = -1;
    int printableErrors = 100;
    bool so_flag = false;
    bool so_coord = false;
    bool so_query = false;
    bool noeof = false;
    bool disableStatistics = false;
    bool verbose = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_PARAMETER("noeof", &noeof)
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_INTPARAMETER("maxErrors", &maxErrors)
        LONG_PARAMETER("verbose", &verbose)
        LONG_INTPARAMETER("printableErrors", &printableErrors)
        LONG_PARAMETER("disableStatistics", &disableStatistics)
        LONG_PARAMETER("params", &params)
        LONG_PARAMETER_GROUP("SortOrder")
        EXCLUSIVE_PARAMETER("so_flag", &so_flag)
        EXCLUSIVE_PARAMETER("so_coord", &so_coord)
        EXCLUSIVE_PARAMETER("so_query", &so_query)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);

    // Determine the sort type for validation based on the parameters.
    SamFile::SortedType sortType = SamFile::UNSORTED;
    if(so_flag)
    {
        sortType = SamFile::FLAG;
    } 
    else if(so_coord)
    {
        sortType = SamFile::COORDINATE;
    }
    else if(so_query)
    {
        sortType = SamFile::QUERY_NAME;
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
        std::cerr << "--in is a mandatory argument for validate, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    // Check to see if the ref file was specified.
    // Open the reference.
    GenomeSequence* refPtr = NULL;
    if(refFile != "")
    {
        refPtr = new GenomeSequence(refFile);
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Since we want to accumulate multiple errors, use RETURN rather
    // than throwing exceptions.
    SamFile samIn(ErrorHandler::RETURN);
    // Open the file for reading.   
    if(!samIn.OpenForRead(inFile))
    {
        fprintf(stderr, "Failed opening the SAM/BAM file, returning: %d (%s)\n",
                samIn.GetStatus(), 
                SamStatus::getStatusString(samIn.GetStatus()));
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    // Set the reference.
    samIn.SetReference(refPtr);

    // Set the sorting validation type.
    samIn.setSortedValidation(sortType);

    // Set that statistics should be generated.
    samIn.GenerateStatistics(!disableStatistics);

    // Read the sam header.
    SamFileHeader samHeader;
    if(!samIn.ReadHeader(samHeader))
    {
        fprintf(stderr, "Failed header validation, returning: %d (%s)\n", 
                samIn.GetStatus(), 
                SamStatus::getStatusString(samIn.GetStatus()));
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    // Read the sam records.
    SamRecord samRecord(ErrorHandler::RETURN);

    // Track the status.
    SamStatus::Status status = SamStatus::SUCCESS;

    // Keep reading records until the end of the file is reached.
    int numValidRecords = 0;
    int numInvalidRecords = 0;
    int numErrorRecords = 0;
    int numRecords = 0;
    int numReportedErrors = 0;
    int totalErrorRecords = 0;

    std::map<SamStatus::Status, uint64_t> errorStats;
    std::map<SamValidationError::Type, uint64_t> invalidStats;

    SamValidationErrors invalidSamErrors;

    // Keep reading records from the file until SamFile::ReadRecord
    // indicates to stop (returns false).
    while( ( (maxErrors < 0) || (totalErrorRecords < maxErrors) ) &&
           ( (samIn.ReadRecord(samHeader, samRecord)) || (SamStatus::isContinuableStatus(samIn.GetStatus())) ) )
    {
        ++numRecords;
        if(samIn.GetStatus() == SamStatus::SUCCESS)
        {
            // Successfully set the record, so check to see if it is valid.
            // Clear any errors in the list.
            invalidSamErrors.clear();
            if(!SamValidator::isValid(samHeader, samRecord, invalidSamErrors))
            {
                // The record is not valid.
                ++numInvalidRecords;
                ++totalErrorRecords;
                if(verbose && (numReportedErrors < printableErrors))
                {
                    std::cerr << "Record " << numRecords << std::endl
                              << invalidSamErrors << std::endl;
                    ++numReportedErrors;
                }
                // Update the statistics for all validation errors found in this record.
                invalidSamErrors.resetErrorIter();
                const SamValidationError* errorPtr = invalidSamErrors.getNextError();
                while(errorPtr != NULL)
                {
                    ++invalidStats[errorPtr->getType()];
                    errorPtr = invalidSamErrors.getNextError();
                }

                // If the status is not yet set, set it.
                if(status == SamStatus::SUCCESS)
                {
                    status = SamStatus::INVALID;
                }
            }
            else
            {
                // Valid record, so increment the counter.
                ++numValidRecords;
            }
        }
        else
        {
            // Error reading the record.
            ++numErrorRecords;
            ++totalErrorRecords;
            if(verbose && (numReportedErrors < printableErrors))
            {
                // report error.
                std::cerr << "Record " << numRecords << std::endl
                          << samIn.GetStatusMessage() << std::endl
                          << std::endl;
                ++numReportedErrors;
            }
            // Increment the statistics
            ++errorStats[samIn.GetStatus()];

            // If the status is not yet set, set it.
            if(status == SamStatus::SUCCESS)
            {
                status = samIn.GetStatus();
            }
        }
    }

    if( (samIn.GetStatus() != SamStatus::NO_MORE_RECS) &&
        ((totalErrorRecords < maxErrors) || (maxErrors < 0)))
    {
        // The last read call had a failure, so report it.
        // If the number of errors is >= ,maxErrors we don't
        // want to print any more failures.
        ++numErrorRecords;
        ++totalErrorRecords;
        if(numReportedErrors < printableErrors)
        {
            std::cerr << "Record " << numRecords << ": ";
            std::cerr << std::endl << samIn.GetStatusMessage() << std::endl;
        }

        // Increment the statistics
        ++errorStats[samIn.GetStatus()];

        if(status == SamStatus::SUCCESS)
        {
            status = samIn.GetStatus();
        }
    }

    if(totalErrorRecords == maxErrors)
    {
        if(maxErrors == 0)
        {
            std::cerr << "WARNING file was not read at all due to maxErrors setting, but returning Success.\n";
        }
        else
        {
            // Print a note that the entire file was not read.
            std::cerr << "File was not completely read due to the number of errors.\n";
            std::cerr << "Statistics only reflect the part of the file that was read.\n";
        }
    }

    fprintf(stderr, "\nNumber of records read = %d\n", numRecords);
    fprintf(stderr, "Number of valid records = %d\n", numValidRecords);

    std::cerr << std::endl;
    if(numRecords != numValidRecords)
    {
        std::cerr << "Error Counts:\n";

        // Loop through the non-validation errors.
        std::map<SamStatus::Status, uint64_t>::iterator statusIter;
        for(statusIter = errorStats.begin(); statusIter != errorStats.end(); statusIter++)
        {
            std::cerr << "\t" << SamStatus::getStatusString(statusIter->first) << ": "
                      << statusIter->second << std::endl;
        }

        std::map<SamValidationError::Type, uint64_t>::iterator invalidIter;
        for(invalidIter = invalidStats.begin(); invalidIter != invalidStats.end(); invalidIter++)
        {
            std::cerr << "\t" << SamValidationError::getTypeString(invalidIter->first) << ": "
                      << invalidIter->second << std::endl;
        }

        std::cerr << std::endl;
    }
    samIn.PrintStatistics();

    fprintf(stderr, "Returning: %d (%s)\n", status, SamStatus::getStatusString(status));
    return(status);
}


