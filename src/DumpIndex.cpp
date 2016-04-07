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
// This file contains the processing for the executable option "dumpIndex"
// which prints a BAM Index File in a readable format.

#include "DumpIndex.h"
#include "SamFile.h"
#include "Parameters.h"
#include <iomanip>

void DumpIndex::printDumpIndexDescription(std::ostream& os)
{
    os << " dumpIndex - Print BAM Index File in English" << std::endl;
}


void DumpIndex::printDescription(std::ostream& os)
{
    printDumpIndexDescription(os);
}


void DumpIndex::printUsage(std::ostream& os)
{
    BamExecutable::printUsage(os);
    os << "\t./bam dumpIndex --bamIndex <bamIndexFile> [--refID <ref#>] [--summary] [--params]" << std::endl;
    os << "\tRequired Parameters:" << std::endl;
    os << "\t\t--bamIndex : the path/name of the bam index file to display" << std::endl;
    os << "\tOptional Parameters:" << std::endl;
    os << "\t\t--refID    : the reference ID to read, defaults to print all\n";
    os << "\t\t--summary  : only print a summary - 1 line per reference.\n";
    os << "\t\t--params   : print the parameter settings" << std::endl;
    os << std::endl;
}


// Dump the specified Bam Index file.
int DumpIndex::execute(int argc, char **argv)
{
    // Extract command line arguments.
    static const int UNSPECIFIED_INT = -1;
    String indexFile = "";
    int refID = UNSPECIFIED_INT;
    bool summary = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
        LONG_INTPARAMETER("refID", &refID)
        LONG_PARAMETER("summary", &summary)
        LONG_PARAMETER("params", &params)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);

    // Check to see if the index file was specified, if not, report an error.
    if(indexFile == "")
    {
        printUsage(std::cerr);
        inputParameters.Status();
        // mandatory argument was not specified.
        std::cerr << "Missing mandatory argument: --bamIndex" << std::endl;
        return(-1);
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Read the index.
    BamIndex bamIndex;
    SamStatus status;
    status = bamIndex.readIndex(indexFile);

    if(status != SamStatus::SUCCESS)
    {
        // Failed to read the index, return.
        fprintf(stderr, "%s\n", status.getStatusMessage());
        return(status.getStatus());
    }

    // Print the index file.
    bamIndex.printIndex(refID, summary);

    return(status.getStatus());
}

