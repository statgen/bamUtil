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
// This file contains the processing for the executable option "dumpAsp"
// which prints a Asp File in a readable format.

#include "DumpAsp.h"
#include <iomanip>
#include "SamFile.h"
#include "Parameters.h"
#include "AspFile.h"
#include "BaseAsciiMap.h"

void DumpAsp::dumpAspDescription()
{
    std::cerr << " dumpAsp - Print Asp File in English" << std::endl;
}


void DumpAsp::description()
{
    dumpAspDescription();
}


void DumpAsp::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam dumpAsp --asp <aspFile> [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--asp : the path/name of the asp file to display" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--statsOnly : Only print the stats on the records."
              << std::endl;
    std::cerr << "\t\t--dataOnly  : Only print the data records." << std::endl;
    std::cerr << "\t\t--params    : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


// Dump the specified Bam Asp file.
int DumpAsp::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String aspFile = "";
    bool statsOnly = false;
    bool dataOnly = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("asp", &aspFile)
        LONG_PARAMETER("statsOnly", &statsOnly)
        LONG_PARAMETER("dataOnly", &dataOnly)
        LONG_PARAMETER("params", &params)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    inputParameters.Read(argc-1, &(argv[1]));

    // Check to see if the asp file was specified, if not, report an error.
    if(aspFile == "")
    {
        usage();
        inputParameters.Status();
        // mandatory argument was not specified.
        std::cerr << "Missing mandatory argument: --asp" << std::endl;
        return(-1);
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Read the asp.
    AspFileReader asp;
    // Throws acception if the open fails.
    asp.open(aspFile);

    AspRecord record;
    int returnVal = 0;
    if(dataOnly)
    {
        while(asp.getNextDataRecord(record))
        {
            // print the chrom/pos of this record.
            std::cout << record.getChromID() << ":" << record.get0BasedPos() << "\t";
            if(record.isRefOnlyType())
            {
                printRefOnly(record);
            }
            else if(record.isDetailedType())
            {
                printDetailed(record);
            }
            else
            {
                std::cout << "UNKNOWN!!!\n";
                returnVal = 1;
            }            
        }


    }
    else if(statsOnly)
    {
        unsigned int numEmpty = 0;
        unsigned int numPos = 0;
        unsigned int numRefOnly = 0;
        unsigned int numDetailed = 0;
        unsigned int numUnknown = 0;
        while(asp.getNextRecord(record))
        {
            if(record.isEmptyType())
            {
                ++numEmpty;
                if(numEmpty == 0)
                {
                    std::cerr << "Empty Wrap Around\n";
                }
            }
            else if(record.isPosType())
            {
                ++numPos;
                if(numPos == 0)
                {
                    std::cerr << "Pos Wrap Around\n";
                }
            }
            else if(record.isRefOnlyType())
            {
                ++numRefOnly;
                if(numRefOnly == 0)
                {
                    std::cerr << "RefOnly Wrap Around\n";
                }
            }
            else if(record.isDetailedType())
            {
                ++numDetailed;
                if(numDetailed == 0)
                {
                    std::cerr << "Detailed Wrap Around\n";
                }
            }
            else
            {
                ++numUnknown;
                returnVal = 1;
            }
        }
        // Output the stats.
        std::cerr << "Number of Position Records = "
                  << numPos << "\n"
                  << "Number of Empty Records = "
                  << numEmpty << "\n"
                  << "Number of Reference Only Records = "
                  << numRefOnly << "\n"
                  << "Number of Detailed Records = "
                  << numDetailed << "\n"
                  << "Number of Unknown Records = "
                  << numUnknown << "\n";
    }
    else
    {
        while(asp.getNextRecord(record))
        {
            // print the chrom/pos of this record.
            std::cout << record.getChromID() << ":" << record.get0BasedPos() << "\t";
            if(record.isEmptyType())
            {
                std::cout << "EMPTY\n";
            }
            else if(record.isPosType())
            {
                std::cout << "POS\n";
            }
            else if(record.isRefOnlyType())
            {
                printRefOnly(record);
            }
            else if(record.isDetailedType())
            {
                printDetailed(record);
            }
            else
            {
                std::cout << "UNKNOWN!!!\n";
                returnVal = 1;
            }
        }
    }

    return(returnVal);
}


void DumpAsp::printRefOnly(AspRecord& record)
{
    std::cout << "REF_ONLY\t" 
              << record.getRefBase() << "\t"
              << record.getNumBases() << "\t"
              << record.getGLH() << "\t" 
              << record.getGLA() << "\n";
}


void DumpAsp::printDetailed(AspRecord& record)
{
    std::cout << "DETAILED\t" 
              << record.getRefBase() << "\t"
              << record.getNumBases() << "\t";
    for(int i = 0; i < record.getNumBases(); i++)
    {
        std::cout << record.getBaseChar(i);
    }
    std::cout << "\t";
    for(int i = 0; i < record.getNumBases(); i++)
    {
        std::cout << record.getCharQual(i);
    }
    std::cout << "\t";
    for(int i = 0; i < record.getNumBases(); i++)
    {
        if(i != 0)
        {
            std::cout << ":";
        }
        std::cout << record.getCycle(i);
    }
    std::cout << "\t";
    for(int i = 0; i < record.getNumBases(); i++)
    {
        std::cout << (bool)record.getStrand(i);
    }
    std::cout << "\t";
    for(int i = 0; i < record.getNumBases(); i++)
    {
        if(i != 0)
        {
            std::cout << ":";
        }
        std::cout << record.getMQ(i);
    }

    std::cout << "\n";
}
