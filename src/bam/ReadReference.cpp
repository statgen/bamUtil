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
#include "ReadReference.h"

#include "Parameters.h"
#include "BgzfFileType.h"
#include "GenomeSequence.h"

void readReferenceDescription()
{
    std::cerr << " readReference - Print the reference string for the specified region" << std::endl;
    std::cerr << "\t./bam readReference --refFile <referenceFilename> --refName <reference Name> --start <0 based start> --end <0 based end>|--numBases <number of bases> [--params]"<< std::endl;
}

void readReferenceUsage()
{
    readReferenceDescription();
    std::cerr << "\tRequired Parameters:\n"
              << "\t\t--refFile  : the reference\n"
              << "\t\t--refName  : the SAM/BAM reference Name to read\n"
              << "\t\t--start    : inclusive 0-based start position (defaults to -1)\n"
              << "\tRequired Length Parameter (one but not both needs to be specified):\n"
              << "\t\t--end      : exclusive 0-based end position (defaults to -1: meaning til the end of the reference)\n"
              << "\t\t--numBases : number of bases from start to display\n"
              << "\t\t--params   : print the parameter settings\n"
              << std::endl;
}



int readReference(int argc, char **argv)
{
    static const int UNSPECIFIED_INT = -1;
    String refFile = "";
    String refName = "";
    int start = UNSPECIFIED_INT;
    int numBases = UNSPECIFIED_INT;
    int end = UNSPECIFIED_INT;
    bool params = false;
    
    // Read in the parameters.    
    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_STRINGPARAMETER("refName", &refName)
        LONG_INTPARAMETER("start", &start)
        LONG_INTPARAMETER("end", &end)
        LONG_INTPARAMETER("numBases", &numBases)
        LONG_PARAMETER("params", &params)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));
    
    inputParameters.Read(argc-1, &(argv[1]));
    
    if((refName == "") || (start == UNSPECIFIED_INT) || 
       ((end == UNSPECIFIED_INT) && (numBases == UNSPECIFIED_INT)))
    {
        readReferenceUsage();
        inputParameters.Status();
        std::cerr << "Missing Required Parameter\n\n";
        return(-1);
    }
    if((end != UNSPECIFIED_INT) && (numBases != UNSPECIFIED_INT))
    {
        readReferenceUsage();
        inputParameters.Status();
        std::cerr << "Only --end or --numBases can be specified\n\n";
        return(-1);
    }
    else if(numBases != UNSPECIFIED_INT)
    {
        end = start + numBases;
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Open the reference.
    GenomeSequence reference(refFile);

    uint32_t refStart = 
        reference.getGenomePosition(refName.c_str());

    std::string refString;
    
    reference.getString(refString, refStart + start, end - start);
    std::cout << refString << std::endl;
    
    return(0);
}
