/*
 *  Copyright (C) 2013  Regents of the University of Michigan
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
// This file contains the processing for the executable option "explainFlags"
// which writes a file with the reads in the specified region.

#include "ExplainFlags.h"
#include "SamFlag.h"
#include "Parameters.h"

ExplainFlags::ExplainFlags()
    : BamExecutable()
{
    
}

void ExplainFlags::explainFlagsDescription()
{
    std::cerr << " explainFlags - Describe flags" << std::endl;
}

void ExplainFlags::description()
{
    explainFlagsDescription();
}

void ExplainFlags::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam explainFlags --hex|dec <flag> [--params]" << std::endl;
    std::cerr << "\tFlagValue Parameters:" << std::endl;
    std::cerr << "\t\t--hex     : Explain this hex flag" << std::endl;
    std::cerr << "\t\t--dec     : Explain this decimal flag" << std::endl;
    std::cerr << "\t\t--params  : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int ExplainFlags::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String hexFlag = "";
    int decFlag = -1;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("FlagValue Parameters")
        LONG_STRINGPARAMETER("hex", &hexFlag)
        LONG_INTPARAMETER("dec", &decFlag)
        LONG_PARAMETER("params", &params)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);

    if(params)
    {
        inputParameters.Status();
    }

    if((hexFlag.Length() != 0) && (decFlag != -1))
    {
        usage();
        inputParameters.Status();
        std::cerr << "Can't specify both --hex & --dec" << std::endl;
        return(-1);
    }
    if((hexFlag.Length() == 0) && (decFlag == -1))
    {
        usage();
        inputParameters.Status();
        std::cerr << "Either --hex or --dec is required" << std::endl;
        return(-1);
    }
    if(hexFlag.Length() != 0)
    {
        if((hexFlag.Length() >= 2) && hexFlag[0] == '0' && 
           ((hexFlag[1] == 'x') || (hexFlag[2] == 'X')))
        {
            // Already has the 0x at the beginning.
            decFlag = hexFlag.AsInteger();
        }
        else
        {
            // Add the 0x before converting to int.
            String newString = "0x";
            newString += hexFlag;
            decFlag = newString.AsInteger();
        }
    }

    std::cout << "0x" << std::hex << decFlag 
              << " (" << std::dec << decFlag << "):\n";

    if(SamFlag::isPaired(decFlag)) { std::cout << "\tpaired\n"; }
    if(SamFlag::isProperPair(decFlag)) { std::cout << "\tproperly paired\n"; }
    if(!SamFlag::isMapped(decFlag)) { std::cout << "\tunmapped\n"; }
    if(!SamFlag::isMateMapped(decFlag)) { std::cout << "\tmate unmapped\n"; }
    if(SamFlag::isReverse(decFlag)) { std::cout << "\treverse strand\n"; }
    if(SamFlag::isMateReverse(decFlag)) { std::cout << "\tmate reverse strand\n"; }
    if(SamFlag::isMidFragment(decFlag)) { std::cout << "\tmiddle fragment\n"; }
    else if(SamFlag::isFirstFragment(decFlag)) { std::cout << "\tfirst fragment\n"; }
    else if(SamFlag::isLastFragment(decFlag)) { std::cout << "\tlast fragment\n"; }
    else { std::cout << "\tunknown fragment\n"; }
    if(SamFlag::isSecondary(decFlag)) { std::cout << "\tsecondary alignment\n"; }
    if(SamFlag::isQCFailure(decFlag)) { std::cout << "\tfails QC checks\n"; }
    if(SamFlag::isDuplicate(decFlag)) { std::cout << "\tduplicate\n"; }

    return(0);
}
