/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan
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

#include <iostream>
#include <string.h>
#include <stdlib.h>

#include "Validate.h"
#include "Convert.h"
#include "DumpHeader.h"
#include "SplitChromosome.h"
#include "WriteRegion.h"
#include "DumpIndex.h"
#include "ReadIndexedBam.h"
#include "DumpRefInfo.h"
#include "Filter.h"
#include "ReadReference.h"
#include "Revert.h"
#include "Diff.h"
#include "Squeeze.h"
#include "FindCigars.h"
#include "Stats.h"
#include "ClipOverlap.h"
#include "SplitBam.h"
#include "GapInfo.h"

void Usage()
{
    BamExecutable::bamExecutableDescription();
    std::cerr << std::endl;
    std::cerr << "Tools to Rewrite SAM/BAM Files: " << std::endl;
    Convert::convertDescription();
    WriteRegion::writeRegionDescription();
    SplitChromosome::splitChromosomeDescription();
    SplitBam::splitBamDescription();
    FindCigars::findCigarsDescription();

    std::cerr << "\nTools to Modify & write SAM/BAM Files: " << std::endl;
    ClipOverlap::clipOverlapDescription();
    Filter::filterDescription();
    Revert::revertDescription();
    Squeeze::squeezeDescription();

    std::cerr << "\nInformational Tools\n";
    Validate::validateDescription();
    Diff::diffDescription();
    Stats::statsDescription();
    GapInfo::gapInfoDescription();

    std::cerr << "\nTools to Print Information In Readable Format\n";
    DumpHeader::dumpHeaderDescription();
    DumpRefInfo::dumpRefInfoDescription();
    DumpIndex::dumpIndexDescription();
    ReadReference::readReferenceDescription();

    std::cerr << "\nDummy/Example Tools\n";
    ReadReference::readReferenceDescription();


    std::cerr << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << "\tbam <tool> [<tool arguments>]" << std::endl;
    std::cerr << "The usage for each tool is described by specifying the tool with no arguments." << std::endl;
}


int main(int argc, char ** argv)
{
    BamExecutable* bamExe = NULL;

    // Verify at least one arg.
    if(argc < 2)
    {
        // Not enough args...
        Usage();
        exit(-1);
    }

    if(strcmp(argv[1], "readIndexedBam") == 0)
    {
        bamExe = new ReadIndexedBam();
    }
    else if(strcmp(argv[1], "dumpHeader") == 0)
    {
        bamExe = new DumpHeader();
    }
    else if(strcmp(argv[1], "dumpIndex") == 0)
    {
        bamExe = new DumpIndex();
    }
    else if(strcmp(argv[1], "writeRegion") == 0)
    {
        bamExe = new WriteRegion();
    }
    else if(strcmp(argv[1], "validate") == 0)
    {
        bamExe = new Validate();
    }
    else if(strcmp(argv[1], "splitChromosome") == 0)
    {
        bamExe = new SplitChromosome();
    }
    else if(strcmp(argv[1], "dumpRefInfo") == 0)
    {
        bamExe = new DumpRefInfo();
    }
    else if(strcmp(argv[1], "filter") == 0)
    {
        bamExe = new Filter();
    }
    else if(strcmp(argv[1], "readReference") == 0)
    {
        bamExe = new ReadReference();
    }
    else if(strcmp(argv[1], "revert") == 0)
    {
        bamExe = new Revert();
    }
    else if(strcmp(argv[1], "diff") == 0)
    {
        bamExe = new Diff();
    }
    else if(strcmp(argv[1], "squeeze") == 0)
    {
        bamExe = new Squeeze();
    }
    else if(strcmp(argv[1], "findCigars") == 0)
    {
        bamExe = new FindCigars();
    }
    else if(strcmp(argv[1], "stats") == 0)
    {
        bamExe = new Stats();
    }
    else if(strcmp(argv[1], "clipOverlap") == 0)
    {
        bamExe = new ClipOverlap();
    }
    else if(strcmp(argv[1], "splitBam") == 0)
    {
        bamExe = new SplitBam();
    }
    else if(strcmp(argv[1], "gapInfo") == 0)
    {
        bamExe = new GapInfo();
    }
    else if(strcmp(argv[1], "convert") == 0)
    {
        bamExe = new Convert();
    }
    else
    {
        // This is the backward compatable version of convert.
        bool noeof = false;
        if(argc != 3)
        {
            if(argc == 4)
            {
                if(strcmp(argv[3], "NOEOF") != 0)
                {
                    Usage();
                    exit(-1);
                }
                else
                {
                    noeof = true;
                }
            }
            else
            {
                Usage();
                exit(-1);
            }
        }
        int numArgs = 6;
        // char ** args;
        char* args[7];
        char arg1[] = "convert";
        char arg2[] = "--in";
        char arg4[] = "--out";
        char arg6[] = "--noeof";
        args[0] = argv[0];
        args[1] = arg1;
        args[2] = arg2;
        args[3] = argv[1];
        args[4] = arg4;
        args[5] = argv[2];
        if(noeof)
        {
            args[6] = arg6;
            ++numArgs;
        }
        bamExe = new Convert();
        int returnVal = bamExe->execute(numArgs, args);
        delete bamExe;
        bamExe = NULL;
        return(returnVal);
    }
    
    if(bamExe != NULL)
    {
        int returnVal = bamExe->execute(argc, argv);
        delete bamExe;
        bamExe = NULL;
        return(returnVal);
    }
    return(-1);
}



