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
#include "ExplainFlags.h"
#include "Filter.h"
#include "ReadReference.h"
#include "Revert.h"
#include "Diff.h"
#include "Squeeze.h"
#include "FindCigars.h"
#include "Stats.h"
#include "ClipOverlap.h"
#include "SplitBam.h"
#include "TrimBam.h"
#include "MergeBam.h"
#include "PolishBam.h"
#include "GapInfo.h"
#include "Dedup.h"
#include "Dedup_LowMem.h"
#include "Recab.h"
#include "Bam2FastQ.h"
#include "PhoneHome.h"

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
    TrimBam::trimBamDescription();
    MergeBam::mergeBamDescription();
    PolishBam::polishBamDescription();
    Dedup::dedupDescription();
    Dedup_LowMem::dedup_LowMemDescription();
    Recab::recabDescription();

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
    ExplainFlags::explainFlagsDescription();

    std::cerr << "\nAdditional Tools\n";
    Bam2FastQ::bam2FastQDescription();

    std::cerr << "\nDummy/Example Tools\n";
    ReadIndexedBam::readIndexedBamDescription();


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

    String cmd = argv[1];
    if(cmd.SlowCompare("readIndexedBam") == 0)
    {
        bamExe = new ReadIndexedBam();
    }
    else if(cmd.SlowCompare("dumpHeader") == 0)
    {
        bamExe = new DumpHeader();
    }
    else if(cmd.SlowCompare("dumpIndex") == 0)
    {
        bamExe = new DumpIndex();
    }
    else if(cmd.SlowCompare("writeRegion") == 0)
    {
        bamExe = new WriteRegion();
    }
    else if(cmd.SlowCompare("validate") == 0)
    {
        bamExe = new Validate();
    }
    else if(cmd.SlowCompare("splitChromosome") == 0)
    {
        bamExe = new SplitChromosome();
    }
    else if(cmd.SlowCompare("dumpRefInfo") == 0)
    {
        bamExe = new DumpRefInfo();
    }
    else if(cmd.SlowCompare("explainFlags") == 0)
    {
        bamExe = new ExplainFlags();
    }
    else if(cmd.SlowCompare("filter") == 0)
    {
        bamExe = new Filter();
    }
    else if(cmd.SlowCompare("readReference") == 0)
    {
        bamExe = new ReadReference();
    }
    else if(cmd.SlowCompare("revert") == 0)
    {
        bamExe = new Revert();
    }
    else if(cmd.SlowCompare("diff") == 0)
    {
        bamExe = new Diff();
    }
    else if(cmd.SlowCompare("squeeze") == 0)
    {
        bamExe = new Squeeze();
    }
    else if(cmd.SlowCompare("findCigars") == 0)
    {
        bamExe = new FindCigars();
    }
    else if(cmd.SlowCompare("stats") == 0)
    {
        bamExe = new Stats();
    }
    else if(cmd.SlowCompare("clipOverlap") == 0)
    {
        bamExe = new ClipOverlap();
    }
    else if(cmd.SlowCompare("splitBam") == 0)
    {
        bamExe = new SplitBam();
    }
    else if(cmd.SlowCompare("trimBam") == 0)
    {
        bamExe = new TrimBam();
    }
    else if((cmd.SlowCompare("mergeBam") == 0) ||
            (cmd.SlowCompare("rgMergeBam") == 0))
    {
        bamExe = new MergeBam();
    }
    else if(cmd.SlowCompare("polishBam") == 0)
    {
        bamExe = new PolishBam();
    }
    else if(cmd.SlowCompare("gapInfo") == 0)
    {
        bamExe = new GapInfo();
    }
    else if(cmd.SlowCompare("dedup") == 0)
    {
        bamExe = new Dedup();
    }
    else if(cmd.SlowCompare("dedup_LowMem") == 0)
    {
        bamExe = new Dedup_LowMem();
    }
    else if(cmd.SlowCompare("recab") == 0)
    {
        bamExe = new Recab();
    }
    else if(cmd.SlowCompare("bam2FastQ") == 0)
    {
        bamExe = new Bam2FastQ();
    }
    else if(cmd.SlowCompare("convert") == 0)
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
        int returnVal = 0;
        String compStatus;
        try
        {
            bamExe = new Convert();
            returnVal = bamExe->execute(numArgs, args);
        }
        catch (std::runtime_error e)
        {
            compStatus = "Exception";
            PhoneHome::completionStatus(compStatus.c_str());
            
            std::string errorMsg = "Exiting due to ERROR:\n\t";
            errorMsg += e.what();
            std::cerr << errorMsg << std::endl;
            return(-1);
        }
        compStatus = returnVal;
        PhoneHome::completionStatus(compStatus.c_str());
        delete bamExe;
        bamExe = NULL;
        return(returnVal);
    }
    
    if(bamExe != NULL)
    {
        int returnVal = 0;
        String compStatus;
        try
        {
            returnVal = bamExe->execute(argc, argv);
        }
        catch (std::runtime_error e)
        {
            compStatus = "Exception";
            PhoneHome::completionStatus(compStatus.c_str());
            std::string errorMsg = "Exiting due to ERROR:\n\t";
            errorMsg += e.what();
            std::cerr << errorMsg << std::endl;
            return(-1);
        }
        compStatus = returnVal;
        PhoneHome::completionStatus(compStatus.c_str());
        delete bamExe;
        bamExe = NULL;
        return(returnVal);
    }
    return(-1);
}



