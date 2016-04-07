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
#include <algorithm>
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

// May add option to print to console in red for errors.
namespace console_color
{
    static const std::string reset = "\033[0m";
    static const std::string red = "\033[1;31m";
}

std::string ToLowerCase(const char* cstr)
{
    std::string str(cstr);
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

void WriteUsage(std::ostream& os)
{
    BamExecutable::printBamExecutableDescription(os);

    os << std::endl;
    os << "Usage: " << std::endl;
    os << "bam <tool> [<tool arguments>]" << std::endl;
    os << "bam (usage|help) - Will print usage and exit." << std::endl;
    os << "bam (usage|help) <tool> - Will print usage for specific tool and exit." << std::endl;

    os << std::endl;
    os << "Tools to Rewrite SAM/BAM Files: " << std::endl;
    Convert::printConvertDescription(os);
    WriteRegion::printWriteRegionDescription(os);
    SplitChromosome::printSplitChromosomeDescription(os);
    SplitBam::printSplitBamDescription(os);
    FindCigars::printFindCigarsDescription(os);

    os << "\nTools to Modify & write SAM/BAM Files: " << std::endl;
    ClipOverlap::printClipOverlapDescription(os);
    Filter::printFilterDescription(os);
    Revert::printRevertDescription(os);
    Squeeze::printSqueezeDescription(os);
    TrimBam::printTrimBamDescription(os);
    MergeBam::printMergeBamDescription(os);
    PolishBam::printPolishBamDescription(os);
    Dedup::printDedupDescription(os);
    Dedup_LowMem::printDedup_LowMemDescription(os);
    Recab::printRecabDescription(os);

    os << "\nInformational Tools\n";
    Validate::printValidateDescription(os);
    Diff::printDiffDescription(os);
    Stats::printStatsDescription(os);
    GapInfo::printGapInfoDescription(os);

    os << "\nTools to Print Information In Readable Format\n";
    DumpHeader::printDumpHeaderDescription(os);
    DumpRefInfo::printDumpRefInfoDescription(os);
    DumpIndex::printDumpIndexDescription(os);
    ReadReference::printReadReferenceDescription(os);
    ExplainFlags::printExplainFlagsDescription(os);

    os << "\nAdditional Tools\n";
    Bam2FastQ::printBam2FastQDescription(os);

    os << "\nDummy/Example Tools\n";
    ReadIndexedBam::printReadIndexedBamDescription(os);
}

BamExecutable* CreateBamExe(const std::string& name)
{
    BamExecutable* ret = NULL;

    if(name == ToLowerCase("readIndexedBam"))
    {
        ret = new ReadIndexedBam();
    }
    else if(name == ToLowerCase("dumpHeader"))
    {
        ret = new DumpHeader();
    }
    else if(name == ToLowerCase("dumpIndex"))
    {
        ret = new DumpIndex();
    }
    else if(name == ToLowerCase("writeRegion"))
    {
        ret = new WriteRegion();
    }
    else if(name == "validate")
    {
        ret = new Validate();
    }
    else if(name == ToLowerCase("splitChromosome"))
    {
        ret = new SplitChromosome();
    }
    else if(name == ToLowerCase("dumpRefInfo"))
    {
        ret = new DumpRefInfo();
    }
    else if(name == ToLowerCase("explainFlags"))
    {
        ret = new ExplainFlags();
    }
    else if(name == "filter")
    {
        ret = new Filter();
    }
    else if(name == ToLowerCase("readReference"))
    {
        ret = new ReadReference();
    }
    else if(name == "revert")
    {
        ret = new Revert();
    }
    else if(name == "diff")
    {
        ret = new Diff();
    }
    else if(name == "squeeze")
    {
        ret = new Squeeze();
    }
    else if(name == ToLowerCase("findCigars"))
    {
        ret = new FindCigars();
    }
    else if(name == "stats")
    {
        ret = new Stats();
    }
    else if(name == ToLowerCase("clipOverlap"))
    {
        ret = new ClipOverlap();
    }
    else if(name == ToLowerCase("splitBam"))
    {
        ret = new SplitBam();
    }
    else if(name == ToLowerCase("trimBam"))
    {
        ret = new TrimBam();
    }
    else if((name == ToLowerCase("mergeBam")) ||
            (name == ToLowerCase("rgMergeBam")))
    {
        ret = new MergeBam();
    }
    else if(name == ToLowerCase("polishBam"))
    {
        ret = new PolishBam();
    }
    else if(name == ToLowerCase("gapInfo"))
    {
        ret = new GapInfo();
    }
    else if(name == "dedup")
    {
        ret = new Dedup();
    }
    else if(name == ToLowerCase("dedup_LowMem"))
    {
        ret = new Dedup_LowMem();
    }
    else if(name == "recab")
    {
        ret = new Recab();
    }
    else if(name == ToLowerCase("bam2FastQ"))
    {
        ret = new Bam2FastQ();
    }
    else if(name == "convert")
    {
        ret = new Convert();
    }

    return ret;
}


int main(int argc, char ** argv)
{
    int ret = 0;


    // Verify at least one arg.
    if (argc<2)
    {
        // Not enough args...
        std::cerr << "Error: Not enough args.\n\n";
        WriteUsage(std::cerr);
        ret = -1;
    }
    else
    {
        std::string cmd = argv[1];
        std::transform(cmd.begin(), cmd.end(), cmd.begin(), ::tolower);

        if(cmd == "usage" || cmd == "help")
        {
            if (argc == 2)
            {
                WriteUsage(std::cout);
            }
            else
            {
                BamExecutable* bamExe = CreateBamExe(argv[2]);

                if (bamExe != NULL)
                {
                    bamExe->printUsage(std::cout);
                }
                else
                {
                    std::cerr << "Error: Invalid tool.\n\n";
                    WriteUsage(std::cerr);
                    ret = -1;
                }
            }
        }
        else
        {
            BamExecutable* bamExe = CreateBamExe(cmd);

            if(bamExe != NULL)
            {
                String compStatus;
                try
                {
                    ret = bamExe->execute(argc, argv);
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
                compStatus = ret;
                PhoneHome::completionStatus(compStatus.c_str());
                delete bamExe;
                bamExe = NULL;
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
                            WriteUsage(std::cerr);
                            exit(-1);
                        }
                        else
                        {
                            noeof = true;
                        }
                    }
                    else
                    {
                        std::cerr << "Error: Invalid command.\n\n";
                        WriteUsage(std::cerr);
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

                String compStatus;
                try
                {
                    bamExe = new Convert();
                    ret = bamExe->execute(numArgs, args);
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
                compStatus = ret;
                PhoneHome::completionStatus(compStatus.c_str());
                delete bamExe;
                bamExe = NULL;
            }
        }
    }

    return ret;
}



