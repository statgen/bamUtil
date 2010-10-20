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

void Usage()
{
    std::cerr << "Usage: " << std::endl;
    std::cerr << std::endl;
    validateDescription();
    convertDescription();
    dumpHeaderDescription();
    splitChromosomeDescription();
    writeRegionDescription();
    dumpRefInfoDescription();
    dumpIndexDescription();
    readIndexedBamDescription();
    filterDescription();
    readReferenceDescription();
}


int main(int argc, char ** argv)
{
    // Verify at least one arg.
    if(argc < 2)
    {
        // Not enough args...
        Usage();
        exit(-1);
    }

    if(strcmp(argv[1], "readIndexedBam") == 0)
    {
        if(argc != 5)
        {
            readIndexedBamUsage();
            exit(-1);
        }
        return(readIndexedBam(argv[2], argv[3], argv[4]));
    }
   
    if(strcmp(argv[1], "dumpHeader") == 0)
    {
        if(argc != 3)
        {
            dumpHeaderUsage();
            exit(-1);
        }
        // Dump the bam index.
        return(dumpHeader(argv[2]));
    }
   
    if(strcmp(argv[1], "dumpIndex") == 0)
    {
        // Dump the bam index.
        return(dumpIndex(argc, argv));
    }
   
    if(strcmp(argv[1], "writeRegion") == 0)
    {
        return(writeRegion(argc, argv));
    }

    if(strcmp(argv[1], "validate") == 0)
    {
        return(validate(argc, argv));
    }

    if(strcmp(argv[1], "splitChromosome") == 0)
    {
        return(splitChromosome(argc, argv));
    }

    if(strcmp(argv[1], "dumpRefInfo") == 0)
    {
        return(dumpRefInfo(argc, argv));
    }

    if(strcmp(argv[1], "filter") == 0)
    {
        return(filter(argc, argv));
    }

    if(strcmp(argv[1], "readReference") == 0)
    {
        return(readReference(argc, argv));
    }

    if(strcmp(argv[1], "readReference") == 0)
    {
        return(readReference(argc, argv));
    }

    // Check usage for the sam/bam converter.
    if(strcmp(argv[1], "convert") == 0)
    {
        return(convert(argc, argv));
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
        return(convert(numArgs, args));
    }
}



