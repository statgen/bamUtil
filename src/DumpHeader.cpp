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
// This file contains the processing for the executable option "dumpHeader"
// which prints the SAM/BAM header file to the screen.

#include "DumpHeader.h"
#include "SamFile.h"
#include "Parameters.h"
#include "PhoneHome.h"

void DumpHeader::dumpHeaderDescription()
{
    std::cerr << " dumpHeader - Print SAM/BAM Header" << std::endl;
}


void DumpHeader::description()
{
    dumpHeaderDescription();
}


void DumpHeader::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam dumpHeader <inputFile>" << std::endl;
}


// Dump the specified Bam Index file.
int DumpHeader::execute(int argc, char **argv)
{
    if(argc != 3)
    {
        String noPhArg = "--noPhoneHome";
        if((argc != 4) || (noPhArg.SlowCompareToStem(argv[3]) != 0))
        {
            usage();
            exit(-1);
        }
    }
    else
    {
        PhoneHome::checkVersion(getProgramName(), VERSION);
    }
    
    // Dump the bam index.
    return(dumpHeader(argv[2]));
}


int DumpHeader::dumpHeader(const char* filename)
{
    // Open the input file for reading.
    SamFile samIn;
    samIn.OpenForRead(filename);

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    std::string header;
    if(samHeader.getHeaderString(header))
    {
        std::cout << header;
    }
    return(SamStatus::SUCCESS);
}

