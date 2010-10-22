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
// This file contains the processing for the executable option "dumpGlf"
// which prints the GLF file to the screen.

#include "GlfFile.h"
#include "GlfRefSection.h"

void dumpGlfDescription()
{
    std::cerr << " dumpGlf - Print GLF File in a readable format:" << std::endl;
    std::cerr << "\t./bam dumpGlf <inputFile>" << std::endl;
}


void dumpGlfUsage()
{
    dumpGlfDescription();
    std::cerr << std::endl;
}


// Dump the specified Bam Index file.
int dumpGlf(const char* filename)
{
    // Open the input file for reading.
    GlfFile glfIn;
    GlfHeader header;
    if(!glfIn.openForRead(filename, header))
    {
        fprintf(stderr, "%s\n", glfIn.getStatusMessage());
        return(glfIn.getStatus());
    }

    std::string text;
    if(header.getHeaderTextString(text))
    {
        std::cerr << "GLF Header Text: " << text << std::endl;
    }
    else
    {
        std::cerr << "Failed to read the GLF Header.\n";
    }

    GlfRefSection refSection;
    GlfRecord glfRecord;
    uint32_t position = 0;
    while(glfIn.getNextRefSection(refSection))
    {
        refSection.print();
        
        // Get the records.
        while(glfIn.getNextRecord(glfRecord))
        {
            position += glfRecord.getOffset();
            std::cerr << "Position = " << position << "\t";
            glfRecord.print();
        }
    }

    return(GlfStatus::SUCCESS);
}
