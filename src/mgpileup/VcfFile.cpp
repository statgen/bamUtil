/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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

#include <stdlib.h>
#include "VcfFile.h"

VcfFile::VcfFile()
    : myFilePtr(NULL)
{
}

VcfFile::~VcfFile()
{
    close();
}


bool VcfFile::openForRead(const char* filename)
{
    // Close if a file is already open.
    close();

    myFilePtr = ifopen(filename, "r");

    if(myFilePtr == NULL) {
        std::cerr << "FATAL error : Cannot open file " << filename << std::endl;
        return(false);
    }
    return(true);
}



bool VcfFile::openForWrite(const char* filename)
{
    // Close if a file is already open.
    close();

    myFilePtr = ifopen(filename, "w");

    if(myFilePtr == NULL) {
        std::cerr << "FATAL error : Cannot open file " << filename << std::endl;
        return(false);
    }

    // Write the header.
    std::string tempStr("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
    tempStr.append(filename);
    tempStr.append("\n");
    return(ifwrite(myFilePtr, tempStr.c_str(), 
                   tempStr.length()) == tempStr.length());
}


void VcfFile::close()
{
    if (myFilePtr != NULL)
    {
        // If we already have an open file, close it.
        ifclose(myFilePtr);
        myFilePtr = NULL;
    }
}

bool VcfFile::isEof()
{
    return(ifeof(myFilePtr));
}


bool VcfFile::getNextDataLine(VcfFile::VcfDataLine& nextLine)
{
    char c = ' ';
    std::string field;
    bool dataLineFound = false;

    field.clear();

    // Loop until the end of the file or
    // until a chromosome/position is found.
    while((!ifeof(myFilePtr)) && (!dataLineFound))
    {
        if((c=ifgetc(myFilePtr))!='#')
        {
            nextLine.chromosome = c;
            
            //read chromosome
            while((c=ifgetc(myFilePtr))!='\t')
            {
                nextLine.chromosome.append(&c, 1);
            }

            // Add the chromosome
            nextLine.wholeLine = nextLine.chromosome;
            // add the tab.
            nextLine.wholeLine += '\t';

            //read position
            while((c=ifgetc(myFilePtr))!='\t')
            {
                field.append(&c, 1);
            }

            // Add the position.
            nextLine.wholeLine += field;

            // add the tab.
            nextLine.wholeLine += '\t';

            //std::cout << " pos:" << field << "\n";
            nextLine.position = atoi(field.c_str());
            dataLineFound = true;
        }

        //read rest of line
        while((c=ifgetc(myFilePtr))!='\n')
        {
            nextLine.wholeLine.append(&c, 1);
        }
        // Append the new line.
        nextLine.wholeLine.append(&c, 1 );
    }

    if(!dataLineFound)
    {
        // Data line not found, so reset the vcfDataLine.
        nextLine.wholeLine.clear();
        nextLine.chromosome.clear();
        nextLine.position = 0;
    }

    return(dataLineFound);
}


bool VcfFile::writeLine(VcfFile::VcfDataLine* line)
{
    if(line != NULL)
    {
        return(writeLine(line->wholeLine));
    }
    // No line, so nothing to write - success.
    return(true);
}


bool VcfFile::writeLine(std::string& line)
{
    unsigned int lineLen = line.length();
    return(ifwrite(myFilePtr, line.c_str(), lineLen) == lineLen);
}
