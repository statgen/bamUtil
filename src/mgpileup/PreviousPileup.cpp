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

#include <stdexcept>
#include "PreviousPileup.h"

PreviousPileup::PreviousPileup()
    : myFile(),
      myCurrentChromosome(""),
      myCurrentChromID(-1),
      myEof(false),
      myMaxStoredLines(1)
{
    myCurrentDataLine = new VcfFile::VcfDataLine();
    if(myCurrentDataLine == NULL)
    {
        throw(std::runtime_error("Failed to allocate a VcfDataLine"));
    }
    myCurrentDataLine->wholeLine.clear();
    myCurrentDataLine->chromosome.clear();
    myCurrentDataLine->position = 0;
}


PreviousPileup::~PreviousPileup()
{
    myFile.close();
    myCurrentDataLine = NULL;
}


void PreviousPileup::open(const std::string& fileName, uint32_t maxStoredLines)
{
    // Previously piledup VCF file
    if(!fileName.empty())
    {
        myFile.openForRead(fileName.c_str());
    }
    myMaxStoredLines = maxStoredLines;
}


bool PreviousPileup::findTarget(int targetChromID, int targetPosition, 
                                SamFileHeader& samHeader)
{
    // Check if current data line is not set.
    if(myCurrentDataLine == NULL)
    {
        throw(std::runtime_error("Code bug, PreviousPileup::findTarget, myCurrentDataLine should have been set."));
    }

    // Loop as long as there is still more to read or until the target is found or
    // past.
    while(!myEof)
    {
        if(myCurrentChromID != -1)
        {
            // The current line is set, so start there.
            if(myCurrentChromID == targetChromID)
            {
                if(myCurrentDataLine->position == targetPosition)
                {
                    // Found the target position.
                    return(true);
                }
                else if(myCurrentDataLine->position > targetPosition)
                {
                    // Past the target position, so return false, the target was not found.
                    return(false);
                }
            }
            else if(myCurrentChromID > targetChromID)
            {
                // Past the target chromsomoe ID, so return false, the target was not found.
                return(false);
            }
        }

        // The current line is still prior to the target or is not set, so get the next line.
        if(!myFile.getNextDataLine(*myCurrentDataLine))
        {
            // no more data lines.
            myEof = true;
            myCurrentChromID = -1;
        }
        // Check if the current chromosome id changed.
        else if((myCurrentChromID == -1) || (myCurrentChromosome != myCurrentDataLine->chromosome))
        {
            // Update the current chromosome info.
            myCurrentChromosome = myCurrentDataLine->chromosome;
            myCurrentChromID = samHeader.getReferenceID(myCurrentChromosome.c_str());
        }
    }

    // Got here, so the target was not found before the end of the file.
    return(false);
}


bool PreviousPileup::writeCurrentLine(VcfFile& outputFile)
{
    return(outputFile.writeLine(myCurrentDataLine));
}


bool PreviousPileup::writeStored(VcfFile& outputFile)
{
    bool status = true;

    while(status && (!myStoredLines.empty()))
    {
        // Keep writing lines until there are no more stored lines
        // or a write fails.
        status = outputFile.writeLine(myStoredLines.front());
        myUnusedLines.push(myStoredLines.front());
        myStoredLines.pop();
    }
    
    return(status);
}

// Write the stored values up until the specified position.
bool PreviousPileup::writeStoredUpToPos(int position, VcfFile& outputFile)
{
    bool status = true;

    while(status && (!myStoredLines.empty()))
    {

        VcfFile::VcfDataLine* firstStored = myStoredLines.front();
        if(firstStored == NULL)
        {
            throw(std::runtime_error("Invalid line stored in PreviousPileup."));
            return(false);
        }

        if(firstStored->position < position)
        {
            // Keep writing lines until there are no more stored lines
            // or a write fails.
            status = outputFile.writeLine(firstStored);
            myUnusedLines.push(firstStored);
            myStoredLines.pop();
        }
        else
        {
            // Gotten to or past the position where no more lines should be written
            // so break out of the write loop.
            break;
        }
    }
    
    return(status);
}


bool PreviousPileup::storeCurrent()
{
    if(myStoredLines.size() >= myMaxStoredLines)
    {
        // At the max number of stored lines, so do not store the
        // current line and return false.
        return(false);
    }

    // Store the current line.
    myStoredLines.push(myCurrentDataLine);

    // Get a new current line either from unused lines or allocate a new one.
    if(!myUnusedLines.empty())
    {
        // Still already allocated, unused lines, so use one of those.
        myCurrentDataLine = myUnusedLines.top();
        myUnusedLines.pop();
    }
    else
    {
        // No more unused lines, so allocate a new one.
        unsigned int currentAllocatedSize = myAllocatedLines.size();
        myAllocatedLines.resize(currentAllocatedSize + 1);
        myCurrentDataLine = &(myAllocatedLines.back());
    }

    myCurrentChromID = -1;
    return(true);
}
