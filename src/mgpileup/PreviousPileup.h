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

#ifndef __PREVIOUS_PILEUP_H__
#define __PREVIOUS_PILEUP_H__

#include <deque>
#include <queue>
#include <stack>
#include "VcfFile.h"
#include "SamFileHeader.h"

class PreviousPileup
{
public:
    PreviousPileup();
    ~PreviousPileup();

    /// Open the file.
    void open(const std::string& fileName, uint32_t maxStoredLines = 1);

    /// Search for the target chromosome id/position in the Previous Pileup.  
    /// Advances the file until the target position is found or past.
    /// \param targetChromID chromID/position to check for in the previous pileup.
    /// \param targetPosition chromID/position to check for in the previous pileup.
    /// \return true if the target position is found, false if not.
    bool findTarget(int chromID, int position, SamFileHeader& samHeader);

    // Write the current line to the specified output file.    
    bool writeCurrentLine(VcfFile& outputFile);

    bool writeStored(VcfFile& outputFile);

    // Write the stored values up until the specified position.
    bool writeStoredUpToPos(int position, VcfFile& outputFile);

    /// Store the current line for later use, resets the currentLine.
    /// \return false if the current line cannot be stored because too many have 
    /// already been stored. 
    bool storeCurrent();
    

private:
    VcfFile myFile;

    VcfFile::VcfDataLine* myCurrentDataLine;
    std::string myCurrentChromosome;
    int myCurrentChromID;
    
    bool myEof;

    uint32_t myMaxStoredLines;
    std::deque<VcfFile::VcfDataLine> myAllocatedLines;
    std::queue<VcfFile::VcfDataLine*> myStoredLines;
    std::stack<VcfFile::VcfDataLine*> myUnusedLines;
};
#endif
