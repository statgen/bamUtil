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

//////////////////////////////////////////////////////////////////////////
// This file contains the processing for the executable option "clipOverlap"
// which clips overlapping read pairs.

#ifndef __CLIP_OVERLAP_H__
#define __CLIP_OVERLAP_H__

#include "BamExecutable.h"
#include "SamFile.h"

class ClipOverlap : public BamExecutable
{
public:
    ClipOverlap();

    static void clipOverlapDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);

private:
    // Clip overlapping reads and strands that extend past the other strand.
    void clip(SamRecord& firstRecord, SamRecord& secondRecord);

    // Clip an entire read.
    void clipEntire(SamRecord& record);

    // Calculate the sum of the qualities at read positions starting at 
    // startPos and ending with endPos (included).  Stops reading if a 0 is
    // found, indicating the end of the string.
    int32_t getSumQual(SamRecord& record, int32_t startPos, int32_t endPos);

    String myStoreOrig;
};

#endif
