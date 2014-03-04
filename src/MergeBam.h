/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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
// This file contains the processing for the executable option "mergeBam"
// which merges multiple BAM files appending ReadGroup IDs if necessary.

#ifndef __MERGE_BAM_H__
#define __MERGE_BAM_H__

#include "BamExecutable.h"

class MergeBam : public BamExecutable
{
public:
    MergeBam();
    static void mergeBamDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);
    virtual const char* getProgramName() {return("bam:mergeBam");}
private:
    bool getNextSection(SamFile *in_bams, uint32_t numBams);

    StringArray myRegionArray;
    int32_t myRegionArrayIndex;
    IFILE myRegionFile;
};

#endif
