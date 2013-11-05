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
// This file contains the processing for the executable option "readIndexedBam"
// which reads an indexed BAM file by chromosome and writes it into a new
// file sorted from reference id -1 to maxRefID.

#ifndef __READ_INDEXED_BAM_H__
#define __READ_INDEXED_BAM_H__

#include "BamExecutable.h"

class ReadIndexedBam : public BamExecutable
{
public:
    static void readIndexedBamDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);
    virtual const char* getProgramName() {return("bam:readIndexedBam");}

private:
    int readIndexedBam(const char* inputFilename,
                       const char* outputFilename,
                       const char* indexFilename);
};

#endif
