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
// This file contains the processing for the executable option "gapInfo"
// which prints information on the gap between read pairs.

#ifndef __GAP_INFO_H__
#define __GAP_INFO_H__

#include "BamExecutable.h"

class GapInfo : public BamExecutable
{
public:
    GapInfo();

    static void gapInfoDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);
    virtual const char* getProgramName() {return("bam:gapInfo");}

private:
    int processFile(const char* inputFileName, const char* outputFileName,
                    const char* refFile, bool detailed,
                    bool checkFirst, bool checkStrand);
};

#endif
