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
// This file contains the processing for the executable option "extractPosition"
// which writes a file with the reads in the specified region.

#ifndef __EXTRACT_POSITION_H__
#define __EXTRACT_POSITION_H__

#include "BamExecutable.h"
#include "SamFile.h"

class ExtractPosition : public BamExecutable
{
public:
    ExtractPosition();
    static void extractPositionDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);

private:
    bool getNextSection();

    static const int UNSPECIFIED_INT = -1;
    static const int UNSET_REF = -2;

    bool myWithinReg;
    bool myWroteReg;

    int myStart;
    int myEnd;
    int myPrevStart;
    int myPrevEnd;

    int myRefID;
    String myRefName;

    String myPrevRefName;
    int myBedRefID;

    IFILE       myBedFile;
    String      myBedBuffer;
    StringArray myBedColumn;

    SamFile mySamIn;
    SamFileHeader mySamHeader;
};

#endif
