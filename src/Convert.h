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
// This file contains the processing for the executable option "convert"
// which reads an SAM/BAM file and writes a SAM/BAM file (it can convert 
// between SAM and BAM formats).

#ifndef __CONVERT_H__
#define __CONVERT_H__

#include "BamExecutable.h"

class Convert : public BamExecutable
{
public:
    static void convertDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);
    virtual const char* getProgramName() {return("bam:convert");}
};

#endif
