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
/// This file contains the processing for the executable option "squeeze"
/// which attempts to reduce the size of SAM/BAM files by removing unneccessary
/// fields and optionally converting matching bases to '='.

#ifndef __SQUEEZE_H__
#define __SQUEEZE_H__

#include <stack>
#include <list>
#include <map>
#include "BamExecutable.h"
#include "SamFile.h"

class Squeeze : public BamExecutable
{
public:
    Squeeze();
    ~Squeeze();
    static void squeezeDescription();
    void description();
    void binningUsageLine();
    void binningUsage();
    void usage();
    int execute(int argc, char **argv);
    virtual const char* getProgramName() {return("bam:squeeze");}

    void addBinningParameters(LongParamContainer& params);
    int processBinningParam();
    int getQualCharFromQemp(uint8_t qemp);

private:
    void binPhredQuals(int binStartPhred, int binEndPhred);
    void bin(SamRecord& samRecord);

    // Non-phred max
    static const int MAX_QUAL_CHAR = 126;
    static const int MAX_PHRED_QUAL = 93;
    static const int QUAL_CONVERT = 33;
    bool myBinMid;
    bool myBinHigh;
    String myBinQualS;
    String myBinQualF;
    // Non-phred indices
    int myQualBinMap[MAX_QUAL_CHAR+1];
};

#endif
