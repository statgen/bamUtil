/*
 *  Copyright (C) 2010-2012  Christian Fuchsberger,
 *                           Regents of the University of Michigan
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
////////////////////////////////////////////////////////////////////////
// Recalibrator
//
// by Christian Fuchsberger
//
// - Last modified on June 9th, 2010
//
////////////////////////////////////////////////////////////////////////

#ifndef __RE_CAP_H__
#define __RE_CAP_H__


#include <stdint.h>
#include <string>
// imports from samtools
#include "SamFile.h"
#include "Generic.h"
#include "GenomeSequence.h"
#include "MemoryMapArray.h"
#include "HashErrorModel.h"
#include "Prediction.h"
#include "BaseAsciiMap.h"
#include "BamExecutable.h"

class Recab : public BamExecutable
{
public:
    Recab();
    ~Recab();

    static void recabDescription();
    void description();
    void usage();
    void recabSpecificUsageLine();
    void recabSpecificUsage();
    int execute(int argc, char **argv);

    bool processReadBuildTable(SamRecord& record);
    bool processReadApplyTable(SamRecord& record);
    void modelFitPrediction(const char* outputBase);

    void addRecabSpecificParameters(LongParamContainer& params);

private:
    static const int DEFAULT_MIN_BASE_QUAL = 5;

    // quality String
    typedef struct {
        std::string oldq;
        std::string newq;
    } quality_t;

    void processParams();

    // So external programs can read recab parameters.
    bool myParamsSetup;
    String myRefFile;
    String myDbsnpFile;
    String myQField;  // Quality TAG
    int myMinBaseQual;
    int myBlendedWeight;

    //stats
    uint64_t myBasecounts;
    uint64_t myMappedCount;
    uint64_t myUnMappedCount;
    uint64_t myMappedCountQ;
    uint64_t myBMismatchCount;
    uint64_t myBMatchCount;
    uint64_t myZeroMapQualCount;
    uint64_t myNumDBSnpSkips;
    uint64_t myDupCount;
    uint64_t myMapQual0Count;
    uint64_t myMapQual255Count;

    GenomeSequence* myReferenceGenome;
    mmapArrayBool_t myDbSNP;
    HashErrorModel hasherrormodel;
    Prediction prediction;

    // Make this member data so it reuses the string structures everytime
    // rather than constructing new ones every time.
    quality_t myQualityStrings;

    std::map<std::string, uint16_t> myRg2Id;
    typedef std::pair<std::map<std::string, uint16_t>::iterator, bool> RgInsertReturn;
    std::vector<std::string> myId2Rg;
};

#endif
