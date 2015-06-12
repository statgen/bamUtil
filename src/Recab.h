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
#include "Squeeze.h"

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
    virtual const char* getProgramName() {return("bam:recab");}

    bool processReadBuildTable(SamRecord& record);
    bool processReadApplyTable(SamRecord& record);
    void modelFitPrediction(const char* outputBase);

    void addRecabSpecificParameters(LongParamContainer& params);
    int processRecabParam();

private:
    static const int DEFAULT_MIN_BASE_QUAL = 5;
    static const int DEFAULT_MAX_BASE_QUAL = 50;

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
    String myStoreQualTag;  // Store old quality into this TAG
    String myBuildExcludeFlags;
    String myApplyExcludeFlags;
    uint16_t myIntBuildExcludeFlags;
    uint16_t myIntApplyExcludeFlags;
    int myMinBaseQual;
    int myMaxBaseQual;
    int myBlendedWeight;
    bool myFitModel;
    bool myFast;
    bool myKeepPrevDbsnp;
    bool myKeepPrevNonAdjacent;
    bool myLogReg;

    // Per read counts
    uint64_t myMappedCount;
    uint64_t myUnMappedCount;
    uint64_t mySecondaryCount;
    uint64_t myDupCount;
    uint64_t myQCFailCount;
    uint64_t myMapQual0Count;
    uint64_t myMapQual255Count;
    uint64_t myNumBuildSkipped;
    uint64_t myNumBuildReads;
    uint64_t myNumApplySkipped;
    uint64_t myNumApplyReads;

    // Couldn't find quality tag, so using current quality.
    uint64_t myNumQualTagErrors;

    // Per base counts.
    uint64_t myNumDBSnpSkips;
    uint64_t mySubMinQual;
    uint64_t myAmbiguous;
    uint64_t myBMatchCount;
    uint64_t myBMismatchCount;

    // should be sum of myBMatchCount & myBMismatchCount
    uint64_t myBasecounts;

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

    // Squeeze
    Squeeze mySqueeze;
};

#endif
