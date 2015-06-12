/*
 *  Copyright (C) 2010-2015  Christian Fuchsberger,
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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <ctime>
#include <unistd.h>
#include <math.h>
//
#include "Error.h"
#include "Recab.h"
#include "Logger.h"
#include "Parameters.h"
#include "SimpleStats.h"
#include "BaseUtilities.h"
#include "SamFlag.h"
#include "BgzfFileType.h"

// STL headers
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>

Recab::Recab()
    : myParamsSetup(false),
      myRefFile(""),
      myDbsnpFile(""),
      myQField(""),
      myStoreQualTag(""),
      myBuildExcludeFlags("0x0704"),
      myApplyExcludeFlags("0x0000"),
      myIntBuildExcludeFlags(0),
      myIntApplyExcludeFlags(0)
{
    myMappedCount = 0;
    myUnMappedCount = 0;
    mySecondaryCount = 0;
    myDupCount = 0;
    myQCFailCount = 0;
    myMapQual0Count = 0;
    myMapQual255Count = 0;
    myNumBuildSkipped = 0;
    myNumBuildReads = 0;
    myNumApplySkipped = 0;
    myNumApplyReads = 0;
    myNumQualTagErrors = 0;
    myNumDBSnpSkips = 0;
    mySubMinQual = 0;
    myAmbiguous = 0;
    myBMatchCount = 0;
    myBMismatchCount = 0;
    myBasecounts = 0;

    myBlendedWeight = 0;
    myFitModel = false;
    myFast = false;
    myKeepPrevDbsnp = false;
    myKeepPrevNonAdjacent = false;
    myLogReg = false;
    myMinBaseQual = DEFAULT_MIN_BASE_QUAL;
    myMaxBaseQual = DEFAULT_MAX_BASE_QUAL;
}


Recab::~Recab()
{
}


void Recab::recabDescription()
{
    std::cerr << " recab - Recalibrate\n";
}


void Recab::description()
{
    recabDescription();
}


void Recab::usage()
{
    std::cerr << "Usage: ./bam recab (options) --in <InputBamFile> --out <OutputFile> [--log <logFile>] [--verbose] [--noeof] [--params] ";
    recabSpecificUsageLine();
    std::cerr << std::endl << std::endl;

    std::cerr << "Required General Parameters :" << std::endl;
    std::cerr << "\t--in <infile>   : input BAM file name" << std::endl;
    std::cerr << "\t--out <outfile> : output recalibration file name" << std::endl;
    std::cerr << "Optional General Parameters : " << std::endl;
    std::cerr << "\t--log <logfile> : log and summary statistics (default: [outfile].log)" << std::endl;
    std::cerr << "\t--verbose       : Turn on verbose mode" << std::endl;
    std::cerr << "\t--noeof         : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t--params        : print the parameter settings" << std::endl;
   recabSpecificUsage();
    std::cerr << "\n" << std::endl;
}


void Recab::recabSpecificUsageLine()
{
    std::cerr << "--refFile <ReferenceFile> [--dbsnp <dbsnpFile>] [--minBaseQual <minBaseQual>] [--maxBaseQual <maxBaseQual>] [--blended <weight>] [--fitModel] [--fast] [--keepPrevDbsnp] [--keepPrevNonAdjacent] [--useLogReg] [--qualField <tag>] [--storeQualTag <tag>] [--buildExcludeFlags <flag>] [--applyExcludeFlags <flag>] ";
    mySqueeze.binningUsageLine();
}

void Recab::recabSpecificUsage()
{
    std::cerr << "\nRecab Specific Required Parameters\n";
    std::cerr << "\t--refFile <reference file>    : reference file name" << std::endl;
    std::cerr << "Recab Specific Optional Parameters : " << std::endl;
    std::cerr << "\t--dbsnp <known variance file> : dbsnp file of positions" << std::endl;
    std::cerr << "\t--minBaseQual <minBaseQual>   : minimum base quality of bases to recalibrate (default: " << DEFAULT_MIN_BASE_QUAL << ")" << std::endl;
    std::cerr << "\t--maxBaseQual <maxBaseQual>   : maximum recalibrated base quality (default: " << DEFAULT_MAX_BASE_QUAL << ")" << std::endl;
    std::cerr << "\t--blended <weight>            : blended model weight" << std::endl;
    std::cerr << "\t--fitModel                    : check if the logistic regression model fits the data" << std::endl;
    std::cerr << "\t                                overriden by fast, but automatically applied by useLogReg" << std::endl;
    std::cerr << "\t--fast                        : use a compact representation that only allows:" << std::endl;
    std::cerr << "\t                                   * at most " << (BaseData::getFastMaxRG() + 1) << " Read Groups" << std::endl;
    std::cerr << "\t                                   * maximum quality " << BaseData::getFastMaxQual() << std::endl;
    std::cerr << "\t                                   * at most " << BaseData::getFastMaxCycles() << " cycles" << std::endl;
    std::cerr << "\t                                overrides fitModel, but is overridden by useLogReg" << std::endl;
    std::cerr << "\t                                uses up to about 2.25G more memory than running without --fast." << std::endl;
    std::cerr << "\t--keepPrevDbsnp               : do not exclude entries where the previous base is in dbsnp when\n";
    std::cerr << "\t                                building the recalibration table" << std::endl;
    std::cerr << "\t                                By default they are excluded from the table." << std::endl;
    std::cerr << "\t--keepPrevNonAdjacent         : do not exclude entries where the previous base is not adjacent\n";
    std::cerr << "\t                                (not a Cigar M/X/=) when building the recalibration table" << std::endl;
    std::cerr << "\t                                By default they are excluded from the table (except the first cycle)." << std::endl;
    std::cerr << "\t--useLogReg                   : use logistic regression calculated quality for the new quality" << std::endl;
    std::cerr << "\t                                automatically applies fitModel and overrides fast." << std::endl;
    std::cerr << "\t--qualField <quality tag>     : tag to get the starting base quality\n";
    std::cerr << "\t                                (default is to get it from the Quality field)" << std::endl;
    std::cerr << "\t--storeQualTag <quality tag>  : tag to store the previous quality into" << std::endl;
    std::cerr << "\t--buildExcludeFlags <flag>    : exclude reads with any of these flags set when building the" << std::endl;
    std::cerr << "\t                                recalibration table" << std::endl;
    std::cerr << "\t--applyExcludeFlags <flag>    : do not apply the recalibration table to any reads with any of these flags set" << std::endl;
    mySqueeze.binningUsage();
}


int Recab::execute(int argc, char *argv[])
{
    bool verboseFlag = false;

    String inFile,outFile,logFile;

    bool noeof = false;
    bool params = false;

    SamFile samIn,samOut;

    ParameterList inputParameters;

    LongParamContainer parameters;

    parameters.addGroup("Required Generic Parameters");
    parameters.addString("in", &inFile);
    parameters.addString("out", &outFile);
    parameters.addGroup("Optional Generic Parameters");
    parameters.addString("log", &logFile);
    parameters.addBool("verbose", &verboseFlag);
    parameters.addBool("noeof", &noeof);
    parameters.addBool("params", &params);
    parameters.addPhoneHome(VERSION);
    addRecabSpecificParameters(parameters);
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            parameters.getLongParameterList()));
    
    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);
    
    // If no eof block is required for a bgzf file, set the bgzf file type to 
    // not look for it.
    if(noeof)
    {
        // Set that the eof block is not required.
        BgzfFileType::setRequireEofBlock(false);
    }

    if(inFile.IsEmpty())
    {
        usage();
        inputParameters.Status();
        std::cerr << "Missing required --in parameter" << std::endl;
        return EXIT_FAILURE;
    }

    if(outFile.IsEmpty())
    {
        usage();
        inputParameters.Status();
        std::cerr << "Missing required --out parameter" << std::endl;
        return EXIT_FAILURE;
    }

    int status = processRecabParam();
    if(status != 0)
    {
        inputParameters.Status();
        return(status);
    }

    if ( logFile.IsEmpty() )
    {
        logFile = outFile + ".log";
    }
  
    if(params)
    {
        inputParameters.Status();
    }
    
    Logger::gLogger = new Logger(logFile.c_str(), verboseFlag);

    ////////////////
    //////  Errormodel
    Logger::gLogger->writeLog("Initialize errormodel structure...");

    ////////////////////////////////////////
    // SAM/BAM file open
    ////////////////////////////////////////
    ////////////////////////////////////////

    // Iterate SAM records
    if(!samIn.OpenForRead(inFile.c_str()))
    {
        Logger::gLogger->error("Failed to open SAM/BAM file %s",inFile.c_str() );
        return EXIT_FAILURE;
    }

    Logger::gLogger->writeLog("Start iterating SAM/BAM file %s",inFile.c_str());

    time_t now = time(0);
    tm* localtm = localtime(&now);

    Logger::gLogger->writeLog("Start: %s", asctime(localtm));
    SamRecord samRecord;
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    srand (time(NULL));

    int numRecs = 0;
    while(samIn.ReadRecord(samHeader, samRecord) == true)
    {
        processReadBuildTable(samRecord);

        //Status info
        numRecs++;
        if(verboseFlag)
        {
            if(numRecs%10000000==0)
                Logger::gLogger->writeLog("%ld records processed", numRecs);
        }
    }

    now = time(0);
    localtm = localtime(&now);
    Logger::gLogger->writeLog("End: %s", asctime(localtm));

    modelFitPrediction(outFile);

    Logger::gLogger->writeLog("Writing recalibrated file %s",outFile.c_str());

    ////////////////////////
    ////////////////////////
    //// Write file
    samIn.OpenForRead(inFile.c_str());
    samOut.OpenForWrite(outFile.c_str());
    samIn.ReadHeader(samHeader);
    samOut.WriteHeader(samHeader);
    
    while(samIn.ReadRecord(samHeader, samRecord) == true)
    {
        // Recalibrate.
        processReadApplyTable(samRecord);
        samOut.WriteRecord(samHeader, samRecord);
    }

    Logger::gLogger->writeLog("Total # Reads recab table not applied to: %ld", myNumApplySkipped);
    Logger::gLogger->writeLog("Total # Reads recab table applied to: %ld", myNumApplyReads);
    Logger::gLogger->writeLog("Recalibration successfully finished");
    return EXIT_SUCCESS;
}


void Recab::addRecabSpecificParameters(LongParamContainer& params)
{
    params.addGroup("Required Recab Parameters");
    params.addString("refFile", &myRefFile);
    params.addGroup("Optional Recab Parameters");
    params.addString("dbsnp", &myDbsnpFile);
    params.addInt("minBaseQual", &myMinBaseQual);
    params.addInt("maxBaseQual", &myMaxBaseQual);
    params.addInt("blended", &myBlendedWeight);
    params.addBool("fitModel", &myFitModel);
    params.addBool("fast", &myFast);
    params.addBool("keepPrevDbsnp", &myKeepPrevDbsnp);
    params.addBool("keepPrevNonAdjacent", &myKeepPrevNonAdjacent);
    params.addBool("useLogReg", &myLogReg);
    params.addString("qualField", &myQField);
    params.addString("storeQualTag", &myStoreQualTag);
    params.addString("buildExcludeFlags", &myBuildExcludeFlags);
    params.addString("applyExcludeFlags", &myApplyExcludeFlags);
    myParamsSetup = false;
    mySqueeze.addBinningParameters(params);

}


int Recab::processRecabParam()
{
    if(myRefFile.IsEmpty())
    {
        std::cerr << "Missing required --refFile parameter" << std::endl;
        return EXIT_FAILURE;
    }
    return(mySqueeze.processBinningParam());
}


bool Recab::processReadBuildTable(SamRecord& samRecord)
{
    static BaseData data;
    static std::string chromosomeName;
    static std::string readGroup;
    static std::string aligTypes;

    int seqLen = samRecord.getReadLength();
    
    // Check if the parameters have been processed.
    if(!myParamsSetup)
    {
        // This throws an exception if the reference cannot be setup.
        processParams();
    }

    uint16_t  flag = samRecord.getFlag();

    if(!SamFlag::isMapped(flag))
    {
        // Unmapped, skip processing
        ++myUnMappedCount;
    }
    else
    {
        // This read is mapped.
        ++myMappedCount;
    }

    if(SamFlag::isSecondary(flag))
    {
        // Secondary read
        ++mySecondaryCount;
    }
    if(SamFlag::isDuplicate(flag))
    {
        ++myDupCount;
    }
    if(SamFlag::isQCFailure(flag))
    {
        ++myQCFailCount;
    }

    // Check if the flag contains an exclude.
    if((flag & myIntBuildExcludeFlags) != 0)
    {
        // Do not use this read for building the recalibration table.
        ++myNumBuildSkipped;
        return(false);
    }

    if(samRecord.getMapQuality() == 0)
    {
        // 0 mapping quality, so skip processing.
        ++myMapQual0Count;
        ++myNumBuildSkipped;
        return(false);
    }
    if(samRecord.getMapQuality() == 255)
    {
        // 255 mapping quality, so skip processing.
        ++myMapQual255Count;
        ++myNumBuildSkipped;
        return(false);
    }
    
    chromosomeName = samRecord.getReferenceName();
    readGroup = samRecord.getString("RG").c_str();

    // Look for the read group in the map.
    // TODO - extra string constructor??
    RgInsertReturn insertRet = 
        myRg2Id.insert(std::pair<std::string, uint16_t>(readGroup, 0));
    if(insertRet.second == true)
    {
        // New element inserted.
        insertRet.first->second = myId2Rg.size();
        myId2Rg.push_back(readGroup);
    }

    data.rgid = insertRet.first->second;


    //reverse
    bool reverse;
    if(SamFlag::isReverse(flag))
        reverse = true;
    else
        reverse = false;

    if(myReferenceGenome == NULL)
    {
        throw std::runtime_error("Failed to setup Reference File.\n");
    }

    genomeIndex_t mapPos = 
        myReferenceGenome->getGenomePosition(chromosomeName.c_str(), 
                                             samRecord.get1BasedPosition());

    if(mapPos==INVALID_GENOME_INDEX)
    {
    	Logger::gLogger->warning("INVALID_GENOME_INDEX (chrom:pos %s:%ld) and record skipped... Reference in BAM is different from the ref used here!", chromosomeName.c_str(), samRecord.get1BasedPosition());

        ++myNumBuildSkipped;
        return false;
    }

    if(!myQField.IsEmpty())
    {
        // Check if there is an old quality.
        const String* oldQPtr = 
            samRecord.getStringTag(myQField.c_str());
        if((oldQPtr != NULL) && (oldQPtr->Length() == seqLen))
        {
            // There is an old quality, so use that.
            myQualityStrings.oldq = oldQPtr->c_str();
        }
        else
        {
            // Tag was not found, so use the current quality.
            ++myNumQualTagErrors;
            if(myNumQualTagErrors == 1)
            {
                Logger::gLogger->warning("Recab: %s tag was not found/invalid, so using the quality field in records without the tag", myQField.c_str());
            }
            myQualityStrings.oldq = samRecord.getQuality();
        }
        //printf("%s\n",samRecord.getQuality());
        //printf("%s:%s\n",myQField.c_str(),temp.c_str());
    }
    else
    {
        myQualityStrings.oldq = samRecord.getQuality();
    }

    if(myQualityStrings.oldq.length() != (unsigned int)seqLen)
    {
        Logger::gLogger->warning("Quality is not the correct length, so skipping recalibration on that record.");
        ++myNumBuildSkipped;
        return(false);
    }

    aligTypes = "";
    Cigar* cigarPtr = samRecord.getCigarInfo();

    if(cigarPtr == NULL)
    {
        Logger::gLogger->warning("Failed to get the cigar");
        ++myNumBuildSkipped;
        return(false);
    }

    // This read will be used for building the recab table.
    ++myNumBuildReads;

    ////////////////
    ////// iterate sequence
    ////////////////
    genomeIndex_t refPos = 0;
    int32_t refOffset = 0;
    int32_t prevRefOffset = Cigar::INDEX_NA;
    int32_t seqPos = 0;
    int seqIncr = 1;
    if(reverse)
    {
        seqPos = seqLen - 1;
        seqIncr = -1;
    }

    // read
    if(!SamFlag::isPaired(flag) || SamFlag::isFirstFragment(flag))
        // Mark as first if it is not paired or if it is the
        // first in the pair.
        data.read = 0;
    else
        data.read = 1;

    // Set unsetbase for curBase.
    // This will be used for the prebase of cycle 0.
    data.curBase = 'K';

    for (data.cycle = 0; data.cycle < seqLen; data.cycle++, seqPos += seqIncr)
    {
        // Store the previous current base in preBase.
        data.preBase = data.curBase;

        // Get the current base before checking if we are going to
        // process this position so it will be set for the next position.
        data.curBase = samRecord.getSequence(seqPos);
        if(reverse)
        {
            // Complement the current base.
            // The prebase is already complemented.
            data.curBase = 
                BaseAsciiMap::base2complement[(unsigned int)(data.curBase)];
        }
        
        // Get the reference offset.
        refOffset = cigarPtr->getRefOffset(seqPos);
        if(refOffset == Cigar::INDEX_NA)
        {
            // Not a match/mismatch, so continue to the next one which will
            // not have a previous match/mismatch.
            // Set previous ref offset to a negative so
            // the next one won't be kept.
            prevRefOffset = -2;
            continue;
        }

        // This one is a match.
        refPos = mapPos + refOffset;

        // Check to see if we should process this position.
        // Do not process if it is cycle 0 and:
        //   1) current base is in dbsnp
        if(data.cycle == 0)
        {
            if(!(myDbsnpFile.IsEmpty()) && myDbSNP[refPos])
            {
                // Save the previous reference offset.
                ++myNumDBSnpSkips;
                prevRefOffset = refOffset;
                continue;
            }
        }
        else
        {
            // Do not process if it is not cycle 0 and:
            //   1) previous reference position not adjacent 
            //      (not a match/mismatch)
            //   2) previous base is in dbsnp
            //   3) current base is in dbsnp
            if((!myKeepPrevNonAdjacent && (refOffset != (prevRefOffset + seqIncr))) ||
               (data.preBase == 'K'))
            {
                // Save the previous reference offset.
                prevRefOffset = refOffset;
                continue;
            }
            if(!(myDbsnpFile.IsEmpty()) && 
               (myDbSNP[refPos] ||
                (!myKeepPrevDbsnp && myDbSNP[refPos - seqIncr])))
            {
                ++myNumDBSnpSkips;
                // Save the previous reference offset.
                prevRefOffset = refOffset;
                continue;
            }
       }
        
        // Save the previous reference offset.
        prevRefOffset = refOffset;

        // Set the reference & read bases in the Covariates
        char refBase = (*myReferenceGenome)[refPos];

        if(BaseUtilities::isAmbiguous(refBase))
        {
            // N reference, so skip it when building the table.
            ++myAmbiguous;
            continue;
        }

        if(reverse)
        {
            refBase = BaseAsciiMap::base2complement[(unsigned int)(refBase)];
        }

        // Get quality char
        data.qual = 
            BaseUtilities::getPhredBaseQuality(myQualityStrings.oldq[seqPos]);

        // skip bases with quality below the minimum set.
        if(data.qual < myMinBaseQual)
        {
            ++mySubMinQual;
            continue;
        }

        if(BaseUtilities::areEqual(refBase, data.curBase)
           && (BaseAsciiMap::base2int[(unsigned int)(data.curBase)] < 4))
            myBMatchCount++;
        else
            myBMismatchCount++;

        hasherrormodel.setCell(data, refBase);
        myBasecounts++;
    }
    return true;
}


bool Recab::processReadApplyTable(SamRecord& samRecord)
{
    static BaseData data;
    static std::string readGroup;
    static std::string aligTypes;

    int seqLen = samRecord.getReadLength();

    uint16_t  flag = samRecord.getFlag();

    // Check if the flag contains an exclude.
    if((flag & myIntApplyExcludeFlags) != 0)
    {
        // Do not apply the recalibration table to this read.
        ++myNumApplySkipped;
        return(false);
    }
    ++myNumApplyReads;
   
    readGroup = samRecord.getString("RG").c_str();

    // Look for the read group in the map.
    // TODO - extra string constructor??
    RgInsertReturn insertRet = 
        myRg2Id.insert(std::pair<std::string, uint16_t>(readGroup, 0));
    if(insertRet.second == true)
    {
        // New element inserted.
        insertRet.first->second = myId2Rg.size();
        myId2Rg.push_back(readGroup);
    }

    data.rgid = insertRet.first->second;

    if(!myQField.IsEmpty())
    {
        // Check if there is an old quality.
        const String* oldQPtr =
            samRecord.getStringTag(myQField.c_str());
        if((oldQPtr != NULL) && (oldQPtr->Length() == seqLen))
        {
            // There is an old quality, so use that.
            myQualityStrings.oldq = oldQPtr->c_str();
        }
        else
        {
            myQualityStrings.oldq = samRecord.getQuality();
        }
    }
    else
    {
        myQualityStrings.oldq = samRecord.getQuality();
    }

    if(myQualityStrings.oldq.length() != (unsigned int)seqLen)
    {
        Logger::gLogger->warning("Quality is not the correct length, so skipping recalibration on that record.");
        return(false);
    }

    myQualityStrings.newq.resize(seqLen);

    ////////////////
    ////// iterate sequence
    ////////////////
    int32_t seqPos = 0;
    int seqIncr = 1;

    bool reverse;
    if(SamFlag::isReverse(flag))
    {
        reverse = true;
        seqPos = seqLen - 1;
        seqIncr = -1;
    }
    else
        reverse = false;

    // Check which read - this will be the same for all positions, so 
    // do this outside of the smaller loop.
    if(!SamFlag::isPaired(flag) || SamFlag::isFirstFragment(flag))
        // Mark as first if it is not paired or if it is the
        // first in the pair.
        data.read = 0;
    else
        data.read = 1;

    // Set unsetbase for curBase.
    // This will be used for the prebase of cycle 0.
    data.curBase = 'K';

    for (data.cycle = 0; data.cycle < seqLen; data.cycle++, seqPos += seqIncr)
    {
        // Set the preBase to the previous cycle's current base.
        // For cycle 0, curBase was set to a default value.
        data.preBase = data.curBase;

        // Get the current base.
        data.curBase = samRecord.getSequence(seqPos);

        if(reverse)
        {
            // Complement the current base.
            data.curBase =
                BaseAsciiMap::base2complement[(unsigned int)(data.curBase)];
        }

        // Get quality
        data.qual = 
            BaseUtilities::getPhredBaseQuality(myQualityStrings.oldq[seqPos]);

        // skip bases with quality below the minimum set.
        if(data.qual < myMinBaseQual)
        {
            myQualityStrings.newq[seqPos] = myQualityStrings.oldq[seqPos];
            continue;
        }

        // Update quality score
        uint8_t qemp = hasherrormodel.getQemp(data);
        if(qemp > myMaxBaseQual)
        {
            qemp = myMaxBaseQual;
        }
        myQualityStrings.newq[seqPos] = mySqueeze.getQualCharFromQemp(qemp);
    }

    if(!myStoreQualTag.IsEmpty())
    {
        samRecord.addTag(myStoreQualTag, 'Z', myQualityStrings.oldq.c_str());
    }
    samRecord.setQuality(myQualityStrings.newq.c_str());

    return true;
}


void Recab::modelFitPrediction(const char* outputBase)
{
    Logger::gLogger->writeLog("# mapped Reads observed: %ld", myMappedCount);
    Logger::gLogger->writeLog("# unmapped Reads observed: %ld", myUnMappedCount);
    Logger::gLogger->writeLog("# Secondary Reads observed: %ld", mySecondaryCount);
    Logger::gLogger->writeLog("# Duplicate Reads observed: %ld", myDupCount);
    Logger::gLogger->writeLog("# QC Failure Reads observed: %ld", myQCFailCount);
    Logger::gLogger->writeLog("# Mapping Quality 0 Reads skipped: %ld", myMapQual0Count);
    Logger::gLogger->writeLog("# Mapping Quality 255 Reads skipped: %ld", myMapQual255Count);
    Logger::gLogger->writeLog("Total # Reads skipped for building recab table: %ld", myNumBuildSkipped);
    Logger::gLogger->writeLog("Total # Reads used for building recab table: %ld", myNumBuildReads);

    Logger::gLogger->writeLog("# Bases observed: %ld - #match: %ld; #mismatch: %ld",
                              myBasecounts, myBMatchCount, myBMismatchCount);
    Logger::gLogger->writeLog("# Bases Skipped for DBSNP: %ld, for BaseQual < %ld: %ld, ref 'N': %ld", 
                              myNumDBSnpSkips, myMinBaseQual, mySubMinQual, myAmbiguous);
    if(myNumQualTagErrors != 0)
    {
        Logger::gLogger->warning("%ld records did not have tag %s or it was invalid, so the quality field was used for those records.", myNumQualTagErrors, myQField.c_str());
    }

    ////////////////////////
    ////////////////////////
    if(myLogReg || myFitModel)
    {
        //// Model fitting + prediction
        std::string modelfile = outputBase;
        modelfile += ".model";
        
        prediction.setErrorModel(&(hasherrormodel));
        
        Logger::gLogger->writeLog("Start model fitting!");
        if(prediction.fitModel(true,modelfile))
            prediction.outModel();
        else
            Logger::gLogger->error("Could not fit model!");
        
        hasherrormodel.addPrediction(prediction.getModel(),myBlendedWeight);

        std::string recabFile = outputBase;
        recabFile += ".recab";
        Logger::gLogger->writeLog("Writing recalibration table %s",recabFile.c_str());
        if(!(hasherrormodel.writeTableQemp(recabFile,
                                           myId2Rg, true)))
            Logger::gLogger->error("Writing errormodel not possible!");
    }

    std::string qempFile = outputBase;
    qempFile += ".qemp";
    Logger::gLogger->writeLog("Writing recalibration table %s",qempFile.c_str());
    if(!(hasherrormodel.writeTableQemp(qempFile,
                                       myId2Rg, false)))
        Logger::gLogger->error("Writing errormodel not possible!");
}


void Recab::processParams()
{
    if(myReferenceGenome == NULL)
    {
        Logger::gLogger->writeLog("Open reference");
        myReferenceGenome = new GenomeSequence(myRefFile);
        if(myReferenceGenome == NULL)
        {
            throw std::runtime_error("Failed to open Reference File.\n");
        }
        Logger::gLogger->writeLog("Done! Sequence length %u",
                                  myReferenceGenome->sequenceLength());
        //dbSNP
        if(myDbsnpFile.IsEmpty())
        {
            Logger::gLogger->writeLog("No dbSNP File");
        }
        else if(myReferenceGenome->loadDBSNP(myDbSNP,myDbsnpFile.c_str()))
        {
            Logger::gLogger->error("Failed to open dbSNP file.");
        }
    }

    if(myLogReg && myFast)
    {
        Logger::gLogger->writeLog("Cannot use both --fast & --useLogReg, so just using --useLogReg");
        myFast = false;
    }

    if(myFast)
    {
        myFitModel = false;
    }

   if(myLogReg)
    {
        // Set fitModel.
        myFitModel = true;
    }

    HashErrorModel::setUseLogReg(myLogReg);
    HashErrorModel::setUseFast(myFast);

    myIntBuildExcludeFlags = myBuildExcludeFlags.AsInteger();
    myIntApplyExcludeFlags = myApplyExcludeFlags.AsInteger();

    myParamsSetup = true;
}
