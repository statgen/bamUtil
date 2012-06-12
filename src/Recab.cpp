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
      myQField("")
{
    myBasecounts = 0;
    myMappedCount = 0;
    myUnMappedCount = 0;
    myMappedCountQ = 0;
    myZeroMapQualCount = 0;
    myBunMappedCount = 0;
    myBMappedCount = 0;
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
    std::cerr << "--refFile <ReferenceFile> [--dbsnp <dbsnpFile>] [--blended <weight>] ";
}

void Recab::recabSpecificUsage()
{
    std::cerr << "\nRecab Specific Required Parameters\n";
    std::cerr << "\t--refFile <reference file>    : reference file name" << std::endl;
    std::cerr << "Recab Specific Optional Parameters : " << std::endl;
    std::cerr << "\t--dbsnp <known variance file> : dbsnp file of positions" << std::endl;
    std::cerr << "\t--blended <weight>            : blended model weight" << std::endl;
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
    addRecabSpecificParameters(parameters);

    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            parameters.getLongParameterList()));
    
    inputParameters.Read(argc-1, &(argv[1]));
    
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

    if(myRefFile.IsEmpty())
    {
        usage();
        inputParameters.Status();
        std::cerr << "Missing required --refFile parameter" << std::endl;
        return EXIT_FAILURE;
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
        processRead(samRecord, ANALYZE);

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
        ////////////
        uint16_t  flag = samRecord.getFlag();
        
        if((SamFlag::isDuplicate(flag)) ||
           (samRecord.getMapQuality() == 0) ||
           (!SamFlag::isMapped(flag)))
        {
            // Do nothing.
        }
        else
        {
            // Recalibrate.
            processRead(samRecord,UPDATE);
        }
        samOut.WriteRecord(samHeader, samRecord);
    }

    Logger::gLogger->writeLog("Recalibration successfully finished");
    return EXIT_SUCCESS;
}


void Recab::addRecabSpecificParameters(LongParamContainer& params)
{
    params.addGroup("Required Recab Parameters");
    params.addString("refFile", &myRefFile);
    params.addGroup("Optional Recab Parameters");
    params.addString("dbsnp", &myDbsnpFile);
    params.addInt("blended", &myBlendedWeight);
    // params.addString("qualTag", &myQField);
    myParamsSetup = false;
}


bool Recab::processRead(SamRecord& samRecord, PROC_TYPE processtype)
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

    // Skip Duplicates.
    if(SamFlag::isDuplicate(flag))
    {
        return(false);
    }
    if(!SamFlag::isMapped(flag))
    {
        // Unmapped, skip processing, but increment unmapped count if PROC_WRITE
        if(processtype == ANALYZE)
        {
            myUnMappedCount++;
        }
        return(false);
    }
    // This read is not a duplicate and is mapped.
    if(processtype == ANALYZE)
    {
        myMappedCount++;
    }
    if(samRecord.getMapQuality() == 0)
    {
        // 0 mapping quality, so skip processing.
        return(false);
    }
    // quality>0 & mapped.
    if(processtype == ANALYZE)
    {
        myMappedCountQ++;
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
    	Logger::gLogger->warning("INVALID_GENOME_INDEX (%u) and record skipped... Reference in BAM is different from the ref used here!", mapPos);
        return false;
    }

    if(myQField.IsEmpty())
    {
        myQualityStrings.oldq = samRecord.getQuality();
    }
    else
    {
        myQualityStrings.oldq = samRecord.getString(myQField.c_str()).c_str();
        //printf("%s\n",samRecord.getQuality());
        //printf("%s:%s\n",myQField.c_str(),temp.c_str());
    }

    aligTypes = "";
    Cigar* cigarPtr = samRecord.getCigarInfo();

    if(cigarPtr == NULL)
    {
        Logger::gLogger->warning("Failed to get the cigar");
        return(false);
    }
    myQualityStrings.newq = samRecord.getQuality();

    ////////////////
    ////// iterate sequence
    ////////////////
    genomeIndex_t refPos = 0;
    int32_t refOffset = 0;
    bool prevMatch = false;
    int32_t prevRefOffset = Cigar::INDEX_NA;
    int32_t seqPos = 0;
    int seqIncr = 1;
    if(reverse)
    {
        seqPos = seqLen - 1;
        seqIncr = -1;
    }

    // TODO - what is prebase??????

    //    for (int cycle = 0; cycle < seqLen; cycle++, seqPos += seqIncr)
    int cycle = 0;
    int cycleIncr = 1;
    seqIncr = 1;
    if(reverse)
    {
        cycle = seqLen - 1;
        cycleIncr = -1;
    }
    for (int seqPos = 0; seqPos < seqLen; seqPos++, cycle += cycleIncr)
    {
        // Get the reference offset.
        refOffset = cigarPtr->getRefOffset(seqPos);
        if(refOffset == Cigar::INDEX_NA)
        {
            // Not a match/mismatch, so continue to the next one which will
            // not have a previous match/mismatch.
            prevMatch = false;
            continue;
        }

        if(refOffset != (prevRefOffset + seqIncr))
        {
            // The previous cigar operation was not a match.
            prevMatch = false;
        }

        if(prevMatch)
        {
            // Since the previous one was a match, carry over
            // the already calculated info.
            data.preBase = data.curBase;
        }
        else
        {
            // The previous position was not a match, so there is no
            // previous.
            data.preBase = 'K';
        }

        // Set the next round values.
        prevMatch = true;
        prevRefOffset = refOffset;
        

        // Get the current base.
        data.curBase = samRecord.getSequence(seqPos);

        // Set the cycle
        data.cycle = cycle;

        if(reverse)
        {
            // Complement the current base before checking for
            // dbsnps so it is done for the next position's prebase.
            data.curBase = BaseAsciiMap::base2complement[(unsigned int)(data.curBase)];
        }

        if(data.curBase == 'N')
        {
            data.curBase = 'K';
        }

        //exclude dbSNP sites
        refPos = mapPos + refOffset;
        if(!(myDbsnpFile.IsEmpty()) && (myDbSNP[refPos]==true))
        {
            continue;
        }
        
        // Set the reference & read bases in the Covariates
        data.refBase = (*myReferenceGenome)[refPos];
        
        if(reverse)
        {
            data.refBase = BaseAsciiMap::base2complement[(unsigned int)(data.refBase)];
        }

        // Get quality string
        data.qual = 
            BaseUtilities::getPhredBaseQuality(myQualityStrings.oldq[seqPos]);

        // TODO - make this configurable.
        // skip bases with quality below <5
        if(data.qual < 5)
        {
            continue;
        }

        if(BaseUtilities::areEqual(data.refBase, data.curBase)
           && (BaseAsciiMap::base2int[(unsigned int)(data.curBase)] < 4))
            myBMappedCount++;
        else
            myBunMappedCount++;

        // read
        if(flag & SamFlag::FIRST_READ)
            data.read = 0;
        else
            data.read = 1;

        //read
        if(processtype==ANALYZE)
        {
            hasherrormodel.setCell(data);
            myBasecounts++;
        }
        //write
        if(processtype==UPDATE)
        {
            // TODO - implement old quality flag.
            bool oldQualityFlag = false;

            // Update quality score
            uint8_t qemp = hasherrormodel.getQemp(data);
            //if(qemp>4)
            myQualityStrings.newq[seqPos] = qemp+33;
            // otherwise keep the old one

            if(oldQualityFlag)
            {
                samRecord.addTag("QQ", 'Z', myQualityStrings.oldq.c_str());
            }
            samRecord.setQuality(myQualityStrings.newq.c_str());
        }
    }
    return true;
}


void Recab::modelFitPrediction(const char* outputBase)
{
    Logger::gLogger->writeLog("# mapped Reads observed: %ld Q>0: %ld",myMappedCount,myMappedCountQ);
    Logger::gLogger->writeLog("# unmapped Reads observed: %ld",myUnMappedCount);
    Logger::gLogger->writeLog("# Bases observed: %ld / ZeroQ: %ld - %ld %ld",myBasecounts,myZeroMapQualCount,
                              myBMappedCount,myBunMappedCount);
    ////////////////////////
    ////////////////////////
    //// Model fitting + prediction
    std::string modelfile = outputBase;
    modelfile += ".model";
    std::string recabFile = outputBase;
    recabFile += ".recab";

    prediction.setErrorModel(&(hasherrormodel));

    Logger::gLogger->writeLog("Start model fitting!");
    if(prediction.fitModel(true,modelfile))
        prediction.outModel();
    else
        Logger::gLogger->error("Could not fit model!");

    hasherrormodel.addPrediction(prediction.getModel(),myBlendedWeight);

    Logger::gLogger->writeLog("Writing recalibration table %s",recabFile.c_str());
    if(!(hasherrormodel.writeAllTableMM(recabFile.c_str(),
                                        myId2Rg)))
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
            Logger::gLogger->writeLog("No dbSNPFile File");
        }
        else if(myReferenceGenome->loadDBSNP(myDbSNP,myDbsnpFile.c_str() ))
        {
            Logger::gLogger->error("Failed to open dbSNP file.");
        }
    }
    myParamsSetup = true;
}
