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
#include "ReCab.h"
#include "Logger.h"
#include "Parameters.h"
#include "SimpleStats.h"
#include "BaseUtilities.h"
#include "SamFlag.h"

// STL headers
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>

#define BASEQC_VERSION "0.31"
#define PROCREAD 0
#define PROCWRITE 1
#define MAXITERATION 10000
#define READSEPARATORS ":_"

void ReCab::recabDescription()
{
    std::cerr << " recab - Recalibrate\n";
}


void ReCab::description()
{
    recabDescription();
}


void ReCab::usage()
{
    std::cerr << "Usage: ReCab (options) --i <InputBamFile> --o <OutputFile>\n" << std::endl;
    std::cerr << "Required parameters :" << std::endl;
    std::cerr << "--i [infile]  : input BAM file name" << std::endl;
    std::cerr << "--o [outfile] : output recalibration file name" << std::endl;
    std::cerr << "--r [reference file]" << std::endl;
    std::cerr << "Optional parameters : " << std::endl;
    std::cerr << "--s [known variance file]" << std::endl;
    std::cerr << "--l [logfile] : log and summary statistics (default: [outfile].log)" << std::endl;
    std::cerr << "--v : Turn on verbose mode" << std::endl;
    std::cerr << "--first [n] : First n reads are processed" << std::endl;
    std::cerr << "--b [weight]: blended model weight" << std::endl;
    std::cerr << "\n" << std::endl;
    std::cerr << "version: " << BASEQC_VERSION << std::endl;
    std::cerr << "\n" << std::endl;
}

ReCab::ReCab()
    : myBaseAsciiMap()
{
    basecounts=0;
    mappedCount = 0;
    unMappedCount = 0;
    mappedCountQ = 0;
    zeroMapQualCount = 0;
    BunMappedCount =0;
    BMappedCount =0;
}

ReCab::~ReCab(){

}

int ReCab::nt2idx(char c) {
    switch(toupper(c)){
        case 'A': return(0);
        case 'C': return(1);
        case 'G': return(2);
        case 'T': return(3);
            //case 'N': return(4);
        default:  return(5);
    }
}

char ReCab::complement(char c)
{
    switch(toupper(c)){
        case 'A': return('T');
        case 'C': return('G');
        case 'G': return('C');
        case 'T': return('A');
        default:  return(toupper(c));
    }
}

int ReCab::nt2idx2[256];

void ReCab::conversionTable(){
    for(int i=0; i<256; i++) 
        nt2idx2[i] = nt2idx(i);
}

// from Hyun
uint32_t ReCab::addTokenizedStrings(const std::string& str, const std::string& delimiters, std::vector<std::string>& tokens)
{
    std::string::size_type delimPos = 0, tokenPos = 0, pos = 0;
    uint32_t numAddedTokens = 0;

    if( str.length() < 1 )  return 0;

    while( true ){
        delimPos = str.find_first_of(delimiters, pos);
        tokenPos = str.find_first_not_of(delimiters, pos);

        if(std::string::npos != delimPos){
            if(std::string::npos != tokenPos){
                if( tokenPos < delimPos ) {
                    tokens.push_back(str.substr(pos,delimPos-pos));
                    ++numAddedTokens;
                }
                else {
                    tokens.push_back("");
                    ++numAddedTokens;
                }
            }
            else {
                tokens.push_back("");
                ++numAddedTokens;
            }
            pos = delimPos+1;
        }
        else {
            if( std::string::npos != tokenPos ){
                tokens.push_back(str.substr(pos));
                ++numAddedTokens;
            }
            else {
                tokens.push_back("");
                ++numAddedTokens;
            }
            break;
        }
    }
    return numAddedTokens;
}


bool ReCab::processRead(SamRecord& samRecord,int processtype, ReCab::quality_t& quality_strings){
    int offset = 0;
    int offset2 = 0;
    int cycleIdx = 0;
    int seqpos;
    int nSeq = samRecord.getReadLength();
    uint16_t cread,cprebase,cnexbase,cqual;
    uint8_t cobs,cref;
    genomeIndex_t refpos;
    bool reverse;
    char refBase, readBase, prebase,nexbase;
    baseCV basecv;
    std::string chromosomeName = samRecord.getReferenceName();
    std::string readGroup = samRecord.getString("RG").c_str();

    uint16_t  flag = samRecord.getFlag();

    //reverse
    if(flag & 0x0010)
        reverse = true;
    else
        reverse = false;

    genomeIndex_t mapPos = referenceGenome.getGenomePosition(chromosomeName.c_str(), samRecord.get1BasedPosition());

    if(mapPos==INVALID_GENOME_INDEX)
    {
    	Logger::gLogger->warning("INVALID_GENOME_INDEX (%u) and record skipped... Reference in BAM is different from the ref used here!", mapPos);
        return false;
    }

    if(qField.empty())
    {
        quality_strings.oldq = samRecord.getQuality();
        //printf("empty\n");
    }
    else
    {
        String temp = samRecord.getString(qField.c_str());
        quality_strings.oldq = temp.c_str();
        //printf("%s\n",samRecord.getQuality());
        //printf("%s:%s\n",qField.c_str(),temp.c_str());
    }

    std::string aligTypes = "";
    Cigar* cigarPtr = samRecord.getCigarInfo();
    if(cigarPtr == NULL)
    {
        Logger::gLogger->warning("Failed to get the cigar");
        return(false);
    }
    cigarPtr->getExpandedString(aligTypes);

    quality_strings.newq = samRecord.getQuality();

    //get tile information
    uint16_t tile;
    if(false){
        std::string readName = samRecord.getReadName();
        std::vector<std::string> readNameFields;
        std::string delimiters = ":_";
        addTokenizedStrings(readName, delimiters, readNameFields);
        // try
        tile = atoi(readNameFields[2].c_str());
    }
    else
    	tile = 0;
    ////////////////
    ////// iterate sequence
    ////////////////
    for (unsigned int i = 0; i < aligTypes.length(); i++)
    {
        switch(aligTypes[i])
        {
            case 'S':
            case 'H':
                offset--;
                continue;
            case 'D':
            case 'N':
                offset2--;
                continue;
            case 'I':
                offset--;
                continue;
        }
        refpos = mapPos + i + offset;
        seqpos = i + offset2;
        
        //exclude dbSNP sites
        if(dbSNP[refpos]==true)
        {
            continue;
        }
        
        refBase = referenceGenome[refpos];
        readBase = samRecord.getSequence(seqpos);
        
        if(reverse ==true)
        {
            cycleIdx = nSeq-seqpos-1;
            refBase = complement(refBase);
            readBase = complement(readBase);
        }
        else
        {
            cycleIdx = seqpos;
        }

        // Get quality string
        // Is other field definied?
        cqual = quality_strings.oldq[seqpos]-33;

        // skip bases with quality below <5
        if(cqual<5)
        {
            continue;
        }

        cobs = nt2idx2[(unsigned int)readBase];
        cref = nt2idx2[(unsigned int)refBase];

        if((cobs==cref) && (cobs<4))
            BMappedCount++;
        else
            BunMappedCount++;

        // read
        if(flag & SamFlag::FIRST_READ)
            cread = 0;
        else
            cread = 1;

        // get cprebase + cnexbase
        // fill in dinucleotide quality
        // previous observed base
        prebase = 'K';
        nexbase = 'K';

        if(cobs < 5)
        {
            if((seqpos != 0) && (aligTypes[i-1] == 'M'))
            {
                // Not the first position, and the cigar indicates the previous
                // position is a match, so there is a prebase.
                // set prebase.
                prebase = samRecord.getSequence(seqpos-1);
            }
            if((seqpos < (nSeq - 1)) && (aligTypes[i+1] == 'M'))
            {
                // Not the last sequence position, and the cigar indicates
                // the previous position is a match, so there is a nexbase.
                // set nexbase.
                nexbase = samRecord.getSequence(seqpos+1);
            }
            if(reverse == true)
            {
                // reverse, so complement the bases.  If they were not set,
                // the default value will be complemented.
                prebase= complement(prebase);
                nexbase = complement(nexbase);
            }
        }

        cprebase = nt2idx2[(unsigned int)prebase];
        cnexbase = nt2idx2[(unsigned int)nexbase];
        //read
        if(processtype==PROCREAD)
        {
            basecv.init(readGroup,cqual,cycleIdx,tile,cread,cprebase,cnexbase,cobs,cobs,cref);
            hasherrormodel.setCell(basecv);
            basecounts++;
        }
        //write
        if(processtype==PROCWRITE)
        {
            // Update quality score
            basecv.init(readGroup,cqual,cycleIdx,tile,cread,cprebase,cnexbase,cobs,0,0);
            uint8_t qemp = hasherrormodel.getQemp(basecv);
            //if(qemp>4)
            quality_strings.newq[seqpos] = qemp+33;
            // otherwise keep the old one
        }
    }
    return true;
}


int ReCab::execute(int argc, char *argv[])
{
    // Shift arguments due to format being ./bam recab and then the args.
    ++argv;
    --argc;

    bool oldQualityFlag = false;
    bool verboseFlag = false;
    bool maxReadFlag = false;
    bool rewriteBam = true;
    int blendedWeight = 0;

    uint64_t skip = 0;
    uint64_t maxi = 0;

    std::string refFile,dbSNPFile,inFile,outFile,logFile,recabFile;

    SamFile samIn,samOut;
    ReCab recab;


    //commandline arguments 
    int arg;
    while((arg = getopt(argc, argv, "i:o:vr:s:l:b:q:")) != -1)
    {
        switch(arg)
        {
            case 'i':
                inFile = optarg;
                break;
            case 'o':
                outFile = optarg;
                break;
            case 'r':
                refFile = optarg;
                break;
            case 's':
                dbSNPFile = optarg;
                break;
            case 'l':
                logFile = optarg;
                break;
            case 'q':
                recab.qField = optarg;
                break;
            case 'b':
                blendedWeight = atoi(optarg);
                break;
            case 'v':
                verboseFlag = true;
                break;
            default:
                return EXIT_FAILURE;
        }
    }

    if (inFile.empty() || outFile.empty() || refFile.empty())
    {
        usage();
        std::cerr << "Missing parameters" << std::endl;
        return EXIT_FAILURE;
    }

    if ( logFile.empty() )
    {
        logFile = outFile + ".log";
    }
    Logger::gLogger = new Logger(logFile.c_str(), verboseFlag);
  
    if ( optind < argc )
    {
        usage();
        Logger::gLogger->error("Argument with no option");
    }

    recab.conversionTable();
    recabFile = outFile + ".recab";

    ////////////////////////////////////////////
    Logger::gLogger->writeLog("Open reference");

    recab.referenceGenome.setReferenceName(refFile.c_str());
    recab.referenceGenome.useMemoryMap();
    if (recab.referenceGenome.open())
    {
        Logger::gLogger->error("Failed to open reference file.");
        return EXIT_FAILURE;
    }

    Logger::gLogger->writeLog("Done! Sequence length %u", recab.referenceGenome.sequenceLength());
    //dbSNP
    bool dbSnpFlag = true;
    if (dbSNPFile == "")
    {
        dbSnpFlag = false;
        Logger::gLogger->writeLog("Ignore dbSNPFile");
    }
    if(recab.referenceGenome.loadDBSNP(recab.dbSNP,dbSNPFile.c_str() ))
    {
    	Logger::gLogger->error("Failed to open dbSNP file.");
    }

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
    bool ok = true;
    ReCab::quality_t quality_strings;

    while ((samIn.ReadRecord(samHeader, samRecord) == true) && (!maxReadFlag || maxi<MAXITERATION))
    {
        //skip records
        if(maxi>=skip)
        {
            uint16_t  flag = samRecord.getFlag();
            if(flag & 0x0400) continue; // duplicates
            if(samRecord.getMapQuality() == 0) continue; // 0 Quality reads
            if(flag & 0x0004)   // unmapped reads
            {
                recab.unMappedCount++;
                continue;
            }
            else
            {
                recab.mappedCount++;
            }

            // quality>0 mapped
            recab.mappedCountQ++;
            ok = recab.processRead(samRecord,PROCREAD, quality_strings);
        }
        //Status info
        if(verboseFlag)
        {
            if(maxi%10000000==0)
                Logger::gLogger->writeLog("%ld records processed", maxi);
        }
        maxi++;
    }

    now = time(0);
    localtm = localtime(&now);
    Logger::gLogger->writeLog("End: %s", asctime(localtm));
    Logger::gLogger->writeLog("# mapped Reads observed: %ld Q>0: %ld",recab.mappedCount,recab.mappedCountQ);
    Logger::gLogger->writeLog("# unmapped Reads observed: %ld",recab.unMappedCount);
    Logger::gLogger->writeLog("# Bases observed: %ld / ZeroQ: %ld - %ld %ld",recab.basecounts,recab.zeroMapQualCount,
                              recab.BMappedCount,recab.BunMappedCount);

    ////////////////////////
    ////////////////////////
    //// Model fitting + prediction
    recab.prediction.setErrorModel(&(recab.hasherrormodel));

    Logger::gLogger->writeLog("Start model fitting!");
    std::string modelfile = outFile + ".model";
    if(recab.prediction.fitModel(true,modelfile))
        recab.prediction.outModel();
    else
        Logger::gLogger->error("Could not fit model!");

    recab.hasherrormodel.addPrediction(recab.prediction.getModel(),blendedWeight);

    Logger::gLogger->writeLog("Writing recalibration table %s",outFile.c_str());
    if(!(recab.hasherrormodel.writeAllTableMM(recabFile.c_str())))
        Logger::gLogger->error("Writing errormodel not possible!");

    Logger::gLogger->writeLog("Writing recalibrated file %s",outFile.c_str());

    ////////////////////////
    ////////////////////////
    //// Write file
    if(rewriteBam)
    {
        samIn.OpenForRead(inFile.c_str());
        samOut.OpenForWrite(outFile.c_str());
        samIn.ReadHeader(samHeader);
        samOut.WriteHeader(samHeader);
        maxi = 0;

        while ((samIn.ReadRecord(samHeader, samRecord) == true) && (!maxReadFlag || maxi<MAXITERATION))
        {
            maxi++;
            ////////////
            uint16_t  flag = samRecord.getFlag();
            if (flag & 0x0400) continue; // duplicates
            if(samRecord.getMapQuality() == 0) continue; // 0 Quality reads
            if(flag & 0x0004) continue;   // unmapped reads
            
            if(recab.processRead(samRecord,PROCWRITE,quality_strings))
            {
                //Logger::gLogger->writeLog("O: %s",quality_strings.oldq.c_str());
                //Logger::gLogger->writeLog("R: %s",quality_strings.newq.c_str());
                //Logger::gLogger->writeLog("---------------------");
                if(oldQualityFlag)
                    samRecord.addTag("QQ",'Z',quality_strings.oldq.c_str());
                samRecord.setQuality(quality_strings.newq.c_str());
                samOut.WriteRecord(samHeader, samRecord);
            }
        }
    }
    Logger::gLogger->writeLog("Recalibration successfully finished");
    return EXIT_SUCCESS;
}
