/*
 *  Copyright (C) 2011-2012  Regents of the University of Michigan
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
// This file contains the processing for the executable option "stats"
// which generates some statistics for SAM/BAM files.
#include "Stats.h"
#include "SamFile.h"
#include "BgzfFileType.h"
#include "Pileup.h"
#include "PileupElementBaseQCStats.h"
#include "SamFlag.h"

void Stats::statsDescription()
{
    std::cerr << " stats - Stats a SAM/BAM File" << std::endl;
}


void Stats::description()
{
    statsDescription();
}


void Stats::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam stats --in <inputFile> [--basic] [--qual] [--phred] [--pBaseQC <outputFileName>] [--cBaseQC <outputFileName>] [--maxNumReads <maxNum>]"
              << "[--unmapped] [--bamIndex <bamIndexFile>] [--regionList <regFileName>] [--requiredFlags <integerRequiredFlags>] [--excludeFlags <integerExcludeFlags>] [--noeof] [--params] [--withinRegion] [--baseSum] [--bufferSize <buffSize>] [--minMapQual <minMapQ>] [--dbsnp <dbsnpFile>]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in : the SAM/BAM file to calculate stats for" << std::endl;
    std::cerr << "\tTypes of Statistics that can be generated:" << std::endl;
    std::cerr << "\t\t--basic         : Turn on basic statistic generation" << std::endl;
    std::cerr << "\t\t--qual          : Generate a count for each quality (displayed as non-phred quality)" << std::endl;
    std::cerr << "\t\t--phred         : Generate a count for each quality (displayed as phred quality)" << std::endl;
    std::cerr << "\t\t--pBaseQC       : Write per base statistics as Percentages to the specified file." << std::endl;
    std::cerr << "\t\t                  pBaseQC & cBaseQC cannot both be specified." << std::endl;
    std::cerr << "\t\t--cBaseQC       : Write per base statistics as Counts to the specified file." << std::endl;
    std::cerr << "\t\t                  pBaseQC & cBaseQC cannot both be specified." << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--maxNumReads   : Maximum number of reads to process" << std::endl;
    std::cerr << "\t\t                  Defaults to -1 to indicate all reads." << std::endl;
    std::cerr << "\t\t--unmapped      : Only process unmapped reads (requires a bamIndex file)" << std::endl;
    std::cerr << "\t\t--bamIndex      : The path/name of the bam index file" << std::endl;
    std::cerr << "\t\t                  (if required and not specified, uses the --in value + \".bai\")" << std::endl;
    std::cerr << "\t\t--regionList    : File containing the regions to be processed chr<tab>start_pos<tab>end_pos." << std::endl;
    std::cerr << "\t\t                  Positions are 0 based and the end_pos is not included in the region." << std::endl;
    std::cerr << "\t\t                  Uses bamIndex." << std::endl;
    std::cerr << "\t\t--excludeFlags  : Skip any records with any of the specified flags set\n";
    std::cerr << "\t\t                  (specify an integer representation of the flags)\n";
    std::cerr << "\t\t--requiredFlags : Only process records with all of the specified flags set\n";
    std::cerr << "\t\t                  (specify an integer representation of the flags)\n";
    std::cerr << "\t\t--noeof         : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params        : Print the parameter settings." << std::endl;
    std::cerr << "\tOptional phred/qual Only Parameters:" << std::endl;
    std::cerr << "\t\t--withinRegion  : Only count qualities if they fall within regions specified.\n";
    std::cerr << "\t\t                  Only applicable if regionList is also specified.\n";
    std::cerr << "\tOptional BaseQC Only Parameters:" << std::endl;
    std::cerr << "\t\t--baseSum       : Print an overall summary of the baseQC for the file to stderr." << std::endl;
    std::cerr << "\t\t--bufferSize    : Size of the pileup buffer for calculating the BaseQC parameters." << std::endl;
    std::cerr << "\t\t                  Default: " << PileupHelper::DEFAULT_WINDOW_SIZE << std::endl;
    std::cerr << "\t\t--minMapQual    : The minimum mapping quality for filtering reads in the baseQC stats." << std::endl;
    std::cerr << "\t\t--dbsnp         : The dbSnp file of positions to exclude from baseQC analysis." << std::endl;
    std::cerr << std::endl;
}


int Stats::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String indexFile = "";
    bool basic = false;
    bool noeof = false;
    bool params = false;
    bool qual = false;
    bool phred = false;
    int maxNumReads = -1;
    bool unmapped = false;
    String pBaseQC = "";
    String cBaseQC = "";
    String regionList = "";
    int excludeFlags = 0;
    int requiredFlags = 0;
    bool withinRegion = false;
    int minMapQual = 0;
    String dbsnp = "";
    PosList *dbsnpListPtr = NULL;
    bool baseSum = false;
    int bufferSize = PileupHelper::DEFAULT_WINDOW_SIZE;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_PARAMETER_GROUP("Types of Statistics")
        LONG_PARAMETER("basic", &basic)
        LONG_PARAMETER("qual", &qual)
        LONG_PARAMETER("phred", &phred)
        LONG_STRINGPARAMETER("pBaseQC", &pBaseQC)
        LONG_STRINGPARAMETER("cBaseQC", &cBaseQC)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_INTPARAMETER("maxNumReads", &maxNumReads)
        LONG_PARAMETER("unmapped", &unmapped)
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
        LONG_STRINGPARAMETER("regionList", &regionList)
        LONG_INTPARAMETER("excludeFlags", &excludeFlags)
        LONG_INTPARAMETER("requiredFlags", &requiredFlags)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PARAMETER_GROUP("Optional phred/qual Only Parameters")
        LONG_PARAMETER("withinRegion", &withinRegion)
        LONG_PARAMETER_GROUP("Optional BaseQC Only Parameters")
        LONG_PARAMETER("baseSum", &baseSum)
        LONG_INTPARAMETER("bufferSize", &bufferSize)
        LONG_INTPARAMETER("minMapQual", &minMapQual)
        LONG_STRINGPARAMETER("dbsnp", &dbsnp)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);

    // If no eof block is required for a bgzf file, set the bgzf file type to 
    // not look for it.
    if(noeof)
    {
        // Set that the eof block is not required.
        BgzfFileType::setRequireEofBlock(false);
    }

    // Check to see if the in file was specified, if not, report an error.
    if(inFile == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--in is a mandatory argument for stats, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    // Use the index file if unmapped or regionList is not empty.
    bool useIndex = (unmapped|| (!regionList.IsEmpty()));

    // IndexFile is required, so check to see if it has been set.
    if(useIndex && (indexFile == ""))
    {
        // In file was not specified, so set it to the in file
        // + ".bai"
        indexFile = inFile + ".bai";
    }
    ////////////////////////////////////////
    // Setup in case pileup is used.
    Pileup<PileupElementBaseQCStats> pileup(bufferSize);
    // Initialize start/end positions.
    myStartPos = 0;
    myEndPos = -1;
    
    // Open the output qc file if applicable.
    IFILE baseQCPtr = NULL;
    if(!pBaseQC.IsEmpty() && !cBaseQC.IsEmpty())
    {
        usage();
        inputParameters.Status();
        // Cannot specify both types of baseQC.
        std::cerr << "Cannot specify both --pBaseQC & --cBaseQC." << std::endl;
        return(-1);
    }
    else if(!pBaseQC.IsEmpty())
    {
        baseQCPtr = ifopen(pBaseQC, "w");
        PileupElementBaseQCStats::setPercentStats(true);
    }
    else if(!cBaseQC.IsEmpty())
    {
        baseQCPtr = ifopen(cBaseQC, "w");
        PileupElementBaseQCStats::setPercentStats(false);
    }

    if(baseQCPtr != NULL)
    {
        PileupElementBaseQCStats::setOutputFile(baseQCPtr);
        PileupElementBaseQCStats::printHeader();
    }
    if((baseQCPtr != NULL) || baseSum)
    {
        PileupElementBaseQCStats::setMapQualFilter(minMapQual);
        PileupElementBaseQCStats::setBaseSum(baseSum);
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Open the file for reading.
    SamFile samIn;
    if(!samIn.OpenForRead(inFile))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    samIn.SetReadFlags(requiredFlags, excludeFlags);

    // Set whether or not basic statistics should be generated.
    samIn.GenerateStatistics(basic);

    // Read the sam header.
    SamFileHeader samHeader;
    if(!samIn.ReadHeader(samHeader))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    // Open the bam index file for reading if we are
    // doing unmapped reads (also set the read section).
    if(useIndex)
    {
        samIn.ReadBamIndex(indexFile);

        if(unmapped)
        {
            samIn.SetReadSection(-1);
        }

        if(!regionList.IsEmpty())
        {
            myRegionList = ifopen(regionList, "r");
        }
    }

    //////////////////////////
    // Read dbsnp if specified and doing baseQC
    if(((baseQCPtr != NULL) || baseSum) && (!dbsnp.IsEmpty()))
    {
        // Read the dbsnp file.
        IFILE fdbSnp;
        fdbSnp = ifopen(dbsnp,"r");
        // Determine how many entries.
        const SamReferenceInfo& refInfo = samHeader.getReferenceInfo();
        int maxRefLen = 0;
        for(int i = 0; i < refInfo.getNumEntries(); i++)
        {
            int refLen = refInfo.getReferenceLength(i);
            if(refLen >= maxRefLen)
            {
                maxRefLen = refLen + 1;
            }
        }
        
        dbsnpListPtr = new PosList(refInfo.getNumEntries(),maxRefLen);

        if(fdbSnp==NULL)
        {
            std::cerr << "Open dbSNP file " << dbsnp.c_str() << " failed!\n";
        }
        else if(dbsnpListPtr == NULL)
        {
            std::cerr << "Failed to init the memory allocation for the dbsnpList.\n";
        }
        else
        {
            // Read the dbsnp file.
            StringArray tokens;
            String buffer;
            int position = 0;
            int refID = 0;

            // Loop til the end of the file.
            while (!ifeof(fdbSnp))
            {
                // Read the next line.
                buffer.ReadLine(fdbSnp);
                // If it does not have at least 2 columns, 
                // continue to the next line.
                if (buffer.IsEmpty() || buffer[0] == '#') continue;
                tokens.AddTokens(buffer);
                if(tokens.Length() < 2) continue;

                if(!tokens[1].AsInteger(position))
                {
                    std::cerr << "Improperly formatted region line, start position "
                              << "(2nd column) is not an integer: "
                              << tokens[1]
                              << "; Skipping to the next line.\n";         
                    continue;
                }

                // Look up the reference name.
                refID = samHeader.getReferenceID(tokens[0]);
                if(refID != SamReferenceInfo::NO_REF_ID)
                {
                    // Reference id was found, so add it to the dbsnp
                    dbsnpListPtr->addPosition(refID, position);
                }
        
                tokens.Clear();
                buffer.Clear();
            }
        }
        ifclose(fdbSnp);
    }

    // Read the sam records.
    SamRecord samRecord;

    int numReads = 0;

    //////////////////////
    // Setup in case doing a quality count.
    // Quality histogram.
    const int MAX_QUAL = 126;
    const int START_QUAL = 33;
    uint64_t qualCount[MAX_QUAL+1];
    for(int i = 0; i <= MAX_QUAL; i++)
    {
        qualCount[i] = 0;
    }
    
    const int START_PHRED = 0;
    const int PHRED_DIFF = START_QUAL - START_PHRED;
    const int MAX_PHRED = MAX_QUAL - PHRED_DIFF;
    uint64_t phredCount[MAX_PHRED+1];
    for(int i = 0; i <= MAX_PHRED; i++)
    {
        phredCount[i] = 0;
    }
    
    int refPos = 0;
    Cigar* cigarPtr = NULL;
    char cigarChar = '?';
    // Exclude clips from the qual/phred counts if unmapped reads are excluded.
    bool qualExcludeClips = excludeFlags & SamFlag::UNMAPPED;

    //////////////////////////////////
    // When not reading by sections, getNextSection returns true
    // the first time, then false the next time.
    while(getNextSection(samIn))
    {
        // Keep reading records from the file until SamFile::ReadRecord
        // indicates to stop (returns false).
        while(((maxNumReads < 0) || (numReads < maxNumReads)) && samIn.ReadRecord(samHeader, samRecord))
        {
            // Another record was read, so increment the number of reads.
            ++numReads;
            // See if the quality histogram should be genereated.
            if(qual || phred)
            {
                // Get the quality.
                const char* qual = samRecord.getQuality();
                // Check for no quality ('*').
                if((qual[0] == '*') && (qual[1] == 0))
                {
                    // This record does not have a quality string, so no 
                    // quality processing is necessary.
                }
                else
                {
                    int index = 0;
                    cigarPtr = samRecord.getCigarInfo();
                    cigarChar = '?';
                    refPos = samRecord.get0BasedPosition();
                    if(!qualExcludeClips && (cigarPtr != NULL))
                    {
                        // Offset the reference position by any soft clips
                        // by subtracting the queryIndex of this start position.
                        // refPos is now the start position of the clips.
                        refPos -= cigarPtr->getQueryIndex(0);
                    }

                    while(qual[index] != 0)
                    {
                        // Skip this quality if it is clipped and we are skipping clips.
                        if(cigarPtr != NULL)
                        {
                            cigarChar = cigarPtr->getCigarCharOpFromQueryIndex(index);
                        }
                        if(qualExcludeClips && Cigar::isClip(cigarChar))
                        {
                            // Skip a clipped quality.
                            ++index;
                            // Increment the position.
                            continue;
                        }

                        if(withinRegion && (myEndPos != -1) && (refPos >= myEndPos))
                        {
                            // We have hit the end of the region, stop processing this
                            // quality string.
                            break;
                        }

                        if(withinRegion && (refPos < myStartPos))
                        {
                            // This position is not in the target.
                            ++index;
                            // Update the position if this is found in the reference or a clip.
                            if(Cigar::foundInReference(cigarChar) || Cigar::isClip(cigarChar))
                            {
                                ++refPos;
                            }
                            continue;
                        }

                        // Check for valid quality.
                        if((qual[index] < START_QUAL) || (qual[index] > MAX_QUAL))
                        {
                            if(qual)
                            {
                                std::cerr << "Invalid Quality found: " << qual[index] 
                                          << ".  Must be between "
                                          << START_QUAL << " and " << MAX_QUAL << ".\n";
                            }
                            if(phred)
                            {
                                std::cerr << "Invalid Phred Quality found: " << qual[index] - PHRED_DIFF
                                          << ".  Must be between "
                                          << START_QUAL << " and " << MAX_QUAL << ".\n";
                            }
                            // Skip an invalid quality.
                            ++index;
                            // Update the position if this is found in the reference or a clip.
                            if(Cigar::foundInReference(cigarChar) || Cigar::isClip(cigarChar))
                            {
                                ++refPos;
                            }
                            continue;
                        }
                        
                        // Increment the count for this quality.
                        ++(qualCount[(int)(qual[index])]);
                        ++(phredCount[(int)(qual[index]) - PHRED_DIFF]);
                        // Update the position if this is found in the reference or a clip.
                        if(Cigar::foundInReference(cigarChar) || Cigar::isClip(cigarChar))
                        {
                            ++refPos;
                        }
                        ++index;
                    }
                }
            }

            // Check the next thing to do for the read.
            if((baseQCPtr != NULL) || baseSum)
            {
                // Pileup the bases for this read.
                pileup.processAlignmentRegion(samRecord, myStartPos, myEndPos, dbsnpListPtr);
            }
        }

        // Done with a section, move on to the next one.

        // New section, so flush the pileup.
        pileup.flushPileup();
    }

    if((baseQCPtr != NULL) || baseSum)
    {
        PileupElementBaseQCStats::printSummary();
        ifclose(baseQCPtr);
    }

    std::cerr << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;

    if(basic)
    {
        std::cerr << std::endl;
        samIn.PrintStatistics();
    }

    // Print the quality stats.
    if(qual)
    {
        std::cerr << std::endl;
        std::cerr << "Quality\tCount\n";
        for(int i = START_QUAL; i <= MAX_QUAL; i++)
        {
            std::cerr << i << "\t" << qualCount[i] << std::endl;
        }
    }
    // Print the phred quality stats.
    if(phred)
    {
        std::cerr << std::endl;
        std::cerr << "Phred\tCount\n";
        for(int i = START_PHRED; i <= MAX_PHRED; i++)
        {
            std::cerr << i << "\t" << phredCount[i] << std::endl;
        }
    }

    SamStatus::Status status = samIn.GetStatus();
    if(status == SamStatus::NO_MORE_RECS)
    {
        // A status of NO_MORE_RECS means that all reads were successful.
        status = SamStatus::SUCCESS;
    }

    return(status);
}


bool Stats::getNextSection(SamFile &samIn)
{
    static bool alreadyRead = false;
    if(myRegionList == NULL)
    {
        // no region list is set, so just read once.
        if(alreadyRead)
        {
            // No regions and it has already been read, so
            // return false, no more to read.
            return(false);
        }
        // Return true that there is more to read, but
        // set the flag that it has already been read
        // so the next call will return false.
        alreadyRead = true;
        return(true);
    }
    else
    {
        // There is a region list, so read process that.
        // Track whether or not a section has been found.
        bool sectionFound = false;
        myStartPos = 0;
        myEndPos = 0;

        // Loop until the end of the file or the end of the file or 
        // a section is found.
        while(!sectionFound && !ifeof(myRegionList))
        {
            myRegBuffer.Clear();
            myRegBuffer.ReadLine(myRegionList);
            if(myRegBuffer.IsEmpty())
            {
                // Nothing read, so continue to the next line.
                continue;
            }
        
            // A line was read, so parse it.
            myRegColumn.ReplaceColumns(myRegBuffer, '\t');
            if(myRegColumn.Length() < 3)
            {
                // Incorrectly formatted line.
                std::cerr << "Improperly formatted reg line: "
                          << myRegBuffer
                          << "; Skipping to the next line.\n";
                continue;
            }
            
            // Check the columns.
            if(!myRegColumn[1].AsInteger(myStartPos))
            {
                // The start position (2nd column) is not an integer.
                std::cerr << "Improperly formatted region line, start position "
                          << "(2nd column) is not an integer: "
                          << myRegColumn[1]
                          << "; Skipping to the next line.\n";         
            }
            else if(!myRegColumn[2].AsInteger(myEndPos))
            {
                // The end position (3rd column) is not an integer.
                std::cerr << "Improperly formatted region line, end position "
                          << "(3rd column) is not an integer: "
                          << myRegColumn[2]
                          << "; Skipping to the next line.\n";         
            }
            else if((myStartPos >= myEndPos) && (myEndPos != -1))
            {
                // The start position is >= the end position
                std::cerr << "Improperly formatted region line, the start position "
                          << "is >= end position: "
                          << myRegColumn[1]
                          << " >= "
                          << myRegColumn[2]
                          << "; Skipping to the next line.\n";         
            }
            else
            {
                sectionFound = true;
                samIn.SetReadSection(myRegColumn[0].c_str(), myStartPos, myEndPos);
            }
        }
        return(sectionFound);
    }
}



