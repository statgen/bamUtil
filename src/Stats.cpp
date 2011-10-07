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
// This file contains the processing for the executable option "stats"
// which generates some statistics for SAM/BAM files.
#include "Stats.h"
#include "SamFile.h"
#include "BgzfFileType.h"
#include "Pileup.h"
#include "PileupElementBaseQCStats.h"

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
    std::cerr << "\t./bam stats --in <inputFile> [--basic] [--qual] [--phred] [--baseQC <outputFileName>] [--maxNumReads <maxNum>]"
              << "[--unmapped] [--bamIndex <bamIndexFile>] [--regionList <regFileName>] [--minMapQual <minMapQ>] [--dbsnp <dbsnpFile>] [--sumStats] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in : the SAM/BAM file to calculate stats for" << std::endl;
    std::cerr << "\tTypes of Statistics that can be generated:" << std::endl;
    std::cerr << "\t\t--basic       : Turn on basic statistic generation" << std::endl;
    std::cerr << "\t\t--qual        : Generate a count for each quality (displayed as non-phred quality)" << std::endl;
    std::cerr << "\t\t--phred       : Generate a count for each quality (displayed as phred quality)" << std::endl;
    std::cerr << "\t\t--baseQC      : Write per base statistics to the specified file." << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--maxNumReads : Maximum number of reads to process" << std::endl;
    std::cerr << "\t\t                Defaults to -1 to indicate all reads." << std::endl;
    std::cerr << "\t\t--unmapped    : Only process unmapped reads (requires a bamIndex file)" << std::endl;
    std::cerr << "\t\t--bamIndex    : The path/name of the bam index file" << std::endl;
    std::cerr << "\t\t                (if required and not specified, uses the --in value + \".bai\")" << std::endl;
    std::cerr << "\t\t--regionList  : File containing the regions to be processed chr<tab>start_pos<tab>end<pos>." << std::endl;
    std::cerr << "\t\t                Positions are 0 based and the end_pos is not included in the region." << std::endl;
    std::cerr << "\t\t                Uses bamIndex." << std::endl;
    std::cerr << "\t\t--minMapQual  : The minimum mapping quality for filtering reads in the baseQC stats." << std::endl;
    std::cerr << "\t\t--dbsnp       : The dbSnp file of positions to exclude from baseQC analysis." << std::endl;
    std::cerr << "\t\t--noeof       : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params      : Print the parameter settings" << std::endl;
    std::cerr << "\tOptional Base QC Only Parameters:" << std::endl;
    std::cerr << "\t\t--sumStats    : Alternate summary output." << std::endl;
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
    bool sumStats = false;
    int maxNumReads = -1;
    bool unmapped = false;
    String baseQC = "";
    String regionList = "";
    int minMapQual = 0;
    String dbsnp = "";
    PosList *dbsnpListPtr = NULL;


    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_PARAMETER_GROUP("Types of Statistics")
        LONG_PARAMETER("basic", &basic)
        LONG_PARAMETER("qual", &qual)
        LONG_PARAMETER("phred", &phred)
        LONG_STRINGPARAMETER("baseQC", &baseQC)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_INTPARAMETER("maxNumReads", &maxNumReads)
        LONG_PARAMETER("unmapped", &unmapped)
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
        LONG_STRINGPARAMETER("regionList", &regionList)
        LONG_INTPARAMETER("minMapQual", &minMapQual)
        LONG_STRINGPARAMETER("dbsnp", &dbsnp)
        LONG_PARAMETER_GROUP("Optional BaseQC Only Parameters")
        LONG_PARAMETER("sumStats", &sumStats)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    inputParameters.Read(argc-1, &(argv[1]));

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
    Pileup<PileupElementBaseQCStats> pileup;
    // Initialize start/end positions.
    myStartPos = 0;
    myEndPos = -1;
    
    // Open the output qc file if applicable.
    IFILE baseQCPtr = NULL;
    if(!baseQC.IsEmpty())
    {
        baseQCPtr = ifopen(baseQC, "w");
        PileupElementBaseQCStats::setOutputFile(baseQCPtr);
        PileupElementBaseQCStats::setMapQualFilter(minMapQual);
        PileupElementBaseQCStats::setSumStats(sumStats);
        PileupElementBaseQCStats::printHeader();
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
    if((baseQCPtr != NULL) && (!dbsnp.IsEmpty()))
    {
        // Read the dbsnp file.
        IFILE fdbSnp;
        fdbSnp = ifopen(dbsnp,"r");
        // Determine how many entries.
        const SamReferenceInfo* refInfo = samHeader.getReferenceInfo();
        if(refInfo != NULL)
        {
            int maxRefLen = 0;
            for(int i = 0; i < refInfo->getNumEntries(); i++)
            {
                int refLen = refInfo->getReferenceLength(i);
                if(refLen >= maxRefLen)
                {
                    maxRefLen = refLen + 1;
                }
            }

            dbsnpListPtr = new PosList(refInfo->getNumEntries(),maxRefLen);
        }

        if(fdbSnp==NULL)
        {
            std::cerr << "Open dbSNP file " << dbsnp.c_str() << " failed!\n";
        }
        else if(dbsnpListPtr == NULL)
        {
            std::cerr << "Failed to init the memory allocation for the dbsnpList.\n";
        }
        else if(refInfo == NULL)
        {
            std::cerr << "Failed to get the reference information from the bam file.\n";
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
    int qualCount[MAX_QUAL+1];
    for(int i = 0; i <= MAX_QUAL; i++)
    {
        qualCount[i] = 0;
    }
    
    const int MAX_PHRED = 93;
    const int START_PHRED = 0;
    const int PHRED_DIFF = START_QUAL - START_PHRED;
    int phredCount[MAX_PHRED+1];
    for(int i = 0; i <= MAX_PHRED; i++)
    {
        phredCount[i] = 0;
    }
    
    
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
                    while(qual[index] != 0)
                    {
                        // Check for valid quality.
                        if(qual && ((qual[index] < START_QUAL) || (qual[index] > MAX_QUAL)))
                        {
                            std::cerr << "Invalid Quality found: " << qual[index] 
                                      << ".  Must be between "
                                      << START_QUAL << " and " << MAX_QUAL << ".\n";
                            ++index;
                            continue;
                        }
                        if(phred && ((qual[index] < START_PHRED) || (qual[index] > MAX_PHRED)))
                        {
                            std::cerr << "Invalid Quality found: " << qual[index] 
                                      << ".  Must be between "
                                      << START_PHRED << " and " << MAX_PHRED << ".\n";
                            ++index;
                            continue;
                        }
                        
                        // Increment the count for this quality.
                        ++(qualCount[(int)(qual[index])]);
                        ++(phredCount[(int)(qual[index]) - PHRED_DIFF]);
                        ++index;
                    }
                }
            }

            // Check the next thing to do for the read.
            if(baseQCPtr != NULL)
            {
                // Pileup the bases for this read.
                pileup.processAlignmentRegion(samRecord, myStartPos, myEndPos, dbsnpListPtr);
            }
        }

        // Done with a section, move on to the next one.

        // New section, so flush the pileup.
        pileup.flushPileup();
    }

    // Flush the rest of the pileup.
    if(baseQCPtr != NULL)
    {
        // Pileup the bases.
        pileup.processAlignmentRegion(samRecord, myStartPos, myEndPos, dbsnpListPtr);
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
            else if(myStartPos >= myEndPos)
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



