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
// This file contains the processing for the executable option "indelDiscordance"
// which generates some statistics for SAM/BAM files.
#include <math.h>
#include "IndelDiscordance.h"
#include "SamFile.h"
#include "BgzfFileType.h"
#include "Pileup.h"
#include "BaseAsciiMap.h"
#include "BaseUtilities.h"

const String IndelDiscordance::DEFAULT_UM_REF_LOC = "/data/local/ref/karma.ref/human.g1k.v37.umfa";
GenomeSequence* IndelDiscordance::PileupElementIndelDiscordance::ourReference = NULL;
int IndelDiscordance::PileupElementIndelDiscordance::ourMinDepth = 2;
bool IndelDiscordance::PileupElementIndelDiscordance::ourPrintPos = false;
uint32_t IndelDiscordance::PileupElementIndelDiscordance::ourTotalMinDepth = 0;
uint32_t IndelDiscordance::PileupElementIndelDiscordance::ourTotalDiscordant = 0;
std::map<uint32_t, uint32_t> IndelDiscordance::PileupElementIndelDiscordance::ourTotalDiscordantRepeats;
std::map<uint32_t, uint32_t> IndelDiscordance::PileupElementIndelDiscordance::ourTotalRepeats;
std::map<uint32_t, uint32_t> IndelDiscordance::PileupElementIndelDiscordance::ourDepthCounts;
std::map<uint32_t, uint32_t> IndelDiscordance::PileupElementIndelDiscordance::ourDepthDiscordantCounts;


IndelDiscordance::IndelDiscordance()
    : BamExecutable()
{
}


void IndelDiscordance::indelDiscordanceDescription()
{
    std::cerr << " indelDiscordance - IndelDiscordance a SAM/BAM File" << std::endl;
}


void IndelDiscordance::description()
{
    indelDiscordanceDescription();
}


void IndelDiscordance::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam indelDiscordance --in <inputFile> [--bamIndex <bamIndexFile] [--refFile <filename>] [--umRef] [--depth minDepth] [--minRepeatLen len] [--sumRepeatLen len] [--printPos] [--chrom <name>] [--start 0basedPos] [--end 0basedPos] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in : the SAM/BAM file to calculate indelDiscordance for" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--bamIndex     : The path/name of the bam index file" << std::endl;
    std::cerr << "\t\t                 (if required and not specified, uses the --in value + \".bai\")" << std::endl;
    std::cerr << "\t\t--refFile      : reference file for determining repeat counts" << std::endl;
    std::cerr << "\t\t--umRef        : use the reference at the default UofM location, "
              << "\t\t                 " << DEFAULT_UM_REF_LOC << std::endl;
    std::cerr << "\t\t--depth        : min depth at which to report indel discordance, DEFAULT >= "
              << DEFAULT_MIN_DEPTH << std::endl;
    std::cerr << "\t\t--minRepeatLen : min repeat length for printing repeat info, DEFAULT = "
              << DEFAULT_MIN_REPEAT << std::endl;
    std::cerr << "\t\t--sumRepeatLen : all repeats this length and longer will be accumulated,\n"
              << "\t\t                 DEFAULT = " << DEFAULT_SUM_REPEAT << std::endl;
    std::cerr << "\t\t--printPos     : print details for each position" << std::endl;
    std::cerr << "\t\t--chrom        : chromosome name other than X" << std::endl;
    std::cerr << "\t\t--start        : use a 0-based inclusive start position other than the default, "
              << DEFAULT_START_POS << std::endl;
    std::cerr << "\t\t--end          : use a 0-based exclusive end position other than the default, "
              << DEFAULT_END_POS << std::endl;
    std::cerr << "\t\t--noeof        : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params       : Print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int IndelDiscordance::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String indexFile = "";
    String refFile = "";
    String chrom = "X";
    bool defaultUMRef = false;
    int minDepth = DEFAULT_MIN_DEPTH;
    int minRepeatLen = DEFAULT_MIN_REPEAT;
    int sumRepeatLen = DEFAULT_SUM_REPEAT;
    bool printPos = false;
    int startPos0Based = DEFAULT_START_POS;
    int endPos0Based = DEFAULT_END_POS;
    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_PARAMETER("umRef", &defaultUMRef)
        LONG_INTPARAMETER("depth", &minDepth)
        LONG_INTPARAMETER("minRepeatLen", &minRepeatLen)
        LONG_INTPARAMETER("sumRepeatLen", &sumRepeatLen)
        LONG_PARAMETER("printPos", &printPos)
        LONG_STRINGPARAMETER("chrom", &chrom)
        LONG_INTPARAMETER("start", &startPos0Based)
        LONG_INTPARAMETER("end", &endPos0Based)
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
        std::cerr << "--in is a mandatory argument for indelDiscordance, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    // Check to see if the IndexFile name has been set.
    if(indexFile == "")
    {
        // Index file was not specified, so set it to the input file + ".bai"
        indexFile = inFile + ".bai";
    }


    if(defaultUMRef)
    {
        refFile = DEFAULT_UM_REF_LOC;
    }

    ////////////////////////////////////////
    // Setup in case pileup is used.
    Pileup<PileupElementIndelDiscordance> pileup;
    
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

    // Read the sam header.
    SamFileHeader samHeader;
    if(!samIn.ReadHeader(samHeader))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    // Open the bam index file for reading and set the read section.
    samIn.ReadBamIndex(indexFile);
    
    samIn.SetReadSection(chrom.c_str(), startPos0Based, endPos0Based);

    GenomeSequence* refPtr = NULL;
    if(refFile.Length() != 0)
    {
        refPtr = new GenomeSequence(refFile.c_str());
        if(refPtr != NULL)
        {
            PileupElementIndelDiscordance::setReference(*refPtr);
        }
    }

    PileupElementIndelDiscordance::setMinDepth(minDepth);
    PileupElementIndelDiscordance::setPrintPos(printPos);


    // Read the sam records.
    SamRecord record;
    
    // Keep reading records from the file until SamFile::ReadRecord
    // indicates to stop (returns false).
    while(samIn.ReadRecord(samHeader, record))
    {
        // Pileup the bases for this read.
        // This expects 0-based positions.
        pileup.processAlignmentRegion(record, startPos0Based, endPos0Based);
    }

    // Flush the rest of the pileup.
    pileup.flushPileup();

    SamStatus::Status status = samIn.GetStatus();
    if(status == SamStatus::NO_MORE_RECS)
    {
        // A status of NO_MORE_RECS means that all reads were successful.
        status = SamStatus::SUCCESS;
    }

    std::cerr << "SUMMARY\n";
    std::cerr << "numDepth >= " << minDepth << ": " 
              << PileupElementIndelDiscordance::ourTotalMinDepth << std::endl;
    std::cerr << "num discordant CIGAR, Depth >= " << minDepth << ": " 
              << PileupElementIndelDiscordance::ourTotalDiscordant << std::endl;
    std::map<uint32_t, uint32_t>::iterator iter;
    for(iter = PileupElementIndelDiscordance::ourTotalDiscordantRepeats.begin(); 
        iter != PileupElementIndelDiscordance::ourTotalDiscordantRepeats.end(); iter++)
    {
        std::cerr << "num discordant CIGAR, Depth >= " << minDepth << " with repeats = " 
                  << (*iter).first << ": "
                  << (*iter).second << std::endl;
    }
    for(iter = PileupElementIndelDiscordance::ourTotalRepeats.begin(); 
        iter != PileupElementIndelDiscordance::ourTotalRepeats.end(); iter++)
    {
        std::cerr << "num Depth >= " << minDepth << " with repeats = " 
                  << (*iter).first << ": "
                  << (*iter).second << std::endl;
    }

    // Calculate the error rate for each depth that has discordance.
    for(iter = PileupElementIndelDiscordance::ourDepthDiscordantCounts.begin(); 
        iter != PileupElementIndelDiscordance::ourDepthDiscordantCounts.end(); iter++)
    {
        // discordant sites with this depth/total sites with this depth
        int n = (*iter).first;
        int count = PileupElementIndelDiscordance::ourDepthCounts[n];
        int numDiscordant = (*iter).second;
        if(count != 0)
        {
            double p_discordant = numDiscordant/(double)count;
            double errorRate = 1. - pow(1. - p_discordant, 1./n);
            std::cerr << "For Depth = " << n << ", total num = " << count
                      << ", num discordant = " << numDiscordant
                      << ", p_discordant = " << p_discordant 
                      << ", and Error Rate = " << errorRate << std::endl;
        }
    }


    delete refPtr;
    return(status);
}


IndelDiscordance::PileupElementIndelDiscordance::PileupElementIndelDiscordance()
    : PileupElement()
{
    initVars();
}


IndelDiscordance::PileupElementIndelDiscordance::~PileupElementIndelDiscordance()
{
}


// Add an entry to this pileup element.  
void IndelDiscordance::PileupElementIndelDiscordance::addEntry(SamRecord& record)
{
    // Call the base class:
    PileupElement::addEntry(record);

    ++myDepth;

    // Analyze the record.
    Cigar* cigar = record.getCigarInfo();
    if(cigar == NULL)
    {
        std::cerr << "Failed to retrieve cigar.\n";
        return;
    }
    if(cigar->size() == 0)
    {
        // Cigar not set, so continue.
        return;
    }

    // Get the query index for this position.
    int queryIndex = 
        cigar->getQueryIndex(getRefPosition(),
                             record.get0BasedPosition());
    if(queryIndex == -1)
    {
        ++myNumDeletion;
    }
    else
    {
        ++myNumMatch;
    }
}


// Perform the alalysis associated with this class.  May be a simple print, 
// a calculation, or something else.  Typically performed when this element
// has been fully populated by all records that cover the reference position.
void IndelDiscordance::PileupElementIndelDiscordance::analyze()
{
    // Check the depth for this position.
    if(myDepth >= ourMinDepth)
    {
        // Update the number of occurrances of this depth.
        ++ourDepthCounts[myDepth];

        ++ourTotalMinDepth;
        int numRepeats = 0;
        if(ourReference != NULL)
        {
            // Add one to ref position since genome position takes 1-based
            uint32_t refPos = 
                ourReference->getGenomePosition(getChromosome(),
                                                getRefPosition()+1);
            uint32_t tempRefPos = refPos - 1;
            char refChar = (*ourReference)[refPos];
            if(!BaseUtilities::isAmbiguous(refChar))
            {
                // Found an actual base, so check for repeats.
                while(tempRefPos > 0)
                {
                    if((*ourReference)[tempRefPos] == refChar)
                    {
                        ++numRepeats;
                    }
                    else
                    {
                        // Not a repeat.
                        break;
                    }
                    --tempRefPos;
                }
                // Loop in the other direction.
                tempRefPos = refPos + 1;
                while(tempRefPos < ourReference->getNumberBases())
                {
                    if((*ourReference)[tempRefPos] == refChar)
                    {
                        ++numRepeats;
                    }
                    else
                    {
                        // Not a repeat.
                        break;
                    }
                    ++tempRefPos;
                }
                if(numRepeats > 0)
                {
                    ++ourTotalRepeats[numRepeats];
                }
            }
        }
        if((myNumDeletion != 0) && (myNumMatch != 0))
        {
            // Discordant, so update counts.
            ++ourTotalDiscordant;
            ++ourDepthDiscordantCounts[myDepth];
            if(ourPrintPos)
            {
                std::cerr << "Position " << getRefPosition() 
                          << " has " << myNumDeletion << " deletions and " 
                          << myNumMatch << " matches.  ";
            }
            if(numRepeats > 0)
            {
                if(ourPrintPos)
                {
                    std:: cerr << numRepeats << " repeats.";
                }

                ++ourTotalDiscordantRepeats[numRepeats];
            }
            if(ourPrintPos)
            {
                std::cerr << "\n";
            }
        }
    }
}


// Resets the entry, setting the new position associated with this element.
void IndelDiscordance::PileupElementIndelDiscordance::reset(int32_t refPosition)
{
    // Call the base class.
    PileupElement::reset(refPosition);

    initVars();
}


void IndelDiscordance::PileupElementIndelDiscordance::initVars()
{
    myDepth = 0;
    myNumDeletion = 0;
    myNumMatch = 0;
    myNumInsertions = 0;
}
