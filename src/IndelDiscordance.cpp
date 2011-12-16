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
#include "IndelDiscordance.h"
#include "SamFile.h"
#include "BgzfFileType.h"
#include "Pileup.h"
#include "BaseAsciiMap.h"


GenomeSequence* IndelDiscordance::PileupElementIndelDiscordance::ourReference = NULL;
uint32_t IndelDiscordance::PileupElementIndelDiscordance::numDepth2Plus = 0;
uint32_t IndelDiscordance::PileupElementIndelDiscordance::numDepth3Plus = 0;
uint32_t IndelDiscordance::PileupElementIndelDiscordance::numDiscordant2Plus = 0;
uint32_t IndelDiscordance::PileupElementIndelDiscordance::numDiscordant3Plus = 0;
std::map<uint32_t, uint32_t> IndelDiscordance::PileupElementIndelDiscordance::numDiscordantRepeats2Plus;
std::map<uint32_t, uint32_t> IndelDiscordance::PileupElementIndelDiscordance::numDiscordantRepeats3Plus;
std::map<uint32_t, uint32_t> IndelDiscordance::PileupElementIndelDiscordance::numRepeats2Plus;
std::map<uint32_t, uint32_t> IndelDiscordance::PileupElementIndelDiscordance::numRepeats3Plus;


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
    std::cerr << "\t./bam indelDiscordance --in <inputFile> [--bamIndex <bamIndexFile] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in : the SAM/BAM file to calculate indelDiscordance for" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--bamIndex    : The path/name of the bam index file" << std::endl;
    std::cerr << "\t\t                (if required and not specified, uses the --in value + \".bai\")" << std::endl;
    std::cerr << "\t\t--noeof       : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params      : Print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int IndelDiscordance::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String indexFile = "";
    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
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
    
    samIn.SetReadSection("X", 2699520, 154931044);

    String refFile = "/data/local/ref/karma.ref/human.g1k.v37.umfa";
    GenomeSequence reference(refFile);

    PileupElementIndelDiscordance::setReference(reference);

    // Read the sam records.
    SamRecord record;
    
    // Keep reading records from the file until SamFile::ReadRecord
    // indicates to stop (returns false).
    while(samIn.ReadRecord(samHeader, record))
    {
        // Pileup the bases for this read.
        // This expects 0-based positions.
        pileup.processAlignmentRegion(record, 2699520, 154931044);
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
    std::cerr << "numDepth >= 2: " << PileupElementIndelDiscordance::numDepth2Plus << std::endl;
    std::cerr << "numDepth >  2: " << PileupElementIndelDiscordance::numDepth3Plus << std::endl;
    std::cerr << "num discordant CIGAR, Depth >= 2: " 
              << PileupElementIndelDiscordance::numDiscordant2Plus << std::endl;
    std::cerr << "num discordant CIGAR, Depth >  2: " 
              << PileupElementIndelDiscordance::numDiscordant3Plus << std::endl;
    std::map<uint32_t, uint32_t>::iterator iter;
    for(iter = PileupElementIndelDiscordance::numDiscordantRepeats2Plus.begin(); 
        iter != PileupElementIndelDiscordance::numDiscordantRepeats2Plus.end(); iter++)
    {
        std::cerr << "num discordant CIGAR, Depth >= 2 with repeats = " 
                  << (*iter).first << ": "
                  << (*iter).second << std::endl;
    }
    for(iter = PileupElementIndelDiscordance::numDiscordantRepeats3Plus.begin(); 
        iter != PileupElementIndelDiscordance::numDiscordantRepeats3Plus.end(); iter++)
    {
        std::cerr << "num discordant CIGAR, Depth >  2 with repeats = " 
                  << (*iter).first << ": "
                  << (*iter).second << std::endl;
    }
    for(iter = PileupElementIndelDiscordance::numRepeats2Plus.begin(); 
        iter != PileupElementIndelDiscordance::numRepeats2Plus.end(); iter++)
    {
        std::cerr << "num Depth >= 2 with repeats = " 
                  << (*iter).first << ": "
                  << (*iter).second << std::endl;
    }
    for(iter = PileupElementIndelDiscordance::numRepeats3Plus.begin(); 
        iter != PileupElementIndelDiscordance::numRepeats3Plus.end(); iter++)
    {
        std::cerr << "num Depth >  2 with repeats = " 
                  << (*iter).first << ": "
                  << (*iter).second << std::endl;
    }
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

    ++depth;
    
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
        ++numDeletion;
    }
    else
    {
        ++numMatch;
    }
}


// Perform the alalysis associated with this class.  May be a simple print, 
// a calculation, or something else.  Typically performed when this element
// has been fully populated by all records that cover the reference position.
void IndelDiscordance::PileupElementIndelDiscordance::analyze()
{
    // Check the size of the records.
    if(depth >= 2)
    {
        ++numDepth2Plus;
        if(depth > 2)
        {
            ++numDepth3Plus;
        }
        if(ourReference != NULL)
        {
            // Add one to ref position since genome position takes 1-based
            uint32_t refPos = 
                ourReference->getGenomePosition(getChromosome(),
                                                getRefPosition()+1);
            uint32_t tempRefPos = refPos - 1;
            char refChar = (*ourReference)[refPos];
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
                ++numRepeats2Plus[numRepeats];

                if(depth > 2)
                {
                    ++numRepeats3Plus[numRepeats];
                }
            }
        }
        
        if((numDeletion != 0) && (numMatch != 0))
        {
            ++numDiscordant2Plus;
            if(depth > 2)
            {
                ++numDiscordant3Plus;
            }
            std::cerr << "Position " << getRefPosition() 
                      << " has " << numDeletion << " deletions and " 
                      << numMatch << " matches.  ";

            if(numRepeats > 0)
            {
                std:: cerr << numRepeats << " repeats.";

                ++numDiscordantRepeats2Plus[numRepeats];

                if(depth > 2)
                {
                    ++numDiscordantRepeats3Plus[numRepeats];
                }
            }
            std::cerr << "\n";
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
    depth = 0;
    numDeletion = 0;
    numMatch = 0;
    numInsertions = 0;
    numRepeats = 0;
}
