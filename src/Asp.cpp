/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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
// This file contains the processing for the executable option "asp"
// which generates some statistics for SAM/BAM files.

#include "Asp.h"
#include "SamFile.h"
#include "BgzfFileType.h"
#include "Pileup.h"
#include "PileupElementAsp.h"
#include "AspReadNameID.h"

Asp::Asp()
    : BamExecutable(),
      myRegionList(NULL),
      myRegColumn()
{
    reset();
}


Asp::~Asp()
{
}


void Asp::aspDescription()
{
    std::cerr << " asp - Pileup on a SAM/BAM File" << std::endl;
}


void Asp::description()
{
    aspDescription();
}


void Asp::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam asp --in <inputFile> --out <outputFile> --refFile <referenceFilename> [--bamIndex <bamIndexFile>] [--regionList <regFileName>] [--gapSize <gapSize>] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in       : the SAM/BAM file to calculate asp for" << std::endl;
    std::cerr << "\t\t--out      : the output file to write" << std::endl;
    std::cerr << "\t\t--refFile  : the reference file" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--bamIndex    : The path/name of the bam index file" << std::endl;
    std::cerr << "\t\t                (if required and not specified uses the --in value + \".bai\")" << std::endl;
    std::cerr << "\t\t--regionList  : File containing the regions to be processed chr<tab>start_pos<tab>end<pos>." << std::endl;
    std::cerr << "\t\t                Positions are 0 based and the end_pos is not included in the region." << std::endl;
    std::cerr << "\t\t                Uses bamIndex." << std::endl;
    std::cerr << "\t\t--gapSize     : Gap Size threshold such that position gaps less than this size have an" << std::endl;
    std::cerr << "\t\t                empty record written, while gaps larger than this size have a new" << std::endl;
    std::cerr << "\t\t                chrom/position header written, Default = " << AspFileWriter::DEFAULT_GAP_SIZE << "." << std::endl;
    std::cerr << "\t\t--noeof       : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params      : Print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int Asp::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    String refFile = "";
    String indexFile = "";
    bool noeof = false;
    bool params = false;
    String regionList = "";
    int gapSize = AspFileWriter::DEFAULT_GAP_SIZE;

    reset();

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
        LONG_STRINGPARAMETER("regionList", &regionList)
        LONG_INTPARAMETER("gapSize", &gapSize)
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
        std::cerr << "--in is a mandatory argument for asp, "
                  << "but was not specified" << std::endl;
        return(-1);
    }
    if(outFile == "")
    {
        usage();
        inputParameters.Status();
        // Out file was not specified but it is mandatory.
        std::cerr << "--out is a mandatory argument for asp, "
                  << "but was not specified" << std::endl;
        return(-1);
    }
    if(refFile == "")
    {
        usage();
        inputParameters.Status();
        // Ref file was not specified but it is mandatory.
        std::cerr << "--refFile is a mandatory argument for asp, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    // Use the index file if unmapped or regionList is not empty.
    bool useIndex = !regionList.IsEmpty();

    // IndexFile is required, so check to see if it has been set.
    if(useIndex && (indexFile == ""))
    {
        // In file was not specified, so set it to the in file
        // + ".bai"
        indexFile = inFile + ".bai";
    }

    GenomeSequence reference(refFile);

    ////////////////////////////////////////
    // Setup pileup.
    Pileup<PileupElementAsp> pileup;
    
    AspFileWriter::setGapSize(gapSize);
    PileupElement::setReference(&reference);

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

    AspFileWriter aspOutputFile;
    aspOutputFile.open(outFile, samHeader);

    PileupElementAsp::setOutputFile(&aspOutputFile);

    // Open the bam index file for reading if we are
    // doing unmapped reads (also set the read section).
    if(useIndex)
    {
        samIn.ReadBamIndex(indexFile);

        if(!regionList.IsEmpty())
        {
            myRegionList = ifopen(regionList, "r");
        }
    }

    // Read the sam records.
    SamRecord samRecord;
    AspReadNameID myReadNameMap;
    while(getNextSection(samIn))
    {
        // Keep reading records from the file until SamFile::ReadRecord
        // indicates to stop (returns false).
        while(samIn.ReadRecord(samHeader, samRecord))
        {
            // Determine the read name id for this record.
            PileupElementAsp::setCurrentReadNameID(myReadNameMap.getReadNameID(samRecord.getReadName()));
            // Pileup the bases for this read.
            pileup.processAlignmentRegion(samRecord, myStartPos, myEndPos);
        }

        // Done with a section, move on to the next one.
        // New section, so flush the pileup.
        pileup.flushPileup();
    }

    // Flush the rest of the pileup.
    pileup.flushPileup();

    std::cerr << "Number of Position Records = "
              << aspOutputFile.getNumPosRecs() << "\n"
              << "Number of Empty Records = "
              << aspOutputFile.getNumEmptyRecs() << "\n"
              << "Number of Reference Only Records = "
              << aspOutputFile.getNumRefOnlyRecs() << "\n"
              << "Number of Detailed Records = "
              << aspOutputFile.getNumDetailedRecs() << "\n";

    SamStatus::Status status = samIn.GetStatus();
    if(status == SamStatus::NO_MORE_RECS)
    {
        // A status of NO_MORE_RECS means that all reads were successful.
        status = SamStatus::SUCCESS;
    }
    return(status);
}


void Asp::reset()
{
    myStartPos = 0;
    myEndPos = -1;
    myRegBuffer.Clear();
    myRegColumn.Clear();
    if(myRegionList != NULL)
    {
        ifclose(myRegionList);
        myRegionList = NULL;
    }
}


bool Asp::getNextSection(SamFile &samIn)
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


