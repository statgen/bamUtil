/*
 *  Copyright (C) 2010-2014  Regents of the University of Michigan
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
// This file contains the processing for the executable option "writeRegion"
// which writes a file with the reads in the specified region.

#include "WriteRegion.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <unordered_set>
typedef std::unordered_set<std::string> ReadNameSet;
#else
#include <set>
typedef std::set<std::string> ReadNameSet;
#endif

WriteRegion::WriteRegion()
    : BamExecutable(),
      myWithinReg(false),
      myWroteReg(false),
      myStart(UNSPECIFIED_INT),
      myEnd(UNSPECIFIED_INT),
      myPrevStart(UNSPECIFIED_INT),
      myPrevEnd(UNSPECIFIED_INT),
      myRefID(UNSPECIFIED_INT),
      myRefName(),
      myPrevRefName(),
      myBedRefID(SamReferenceInfo::NO_REF_ID),
      myBedFile(NULL)
{
    
}

void WriteRegion::writeRegionDescription()
{
    std::cerr << " writeRegion - Write a file with reads in the specified region and/or have the specified read name" << std::endl;
}

void WriteRegion::description()
{
    writeRegionDescription();
}

void WriteRegion::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam writeRegion --in <inputFilename>  --out <outputFilename> [--bamIndex <bamIndexFile>] "
              << "[--refName <reference Name> | --refID <reference ID>] [--start <0-based start pos>] "
              << "[--end <0-based end psoition>] [--bed <bed filename>] [--withinRegion] [--readName <readName>] [--rnFile <readNameFileName>] "
              << "[--lshift] [--params] [--noeof]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in        : the BAM file to be read" << std::endl;
    std::cerr << "\t\t--out       : the SAM/BAM file to write to" << std::endl;
    std::cerr << "\tOptional Parameters for Specifying a Region:" << std::endl;
    std::cerr << "\t\t--bamIndex  : the path/name of the bam index file" << std::endl;
    std::cerr << "\t\t              (if not specified, uses the --in value + \".bai\")" << std::endl;
    std::cerr << "\t\t--refName   : the BAM reference Name to read" << std::endl;
    std::cerr << "\t\t              Either this or refID can be specified." << std::endl;
    std::cerr << "\t\t              Defaults to all references." << std::endl;
    std::cerr << "\t\t--refID     : the BAM reference ID to read." << std::endl;
    std::cerr << "\t\t              Either this or refName can be specified." << std::endl;
    std::cerr << "\t\t              Defaults to all references." << std::endl;
    std::cerr << "\t\t              Specify -1 for unmapped" << std::endl;
    std::cerr << "\t\t--start     : inclusive 0-based start position." << std::endl;
    std::cerr << "\t\t              Defaults to -1: meaning from the start of the reference." << std::endl;
    std::cerr << "\t\t              Only applicable if refName/refID is set." << std::endl;
    std::cerr << "\t\t--end       : exclusive 0-based end position." << std::endl;
    std::cerr << "\t\t              Defaults to -1: meaning til the end of the reference." << std::endl;
    std::cerr << "\t\t              Only applicable if refName/refID is set." << std::endl;
    std::cerr << "\t\t--bed       : use the specified bed file for regions." << std::endl;
    std::cerr << "\t\t--withinReg : only print reads fully enclosed within the region." << std::endl;
    std::cerr << "\t\t--readName  : only print reads with this read name." << std::endl;
    std::cerr << "\t\t--rnFile    : only print reads with read names found in the specified file,\n";
    std::cerr << "\t\t              delimited by comma, space, tab, or new line (',', ' ', '\\t', or '\\n')." << std::endl;
    std::cerr << "\tOptional Parameters For Other Operations:\n";
    std::cerr << "\t\t--lshift        : left shift indels when writing records\n";
    std::cerr << "\t\t--excludeFlags  : Skip any records with any of the specified flags set\n";
    std::cerr << "\t\t                  (specify an integer representation of the flags)\n";
    std::cerr << "\t\t--requiredFlags : Only process records with all of the specified flags set\n";
    std::cerr << "\t\t                  (specify an integer representation of the flags)\n";
    std::cerr << "\t\t--params        : print the parameter settings" << std::endl;
    std::cerr << "\t\t--noeof         : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << std::endl;
}


int WriteRegion::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    String indexFile = "";
    String readName = "";
    String rnFile = "";
    String bed = "";
    ReadNameSet rnSet;
    myStart = UNSPECIFIED_INT;
    myEnd = UNSPECIFIED_INT;
    myPrevStart = UNSPECIFIED_INT;
    myPrevEnd = UNSPECIFIED_INT;
    myRefID = UNSET_REF;
    myRefName.Clear();
    myPrevRefName.Clear();
    myBedRefID = SamReferenceInfo::NO_REF_ID;
    bool lshift = false;
    bool noeof = false;
    bool params = false;
    String excludeFlags = "";
    String requiredFlags = "";
    myWithinReg = false;
    myWroteReg = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_PARAMETER_GROUP("Optional Region Parameters")        
        LONG_STRINGPARAMETER("bamIndex", &indexFile)
        LONG_STRINGPARAMETER("refName", &myRefName)
        LONG_INTPARAMETER("refID", &myRefID)
        LONG_INTPARAMETER("start", &myStart)
        LONG_INTPARAMETER("end", &myEnd)
        LONG_STRINGPARAMETER("bed", &bed)
        LONG_PARAMETER("withinReg", &myWithinReg)
        LONG_STRINGPARAMETER("readName", &readName)
        LONG_STRINGPARAMETER("rnFile", &rnFile)
        LONG_PARAMETER_GROUP("Optional Other Parameters")
        LONG_PARAMETER("lshift", &lshift)
        LONG_STRINGPARAMETER("excludeFlags", &excludeFlags)
        LONG_STRINGPARAMETER("requiredFlags", &requiredFlags)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    // bam writeRegion, parameters start at index 2 rather than 1.
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
        // mandatory argument was not specified.
        inputParameters.Status();
        std::cerr << "Missing mandatory argument: --in" << std::endl;
        return(-1);
    }
    if(outFile == "")
    {
        usage();
        // mandatory argument was not specified.
        inputParameters.Status();
        std::cerr << "Missing mandatory argument: --out" << std::endl;
        return(-1);
    }
    
    if(indexFile == "")
    {
        // In file was not specified, so set it to the in file
        // + ".bai"
        indexFile = inFile + ".bai";
    }

    if(myRefID != UNSET_REF && myRefName.Length() != 0)
    {
        std::cerr << "Can't specify both refID and refName" << std::endl;
        inputParameters.Status();
        return(-1);
    }
    if(myRefID != UNSET_REF && bed.Length() != 0)
    {
        std::cerr << "Can't specify both refID and bed" << std::endl;
        inputParameters.Status();
        return(-1);
    }
    if(myRefName.Length() != 0 && bed.Length() != 0)
    {
        std::cerr << "Can't specify both refName and bed" << std::endl;
        inputParameters.Status();
        return(-1);
    }

    if(!bed.IsEmpty())
    {
        myBedFile = ifopen(bed, "r");
    }

    if(!rnFile.IsEmpty())
    {
        // Open the read name file.
        IFILE rnFileRdr = ifopen(rnFile, "r");
        if(rnFileRdr == NULL)
        {
            std::cerr << "ERROR: Could not open read name file: " 
                      << rnFile << std::endl;
            inputParameters.Status();
            return(-1);
        }
        // Read the read names into a map.
        std::string rn;
        while(rnFileRdr->readTilChar(",\n\t ", rn) != -1)
        {
            rnSet.insert(rn);
            rn.clear();
        }
        // Check if rn has a value so it works if there is no final new line.
        if(!rn.empty())
        {
            rnSet.insert(rn);
            rn.clear();
        }
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Open the file for reading.   
    mySamIn.OpenForRead(inFile);

    mySamIn.SetReadFlags(requiredFlags.AsInteger(), excludeFlags.AsInteger());

    // Open the output file for writing.
    SamFile samOut;
    samOut.OpenForWrite(outFile);

    // Open the bam index file for reading if a region was specified.
    if((myRefName.Length() != 0) || (myRefID != UNSET_REF) || (myBedFile != NULL))
    {
        mySamIn.ReadBamIndex(indexFile);
    }

    // Read & write the sam header.
    mySamIn.ReadHeader(mySamHeader);
    samOut.WriteHeader(mySamHeader);

    // Read the sam records.
    SamRecord samRecord;
    // Track the status.
    int numSectionRecords = 0;

    // Set returnStatus to success.  It will be changed
    // to the failure reason if any of the writes fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;
        
    while(getNextSection())
    {
        // Keep reading records until they aren't anymore.
        while(mySamIn.ReadRecord(mySamHeader, samRecord))
        {
            if(!readName.IsEmpty())
            {
                // Check for readname.
                if(strcmp(samRecord.getReadName(), readName.c_str()) != 0)
                {
                    // not a matching read name, so continue to the next record.
                    continue;
                }
            }
            if(!rnSet.empty())
            {
                if(rnSet.count(samRecord.getReadName()) == 0)
                {
                    // Not in the read name set, so continue to the next record.
                    continue;
                }
            }

            // Check to see if the read has already been processed.
            if(myPrevEnd != UNSPECIFIED_INT)
            {
                // Because we already know that the bed was sorted, 
                // we know that the previous section started before
                // this one, so if the previous end is greater than
                // this record's start position we know that it
                // was already written in the previous section.
                // Note: can't be equal to the previous end since
                // the end range was exclusive, while
                // get0BasedPosition is inclusive.
                // myPrevEnd is reset by getNextSection when a new
                // chromosome is hit.
                if(samRecord.get0BasedPosition() < myPrevEnd)
                {
                    // This record was already written.
                    continue;
                }
            }

            // Shift left if applicable.
            if(lshift)
            {
                samRecord.shiftIndelsLeft();
            }

            // Successfully read a record from the file, so write it.
            samOut.WriteRecord(mySamHeader, samRecord);
            ++numSectionRecords;
        }
        myWroteReg = true;
    }

    if(myBedFile != NULL)
    {
        ifclose(myBedFile);
    }
    std::cerr << "Wrote " << outFile << " with " << numSectionRecords
              << " records.\n";
    return(returnStatus);
}


bool WriteRegion::getNextSection()
{
    bool anotherSection = false;
    // If refName is set, use that.
    if(myRefName.Length() != 0)
    {
        // Use Reference Name for the next section.
        anotherSection = true;
        mySamIn.SetReadSection(myRefName.c_str(), myStart, myEnd, !myWithinReg);
        // Already processed this section, so clear the reference name
        // so it will not be used again.
        myRefName.Clear();
        myStart = UNSPECIFIED_INT;
        myEnd = UNSPECIFIED_INT;
    }
    else if(myRefID != UNSET_REF)
    {
        // Use Reference ID for the next section.
        anotherSection = true;
        mySamIn.SetReadSection(myRefID, myStart, myEnd, !myWithinReg);
        // Already processed this section, so clear the reference id
        // so it will not be used again.
        myRefID = UNSET_REF;
        myStart = UNSPECIFIED_INT;
        myEnd = UNSPECIFIED_INT;
    }
    else if(myBedFile != NULL)
    {
        // There is a bed file, so read the next line.
        while(!anotherSection)
        {
            myBedBuffer.Clear();
            myBedBuffer.ReadLine(myBedFile);
            if(ifeof(myBedFile) && myBedBuffer.IsEmpty())
            {
                // End of the file, so break.
                break;
            }
            // Not the end of the file, so parse the line.
            myBedColumn.ReplaceColumns(myBedBuffer, '\t');
            if(myBedColumn.Length() != 3)
            {
                // Incorrectly formatted line.
                std::cerr << "Improperly formatted bed line: "
                          << myBedBuffer
                          << "; Skipping to the next line.\n";
            }
            else
            {
                // Check the reference name.
                if(myPrevRefName != myBedColumn[0])
                {
                    // New reference name (chromosome), so clear the previous
                    // start/end.
                    myPrevStart = UNSPECIFIED_INT;
                    myPrevEnd = UNSPECIFIED_INT;
                    myPrevRefName = myBedColumn[0];

                    // Get the reference ID for the reference name.
                    myBedRefID = mySamHeader.getReferenceID(myPrevRefName);
                    
                    // Check to see if the reference ID is found.
                    if(myBedRefID == SamReferenceInfo::NO_REF_ID)
                    {
                        // The specified Reference ID is not in the file,
                        // so check to see if it has chr.
                        // Check to see if it is the same except for 'chr' appended.
                        if((myPrevRefName[0] == 'c') && 
                           (myPrevRefName[1] == 'h') && 
                           (myPrevRefName[2] == 'r'))
                        {
                            // It starts with chr, so look up with out the chr
                            myBedRefID = mySamHeader.getReferenceID(myPrevRefName.c_str() + 3);
                        }
                    }
                }
                else
                {
                    // Not a new reference name.
                    // Store the previous positions before overwriting them.
                    myPrevStart = myStart;
                    if(myPrevEnd < myEnd)
                    {
                        // The last section ends later than the previous one,
                        // So update the previous latest end.
                        myPrevEnd = myEnd;
                    }
                }

                // If the refID is still NO_REF_ID, just continue to the next bed line.
                if(myBedRefID == SamReferenceInfo::NO_REF_ID)
                {
                    continue;
                }

                // Correct number of columns, check the columns.
                if(!myBedColumn[1].AsInteger(myStart))
                {
                    // The start position (2nd column) is not an integer.
                    std::cerr << "Improperly formatted bed line, start position (2nd column) is not an integer: "
                              << myBedColumn[1]
                              << "; Skipping to the next line.\n";         
                }
                else if(!myBedColumn[2].AsInteger(myEnd))
                {
                    // The end position (3rd column) is not an integer.
                    std::cerr << "Improperly formatted bed line, end position (3rd column) is not an integer: "
                              << myBedColumn[2]
                              << "; Skipping to the next line.\n";         
                }
                else if(myStart >= myEnd)
                {
                    // The start position is >= the end
                    std::cerr << "Improperly formatted bed line, the start position is >= end position: "
                              << myBedColumn[1]
                              << " >= "
                              << myBedColumn[2]
                              << "; Skipping to the next line.\n";         
                }
                else if(myPrevStart > myStart)
                {
                    // Same reference name, but the position goes backwards.
                    // This is against the assumption that the bed is sorted.
                    std::cerr << "Improperly formatted bed, the start position is < the previous start (bed is assumed to be sorted): "
                              << myStart
                              << " < "
                              << myPrevStart
                              << "; Skipping to the next line.\n";
                }
                else
                {
                    anotherSection = true;
                    mySamIn.SetReadSection(myBedRefID, myStart, myEnd, !myWithinReg);
                }
            }
        }
    }
    else
    {
        // If we have no bed file, then we only have another section
        // if we have not already written a region.
        anotherSection = !myWroteReg;
    }
    
    return(anotherSection);
}
