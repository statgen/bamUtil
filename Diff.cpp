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
// This file contains the processing for the executable option "diff"
// which reads an SAM/BAM file and writes a SAM/BAM file with the 
// specified previous values restored if the values are known.

#include "Diff.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"
#include "SamTags.h"
#include "SamFlag.h"

Diff::Diff()
    : myFreeSamRecords(),
      myFile1Unmatched(),
      myFile2Unmatched(),
      myCompCigar(true),
      myCompPos(true),
      myCompBaseQual(false),
      myMaxAllowedRecs(100),
      myAllocatedRecs(0),
      myFile1(),
      myFile2(),
      myDiffFileName("-"),
      myDiffFile(NULL),
      myDiff1(""),
      myDiff2("")
{
}

void Diff::diffDescription()
{
    std::cerr << " diff - Diff 2 SAM/BAM files." << std::endl;
}


void Diff::description()
{
    diffDescription();
}


void Diff::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam diff --in1 <inputFile> --in2 <inputFile> [--cigar] [--qual] [--keepTags] [--rmBQ] [--rmTags <Tag:Type[;Tag:Type]*>] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in1         : first SAM/BAM file to be diffed" << std::endl;
    std::cerr << "\t\t--in2         : second SAM/BAM file to be diffed" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--baseQual    : diff the base qualities." << std::endl;
    std::cerr << "\t\t--tags        : diff the specified Tags formatted as Tag:Type;Tag:Type;Tag:Type..." << std::endl;
    std::cerr << "\t\t--noCigar     : do not diff the the cigars." << std::endl;
    std::cerr << "\t\t--noPos       : do not diff the positions." << std::endl;
    std::cerr << "\t\t--recPoolSize : number of records to allow to be stored at a time" << std::endl;
    std::cerr << "\t\t--noeof       : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params      : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int Diff::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile1 = "";
    String inFile2 = "";
    String tags = "";
    bool noCigar = false;
    bool noPos = false;
    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in1", &inFile1)
        LONG_STRINGPARAMETER("in2", &inFile2)
        LONG_STRINGPARAMETER("out",&myDiffFileName)
        LONG_PARAMETER("baseQual", &myCompBaseQual)
        LONG_STRINGPARAMETER("tags", &tags)
        LONG_PARAMETER("noCigar", &noCigar)
        LONG_PARAMETER("noPos", &noPos)
        LONG_INTPARAMETER("recPoolSize", &myMaxAllowedRecs)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));
    
    myCompCigar = !noCigar;
    myCompPos = !noPos;

    inputParameters.Read(argc-1, &(argv[1]));
    

    // If no eof block is required for a bgzf file, set the bgzf file type to 
    // not look for it.
    if(noeof)
    {
        // Set that the eof block is not required.
        BgzfFileType::setRequireEofBlock(false);
    }
    
    // Check to see if the in file was specified, if not, report an error.
    if(inFile1 == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--in1 is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(inFile2 == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--in2 is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Open the input files for reading.
    myFile1.file.OpenForRead(inFile1);
    myFile2.file.OpenForRead(inFile2);

    // Read the sam headers.
    myFile1.file.ReadHeader(myFile1.header);
    myFile2.file.ReadHeader(myFile2.header);

    myFile1.lastRecord = getSamRecord();
    myFile2.lastRecord = getSamRecord();

    if((myFile1.lastRecord == NULL) || (myFile2.lastRecord == NULL))
    {
        fprintf(stderr, "Failed to allocate initial records, exiting!");
        return(-1);
    }

    while(myFile1.file.ReadRecord(myFile1.header, *myFile1.lastRecord))
    {
        // Read from the 2nd file.
        if(!myFile2.file.ReadRecord(myFile2.header, *myFile2.lastRecord))
        {
            // Failed to read from the 2nd file, so record the first one
            // as a difference.
            writeDiffs(myFile1.lastRecord, NULL);
        }
        else
        {
            // In both files, so compare them.
            if((SamFlag::getFragmentType(myFile1.lastRecord->getFlag()) == 
                SamFlag::getFragmentType(myFile2.lastRecord->getFlag())) &&
               (strcmp(myFile1.lastRecord->getReadName(),
                       myFile2.lastRecord->getReadName()) == 0))
            {
                // Same fragment and read name.
                writeDiffs(myFile1.lastRecord, myFile2.lastRecord);
            }
            else
            {
                // Different fragment and read name.
                writeDiffs(myFile1.lastRecord, NULL);
                writeDiffs(NULL, myFile2.lastRecord);
            }
        }
    }
    // Check for any extra file 2 records.
    while(myFile2.file.ReadRecord(myFile2.header, *myFile2.lastRecord))
    {
        // Not in file1.
        writeDiffs(NULL, myFile2.lastRecord);
    }
    if(myFile1.file.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        // Error.
        return(myFile1.file.GetStatus());
    }
    else if(myFile2.file.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        // Error.
        return(myFile2.file.GetStatus());
    }

    if(myDiffFile != NULL)
    {
        ifclose(myDiffFile);
    }

    return(0);
}
//     int currentRecordChromID = -2;
//     int prevRecordChromID = -2;
//     uint32_t currentPosRec1 = -2;
//     uint32_t currentPosRec2 = -2;

//     // Set returnStatus to success.  It will be changed to the
//     // failure reason if any of the writes or updates fail.
//     SamStatus::Status returnStatus = SamStatus::SUCCESS;

//     // Keep reading records until ReadRecord returns false - no
//     // more records to read.
//     while(samIn1.ReadRecord(samHeader1, myFile1.lastRecord))
//     {
//         bool matchFound = false;

//         // Prune based on positions, writing to the file that
//         // contains mismatches that are just in file2.
//         myFile2Unmatched.pruneUnmatchedByPos(myFile1.lastRecord, mismatchFile2);

//         SamRecord* matchRecord = 
//             myFile2Unmatched.removeFragmentMatch(myFile1.lastRecord);
//         if(matchRecord != NULL)
//         {
//             // Print the diffs for the matching records.
//             printDiffs(record1, returnRecord);
//         }
//         else
//         {
//             // Did not find a match, so read from file 2 until a match
//             // is found or until 

//         }
//         if(!lookupInFile2(myFile1.lastRecord))
//         {
//             // Not yet found in file2, so add it to the file1 list.
//             file1UnMatched.add(myFile1.lastRecord);
//             myFile1.lastRecord = getSamRecord();
//         }
//         // Get the chromosome id of this record.
//         chromIDFile1 = myFile1.lastRecord.getReferenceID();
//         if(chromIDFile1 != chromIDFile2)
//         {
//             flushCurrent();

//         }

        
//         // Get the position of this record.
//         currentPosRec1 = get0BasedPosition();
//         while(currentPosRec1 >= currentPosRec2)
//         {
//             // Read File2.
//             if(!samIn2.ReadRecord(samHeader2, myFile2.lastRecord))
//             {
//                 if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
//                 {
//                     // FAILED to read from file 2.
//                     fprintf(stderr, "%s\n", samIn.GetStatusMessage());
//                     returnStatus = samIn.GetStatus();
//                     // Try reading the next record from file 2.
//                     continue;
//                 }
//                 else
//                 {
//                     // Done reading from file 2, so just flush file 1.
//                     flushFile1();
//                     // Stop reading file 2.
//                     break;
//                 }
//             }
//             // Successfully read a record.
//             // Get the chromsome id & position
//             chromIDFile2 = myFile2.lastRecord.getReferenceID();
//             currentPosRec2 = 
//             if(chromIDFile2 != chromIDFile1)
//             {
//                 flushCurrentChrom();
//             }
//             else
//             {
//                 ...matchFound = true;
//             }
//         }
//         if(!matchFound)
//         {
//             // Didn't find a match for this record read from file 1, so
//             // store it.
//             myFile1Unmatched[samRecord.getReadName()] = myFile1.lastRecord;
//         }
//     }
    
//     // Finished reading file 1, so flush the rest of file 2.
//     flushFile2();


//     return(returnStatus);
// }

// {
//     while(samCurrent->file.ReadRecord(samCurrent->header, currentRecord))
//     {
//         // Read a record from the currently being processed file.
//         // Compare it with the last record from the other file.
//         // Compare the flags and read names
//         if((SamFlag::getFragmentType(lastRecord.getFlag()) == 
//             SamFlag::getFragmentType(currentRecord.getFlag())) &&
//            (samPrev->lastRecord.getReadName() == currentRecord.getReadName()))
//         {
//             // Found a match.
//             processMatch(samPrev->lastRecord, currentRecord);
//         }
//         else
//         {
//             // Does not match.
//             // First check to see if it is in the list for the current file.
//             SamRecord* match = samPrev->unmatchedRecs.removeFragmentMatch(currentRecord);
//             // Check to see if the match was found.
//             if(match != NULL)
//             {
//                 // Found a match.
//                 processMatch();
//             }
//         }
//     }
// }

// void Diff::writeDiff(IFILE file, SamRecord& record)
// {
//     myTempBuffer = '<';
//     if(myCompCigar)
//     {
//         myTempBuffer += "\tCIGAR:";
//         myTempBuffer += record.getCigar();
//     }
//     if(myCompPos)
//     {
//         myTempBuffer += "\tPos:";
//         myTempBuffer += record.get1BasedPosition();
//     }
//     if(myCompBaseQual)
//     {
//         myTempBuffer += "\tBaseQual:";
//         myTempBuffer += record.getQuality();
//     }
//     // TODO Loop through user specified tags.
    
//     if(checkDiffFile())
//     {
//         writeReadNameFrag(file, record);
        
//         ifwrite(file, "<", 1);
//         ifwrite(file, myTempBuffer, myTempBuffer.Length());
//     }
// }


bool Diff::writeDiffs(SamRecord* rec1, SamRecord* rec2)
{
    bool writeStatus = true;
    if((rec1 == NULL) && (rec2 == NULL))
    {
        // No records to compare, so no diffs
        return(false);
    }

    // Check to see if there are any diffs:
    myDiff1.Clear();
    myDiff2.Clear();

    bool only1 = ((rec2 == NULL) && (rec1 != NULL));
    bool only2 = ((rec1 == NULL) && (rec2 != NULL));
    bool both = ((rec1 != NULL) && (rec2 != NULL));

    if(myCompCigar)
    {
        bool cigarDiff = false;
        if(both && (strcmp(rec1->getCigar(), rec2->getCigar()) != 0))
        {
            cigarDiff = true;
        }

        if(only1 || cigarDiff)
        {
            // Get the info for rec1's cigar.
            myDiff1 += "\tCIGAR:";
            myDiff1 += rec1->getCigar();
        }
        if(only2 || cigarDiff)
        {
            // Get the infor for rec2's cigar
            myDiff2 += "\tCIGAR:";
            myDiff2 += rec2->getCigar();
        }
    }
    if(myCompPos)
    {
        bool posDiff = false;
        if(both &&
           (rec1->get1BasedPosition() != rec2->get1BasedPosition()))
        {
            posDiff = true;
        }
        if(only1 || posDiff)
        {
            myDiff1 += "\tPos:";
            myDiff1 += rec1->get1BasedPosition();
        }
        if(only2 || posDiff)
        {
            myDiff2 += "\tPos:";
            myDiff2 += rec2->get1BasedPosition();
        }
    }
    if(myCompBaseQual)
    {
        bool qualDiff = false;
        if(both && (strcmp(rec1->getQuality(), rec2->getQuality()) != 0))
        {
            qualDiff = true;
        }
        if(only1 || qualDiff)
        {
            myDiff1 += "\tBaseQual:";
            myDiff1 += rec1->getQuality();
        }
        if(only2 || qualDiff)
        {
            myDiff2 += "\tBaseQual:";
            myDiff2 += rec2->getQuality();
        }
    }
    // TODO Loop through user specified tags.
    
    if(!myDiff1.IsEmpty())
    {
        if(checkDiffFile())
        {
            writeStatus &= writeReadNameFrag(*rec1);
            myDiff1 += '\n';
            if((ifwrite(myDiffFile, "<", 1) != 1) || 
               (ifwrite(myDiffFile, myDiff1, myDiff1.Length()) != 
                (unsigned int)myDiff1.Length()))
            {
                writeStatus = false;
            }
        }
    }
    if(!myDiff2.IsEmpty())
    {
        if(checkDiffFile())
        {
            if(myDiff1.IsEmpty())
            {
                // Only diff2, so write the readname/fragment.
                writeStatus &= writeReadNameFrag(*rec2);
            }
            myDiff2 += '\n';
            if((ifwrite(myDiffFile, ">", 1) != 1) || 
               (ifwrite(myDiffFile, myDiff2, myDiff2.Length()) != 
                        (unsigned int)myDiff2.Length()))
            {
                writeStatus = false;
            }
        }
    }
    return(writeStatus);
}


bool Diff::writeReadNameFrag(SamRecord& record)
{
    uint8_t nameLen = record.getReadNameLength();
    if(nameLen > 0)
    {
        // it has a read name, so write it.
        // Subtract 1 since the length includes the null.
        --nameLen;
        if(ifwrite(myDiffFile, record.getReadName(), nameLen) != nameLen)
        {
            // Failed to write the entire read name.
            return(false);
        }
    }
    uint16_t flag = record.getFlag();
    if(SamFlag::isFirstFragment(flag))
    {
        myTempBuffer = "\tFirst Fragment\n";
    }
    else if(SamFlag::isLastFragment(flag))
    {
        myTempBuffer = "\tLast Fragment\n";
    }
    else if(SamFlag::isMidFragment(flag))
    {
        myTempBuffer = "\tNot First/Last Fragment\n";
    }
    else
    {
        myTempBuffer = "\tUnknown Fragment\n";
    }
    return(ifwrite(myDiffFile, myTempBuffer.c_str(), myTempBuffer.Length()) == 
           (unsigned int)myTempBuffer.Length());
}


bool Diff::checkDiffFile()
{
    if(myDiffFile == NULL)
    {
        myDiffFile = ifopen(myDiffFileName, "w");

        return(myDiffFile != NULL);
    }
    else
    {
        // File is opened.
        return(true);
    }
}

SamRecord* Diff::getSamRecord()
{
    // Get new samRecord.
    SamRecord* returnSam = NULL;
    if(!myFreeSamRecords.empty())
    {
        // have free already allocated records, so get one of those.
        returnSam = myFreeSamRecords.top();
        myFreeSamRecords.pop();
    }
    else if((myMaxAllowedRecs == -1) || (myAllocatedRecs < myMaxAllowedRecs))
    {
        // There were no free records, but either there is no max or
        // there is still room to allocate more.
        returnSam = new SamRecord();
        ++myAllocatedRecs;
        if(returnSam == NULL)
        {
            // Failed allocation.
            throw(std::runtime_error("Failed to allocate SamRecord"));
        }
    }
    else
    {
        // There are no more free ones and we have already hit the
        // max number allowed to be allocated, so flush the oldest one.
        //TODO  free somehow determine if free from file1 or file2...
    }
    return(returnSam);
}




// SamRecord* Diff::UnmatchedRecords::removeFragmentMatch(SamRecord& record)
// {
//     // Lookup which read this is, first, last, or intermediate
//     uint16_t flag = record.getFlag();

//     UnmatchedRecords::mapType* mapPtr =
//         &(myFragmentMaps[SamFlag::getFragmentType(flag)]);

//     SamRecord* returnRecord = NULL;

//     // Check to see if it was found.
//     myUnmatchedFileIter = mapPtr->find(record.getReadName());
//     if(myUnmatchedFileIter != mapPtr->end())
//     {
//         // Found a match.  Iter->second is an iterator into the list.
//         returnRecord = *(myUnmatchedFileIter->second);
//         // Remove it from the list and map.
//         myListUnmatched.erase(myUnmatchedFileIter->second);
//         mapPtr->erase(myUnmatchedFileIter);
//     }

//     return(returnRecord);
// }

//  so read from file 2 until we get past this
//         // position.
//         // Get the position of this record.
//         currentPosRec1 = get0BasedPosition();
//         while(currentPosRec1 >= currentPosRec2)
//         {
//             // Read File2.
//             if(!samIn2.ReadRecord(samHeader2, samRecord2))
//             {
//                 if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
//                 {
//                     // FAILED to read from file 2.
//                     fprintf(stderr, "%s\n", samIn.GetStatusMessage());
//                     returnStatus = samIn.GetStatus();
//                     // Try reading the next record from file 2.
//                     continue;
//                 }
//                 else
//                 {
//                     // Done reading from file 2, so just flush file 1.
//                     flushFile1();
//                     // Stop reading file 2.
//                     break;
//                 }
//             }
//             // Successfully read a record.
//             // Get the chromsome id & position
//             chromIDFile2 = samRecord2.getReferenceID();
//             currentPosRec2 = 
//             if(chromIDFile2 != chromIDFile1)
//             {
//                 flushCurrentChrom();
//             }
//             else
//             {
//                 // Check the read name.
//                 if(readName Match)
//                 {
//                     ...matchFound = true;
//                 }
//                 else
//                 {
//                     // Store 
//                 }
//             }
//         }
//         if(!matchFound)
//         {
//             // Didn't find a match for this record read from file 1, so
//             // store it.
//             myFile1Unmatched[samRecord.getReadName()] = samRecord1;
//         }
//     }
     
//     }
// }


