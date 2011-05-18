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
      myCompSeq(false),
      myTags(""),
      myOnlyDiffs(false),
      myMaxAllowedRecs(1000000),
      myAllocatedRecs(0),
      myThreshold(100000),
      myFile1(),
      myFile2(),
      myDiffFileName("-"),
      myDiffFile(NULL),
      myDiff1(""),
      myDiff2(""),
      myTags1(""),
      myTags2("")
{
}


Diff::~Diff()
{
    if(myDiffFile != NULL)
    {
        ifclose(myDiffFile);
    }
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
    std::cerr << "\t./bam diff --in1 <inputFile> --in2 <inputFile> [--out <outputFile>] [--baseQual] [--tags <Tag:Type[;Tag:Type]*>] [--noCigar] [--noPos] [--onlyDiffs] [--recPoolSize <int>] [--posDiff <int>] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in1         : first SAM/BAM file to be diffed" << std::endl;
    std::cerr << "\t\t--in2         : second SAM/BAM file to be diffed" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--seq         : diff the sequence bases." << std::endl;
    std::cerr << "\t\t--baseQual    : diff the base qualities." << std::endl;
    std::cerr << "\t\t--tags        : diff the specified Tags formatted as Tag:Type;Tag:Type;Tag:Type..." << std::endl;
    std::cerr << "\t\t--noCigar     : do not diff the the cigars." << std::endl;
    std::cerr << "\t\t--noPos       : do not diff the positions." << std::endl;
    std::cerr << "\t\t--onlyDiffs   : only print the fields that are different, otherwise for any diff all the fields that are compared are printed." << std::endl;
    std::cerr << "\t\t--recPoolSize : number of records to allow to be stored at a time, default value: " << myMaxAllowedRecs << std::endl;
    std::cerr << "\t\t--posDiff     : max base pair difference between possibly matching records" << myThreshold << std::endl;
    std::cerr << "\t\t--noeof       : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params      : print the parameter settings" << std::endl;
    std::cerr << std::endl;
}


int Diff::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile1 = "";
    String inFile2 = "";
    bool noCigar = false;
    bool noPos = false;
    bool noeof = false;
    bool params = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in1", &inFile1)
        LONG_STRINGPARAMETER("in2", &inFile2)
        LONG_STRINGPARAMETER("out",&myDiffFileName)
        LONG_PARAMETER("seq", &myCompSeq)
        LONG_PARAMETER("baseQual", &myCompBaseQual)
        LONG_STRINGPARAMETER("tags", &myTags)
        LONG_PARAMETER("noCigar", &noCigar)
        LONG_PARAMETER("noPos", &noPos)
        LONG_PARAMETER("onlyDiffs", &myOnlyDiffs)
        LONG_INTPARAMETER("recPoolSize", &myMaxAllowedRecs)
        LONG_INTPARAMETER("posDiff", &myThreshold)
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

    SamRecord* tempRecord = NULL;

    bool need1 = true;
    bool need2 = true;

    SamRecord* rec1 = getSamRecord();
    SamRecord* rec2 = getSamRecord();

    if((rec1 == NULL) || (rec2 == NULL))
    {
        fprintf(stderr, "Failed to allocate initial records, exiting!");
        return(-1);
    }

    // While one of the files still has records.
    while((rec1 != NULL) || (rec2 != NULL))
    {
        if(need1 && (rec1 != NULL))
        {
            // Read from the 1st file.
            if(!myFile1.file.ReadRecord(myFile1.header, *rec1))
            {
                // Failed to read from file1, so at the end of the file.
                releaseSamRecord(rec1);
                rec1 = NULL;

                // Don't need file2's unmatched list anymore since there will
                // be no more records from file 1, so flush list.
                tempRecord = myFile2Unmatched.removeFirst();
                while(tempRecord != NULL)
                {
                    writeDiffs(NULL, tempRecord);
                    releaseSamRecord(tempRecord);
                    tempRecord = myFile2Unmatched.removeFirst();
                }
            }
            need1 = false;
        }
        if(need2 && (rec2 != NULL))
        {
            // Read from the 2nd file.
            if(!myFile2.file.ReadRecord(myFile2.header, *rec2))
            {
                // Failed to read from the 2nd file, so at the end of the file.
                releaseSamRecord(rec2);
                rec2 = NULL;
                // Don't need file1's unmatched list anymore since there will
                // be no more records from file 2, so flush list.
                tempRecord = myFile1Unmatched.removeFirst();
                while(tempRecord != NULL)
                {
                    writeDiffs(tempRecord, NULL);
                    releaseSamRecord(tempRecord);
                    tempRecord = myFile1Unmatched.removeFirst();
                }
            }
            need2 = false;
        }
        if((rec1 == NULL) && (rec2 == NULL))
        {
            // Both records are null, so break.
            break;
        }

        // Prune the unmatched lists for any records that are further from
        // the current record than the threshold.
        tempRecord = myFile2Unmatched.getFirst();
        // Do not need to check for null, because less than returns false if tempRecord is null
        while(lessThan(tempRecord, rec1, myThreshold))
        {
            // this record is more than the specified distance, so remove it
            // and write it.
            tempRecord = myFile2Unmatched.removeFirst();
            writeDiffs(NULL, tempRecord);
            releaseSamRecord(tempRecord);
            tempRecord = myFile2Unmatched.getFirst();
        }
        
        // Prune the unmatched lists for any records that are further from
        // the current record than the threshold.
        tempRecord = myFile1Unmatched.getFirst();
        // Do not need to check for null, because less than returns false if tempRecord is null
        while(lessThan(tempRecord, rec2, myThreshold))
        {
            // this record is more than the specified distance, so remove it
            // and write it.
            tempRecord = myFile1Unmatched.removeFirst();
            writeDiffs(tempRecord, NULL);
            releaseSamRecord(tempRecord);
            tempRecord = myFile1Unmatched.getFirst();
       }

        if(matchingRecs(rec1, rec2))
        {
            // Same fragment and read name.
            writeDiffs(rec1, rec2);
            need1 = true;
            need2 = true;
        }
        else if(lessThan(rec1, rec2))
        {
            // Rec1 is smaller or equal to rec2, so process rec1.
            // Look up rec1 in file2's unmatched list.            
            tempRecord = 
                myFile2Unmatched.removeFragmentMatch(*rec1);

            if((tempRecord == NULL) && (rec2 != NULL))
            {
                // Didn't find the match in the previous file2 records, 
                // and there are more records in file2 that will need
                // to compare against this one, so store this record.
                myFile1Unmatched.addUnmatchedRecord(*rec1);
                
                // Need to get a new record since this one is stored for later use.
                rec1 = getSamRecord();
            }
            else
            {
                // Either a match was found or there was no need to store the record,
                // so write out the diffs.
                writeDiffs(rec1, tempRecord);
                // Release SamRecord ignores a call with NULL.
                releaseSamRecord(tempRecord);
            }
            need1 = true;
            need2 = false;
        }
        else
        {
            // Rec2 is smaller than rec1, so process rec2.
             // Look up rec2 in file1's unmatched list.            
            tempRecord = 
                myFile1Unmatched.removeFragmentMatch(*rec2);

            if((tempRecord == NULL) && (rec1 != NULL))
            {
                // Didn't find the match in the previous file1 records, 
                // and there are more records in file1 that will need
                // to compare against this one, so store this record.
                myFile2Unmatched.addUnmatchedRecord(*rec2);
                
                // Need to get a new record.
                rec2 = getSamRecord();
            }
            else
            {
                // Either a match was found or there was no need to store the record,
                // so write out the diffs.
                writeDiffs(tempRecord, rec2);
                // Release SamRecord ignores a call with NULL.
                releaseSamRecord(tempRecord);
            }
            need1 = false;
            need2 = true;
        }
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

    return(0);
}


bool Diff::matchingRecs(SamRecord* rec1, SamRecord* rec2)
{
    if((rec1 == NULL) || (rec2 == NULL))
    {
        // one or both of the records is NULL, so return false.
        return(false);
    }

    // Have 2 records, so compare them.
    if((SamFlag::getFragmentType(rec1->getFlag()) == 
        SamFlag::getFragmentType(rec2->getFlag())) &&
       (strcmp(rec1->getReadName(),
               rec2->getReadName()) == 0))
    {
        // Same fragment and read name.
        return(true);
    }
    return(false);
}
 
bool Diff::lessThan(SamRecord* rec1, SamRecord* rec2, int threshold)
{
    if(rec1 == NULL)
    {
        return(false);
    }
    if(rec2 == NULL)
    {
        return(true);
    }

    int32_t rec1ChromID = rec1->getReferenceID();
    int32_t rec2ChromID = rec2->getReferenceID();
    // The records are different reads, so check which record has the lower
    // chrom id/position to determine what to return.
    if(rec1ChromID < rec2ChromID)
    {
        // The record from file1 has a lower chromosome id, so return true.
        return(true);
    }
    else if(rec1ChromID > rec2ChromID)
    {
        // The record from file2 has a lower chromosome id, so return false.
        return(false);
    }
    else
    {
        // Same chromosome, so see which record has a lower position.
        int32_t rec1Pos = rec1->get1BasedPosition();
        int32_t rec2Pos = rec2->get1BasedPosition();
        if((rec1Pos + threshold) <= rec2Pos)
        {
            // The record from file1 has a lower or equal position, so return true.
            return(true);
        }
        else
        {
            // The record from file2 has a lower position, so return false.
            return(false);
        }
    }
}


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
    myTags1.Clear();
    myTags2.Clear();

    bool only1 = ((rec2 == NULL) && (rec1 != NULL));
    bool only2 = ((rec1 == NULL) && (rec2 != NULL));
    bool both = ((rec1 != NULL) && (rec2 != NULL));

    bool posDiff = false;
    bool cigarDiff = false;
    bool seqDiff = false;
    bool qualDiff = false;
    bool tagsDiff = false;

    if(!myTags.IsEmpty())
    {
        // Compare the tags.
        if(((rec1 != NULL) && !rec1->getTagsString(myTags.c_str(), myTags1)) ||
           ((rec2 != NULL) && !rec2->getTagsString(myTags.c_str(), myTags2)))
        {
            // Failure reading tags.
            std::cerr << "Failed to search for tags: " << myTags.c_str() << std::endl;
        }
    }
    
    // Check if there are any differences.
    bool bothDiff = false;
    if(both)
    {
        // Have both files, so need to compare.
        // Check if different chromosomes or positions.
        if(myCompPos && 
           ((rec1->getReferenceID() != rec2->getReferenceID()) ||
            (rec1->get1BasedPosition() != rec2->get1BasedPosition())))
        {
            // It moved chromosome/position.
            posDiff = true;
        }
        if(myCompCigar && (strcmp(rec1->getCigar(), rec2->getCigar()) != 0))
        {
            cigarDiff = true;
        }
        if(myCompSeq && (strcmp(rec1->getSequence(), rec2->getSequence()) != 0))
        {
            seqDiff = true;
        }
        if(myCompBaseQual && (strcmp(rec1->getQuality(), rec2->getQuality()) != 0))
        {
            qualDiff = true;
        }
        if(myTags1 != myTags2)
        {
            tagsDiff = true;
        }
        bothDiff = posDiff || cigarDiff || seqDiff || qualDiff || tagsDiff;
    }
    else if(only1 || only2)
    {
        // Only 1 or 2, so set which ones are diff.
        posDiff = myCompPos;
        cigarDiff = myCompCigar;
        seqDiff = myCompSeq;
        qualDiff = myCompBaseQual;
        tagsDiff = !myTags.IsEmpty();
    }
    
    if(only1 || bothDiff)
    {
        // If print all values on any diff or if there is a pos diff
        if(!myOnlyDiffs || posDiff)
        {
            myDiff1 += '\t';
            myDiff1 += rec1->getReferenceName();
            myDiff1 += ':';
            myDiff1 += rec1->get1BasedPosition();
        }
        if(!myOnlyDiffs || cigarDiff)
        {
            // Get the info for rec1's cigar.
             myDiff1 += "\t";
             myDiff1 += rec1->getCigar();
        }
        if(!myOnlyDiffs || seqDiff)
        {
            myDiff1 += "\t";
            myDiff1 += rec1->getSequence();
        }
        if(!myOnlyDiffs || qualDiff)
        {
            myDiff1 += "\t";
            myDiff1 += rec1->getQuality();
        }
        if((!myOnlyDiffs || tagsDiff) && !myTags1.IsEmpty())
        {
            myDiff1 += "\t";
            myDiff1 += myTags1;
        }
    }
    if(only2 || bothDiff)
    {
        // If print all values on any diff or if there is a pos diff
        if(!myOnlyDiffs || posDiff)
        {
            myDiff2 += '\t';
            myDiff2 += rec2->getReferenceName();
            myDiff2 += ':';
            myDiff2 += rec2->get1BasedPosition();
        }
        if(!myOnlyDiffs || cigarDiff)
        {
            // Get the info for rec2's cigar.
             myDiff2 += "\t";
             myDiff2 += rec2->getCigar();
        }
        if(!myOnlyDiffs || seqDiff)
        {
            myDiff2 += "\t";
            myDiff2 += rec2->getSequence();
        }
        if(!myOnlyDiffs || qualDiff)
        {
            myDiff2 += "\t";
            myDiff2 += rec2->getQuality();
        }
        if((!myOnlyDiffs || tagsDiff) && !myTags2.IsEmpty())
        {
            myDiff2 += "\t";
            myDiff2 += myTags2;
        }
    }

    if(!myDiff1.IsEmpty())
    {
        if(checkDiffFile())
        {
            writeStatus &= writeReadName(*rec1);
            myDiff1 += '\n';
            if((ifprintf(myDiffFile, "\n<\t%x", rec1->getFlag()) < 2) || 
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
                writeStatus &= writeReadName(*rec2);
                if(ifwrite(myDiffFile, "\n", 1) != 1)
                {
                    writeStatus = false;
                }
            }
            myDiff2 += '\n';
            if((ifprintf(myDiffFile, ">\t%x", rec2->getFlag()) < 2) || 
               (ifwrite(myDiffFile, myDiff2, myDiff2.Length()) != 
                        (unsigned int)myDiff2.Length()))
            {
                writeStatus = false;
            }
        }
    }
    return(writeStatus);
}


bool Diff::writeReadName(SamRecord& record)
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
    return(true);
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
        // max number allowed to be allocated, so flush the first record from
        // the larger list.
        if(myFile1Unmatched.size() <= myFile2Unmatched.size())
        {
            // remove from file1.
            returnSam = myFile1Unmatched.removeFirst();
            // write out the record as unmatched.
            writeDiffs(returnSam, NULL);
        }
        else
        {
            // remove from file2.
            returnSam = myFile2Unmatched.removeFirst();
            // write out the record as unmatched.
            writeDiffs(NULL, returnSam);
        }
    }
    return(returnSam);
}


void Diff::releaseSamRecord(SamRecord* record)
{
    if(record == NULL)
    {
        // Nothing to release, so just return.
        return;
    }

    // Release the samRecord to be reused.
    myFreeSamRecords.push(record);
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


void Diff::UnmatchedRecords::addUnmatchedRecord(SamRecord& record)
{
    // Lookup which read this is, first, last, or intermediate
    uint16_t flag = record.getFlag();

    UnmatchedRecords::mapType* mapPtr =
        &(myFragmentMaps[SamFlag::getFragmentType(flag)]);

    // Add the record to the unmatched list, and add the position of this record
    // in the list to the map.
    (*mapPtr)[record.getReadName()] = myListUnmatched.insert(myListUnmatched.end(), &record);
}


SamRecord* Diff::UnmatchedRecords::removeFragmentMatch(SamRecord& record)
{
    // Lookup which read this is, first, last, or intermediate
    uint16_t flag = record.getFlag();
    
    UnmatchedRecords::mapType* mapPtr =
        &(myFragmentMaps[SamFlag::getFragmentType(flag)]);
    
    SamRecord* returnRecord = NULL;
    
    // Check to see if it was found.
    myUnmatchedFileIter = mapPtr->find(record.getReadName());
    if(myUnmatchedFileIter != mapPtr->end())
    {
        // Found a match.  Iter->second is an iterator into the list.
        returnRecord = *(myUnmatchedFileIter->second);
        // Remove it from the list and map.
        myListUnmatched.erase(myUnmatchedFileIter->second);
        mapPtr->erase(myUnmatchedFileIter);
    }

    return(returnRecord);
}


SamRecord* Diff::UnmatchedRecords::removeFirst()
{
    if(myListUnmatched.empty())
    {
        // The list is empty, so return NULL.
        return(NULL);
    }

    // Get the first entry from the list.
    SamRecord* returnRec = myListUnmatched.front();
    
    // Remove the record from the map, by first looking up which read this is, 
    // first, last, or intermediate
    uint16_t flag = returnRec->getFlag();

    UnmatchedRecords::mapType* mapPtr =
        &(myFragmentMaps[SamFlag::getFragmentType(flag)]);

    // Remove it from the map.
    mapPtr->erase(returnRec->getReadName());

    // Remove it from the list.
    myListUnmatched.pop_front();
    
    return(returnRec);
}


SamRecord* Diff::UnmatchedRecords::getFirst()
{
    if(myListUnmatched.empty())
    {
        // The list is empty, so return NULL.
        return(NULL);
    }

    // Get the first entry from the list.
    SamRecord* returnRec = myListUnmatched.front();
    
    return(returnRec);
}
