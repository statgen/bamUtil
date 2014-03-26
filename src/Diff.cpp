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
#include "SamRecordHelper.h"

const char* Diff::FLAG_DIFF_TAG = "ZF";
const char* Diff::POS_DIFF_TAG = "ZP";
const char* Diff::CIGAR_DIFF_TAG = "ZC";
const char* Diff::MAPQ_DIFF_TAG = "ZM";
const char* Diff::MATE_DIFF_TAG = "ZN";
const char* Diff::ISIZE_DIFF_TAG = "ZI";
const char* Diff::SEQ_DIFF_TAG = "ZS";
const char* Diff::QUAL_DIFF_TAG = "ZQ";
const char* Diff::TAGS_DIFF_TAG = "ZT";

Diff::Diff()
    : myFreeSamRecords(),
      myFile1Unmatched(),
      myFile2Unmatched(),
      myCompAll(false),
      myCompCigar(true),
      myCompPos(true),
      myCompBaseQual(false),
      myCompSeq(false),
      myCompFlag(false),
      myCompMapQ(false),
      myCompMate(false),
      myCompISize(false),
      myTags(""),
      myOnlyDiffs(false),
      myBamOut(false),
      myMaxAllowedRecs(1000000),
      myAllocatedRecs(0),
      myThreshold(100000),
      myNumPoolOverflows(0),
      myFile1(),
      myFile2(),
      myDiffFileName("-"),
      myBamOnly1Name(""),
      myBamOnly2Name(""),
      myBamDiffName(""),
      myDiffFile(NULL),
      myBamOnly1(),
      myBamOnly2(),
      myBamDiff(),
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
    std::cerr << " diff - Diff 2 coordinate sorted SAM/BAM files." << std::endl;
}


void Diff::description()
{
    diffDescription();
}


void Diff::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam diff --in1 <inputFile> --in2 <inputFile> [--out <outputFile>] [--all] [--flag] [--mapQual] [--mate] [--isize] [--seq] [--baseQual] [--tags <Tag:Type[,Tag:Type]*>] [--everyTag] [--noCigar] [--noPos] [--onlyDiffs] [--recPoolSize <int>] [--posDiff <int>] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in1         : first coordinate sorted SAM/BAM file to be diffed" << std::endl;
    std::cerr << "\t\t--in2         : second coordinate sorted SAM/BAM file to be diffed" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--out         : output filename, use .bam extension to output in SAM/BAM format instead of diff format." << std::endl;
    std::cerr << "\t\t                In SAM/BAM format there will be 3 output files:" << std::endl;
    std::cerr << "\t\t                    1) the specified name with record diffs" << std::endl;
    std::cerr << "\t\t                    2) specified name with _only_<in1>.sam/bam with records only in the in1 file" << std::endl;
    std::cerr << "\t\t                    3) specified name with _only_<in2>.sam/bam with records only in the in2 file" << std::endl;
    std::cerr << "\t\t--all         : diff all the SAM/BAM fields." << std::endl;
    std::cerr << "\t\t--flag        : diff the flags." << std::endl;
    std::cerr << "\t\t--mapQual     : diff the mapping qualities." << std::endl;
    std::cerr << "\t\t--mate        : diff the mate chrom/pos." << std::endl;
    std::cerr << "\t\t--isize       : diff the insert sizes." << std::endl;
    std::cerr << "\t\t--seq         : diff the sequence bases." << std::endl;
    std::cerr << "\t\t--baseQual    : diff the base qualities." << std::endl;
    std::cerr << "\t\t--tags        : diff the specified Tags formatted as Tag:Type,Tag:Type,Tag:Type..." << std::endl;
    std::cerr << "\t\t--everyTag    : diff all the Tags" << std::endl;
    std::cerr << "\t\t--noCigar     : do not diff the the cigars." << std::endl;
    std::cerr << "\t\t--noPos       : do not diff the positions." << std::endl;
    std::cerr << "\t\t--onlyDiffs   : only print the fields that are different, otherwise for any diff all the fields that are compared are printed." << std::endl;
    std::cerr << "\t\t--recPoolSize : number of records to allow to be stored at a time, default value: " << myMaxAllowedRecs << std::endl;
    std::cerr << "\t\t                Set to -1 for unlimited number of records" << std::endl;
    std::cerr << "\t\t--posDiff     : max base pair difference between possibly matching records, default value: " << myThreshold << std::endl;
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
    myNumPoolOverflows = 0;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in1", &inFile1)
        LONG_STRINGPARAMETER("in2", &inFile2)
        LONG_STRINGPARAMETER("out",&myDiffFileName)
        LONG_PARAMETER("all", &myCompAll)
        LONG_PARAMETER("flag", &myCompFlag)
        LONG_PARAMETER("mapQual", &myCompMapQ)
        LONG_PARAMETER("mate", &myCompMate)
        LONG_PARAMETER("isize", &myCompISize)
        LONG_PARAMETER("seq", &myCompSeq)
        LONG_PARAMETER("baseQual", &myCompBaseQual)
        LONG_STRINGPARAMETER("tags", &myTags)
        LONG_PARAMETER("everyTags", &myEveryTag)
        LONG_PARAMETER("noCigar", &noCigar)
        LONG_PARAMETER("noPos", &noPos)
        LONG_PARAMETER("onlyDiffs", &myOnlyDiffs)
        LONG_INTPARAMETER("recPoolSize", &myMaxAllowedRecs)
        LONG_INTPARAMETER("posDiff", &myThreshold)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));
    
    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);
    
    myCompCigar = !noCigar;
    myCompPos = !noPos;

    // If all is specified, turn all comparisons on.
    if(myCompAll)
    {
        myCompCigar = true;
        myCompPos = true;
        myCompBaseQual = true;
        myCompSeq = true;
        myCompFlag = true;
        myCompMapQ = true;
        myCompMate = true;
        myCompISize = true;
        myEveryTag = true;
    }

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

    // Check the type of the diff by looking at the filename.
    int extIndex = myDiffFileName.FindLastChar('.');
    String extension = myDiffFileName.SubStr(extIndex);
    if(extension == ".sam" || extension == ".bam" ||
       extension == ".ubam")
    {
        myBamOut = true;

        // Set the diff filenames.
        myBamDiffName = myDiffFileName;
        
        String baseName = myDiffFileName.Left(extIndex);
        // Start after the last '/' and go to the last '.'
        // Note: if not found, -1 is returned.
        int startBasePos = inFile1.FindLastChar('/') + 1;
        int dotPos = inFile1.FindLastChar('.');
        if(dotPos == -1)
        {
            // no dot, so return the len.
            dotPos = inFile1.Length();
        }
        String inBase1 = inFile1.SubStr(startBasePos,
                                        dotPos - startBasePos);
        startBasePos = inFile2.FindLastChar('/') + 1;
        dotPos = inFile2.FindLastChar('.');
        if(dotPos == -1)
        {
            // no dot, so return the len.
            dotPos = inFile2.Length();
        }
        String inBase2 = inFile2.SubStr(startBasePos,
                                        dotPos - startBasePos);

        myBamOnly1Name = baseName + "_only1_" + inBase1 + extension;
        myBamOnly2Name = baseName + "_only2_" + inBase2 + extension;
    }
    else
    {
        myBamOut = false;
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
        fprintf(stderr, "Failed to allocate initial records, exiting!\n");
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
        // Do not need to check for null, because less than returns false if
        // tempRecord is null
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

    if(myNumPoolOverflows != 0)
    {
        std::cerr << "WARNING: Matching records may incorrectly be reported as "
                  << "mismatches due to running out of available records " 
                  << myNumPoolOverflows
                  << " time";
        if(myNumPoolOverflows != 1)
        {
            std::cerr << "s";
        }
        std::cerr << ".\nTry increasing --recPoolSize from "
                  << myMaxAllowedRecs << " or setting it "
                  << "to -1 (unlimited).\n";
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


void Diff::writeBamDiffs(SamRecord* rec1, SamRecord* rec2)
{
    static String tempString;

    if((rec1 == NULL) && (rec2 != NULL))
    {
        // If myBamOnly2 file has not open yet, initialize it.
        if(!myBamOnly2.IsOpen())
        {
            // not yet open.
            myBamOnly2.OpenForWrite(myBamOnly2Name.c_str(), &myFile2.header);
        }
        myBamOnly2.WriteRecord(myFile2.header, *rec2);
    }
    else if((rec1 != NULL) && (rec2 == NULL))
    {
        // If myBamOnly1 file has not open yet, initialize it.
        if(!myBamOnly1.IsOpen())
        {
            // not yet open.
            myBamOnly1.OpenForWrite(myBamOnly1Name.c_str(), &myFile1.header);
        }
        myBamOnly1.WriteRecord(myFile1.header, *rec1);
    }
    else if((rec1 != NULL) && (rec2 != NULL))
    {
        if(!getDiffs(rec1, rec2))
        {
            // No Diffs.
            return;
        }
        // If myDiffOnly2 file has not open yet, initialize it.
        if(!myBamDiff.IsOpen())
        {
            // not yet open.
            myBamDiff.OpenForWrite(myBamDiffName.c_str(), &myFile1.header);
        }

        //  Add the fields from rec2.
        if(myCompPos && (!myOnlyDiffs || myDiffStruct.posDiff))
        {
            tempString = rec2->getReferenceName();
            tempString += ':';
            tempString += rec2->get1BasedPosition();
            rec1->addTag(POS_DIFF_TAG, POS_DIFF_TYPE, tempString.c_str());
        }
        if(myCompCigar && (!myOnlyDiffs || myDiffStruct.cigarDiff))
        {
            rec1->addTag(CIGAR_DIFF_TAG, CIGAR_DIFF_TYPE, rec2->getCigar());
        }
        if(myCompFlag && (!myOnlyDiffs || myDiffStruct.flagDiff || (rec1->getFlag() != rec2->getFlag())))
        {
            rec1->addIntTag(FLAG_DIFF_TAG, rec2->getFlag());
        }
        if(myCompMapQ && (!myOnlyDiffs || myDiffStruct.mapqDiff))
        {
            rec1->addIntTag(MAPQ_DIFF_TAG, rec2->getMapQuality());
        }
        if(myCompMate && (!myOnlyDiffs || myDiffStruct.mateDiff))
        {
            tempString = rec2->getMateReferenceName();
            tempString += ':';
            tempString += rec2->get1BasedMatePosition();
            rec1->addTag(MATE_DIFF_TAG, MATE_DIFF_TYPE, tempString.c_str());
        }
        if(myCompISize && (!myOnlyDiffs || myDiffStruct.isizeDiff))
        {
            rec1->addIntTag(ISIZE_DIFF_TAG, rec2->getInsertSize());
        }
        if(myCompSeq && (!myOnlyDiffs || myDiffStruct.seqDiff))
        {
            rec1->addTag(SEQ_DIFF_TAG, SEQ_DIFF_TYPE, rec2->getSequence());
        }
        if(myCompBaseQual && (!myOnlyDiffs || myDiffStruct.qualDiff))
        {
            rec1->addTag(QUAL_DIFF_TAG, QUAL_DIFF_TYPE, rec2->getQuality());
        }
        if(!myTags2.IsEmpty() && (!myOnlyDiffs || myDiffStruct.tagsDiff))
        {
            rec1->addTag(TAGS_DIFF_TAG, TAGS_DIFF_TYPE, myTags2);
        }

        myBamDiff.WriteRecord(myFile1.header, *rec1);
    }
}


void Diff::writeDiffDiffs(SamRecord* rec1, SamRecord* rec2)
{
    if(!getDiffs(rec1, rec2))
    {
        // No diffs, so just return.
        return;
    }

    // We know there is some diff.
    if(rec1 != NULL)
    {
        // Flag is always automatically written so don't add it here.
        myDiff1.Clear();
        // If print all values on any diff or if there is a pos diff
        if(myCompPos && (!myOnlyDiffs ||  myDiffStruct.posDiff))
        {
            myDiff1 += '\t';
            myDiff1 += rec1->getReferenceName();
            myDiff1 += ':';
            myDiff1 += rec1->get1BasedPosition();
        }
        if(myCompCigar && (!myOnlyDiffs ||  myDiffStruct.cigarDiff))
        {
            // Get the info for rec1's cigar.
             myDiff1 += "\t";
             myDiff1 += rec1->getCigar();
        }
        if(myCompMapQ && (!myOnlyDiffs ||  myDiffStruct.mapqDiff))
        {
            myDiff1 += "\t";
            myDiff1 += rec1->getMapQuality();
        }
        if(myCompMate && (!myOnlyDiffs ||  myDiffStruct.mateDiff))
        {
            myDiff1 += "\t";
            myDiff1 += rec1->getMateReferenceName();
            myDiff1 += ":";
            myDiff1 += rec1->get1BasedMatePosition();
        }
        if(myCompISize && (!myOnlyDiffs ||  myDiffStruct.isizeDiff))
        {
            myDiff1 += "\t";
            myDiff1 += rec1->getInsertSize();
        }
        if(myCompSeq && (!myOnlyDiffs ||  myDiffStruct.seqDiff))
        {
            myDiff1 += "\t";
            myDiff1 += rec1->getSequence();
        }
        if(myCompBaseQual && (!myOnlyDiffs ||  myDiffStruct.qualDiff))
        {
            myDiff1 += "\t";
            myDiff1 += rec1->getQuality();
        }
        if((!myOnlyDiffs ||  myDiffStruct.tagsDiff) && !myTags1.IsEmpty())
        {
            myDiff1 += "\t";
            myDiff1 += myTags1;
        }
    }
    if(rec2 != NULL)
    {
        // Flag is always automatically written so don't add it here.
        myDiff2.Clear();
        // If print all values on any diff or if there is a pos diff
        if(myCompPos && (!myOnlyDiffs ||  myDiffStruct.posDiff))
        {
            myDiff2 += '\t';
            myDiff2 += rec2->getReferenceName();
            myDiff2 += ':';
            myDiff2 += rec2->get1BasedPosition();
        }
        if(myCompCigar && (!myOnlyDiffs ||  myDiffStruct.cigarDiff))
        {
            // Get the info for rec2's cigar.
             myDiff2 += "\t";
             myDiff2 += rec2->getCigar();
        }
        if(myCompMapQ && (!myOnlyDiffs ||  myDiffStruct.mapqDiff))
        {
            myDiff2 += "\t";
            myDiff2 += rec2->getMapQuality();
        }
        if(myCompMate && (!myOnlyDiffs ||  myDiffStruct.mateDiff))
        {
            myDiff2 += "\t";
            myDiff2 += rec2->getMateReferenceName();
            myDiff2 += ":";
            myDiff2 += rec2->get1BasedMatePosition();
        }
        if(myCompISize && (!myOnlyDiffs ||  myDiffStruct.isizeDiff))
        {
            myDiff2 += "\t";
            myDiff2 += rec2->getInsertSize();
        }
        if(myCompSeq && (!myOnlyDiffs ||  myDiffStruct.seqDiff))
        {
            myDiff2 += "\t";
            myDiff2 += rec2->getSequence();
        }
        if(myCompBaseQual && (!myOnlyDiffs ||  myDiffStruct.qualDiff))
        {
            myDiff2 += "\t";
            myDiff2 += rec2->getQuality();
        }
        if((!myOnlyDiffs ||  myDiffStruct.tagsDiff) && !myTags2.IsEmpty())
        {
            myDiff2 += "\t";
            myDiff2 += myTags2;
        }
    }

    if(rec1 != NULL)
    {
        if(checkDiffFile())
        {
            writeReadName(*rec1);
            myDiff1 += '\n';
            ifprintf(myDiffFile, "\n<\t%x", rec1->getFlag());
            ifwrite(myDiffFile, myDiff1, myDiff1.Length());
        }
    }
    if(rec2 != NULL)
    {
        if(checkDiffFile())
        {
            if(rec1 == NULL)
            {
                // Only diff2, so write the readname/fragment.
                writeReadName(*rec2);
                ifwrite(myDiffFile, "\n", 1);
            }
            myDiff2 += '\n';
            ifprintf(myDiffFile, ">\t%x", rec2->getFlag());
            ifwrite(myDiffFile, myDiff2, myDiff2.Length());
        }
    }
}


void Diff::writeDiffs(SamRecord* rec1, SamRecord* rec2)
{
    if(myBamOut)
    {
        writeBamDiffs(rec1, rec2);
    }
    else
    {
        writeDiffDiffs(rec1, rec2);
    }
}


bool Diff::getDiffs(SamRecord* rec1, SamRecord* rec2)
{
    if((rec1 == NULL) && (rec2 == NULL))
    {
        // Neither is set, so no diffs.
        return(false);
    }
    
    myTags1.Clear();
    myTags2.Clear();
    char tagDelim = '\t';
    if(myBamOut)
    {
        tagDelim = ';';
    }

    // Read the tags.
    if(myEveryTag)
    {
        if(rec1 != NULL)
        {
            // Get all the tags from the records.
            if(!SamRecordHelper::genSamTagsString(*rec1, myTags1, tagDelim))
            {
                std::cerr << "Failed to read tags\n";
            }
        }
        if(rec2 != NULL)
        {
            // Get all the tags from the records.
            if(!SamRecordHelper::genSamTagsString(*rec2, myTags2, tagDelim))
            {
                    std::cerr << "Failed to read tags\n";
            }
        }
    }
    else if(!myTags.IsEmpty())
    {
        // Compare the tags.
        if(((rec1 != NULL) &&
            (!rec1->getTagsString(myTags.c_str(), myTags1, tagDelim))) ||
           ((rec2 != NULL) && 
            (!rec2->getTagsString(myTags.c_str(), myTags2, tagDelim))))
        {
            // Failure reading tags.
            std::cerr << "Failed to search for tags: "
                      << myTags.c_str() << std::endl;
        }
    }

    if(((rec1 == NULL) && (rec2 != NULL)) || 
       ((rec1 != NULL) && (rec2 == NULL)))
    {
        // Only one of the records is set.
        myDiffStruct.posDiff = myCompPos;
        myDiffStruct.cigarDiff = myCompCigar;
        myDiffStruct.flagDiff = myCompFlag;
        myDiffStruct.mapqDiff = myCompMapQ;
        myDiffStruct.mateDiff = myCompMate;
        myDiffStruct.isizeDiff = myCompISize;
        myDiffStruct.seqDiff = myCompSeq;
        myDiffStruct.qualDiff = myCompBaseQual;
        myDiffStruct.tagsDiff = (!myTags.IsEmpty()) || myEveryTag;
        return(true);
    }
    else
    {
        // Check to see if there are any diffs:
        myDiffStruct.posDiff = false;
        myDiffStruct.cigarDiff = false;
        myDiffStruct.flagDiff = false;
        myDiffStruct.mapqDiff = false;
        myDiffStruct.mateDiff = false;
        myDiffStruct.isizeDiff = false;
        myDiffStruct.seqDiff = false;
        myDiffStruct.qualDiff = false;
        myDiffStruct.tagsDiff = false;
   
        // Check if there are any differences.
        
        // Check if different chromosomes or positions.
        if(myCompPos && 
           ((rec1->getReferenceID() != rec2->getReferenceID()) ||
            (rec1->get1BasedPosition() != rec2->get1BasedPosition())))
        {
            // It moved chromosome/position.
            myDiffStruct.posDiff = true;
        }
        if(myCompCigar && (strcmp(rec1->getCigar(), rec2->getCigar()) != 0))
        {
            myDiffStruct.cigarDiff = true;
        }
        if(myCompFlag && (rec1->getFlag() != rec2->getFlag()))
        {
            myDiffStruct.flagDiff = true;
        }
        if(myCompMapQ && (rec1->getMapQuality() != rec2->getMapQuality()))
        {
            myDiffStruct.mapqDiff = true;
        }
        if(myCompMate && 
           ((rec1->getMateReferenceID() != rec2->getMateReferenceID()) ||
            (rec1->get1BasedMatePosition() != rec2->get1BasedMatePosition())))
        {
            myDiffStruct.mateDiff = true;
        }
        if(myCompISize && (rec1->getInsertSize() != rec2->getInsertSize()))
        {
            myDiffStruct.isizeDiff = true;
        }
        if(myCompSeq && (strcmp(rec1->getSequence(), rec2->getSequence()) != 0))
        {
            myDiffStruct.seqDiff = true;
        }
        if(myCompBaseQual && (strcmp(rec1->getQuality(), rec2->getQuality()) != 0))
        {
            myDiffStruct.qualDiff = true;
        }
        if(myTags1 != myTags2)
        {
            myDiffStruct.tagsDiff = true;
        }
    }
    return(myDiffStruct.posDiff  || myDiffStruct.cigarDiff || 
           myDiffStruct.flagDiff || myDiffStruct.mapqDiff || 
           myDiffStruct.mateDiff || myDiffStruct.isizeDiff || 
           myDiffStruct.seqDiff  || myDiffStruct.qualDiff || 
           myDiffStruct.tagsDiff);
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
        ++myNumPoolOverflows;

        if(myFile1Unmatched.size() >= myFile2Unmatched.size())
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
