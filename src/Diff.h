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
// which reads two SAM/BAM files and writes the differences.

#ifndef __DIFF_H__
#define __DIFF_H__

#include <stack>
#include <list>
#include <map>
#include "BamExecutable.h"
#include "SamFile.h"

class Diff : public BamExecutable
{
public:
    Diff();
    ~Diff();
    static void diffDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);
    virtual const char* getProgramName() {return("bam:diff");}

private:
    struct diffStruct
    {
        bool posDiff;
        bool cigarDiff;
        bool flagDiff;
        bool mapqDiff;
        bool mateDiff;
        bool isizeDiff;
        bool seqDiff;
        bool qualDiff;
        bool tagsDiff;
    } myDiffStruct;
                   
    class FileInfo
    {
    public:
        SamFile file;
        SamFileHeader header;
    };
    
    
    class UnmatchedRecords
    {
    public:
        UnmatchedRecords()
            : myListUnmatched(),
              myFragmentMaps(4),
              myUnmatchedFileIter()
        {
        }

        // Return the number of elements in this unmatched record container.
        inline int size() { return(myListUnmatched.size()); }
        
        // Add the specified record to the list of unmatched records.
        void addUnmatchedRecord(SamRecord& record);

        // Get and remove the record that matches the specified record's
        // query(read) name and fragment (first/last/mid/unknown).
        // If no match is found, return NULL.
        SamRecord* removeFragmentMatch(SamRecord& record);
        
        // Remove the first entry from this unmatched file container, returning a pointer
        // to the record.
        SamRecord* removeFirst();
        // Get the first entry from this unmatched file container, returning a pointer
        // to the record without removing it.
        SamRecord* getFirst();

    private:

        typedef std::map<std::string, std::list<SamRecord*>::iterator> mapType;
        std::list<SamRecord*> myListUnmatched;
        std::vector<mapType> myFragmentMaps;
        std::map<std::string,std::list<SamRecord*>::iterator>::iterator myUnmatchedFileIter;

    };

    // Check to see if the two records are a match - same read name & fragment.
    // Return true if they match, false if not.
    bool matchingRecs(SamRecord* rec1, SamRecord* rec2);

    // Compare two records.  Returns true if the first record chrom/pos are less than or
    // equal to record2's chrom/pos, false if rec2 is less than rec1.
    // Threshold is a value to be added to rec1's position before comparing to rec2.  It is
    // a way to determine if rec1 is more than a certain number of positions less than rec2.
    // Threshold is defaulted to 0, so the excact positions are compared.
    // If they are on different chromosomes, threshold is not used.
    bool lessThan(SamRecord* rec1, SamRecord* rec2, int threshold = 0);

    void writeBamDiffs(SamRecord* rec1, SamRecord* rec2);
    void writeDiffDiffs(SamRecord* rec1, SamRecord* rec2);
    void writeDiffs(SamRecord* rec1, SamRecord* rec2);
    bool getDiffs(SamRecord* rec1, SamRecord* rec2);
    bool writeReadName(SamRecord& record);
    SamRecord* getSamRecord();
    // If record is not NULL, adds it back to the free list.  If record is NULL, nothing is done.
    void releaseSamRecord(SamRecord* record);

    bool checkDiffFile();
    
    static const char* FLAG_DIFF_TAG;
    static const char* POS_DIFF_TAG;
    static const char* CIGAR_DIFF_TAG;
    static const char* MAPQ_DIFF_TAG;
    static const char* MATE_DIFF_TAG;
    static const char* ISIZE_DIFF_TAG;
    static const char* SEQ_DIFF_TAG;
    static const char* QUAL_DIFF_TAG;
    static const char* TAGS_DIFF_TAG;
    static const char POS_DIFF_TYPE = 'Z';
    static const char MATE_DIFF_TYPE = 'Z';
    static const char CIGAR_DIFF_TYPE = 'Z';
    static const char SEQ_DIFF_TYPE = 'Z';
    static const char QUAL_DIFF_TYPE = 'Z';
    static const char TAGS_DIFF_TYPE = 'Z';

    std::stack<SamRecord*> myFreeSamRecords;

    UnmatchedRecords myFile1Unmatched;
    UnmatchedRecords myFile2Unmatched;

    bool myCompAll;
    bool myCompCigar;
    bool myCompPos;
    bool myCompBaseQual;
    bool myCompSeq;
    bool myCompFlag;
    bool myCompMapQ;
    bool myCompMate;
    bool myCompISize;
    String myTags;
    bool myEveryTag;
    bool myOnlyDiffs;
    bool myBamOut;

    int myMaxAllowedRecs;
    int myAllocatedRecs;
    int myThreshold;
    int myNumPoolOverflows;

    FileInfo myFile1;
    FileInfo myFile2;

    String myDiffFileName;
    String myBamOnly1Name;
    String myBamOnly2Name;
    String myBamDiffName;
    IFILE myDiffFile;
    SamFile myBamOnly1;
    SamFile myBamOnly2;
    SamFile myBamDiff;
    String myDiff1;
    String myDiff2;
    String myTags1;
    String myTags2;
    String myTempBuffer;
};

#endif
