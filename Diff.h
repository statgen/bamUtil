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
    static void diffDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);

private:
    class FileInfo
    {
    public:
        SamFile file;
        SamFileHeader  header;
        SamRecord* lastRecord;
        
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
        // Get and remove the record that matches the specified record's
        // query(read) name and fragment (first/last/mid/unknown).
        // If no match is found, return NULL.
        SamRecord* removeFragmentMatch(SamRecord& record);
        
    private:
        typedef std::map<std::string, std::list<SamRecord*>::reverse_iterator> mapType;
        std::list<SamRecord*> myListUnmatched;
        std::vector<mapType> myFragmentMaps;
        std::map<std::string,std::list<SamRecord*>::reverse_iterator>::iterator myUnmatchedFileIter;
    };

    //    void writeDiff(IFILE file, SamRecord& record);
    bool writeDiffs(SamRecord* rec1, SamRecord* rec2);
    bool writeReadNameFrag(SamRecord& record);
    bool checkDiffFile();
    SamRecord* getSamRecord();

    std::stack<SamRecord*> myFreeSamRecords;

    UnmatchedRecords myFile1Unmatched;
    UnmatchedRecords myFile2Unmatched;

    bool myCompCigar;
    bool myCompPos;
    bool myCompBaseQual;

    int myMaxAllowedRecs;
    int myAllocatedRecs;

    FileInfo myFile1;
    FileInfo myFile2;

    String myDiffFileName;
    IFILE myDiffFile;
    String myDiff1;
    String myDiff2;
    String myTempBuffer;
};





#endif
