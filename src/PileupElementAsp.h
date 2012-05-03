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
// This file contains the processing for the executable option "stats"
// which generates some statistics for SAM/BAM files.

#ifndef __PILEUP_ELEMENT_ASP_H__
#define __PILEUP_ELEMENT_ASP_H__

#include "PileupElement.h"
#include "AspFile.h"


class PileupElementAsp : public PileupElement
{
public:

    /// Set a pointer to the output file to write to.
    /// It is assumed that the file is already opened and the
    /// header already written to it.  This class will not cleanup the
    /// output file, but expects it to exist.
    static void setOutputFile(AspFileWriter* outputFile);

    /// Set whether or not a deletion should be ignored.  Ignored deletions
    /// do not show up in the output. 
    /// If the deletion is not ignored (default), false, the deletion will be
    /// represented as a 'D' in the AspRecord for the base and 0 for the
    /// quality.
    static void setIgnoreDeletion(bool ignoreDeletion);
    
    /// Set the current read name ID for the record that will be processed.
    /// Should be reset for each record.
    static void setCurrentReadNameID(uint16_t readNameID)
    { ourCurrentReadNameID = readNameID; }

    PileupElementAsp();

    /// Constructor that resets the pileup element, does not copy, just resets.
    PileupElementAsp(const PileupElementAsp& q);

    virtual ~PileupElementAsp();

    // Add an entry to this pileup element.  
    virtual void addEntry(SamRecord& record);

    // Perform the analysis associated with this class.
    virtual void analyze();

    // Resets the entry, setting the new position associated with this element.
    virtual void reset(int32_t refPosition);

private:
    void initVars();

    static AspFileWriter* ourOutputFile;
    static bool ourIgnoreDeletion;

    static bool ourReportOverMax;
    static const char UNKNOWN_REF_BASE = 'U';

    static uint16_t ourCurrentReadNameID;

    String myOutputString;

    char myRefBase;
    bool myAllRef;

    AspRecord myAspRecord;
};

#endif
