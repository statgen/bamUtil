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
// This file contains the processing for the executable option "stats"
// which generates some statistics for SAM/BAM files.

#ifndef __INDEL_DISCORDANCE_H__
#define __INDEL_DISCORDANCE_H__

#include "BamExecutable.h"
#include "SamRecordPool.h"
#include "PileupElement.h"

class IndelDiscordance : public BamExecutable
{
public:
    IndelDiscordance();

    static void indelDiscordanceDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);

private:
    class PileupElementIndelDiscordance : public PileupElement
    {
    public:
        /// Set record pool.
        static void setRecordPool(SamRecordPool& pool)
        {
            ourPool = &pool;
        }

        PileupElementIndelDiscordance();
        
        virtual ~PileupElementIndelDiscordance();
        
        // Add an entry to this pileup element.  
        virtual void addEntry(SamRecord& record);
        
        // Perform the analysis associated with this class.
        virtual void analyze();
        
        // Resets the entry, setting the new position associated with this
        // element.
        virtual void reset(int32_t refPosition);
        
    private:
        PileupElementIndelDiscordance(const PileupElement& q);
        
        void initVars();
        void releaseRecords();

        int depth;

        static SamRecordPool* ourPool;
        std::vector<SamRecord*> myRecords;
    };
};

#endif
