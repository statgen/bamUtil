/*
 *  Copyright (C) 2010-2012  Christian Fuchsberger,
 *                           Regents of the University of Michigan
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

#ifndef __COVARIATES_H__
#define __COVARIATES_H__

#include <stdint.h>
#include <string>
#include <vector>
#include "BaseAsciiMap.h"

// TODO - move logic for looping through a read and extracting each covariate 
// to here.

struct BaseData
{
public:
    uint8_t qual;
    int16_t cycle;
    bool read;
    char preBase;
    char curBase;

    int rgid;

    char refBase; // Not part of key.

    BaseData()
        : qual(0), cycle(0), read(false), preBase('?'), 
          curBase('?'), rgid(0), refBase('?')
    {}

    BaseData(const BaseData& other)
    {
        qual = other.qual;
        cycle = other.cycle;
        read = other.read;
        preBase = other.preBase;
        curBase = other.curBase;
        rgid = other.rgid;
        refBase = other.refBase;
    }

    inline bool operator <(const BaseData& other) const
    {
        if(rgid == other.rgid)
        {
            // Same RGIDs, so check others
            if(qual == other.qual)
            {
                // Qualities are the same, so check other diffs.
                if(cycle == other.cycle)
                {
                    // cycles are the same, so check read.
                    if(read == other.read)
                    {
                        // reads are the same, so check the bases.
                        if(preBase == other.preBase)
                        {
                            // Same prebase, so check the curBase, return
                            // the result since this is the last field in the key.
                            return(curBase < other.curBase);
                        }
                        // Different prebases, so return the diff.
                        return(preBase < other.preBase);
                    }
                    // Different reads, return less if read is true.
                    return(read);
                }
                // Different cycles, so return the result.
                return(cycle < other.cycle);
            }
            // Different qualities, so return the result.
            return(qual < other.qual);
        }
        // Different rgids
        return(rgid < other.rgid);
    }
private:
};

struct Covariates
{
public:
    std::vector<uint16_t> covariates;

    Covariates() : covariates() {}

    void setCovariates(const BaseData& baseData)
    {
        covariates.clear();
    	covariates.push_back(baseData.qual);
        covariates.push_back(baseData.cycle);
        covariates.push_back(baseData.read);
        covariates.push_back(BaseAsciiMap::base2int[(int)baseData.preBase]);
        covariates.push_back(BaseAsciiMap::base2int[(int)baseData.curBase]);
    }

private:
    Covariates(const Covariates& other);
};

#endif
