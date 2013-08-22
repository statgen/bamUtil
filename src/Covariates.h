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
#include <stdexcept>
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

    int32_t rgid;

    BaseData()
        : qual(0), cycle(0), read(false), preBase('?'), 
          curBase('?'), rgid(0)
    {}

    BaseData(const BaseData& other)
    {
        qual = other.qual;
        cycle = other.cycle;
        read = other.read;
        preBase = other.preBase;
        curBase = other.curBase;
        rgid = other.rgid;
    }

    inline uint64_t getKey() const
    {
        return( ((uint64_t)(read & 0x1) << 63) |
                ((uint64_t)(qual & 0x7F) << 56) | 
                ((uint64_t)cycle << 40) | 
                ((uint64_t)(BaseAsciiMap::base2int[(int)preBase]) << 36) | 
                ((uint64_t)(BaseAsciiMap::base2int[(int)curBase]) << 32) | 
                (rgid) );
    }

    inline void parseKey(uint64_t key)
    {
        read = (key >> 63);
        qual = (key >> 56) & 0x7F;
        cycle = (key >> 40) & 0xFF;
        preBase = BaseAsciiMap::int2base[(key >> 36) & 0xF];
        curBase = BaseAsciiMap::int2base[(key >> 32) & 0xF];
        rgid = key & 0xFFFFFFFF;
           
    }

    inline static uint32_t getFastKeySize()
    {
        return(0x0FFFFFFF);
    }


    inline static int32_t getFastMaxRG()
    {
        return(0xFF);
    }


    inline static uint32_t getFastMaxQual()
    {
        return(0x3F);
    }


    inline static int32_t getFastMaxCycles()
    {
        return(0x7F);
    }


    inline uint32_t getFastKey() const
    {
        if(rgid > getFastMaxRG())
        {
            String errorMsg = "RECAB: Cannot use --fast option with more than ";
            errorMsg += (getFastMaxRG() + 1);
            errorMsg += " read groups";
            throw(std::runtime_error(errorMsg.c_str()));
        }
        if(qual > getFastMaxQual())
        {
            String errorMsg = "RECAB: Cannot use --fast option with qualities > ";
            errorMsg += getFastMaxQual();
            throw(std::runtime_error(errorMsg.c_str()));
        }
        if(cycle > getFastMaxCycles())
        {
            String errorMsg = "RECAB: Cannot use --fast option with cycles > ";
            errorMsg += getFastMaxCycles();
            throw(std::runtime_error(errorMsg.c_str()));
        }
        return( ((uint32_t)(read & 0x1) << 27) |
                ((uint32_t)(cycle & 0x7F) << 20) | 
                ((uint32_t)(qual & 0x3F) << 14) | 
                ((uint32_t)(BaseAsciiMap::base2int[(int)preBase] & 0x7) << 11) | 
                ((uint32_t)(BaseAsciiMap::base2int[(int)curBase] & 0x7) << 8) | 
                rgid );
    }

    inline void parseFastKey(uint32_t key)
    {
        read = (key >> 27);
        cycle = (key >> 20) & 0x7F;
        qual = (key >> 14) & 0x3F;
        preBase = BaseAsciiMap::int2base[(key >> 11) & 0x7];
        curBase = BaseAsciiMap::int2base[(key >> 8) & 0x7];
        rgid = key & 0xFF;
           
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
