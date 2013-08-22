/*
 *  Copyright (C) 2010-2012  Christian Fuchsberger,
 *                      Regents of the University of Michigan
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
#ifndef __HASH_ERROR_MODEL_H__
#define __HASH_ERROR_MODEL_H__


#include <stdint.h>
#include <vector>

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <unordered_map>
#else
#include <map>
#endif

#include "StringArray.h"
#include "Covariates.h"
#include "MathMatrix.h"
#include "MathVector.h"

#define MINSIZE 1048576


class HashErrorModel {
    
public:

    static void setUseLogReg(bool useLogReg) { ourUseLogReg = useLogReg; }
    static void setUseFast(bool useFast) { ourUseFast = useFast; }
    
    
    typedef std::vector<double> Model;
    struct SMatches{
        SMatches() : m(0), mm(0), qempSimple(0), qempLogReg(0) {}
        uint32_t m;
        uint32_t mm;
        uint8_t qempSimple;
        uint8_t qempLogReg;
    };
    struct SMatchesFast{
        SMatchesFast() : m(0), mm(0), qempSimple(0) {}
        uint32_t m;
        uint32_t mm;
        uint8_t qempSimple;
    };
    

    #ifdef __GXX_EXPERIMENTAL_CXX0X__
    typedef std::unordered_map<uint64_t, HashErrorModel::SMatches> HashMatch;
    #else
    typedef std::map<uint64_t, HashErrorModel::SMatches> HashMatch;
    #endif

    HashMatch mismatchTable;
    std::vector<SMatchesFast> mismatchTableFast;
    uint16_t lastElement;
    
    HashErrorModel();
    ~HashErrorModel();
    
    void setCell(const BaseData& data, char refBase);
    uint8_t getQemp(BaseData& data);
    uint8_t getQempSimple(uint32_t matches, uint32_t mismatches);
    int writeTableQemp(std::string& filename, 
                       const std::vector<std::string>& id2rg,
                       bool logReg);

    void setDataforPrediction(Matrix & X, Vector & succ, Vector& total,bool binarizeFlag);
    void addPrediction(Model model, int blendedWeight);
    
private:
    static bool ourUseLogReg;
    static bool ourUseFast;
};

#endif
