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

#include <sstream>
#include <iostream>
#include <stdio.h>
#include <stdexcept>
#include "HashErrorModel.h"
#include "BaseUtilities.h"
#include <math.h>

bool HashErrorModel::ourUseLogReg = true;
bool HashErrorModel::ourUseFast = false;

HashErrorModel::HashErrorModel()
{
	lastElement = 0;
}

HashErrorModel::~HashErrorModel() 
{
}


void HashErrorModel::setCell(const BaseData& data, char refBase)
{
    if(ourUseFast)
    {
        if(mismatchTableFast.empty())
        {
            mismatchTableFast.resize(BaseData::getFastKeySize());
        }
        SMatchesFast& matchInfo = mismatchTableFast[data.getFastKey()];
        if(BaseUtilities::areEqual(refBase, data.curBase))
        {
            ++(matchInfo.m);
        }
        else
        {
            ++(matchInfo.mm);
        }
        matchInfo.qempSimple = 255;
        return;
    }
    SMatches& matchInfo = mismatchTable[data.getKey()];

    if(BaseUtilities::areEqual(refBase, data.curBase))
    {
        ++(matchInfo.m);
    }
    else
    {
        ++(matchInfo.mm);
    }
    matchInfo.qempSimple = 255;
}

uint8_t HashErrorModel::getQemp(BaseData& data)
{
    if(ourUseFast)
    {
        SMatchesFast& matchInfo = mismatchTableFast[data.getFastKey()];

        // If no matches or mismatches, return the original quality.
        if((matchInfo.m == 0) && (matchInfo.mm == 0))
        {
            return(data.qual);
        }

        if(matchInfo.qempSimple == 255)
        {
            matchInfo.qempSimple = 
                getQempSimple(matchInfo.m, matchInfo.mm);
        }
        return(matchInfo.qempSimple);
    }
    
    // If it is not in the table, return the original quality.
    HashMatch::iterator iter = mismatchTable.find(data.getKey());
    if(iter == mismatchTable.end())
    {
        // Not in the table, so just return the original quality.
        return(data.qual);
    }

    // in the table, so get the qemp
    if(ourUseLogReg)
    {
        return(iter->second.qempLogReg);
    }
    if(iter->second.qempSimple == 255)
    {
        iter->second.qempSimple = 
            getQempSimple(iter->second.m, iter->second.mm);
    }
    return(iter->second.qempSimple);
}


uint8_t HashErrorModel::getQempSimple(uint32_t matches, uint32_t mismatches)
{
    // TODO, if this is used, then make it take the configuration from Recab.cpp
    double qs = 40;
    //   if(mismatches != 0)
    {
        qs = log10((double)(mismatches+1)/
                   (double)(matches + mismatches+1)) * (-10.0);
//         qs = log10((double)(mismatches)/
//                    (double)(matches + mismatches)) * (-10.0);
    }
    return(floor(qs + 0.5));
}


int HashErrorModel::writeTableQemp(std::string& filename, 
                                   const std::vector<std::string>& id2rg,
                                   bool logReg)
{
    FILE *pFile;
    pFile = fopen(filename.c_str(), "w");
    if(!pFile) return 0;
    
    int maxId = id2rg.size() - 1;

    BaseData data;

    if(ourUseFast)
    {
        for(unsigned int i = 0; i < mismatchTableFast.size(); i++)
        {
            SMatchesFast& matchInfo = mismatchTableFast[i];
            if((matchInfo.m == 0) && (matchInfo.mm == 0))
            {
                // No data, so continue to the next entry.
                continue;
            }
            data.parseFastKey(i);
            
            if(data.rgid <= maxId)
            {
                int16_t cycle = data.cycle + 1;
                if(data.read)
                {
                    // 2nd.
                    cycle = -cycle;
                }
                
                if(matchInfo.qempSimple == 255)
                {
                    matchInfo.qempSimple = 
                        getQempSimple(matchInfo.m, matchInfo.mm);
                }
                fprintf(pFile,"%s,%d,%d,%c%c,%d,%d,%d\n",
                        id2rg[data.rgid].c_str(), data.qual, cycle, 
                        data.preBase, data.curBase,
                        matchInfo.m + matchInfo.mm, matchInfo.mm, 
                        matchInfo.qempSimple);
            }
        }
        fclose(pFile);
        return(1);
    }
    for(HashMatch::iterator it = mismatchTable.begin();
        it != mismatchTable.end();
        ++it)
    {
        data.parseKey(it->first);
        if(data.rgid <= maxId)
        {
            int16_t cycle = data.cycle + 1;
            if(data.read)
            {
                // 2nd.
                cycle = -cycle;
            }

            if(it->second.qempSimple == 255)
            {
                it->second.qempSimple = 
                    getQempSimple(it->second.m, it->second.mm);
            }
            uint8_t qemp = it->second.qempSimple;
            if(logReg)
            {
                qemp = it->second.qempLogReg;
            }
            
            fprintf(pFile,"%s,%d,%d,%c%c,%d,%d,%d\n",
                    id2rg[data.rgid].c_str(), data.qual, cycle, 
                    data.preBase, data.curBase,
                    it->second.m + it->second.mm, it->second.mm, qemp);
        }
    }
    fclose(pFile);
    return 1;
}


void HashErrorModel::addPrediction(Model model,int blendedWeight)
{
    BaseData data;
    if(ourUseFast)
    {
        return;
    }

    for(HashMatch::iterator it = mismatchTable.begin();
        it != mismatchTable.end();
        ++it)
    {
        Covariates cov;
        data.parseKey(it->first);
        cov.setCovariates(data);
        int j = 1;
        double qemp = model[0]; //slope
        for(std::vector<uint16_t>::const_iterator itv = cov.covariates.begin();
            itv != cov.covariates.end(); ++itv)
        {
            qemp += model[j]*(double)(*itv);
            j++;
        }
        //phred-score transformation
        //conservative (otherwise +0.5)
        //blended Model?
        int phred = 0;
        if(blendedWeight==9999)
        {
            uint32_t m = it->second.m;
            uint32_t mm = it->second.mm;
            qemp = (mm-blendedWeight)/(m+mm-blendedWeight);
        }
        else
        {
            if(blendedWeight>0)
            {
                uint32_t m = it->second.m;
                uint32_t mm = it->second.mm;
                //qemp = (mm+qemp*blendedWeight)/(m+mm+blendedWeight);
                qemp = (mm+qemp*mm)/(2.0*(m+mm));
            }
        }

        phred = trunc((-10.0*log10(1.0-(1.0/(1.0+exp(-qemp)))))+0.5);
        it->second.qempLogReg = phred;
    }
};


void HashErrorModel::setDataforPrediction(Matrix & X, Vector & succ, Vector & total,bool binarizeFlag)
{
    int i = 0;
    BaseData data;

    if(ourUseFast)
    {
        return;
    }

    for(HashMatch::const_iterator it = mismatchTable.begin();
        it != mismatchTable.end();
        ++it)
    {
        data.parseKey(it->first);
        Covariates cov;
        cov.setCovariates(data);
        if(i == 0)
        {
            uint32_t rows = mismatchTable.size();

            succ.Dimension(rows);
            total.Dimension(rows);
            X.Dimension(rows, cov.covariates.size() + 1);
            X.Zero();
        }

        // The first column of the design matrix is constant one, for the slope
        X[i][0] = 1.0;

        int j = 0;
                
        //binarize a couple of co-variates
        for(std::vector<uint16_t>::const_iterator itv = cov.covariates.begin();
            itv != cov.covariates.end();
            ++itv)
        {
            if(binarizeFlag)
            {
                //hardcoded pos is 2
                if(j==1){
                    uint16_t pos = (uint16_t)*itv;
                    // hardcoded
                    pos += 7;
                    X[i][pos] = 1;
                }
                else
                {
                    if(j>1)
                        X[i][j] = (uint16_t)*itv;
                    else
                        X[i][j+1] = (uint16_t)*itv;
                }
                j++;
            }
            else
            {
                X[i][j+1] = (uint16_t)*itv;
                j++;
            }
        }
        total[i] = it->second.mm + it->second.m;
        succ[i] = it->second.m;
        i++;
    }
}

