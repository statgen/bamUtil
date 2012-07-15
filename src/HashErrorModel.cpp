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
#include "HashErrorModel.h"
#include "BaseUtilities.h"
#include <math.h>

bool HashErrorModel::ourUseLogReg = true;

HashErrorModel::HashErrorModel()
{
	lastElement = 0;
}

HashErrorModel::~HashErrorModel() 
{
}


void HashErrorModel::setCell(const BaseData& data, char refBase)
{
    SMatches& matchInfo = mismatchTable[data.getKey()];

    if(BaseUtilities::areEqual(refBase, data.curBase))
    {
        ++(matchInfo.m);
    }
    else
    {
        ++(matchInfo.mm);
    }
}

uint8_t HashErrorModel::getQemp(BaseData& data)
{
    double qs;

    SMatches& matchMismatch = mismatchTable[data.getKey()];
    if(ourUseLogReg)
    {
        return(matchMismatch.qemp);
    }
    else
    {
        qs = log10((double)(matchMismatch.mm+1)/
                   (double)(matchMismatch.mm+matchMismatch.m+1)) * (-10.0);
        return floor(qs + 0.5);
    }
}


int HashErrorModel::writeAllTableMM(String filename, 
                                    const std::vector<std::string>& id2rg)
{
    FILE *pFile;
    pFile = fopen(filename.c_str(), "w");
    if(!pFile) return 0;
    
    int maxId = id2rg.size() - 1;

    BaseData data;

    for(HashMatch::const_iterator it = mismatchTable.begin();
        it != mismatchTable.end();
        ++it)
    {
        data.parseKey(it->first);
        if(data.rgid <= maxId)
        {
//TODO             fprintf(pFile,"%s %d %d %d %d %d %d %d - %ld %ld %d\n",
//                     id2rg[it->first.rgid].c_str(), it->first.qual, it->first.cycle,
//                     it->first.read, 0, BaseAsciiMap::base2int[(int)(it->first.preBase)], 0, BaseAsciiMap::base2int[(int)(it->first.curBase)],
//                     it->second.m, it->second.mm, it->second.qemp);
            int16_t cycle = data.cycle + 1;
            if(data.read)
            {
                // 2nd.
                cycle = -cycle;
            }
            fprintf(pFile,"%s,%d,%d,%c%c,%ld,%ld,%d\n",
                    id2rg[data.rgid].c_str(), data.qual, cycle, 
                    data.preBase, data.curBase,
                    it->second.m + it->second.mm, it->second.mm, it->second.qemp);
        }
    }
    fclose(pFile);
    return 1;
}


uint32_t HashErrorModel::getSize(){
	return mismatchTable.size();
};


void HashErrorModel::addPrediction(Model model,int blendedWeight)
{
    BaseData data;
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
            uint64_t m = it->second.m;
            uint64_t mm = it->second.mm;
            qemp = (mm-blendedWeight)/(m+mm-blendedWeight);
        }
        else
        {
            if(blendedWeight>0)
            {
                uint64_t m = it->second.m;
                uint64_t mm = it->second.mm;
                //qemp = (mm+qemp*blendedWeight)/(m+mm+blendedWeight);
                qemp = (mm+qemp*mm)/(2.0*(m+mm));
            }
        }

        phred = trunc((-10.0*log10(1.0-(1.0/(1.0+exp(-qemp)))))+0.5);
        it->second.qemp = phred;
    }
};


void HashErrorModel::setDataforPrediction(Matrix & X, Vector & succ, Vector & total,bool binarizeFlag)
{
    int i = 0;
    BaseData data;
    for(HashMatch::const_iterator it = mismatchTable.begin();
        it != mismatchTable.end();
        ++it)
    {
        data.parseKey(it->first);
        Covariates cov;
        cov.setCovariates(data);
        if(i == 0)
        {
            uint32_t rows = getSize();
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

