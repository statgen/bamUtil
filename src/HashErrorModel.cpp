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
#include <math.h>


HashErrorModel::HashErrorModel(){
	lastElement = 0;
	//set minimum size
	mismatchTable.rehash((std::size_t) MINSIZE);
}

HashErrorModel::~HashErrorModel() {
};

//ReadGroup handling
uint8_t HashErrorModel::mapReadGroup(std::string readgroup) {

    if(readgroupDic.find(readgroup) != readgroupDic.end())
        return 	readgroupDic[readgroup];
    else
    {
    	lastElement++;
    	readgroupDic[readgroup] = lastElement;
    	rereadgroupDic[lastElement] = readgroup;
    	return lastElement;
    }
}

uint64_t HashErrorModel::calcKey(baseCV basecv) {

	//Key construction
    uint16_t rg = mapReadGroup(basecv.rg);
    return (
       	 (0xffff000000000000 & ( static_cast<uint64_t>(basecv.tile)    << 48 ) ) |
    	 (0x0000ffff00000000 & ( static_cast<uint64_t>(rg)             << 32 ) ) |
	     (0x00000000ff000000 & ( static_cast<uint64_t>(basecv.q)       << 24 ) ) |
	     (0x0000000000ff0000 & ( static_cast<uint64_t>(basecv.pos)     << 16 ) ) |
	   	 (0x000000000000f000 & ( static_cast<uint64_t>(basecv.read)    << 12 ) ) |
         (0x0000000000000f00 & ( static_cast<uint64_t>(basecv.prebase) << 8  ) ) |
	     //(0x00000000000000f0 & ( static_cast<uint64_t>(basecv.nexbase) << 4  ) ) |
	     (0x000000000000000f & ( static_cast<uint64_t>(basecv.curbase)     ) )
	     );
}

baseCV HashErrorModel::revKey(uint64_t key) {

	baseCV basecv;
	//Key construction
    basecv.tile =    static_cast<uint8_t> ((key & 0xffff000000000000) >> 48);
    uint16_t rg =    static_cast<uint16_t>((key & 0x0000ffff00000000) >> 32);
    basecv.rg = rereadgroupDic[rg];
    basecv.q =       static_cast<uint8_t>((key &  0x00000000ff000000) >> 24);
    basecv.pos =     static_cast<uint16_t>((key & 0x0000000000ff0000) >> 16);
    basecv.read =    static_cast<uint8_t>((key &  0x000000000000f000) >> 12);
    basecv.prebase = static_cast<uint8_t>((key &  0x0000000000000f00) >> 8);
    //basecv.nexbase = static_cast<uint8_t>((key &  0x00000000000000f0) >> 4);
    basecv.curbase = static_cast<uint8_t>(key &   0x000000000000000f);
    return basecv;
}

void HashErrorModel::setCell(baseCV basecv){
	uint64_t key = calcKey(basecv);

	SMatches smatch = {0,0,0.0};
	HashMatch::iterator it = mismatchTable.insert(HashMatch::value_type(key, smatch)).first;
	if(basecv.obs==basecv.ref)
	  it->second.m++;
	else{
	  it->second.mm++;
	}
}

uint8_t HashErrorModel::getQemp(baseCV basecv){
	double qs;
	uint64_t key = calcKey(basecv);
    uint64_t match =  mismatchTable[key].m;
    uint64_t mismatch = mismatchTable[key].mm;
    qs = log10((double)(mismatch+1)/(double)(mismatch+match+1)) * (-10.0);
    return floor(qs + 0.5);
}


int HashErrorModel::writeAllTableMM(String filename){

    FILE *pFile;
	pFile = fopen(filename.c_str(), "w");
	if(!pFile) return 0;

	for(HashMatch::const_iterator it = mismatchTable.begin();
	    it != mismatchTable.end();
	    ++it)
	{
		baseCV basecv = revKey(it->first);
	    uint64_t m = it->second.m;
	    uint64_t mm = it->second.mm;
	    uint8_t qemp = it->second.qemp;
		fprintf(pFile,"%s %d %d %d %d %d %d %d - %ld %ld %d\n",basecv.rg.c_str(),basecv.q,basecv.pos,basecv.read,basecv.tile,basecv.prebase,basecv.nexbase,basecv.curbase,m,mm,qemp);
	}
   fclose(pFile);
   return 1;
}


int HashErrorModel::getKeyLength(){
	// pos 0..119
	//return 7-1+120
	//return 6;
	return 5; // next information
};

uint32_t HashErrorModel::getSize(){
	return mismatchTable.size();
};

void HashErrorModel::setModel(Model model){
  this->model = model;
};

void HashErrorModel::addPrediction(Model model,int blendedWeight){
	this->model = model;
	for(HashMatch::iterator it = mismatchTable.begin();
	    it != mismatchTable.end();
	    ++it)
	{
		baseCV basecv = revKey(it->first);
		basecv.setCovariates();
		int j = 1;
		double qemp = model[0]; //slope
		for(std::vector<uint16_t>::const_iterator itv = basecv.covariates.begin();
		    itv != basecv.covariates.end();
		    ++itv)
		{
			qemp += model[j]*(double)(*itv);
			j++;
		}
		//phred-score transformation
		//conservative (otherwise +0.5)
		//blended Model?
		int phred = 0;
		if(blendedWeight==9999){
			uint64_t m = it->second.m;
			uint64_t mm = it->second.mm;
			qemp = (mm-blendedWeight)/(m+mm-blendedWeight);
		}
		else{
			if(blendedWeight>0){
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


void HashErrorModel::setDataforPrediction(Matrix & X, Vector & succ, Vector & total,bool binarizeFlag){

	int cols = getKeyLength()+1;
	uint32_t rows = getSize();

	succ.Dimension(rows);
	total.Dimension(rows);
	X.Dimension(rows, cols);
	X.Zero();

	int i = 0;
	for(HashMatch::const_iterator it = mismatchTable.begin();
	    it != mismatchTable.end();
	    ++it)
	{
		// The first column of the design matrix is constant one, for the slope
		X[i][0] = 1.0;

		baseCV basecv = revKey(it->first);
		basecv.setCovariates();
		int j = 0;

		//binarize a couple of co-variates
		for(std::vector<uint16_t>::const_iterator itv = basecv.covariates.begin();
		    itv != basecv.covariates.end();
		    ++itv)
		{
            if(binarizeFlag){
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

