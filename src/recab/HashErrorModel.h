/*
 *  Copyright (C) 2010  Regents of the University of Michigan
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
#include <tr1/unordered_map>
#include "StringArray.h"
#include "BaseCV.h"
#include "MathMatrix.h"
#include "MathVector.h"

#define MINSIZE 1048576


class HashErrorModel {
    
  public:
    
   typedef std::vector<double> Model;
   typedef std::vector<uint64_t> Matches;
   typedef struct{
	   uint64_t m;
	   uint64_t mm;
	   uint8_t qemp;
   } SMatches;


   typedef std::tr1::unordered_map<uint64_t,HashErrorModel::SMatches> HashMatch;
   typedef std::tr1::unordered_map<std::string,uint16_t > RGdict;
   typedef std::tr1::unordered_map<uint16_t,std::string > ReRGdict;


   Model model;
   RGdict readgroupDic;
   ReRGdict rereadgroupDic;
   HashMatch mismatchTable;
   uint16_t lastElement;

   HashErrorModel();
   ~HashErrorModel();

   void setCell(baseCV basecv);
   uint8_t mapReadGroup(std::string readgroup);
   uint64_t calcKey(baseCV basecv);
   uint8_t getQemp(baseCV basecv);
   int writeAllTableMM(String filename);
   baseCV revKey(uint64_t key);
   int getKeyLength();
   uint32_t getSize();
   void setModel(Model model);
   void setDataforPrediction(Matrix & X, Vector & succ, Vector& total,bool binarizeFlag);
   void addPrediction(Model model, int blendedWeight);
};

#endif
