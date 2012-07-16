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
#include <map>
#include "StringArray.h"
#include "Covariates.h"
#include "MathMatrix.h"
#include "MathVector.h"

#define MINSIZE 1048576


class HashErrorModel {
    
  public:

    static void setUseLogReg(bool useLogReg) { ourUseLogReg = useLogReg; }

    
   typedef std::vector<double> Model;
   typedef std::vector<uint64_t> Matches;
   typedef struct{
	   uint64_t m;
	   uint64_t mm;
	   uint8_t qemp;
   } SMatches;


   typedef std::map<BaseData,HashErrorModel::SMatches> HashMatch;

   HashMatch mismatchTable;
   uint16_t lastElement;

   HashErrorModel();
   ~HashErrorModel();

   void setCell(const BaseData& data, char refBase);
   uint8_t getQemp(BaseData& data);
    int writeAllTableMM(String filename, const std::vector<std::string>& id2rg);
   uint32_t getSize();
   void setDataforPrediction(Matrix & X, Vector & succ, Vector& total,bool binarizeFlag);
   void addPrediction(Model model, int blendedWeight);

private:
    static bool ourUseLogReg;
};

#endif
