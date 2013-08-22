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
#ifndef __PREDICTION_H__
#define __PREDICTION_H__

#include "LogisticRegression.h"
#include "IntArray.h"
#include "Error.h"
#include "StringHash.h"
#include "Parameters.h"
#include "MemoryAllocators.h"
#include "HashErrorModel.h"
#include <vector>



class Prediction {
  public:
   LogisticRegression lrengine; // our logistic regression engine
   bool fullmodel;
   HashErrorModel *phasherrormodel;

   // These are the count vectors
   Vector succ;
   Vector total;
   // This is the full design matrix
   // # parameters: slope + # co-variates;
   Matrix X;

   Prediction(HashErrorModel *phasherrormodel);
   Prediction();
   ~Prediction();
   int fitModel(bool writeModelFlag, std::string& filename);
   void outModel();
   void setErrorModel(HashErrorModel *phasherrormodel);
   std::vector<double> getModel();
   int writeLogRegdata(std::string& filename);
};

#endif
