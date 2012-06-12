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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <ctime>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <math.h>

#include "Prediction.h"


Prediction::Prediction(HashErrorModel *phasherrormodel){
   this->phasherrormodel = phasherrormodel;
   Prediction();
}

Prediction::Prediction(){
   fullmodel = false;
}

Prediction::~Prediction(){
}

int Prediction::fitModel(bool writeModelFlag, std::string& filename)
{
  //double rsqcutoff = 0.0001;
  int nrrounds = 30;

  //set Matrix, succ, total
  // !!! all co-variates not 0?
  phasherrormodel->setDataforPrediction(X, succ, total, false);
  if(!(filename.empty()))
      writeLogRegdata(filename);


  fullmodel = lrengine.FitLogisticModel(X, succ,total, nrrounds);

  if(fullmodel)
	return 1;
  else
	return 0;
}

void Prediction::outModel()
{
    if(fullmodel)
        printf("Model\n");
    printf("Effect Variance SE\n");
    for(int i=0;i<lrengine.B.Length();i++)
    {
        double effect = lrengine.B[i];
        double variance = lrengine.covB[i][i] + 1e-30;
        double serror = sqrt(variance);
        printf("%f %f %f \n",
               effect, variance, serror);
    }
}

void Prediction::setErrorModel(HashErrorModel *phasherrormodel)
{
   this->phasherrormodel = phasherrormodel;
}

std::vector<double> Prediction::getModel()
{
  std::vector<double> model;

  if(fullmodel)
    for(int i=0;i<lrengine.B.Length();i++)
      model.push_back(lrengine.B[i]);
  return model;
};

int Prediction::writeLogRegdata(std::string& filename){

	FILE *pFile;
	pFile = fopen(filename.c_str(), "w");
	for (int i = 0; i < X.rows; i++)
	{
		for (int j = 0; j < X.cols; j++)
			fprintf(pFile,"%1.0f ",X[i][j]);
		fprintf(pFile,"- %0.0f",succ[i]);
		fprintf(pFile," %0.0f\n",total[i]);
	}
	fclose(pFile);
	return 1;
}
