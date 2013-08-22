////////////////////////////////////////////////////////////////////// 
// mach2dat/LogisticRegression.h
// (c) 2008 Yun Li
//
// March 15, 2008
//

#ifndef __LOGISTIC_REGRESSION_H__
#define __LOGISTIC_REGRESSION_H__

#include "MathMatrix.h"
#include "MathCholesky.h"
#include "StringHash.h"
#include "StringArray.h"

class Pedigree;

class LogisticRegression
{
public:
    LogisticRegression();
    ~LogisticRegression();

    bool FitLogisticModel(Matrix & X, Vector & y, int rnrounds); // return false if not converging
    bool FitLogisticModel(Matrix & X, Vector & succ, Vector& total, int nrrounds);
    double GetDeviance(Matrix & X, Vector & y);
    double GetDeviance(Matrix & X, Vector & succ, Vector& total);
    void reset(Matrix& X); // get everything cleared

    Vector B; // coefficient vector
    Matrix covB; // coefficient covariance matrix
private:
	Vector p, V, W;
    Vector residuals;
    Vector deltaB;
    Matrix D;
    Matrix Dinv;
    Cholesky chol;
    Matrix Dtwo;
	Matrix XtV;
};

#endif



