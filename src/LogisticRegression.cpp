//////////////////////////////////////////////////////////////////////
// mach2dat/LogisticRegression.h
// (c) 2008 Yun Li
//
// March 15, 2008
//

#include "LogisticRegression.h"
#include "StringHash.h"
#include <math.h>

LogisticRegression::LogisticRegression()
{
}

LogisticRegression::~LogisticRegression()
{
}

double LogisticRegression::GetDeviance(Matrix & X, Vector & y)
{
    double ll = 0.0;
    for (int i = 0; i < X.rows; i++)
    {
        double t = 0.0;
        for (int j = 0; j < X.cols; j++)
            t += B[j] * X[i][j];
        double yhat = 1 / (1 + exp(-t));

        ll += y[i] == 1 ? log (yhat) : log(1 - yhat);
    }

    double deviance = -2.0 * ll;
    return deviance;
}

double LogisticRegression::GetDeviance(Matrix & X, Vector & succ, Vector& total)
{
    double ll = 0.0;
    for (int i = 0; i < X.rows; i++)
    {
        double t = 0.0;
        for (int j = 0; j < X.cols; j++)
            t += B[j] * X[i][j];
        double yhat = 1 / (1 + exp(-t));

        ll += succ[i] * log(yhat) + (total[i] - succ[i]) * log(1-yhat);
    }

    double deviance = -2.0 * ll;
    return deviance;
}

void LogisticRegression::reset(Matrix& X){
    B.Dimension(X.cols);
    B.Zero();

    covB.Dimension(X.cols, X.cols);
    covB.Zero();

    p.Dimension(X.rows);
    V.Dimension(X.rows);
    W.Dimension(X.rows);
    p.Zero();
    V.Zero();
    W.Zero();

    residuals.Dimension(X.rows);

    deltaB.Dimension(X.cols);

    D.Dimension(X.cols, X.cols);

    Dinv.Dimension(X.cols, X.cols);
    Dtwo.Dimension(X.cols, X.rows);
    XtV.Dimension(X.cols, X.rows);

}

bool LogisticRegression::FitLogisticModel(Matrix & X, Vector & succ, Vector& total, int nrrounds)
{
    this-> reset(X);
    int rounds = 0;

    // Newton-Raphson
    while (rounds < nrrounds)
    {
        // beta = beta + solve( t(X)%*%diag(p*(1-p)) %*%X) %*% t(X) %*% (Y-p);
        for (int i = 0; i < X.rows; i++)
        {
            double d = 0;
            for (int k = 0; k < X.cols; k ++)
            {
                d += B[k] * X[i][k]; // \eta,
            } // for parameter beta-k

            p[i] = 1.0 / (1.0 + exp(-d)); // \mu = E(prob)
            V[i] = p[i] * (1 - p[i]); 
            W[i] = total[i] * V[i]; // weight
        } // for observation i

        // The first part: solve / inverse
        D.Zero();

        for (int k = 0; k < X.cols; k++)
        {
            for (int l = k; l < X.cols; l++)
            {
                double Dentry = 0.0;
                for (int i = 0; i < X.rows; i++)
                    Dentry += X[i][k] * W[i] * X[i][l];
                D[k][l] = D[l][k] = Dentry;
            }
        }

        Dinv.Zero();
        Dinv = D;
        // SVD svd;
        // svd.SVDInverse(Dinv);
        if (!chol.TryDecompose(D))
            return false;
        chol.Decompose(D);
        chol.Invert();
        Dinv = chol.inv; // (X' W X)^{-1}

        // Accumulate up to the second part: multiply by t(X)
        Dtwo.Zero();

        for (int k = 0; k < X.cols; k++)
            for (int i = 0; i < X.rows; i++)
                for (int l = 0; l < X.cols; l++)
                    Dtwo[k][i] += Dinv[k][l] * X[i][l];

        // The last part, residuals
        residuals.Zero();
        for (int i = 0; i < X.rows; i++)
            residuals[i] = succ[i] - total[i] * p[i];

        deltaB.Zero();
        for (int k = 0; k < X.cols; k++)
            for (int i = 0; i < X.rows; i++)
                deltaB[k] += Dtwo[k][i] * residuals[i];

        // update beta's and check for convergence
        double delta = 0.0;
        for (int k = 0; k < X.cols; k++)
        {
            delta += fabs(deltaB[k]);
            B[k] += deltaB[k];
        }

        //printf ("%d %-11.4g\n", rounds, delta);

        if (delta < 1e-6)
        {
            rounds = 0;
            break;
        }

        rounds ++;
    } // Newton-Raphson iterations

    if (rounds == nrrounds)
        return false;
    // add 10 more rounds before givin up

    // obtain covariance matrix to perform Wald test
    // covB = solve(t(X)%*%V%*%X)

    // Transpose X and multiply by diagonal V
    //XtV.Zero();
    for (int i = 0; i < X.rows; i++)
        for (int k = 0; k < X.cols; k++)
            XtV[k][i] = X[i][k] * W[i];

    covB.Product(XtV, X);

    if (!chol.TryDecompose(covB))
        return false;
    chol.Decompose(covB);
    chol.Invert();
    covB = chol.inv;
    return true;
}
