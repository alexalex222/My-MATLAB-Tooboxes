#include "mex.h"
#include "mxutils.c"
#include <math.h>

#define ELPX prhs[0]
#define PREC prhs[1]
#define X_IN prhs[2]
#define MU_IN prhs[3]
#define LOG_PAR prhs[4]

void mexFunction( int nlhs, mxArray *plhs[], 
                          int nrhs, const mxArray*prhs[] )
     
{
    double *E_log_p_of_x, *P, *X, *mu;
    double log_par,temp;
        
    int mrows, ncols, ii, jj;
    
    mrows = mxGetM(X_IN);
    ncols = mxGetN(X_IN);
       
    /*E_log_p_of_x is altered in-place (pass by reference)*/
    E_log_p_of_x = mxGetPr(ELPX);
    X            = mxGetPr(X_IN);
    mu           = mxGetPr(MU_IN);
    P            = mxGetPr(PREC);
    log_par      = mxGetScalar(LOG_PAR);
    
    for (jj = 0; jj < ncols; jj++)
    {
        E_log_p_of_x[jj] = log_par;
        for (ii = 0; ii < mrows; ii++)
        {
            temp = (X[jj*mrows+ii] - mu[ii]);
            E_log_p_of_x[jj] -= P[ii]*temp*temp;
        }
    }
    
}

