#include "mex.h"
#include "mxutils.c"
#include <math.h>

#define X_IN prhs[0]
#define Q_OF_Z prhs[1]
#define B_OUT plhs[0]


void mexFunction( int nlhs, mxArray *plhs[], 
                          int nrhs, const mxArray*prhs[] )
     
{
    double *B, *X, *q_of_z;
    double temp;
    int mrows, ncols, ii, jj, kk;
    
    mrows = mxGetM(X_IN);
    ncols = mxGetN(X_IN);
       
    X = mxGetPr(X_IN);
    
    q_of_z = mxGetPr(Q_OF_Z);
            
    B_OUT = mxCreateDoubleMatrix(mrows,mrows,mxREAL);    
    B  = mxGetPr(B_OUT);
      
    for (jj = 0; jj < mrows; jj++)
    {
        for (ii = 0; ii < mrows; ii++)
        {
            B[jj*mrows+ii] = 0;
        }
    }
    
    /*traverse columnwise to prevent excessive cache hits*/
    for (kk = 0; kk < ncols; kk++)
    {
        for (jj = 0; jj < mrows; jj++)
        {
            for (ii = 0; ii < mrows; ii++)
            {
                B[jj*mrows+ii] += q_of_z[kk]*X[kk*mrows+ii]*X[kk*mrows+jj];
                
            }
        }
    }
}

