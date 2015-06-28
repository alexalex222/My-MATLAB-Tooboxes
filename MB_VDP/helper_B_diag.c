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
    int mrows, ncols, ii, kk;
    
    mrows = mxGetM(X_IN);
    ncols = mxGetN(X_IN);
       
    X = mxGetPr(X_IN);
    
    q_of_z = mxGetPr(Q_OF_Z);
            
    B_OUT = mxCreateDoubleMatrix(mrows,1,mxREAL);    
    B  = mxGetPr(B_OUT);
      
    for (ii = 0; ii < mrows; ii++)
    {
        B[ii] = 0;
    }
    
    /*traverse columnwise to prevent excessive cache hits*/
    for (kk = 0; kk < ncols; kk++)
    {
        for (ii = 0; ii < mrows; ii++)
        {
            temp = (X[kk*mrows+ii]);
            B[ii] += q_of_z[kk]*temp*temp;
        }
    }
    
}

