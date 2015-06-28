#include "mex.h"
#include "mxutils.c"
#include <math.h>

#define X_IN prhs[0]
#define X_MAX prhs[1]
#define VAL_OUT plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], 
                          int nrhs, const mxArray*prhs[] )
     
{
    double *x, *x_max, *val;
    double temp;
    int mrows, ncols, ii, jj;
    
    mrows = mxGetM(X_IN);
    ncols = mxGetN(X_IN);
       
    x = mxGetPr(X_IN);
    x_max = mxGetPr(X_MAX);
    
    VAL_OUT = mxCreateDoubleMatrix(mrows,1,mxREAL);
    val = mxGetPr(VAL_OUT);
         
    for (ii = 0; ii < mrows; ii++)
    {
        val[ii] = 0;
    }
    
    /*traverse columnwise to prevent excessive cache hits*/
    for (jj = 0; jj < ncols; jj++)
    {
        for (ii = 0; ii < mrows; ii++)
        {
            val[ii] += exp(x[jj*mrows+ii] - x_max[ii]);
        }
    }
    
    for (ii = 0; ii < mrows; ii++)
    {
        val[ii] = log(val[ii]) + x_max[ii];
    }
}

