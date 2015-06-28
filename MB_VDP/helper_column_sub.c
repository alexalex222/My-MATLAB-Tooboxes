#include "mex.h"
#include "mxutils.c"
#include <math.h>

#define X_IN prhs[0]
#define LSE prhs[1]

void mexFunction( int nlhs, mxArray *plhs[], 
                          int nrhs, const mxArray*prhs[] )
     
{
    double *x, *lse;
    int mrows, ncols, ii, jj;
    
    mrows = mxGetM(X_IN);
    ncols = mxGetN(X_IN);
       
    x = mxGetPr(X_IN);
    lse = mxGetPr(LSE);
    /*mexPrintf("ncols=%d mrows=%d",ncols,mrows); */
    /*traverse columnwise to prevent excessive cache hits*/
    for (jj = 0; jj < ncols; jj++)
    {
        for (ii = 0; ii < mrows; ii++)
        {
            /*mexPrintf("x=%f lse=%f\n",x[jj*ncols+ii],lse[ii]);*/
            x[jj*mrows+ii] -= lse[ii];
            /*mexPrintf("x_after=%f\n",x[jj*ncols+ii]);*/
        }
    }
    
}

