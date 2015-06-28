#include "mex.h"
#include "mxutils.c"
#include <math.h>

#define ELPX prhs[0]
#define SUM_XX prhs[1]
#define SUM_X  prhs[2]
#define PREC prhs[3]
#define MU_IN prhs[4]
#define LOG_PAR prhs[5]
#define NC prhs[6]

void mexFunction( int nlhs, mxArray *plhs[], 
                          int nrhs, const mxArray*prhs[] )
     
{
    double *E_log_p_of_x, *P, *sum_xx, *sum_x, *mu, *Nc;
    double log_par;
        
    int mrows, ncols, ii, jj;
    
    mrows = mxGetM(SUM_XX);
    ncols = mxGetN(SUM_XX);
       
    /*E_log_p_of_x is altered in-place (pass by reference)*/
    E_log_p_of_x = mxGetPr(ELPX);
    sum_xx       = mxGetPr(SUM_XX);
    sum_x        = mxGetPr(SUM_X);
    mu           = mxGetPr(MU_IN);
    P            = mxGetPr(PREC);
    log_par      = mxGetScalar(LOG_PAR);
    Nc           = mxGetPr(NC);
    
    for (jj = 0; jj < ncols; jj++)
    {
        E_log_p_of_x[jj] = log_par;
        for (ii = 0; ii < mrows; ii++)
        {
            
            E_log_p_of_x[jj] -= P[ii]
                              *((sum_xx[jj*mrows+ii] - 2*mu[ii]*sum_x[jj*mrows+ii])/Nc[jj] 
                                + mu[ii]*mu[ii]);
        }
    }
    
}

