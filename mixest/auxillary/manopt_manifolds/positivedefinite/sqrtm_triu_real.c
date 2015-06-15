/* Function for calculating the square root of upper triangular matrix
 * when matrix is real with positive diagonals
 * Written by Reshad Hosseini*/
#include "mex.h"
#include <math.h> //needed for sqrt
//#include <string.h> //needed for memcpy

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    int m, n, i, j, k;
    double *Tr, *Sr;
    double num, denum, s;
    
    if(nrhs != 1 || nlhs > 1)
        mexErrMsgTxt("Usage: Ts = sqrtm_triu(T)");
    
    /* prhs[0] is first argument.
     * mxGetPr returns double*  (data, col-major)
     * mxGetM returns int  (rows)
     * mxGetN returns int  (cols)
     */
    /* m = rows(T) */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if(m != n) mexErrMsgTxt("matrix must be square");

    //if(mxIsSparse(prhs[1])) {
    //    mexErrMsgTxt("Can not handle sparse matrices yet.");
    //}
    
    if(mxGetNumberOfDimensions(prhs[0]) != 2) {
        mexErrMsgTxt("Arguments must be matrices.");
    }
    
    Tr = mxGetPr(prhs[0]);
    
    /* Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    
     /* Create a C pointer to a copy of the output matrix. */
    Sr = mxGetPr(plhs[0]);
    
    /* copy T into Ts to speed up memory access */
    //memcpy(Ts, T, m*n*sizeof(double));
    for(j=0;j<n;j++) {
        for(i=j;i>=0;i--) {
            *(Sr + i + m*j) = *(Tr + i + m*j);
        }
    }

    /* Upper triangular */
    for(j=0;j<n;j++) {
         *(Sr + j + m*j) = sqrt(*(Tr + j + m*j));
    }
    
    for(j=1;j<n;j++) {
        for(i=j-1;i>=0;i--) {
            s = 0;
            for(k=i+1;k<j;k++) {
                s += *(Sr + i + m*k) * *(Sr + k + m*j);
            }
            
            num = (*(Tr+i+m*j) - s);
            denum = (*(Sr+i+m*i) + *(Sr+j+m*j));
            if(num == 0){
                *(Sr+i+m*j) =  0;
            }
            else{
                *(Sr+i+m*j) = num / denum;
            }
           //mexPrintf("DONE3 i=%i,j=%i,num=%2.4f,denum=%2.4f \n",i,j,s,denum);
        }
    }
}