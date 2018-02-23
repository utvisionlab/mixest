/* Function for calculating the square root of upper triangular matrix
 * Written by Reshad Hosseini*/
#include "mex.h"
#include <math.h> //needed for sqrt
//#include <string.h> //needed for memcpy

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    int m, n, i, j, k;
    double *Tr, *Ti, *Si, *Sr;
    double numreal, denumreal, numimag, denumimag, sreal, simag, absval;
    //double complex num;
    
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
    Ti = mxGetPi(prhs[0]);
    
    /* Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    
     /* Create a C pointer to a copy of the output matrix. */
    Sr = mxGetPr(plhs[0]);
    Si = mxGetPi(plhs[0]);
    
    /* copy T into Ts to speed up memory access */
    //memcpy(Ts, T, m*n*sizeof(double));
    for(j=0;j<n;j++) {
        for(i=j-1;i>=0;i--) {
            *(Sr + i + m*j) = *(Tr + i + m*j);
            *(Si + i + m*j) = *(Ti + i + m*j);
        }
    }

    /* Upper triangular */
    for(j=0;j<n;j++) {
        absval = *(Tr + j + m*j) * *(Tr + j + m*j) + 
                *(Ti + j + m*j) * *(Ti + j + m*j);
        absval = sqrt(absval);
        *(Sr + j + m*j) = sqrt( (*(Tr + j + m*j) + absval)/2 );
        if(*(Ti + j + m*j) > 0)
            *(Si + j + m*j) = sqrt( (absval - *(Tr + j + m*j))/2 );
        else
            *(Si + j + m*j) = - sqrt( (absval - *(Tr + j + m*j))/2 );
        // num = csqrt(*(Tr + j + m*j) + *(Ti + j + m*j) *I);
        // *(Sr + j + m*j) = creal(num);
        // *(Si + j + m*j) = cimag(num);
    }
    
    for(j=1;j<n;j++) {
        for(i=j-1;i>=0;i--) {
            sreal = 0;
            simag = 0;
            for(k=i+1;k<j;k++) {
                sreal += *(Sr + i + m*k) * *(Sr + k + m*j);
                sreal -= *(Si + i + m*k) * *(Si + k + m*j);
                simag += *(Sr + i + m*k) * *(Si + k + m*j);
                simag += *(Si + i + m*k) * *(Sr + k + m*j);
            }
            
            numreal = *(Tr+i+m*j) - sreal;
            numimag = *(Ti+i+m*j) - simag;
            denumreal = *(Sr+i+m*i) + *(Sr+j+m*j);
            denumimag = *(Si+i+m*i) + *(Si+j+m*j);
            
            absval = denumreal*denumreal + denumimag*denumimag;
            
            //*(Sr+i+m*j) = (numreal*denumreal + numimag*denumimag)/absval;
            //*(Si+i+m*j) = (numimag*denumreal - numreal*denumimag)/absval;
            
            // Instead of the above line, use following to solve 0/0
            sreal = numreal*denumreal + numimag*denumimag;
            simag = numimag*denumreal - numreal*denumimag;
            if(sreal == 0)
                *(Sr+i+m*j) =  0;
            else
                *(Sr+i+m*j) = sreal / absval;
            if(simag == 0)
                *(Si+i+m*j) =  0;
            else
                *(Si+i+m*j) = simag / absval;
            
        }
    }
}