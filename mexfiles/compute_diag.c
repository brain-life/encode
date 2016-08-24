/* compute_diag.c

   ----------------------------------------------------------------------
   This file is part of LiFE toolbox

   Copyright (C) 2016 Cesar Caiafa & Franco Pestilli
   ----------------------------------------------------------------------
*/
#include "compute_diag.h"
#include "mex.h"


/* ----------------------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* ----------------------------------------------------------------------- */
{  const mxArray *atoms;
   const mxArray *fibers;
   const mxArray *values;
   const mxArray *D;
   
   mxArray *d;
   
   int
        nFibers, nTheta, nCoeffs, k, i;
   
   /*unsigned int   dims[2]; */
     
   /* Free memory and exit if no parameters are given */
   if (nrhs == 0)
   {  if (nlhs != 0) { mexErrMsgTxt("No output arguments expected."); }
      return ;
   }

   /* Check for proper number of arguments */
   if (nrhs != 5)  { mexErrMsgTxt("Five input arguments required."); }
   if ((nlhs  > 1)               ) { mexErrMsgTxt("Too many output arguments.");             }

    /* Extract the arguments */
   atoms = (mxArray *)prhs[0];
   fibers = (mxArray *)prhs[1];
   values = (mxArray *)prhs[2];
   D = (mxArray *)prhs[3];
   
   nFibers = (int) mxGetScalar(prhs[4]);
   
   d = mxCreateNumericMatrix(nFibers, 1, mxDOUBLE_CLASS, mxREAL);
   
    /* Verify validity of input arguments */
    if (mxIsEmpty(atoms) || mxIsEmpty(fibers) || mxIsEmpty(values) || mxIsEmpty(D))
    {   /* If arguments are empty, simply return an empty weights */
        plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        mexWarnMsgTxt("Returning empty d.");
        return ;       
    }

    /* Verify validity of input argument 'atoms' */
    if (atoms != NULL)
    {  if (!mxIsDouble(atoms)  || ((mxGetM(atoms) > 1) &&
          (mxGetN(atoms) > 1)) || (mxGetNumberOfDimensions(atoms) != 2))
       {   mexErrMsgTxt("Parameter 'atoms' has to be a double vector.");
       }
    }
   
    /* Verify validity of input argument 'fibers' */
    if (fibers != NULL)
    {  if (!mxIsDouble(fibers)  || ((mxGetM(fibers) > 1) &&
          (mxGetN(fibers) > 1)) || (mxGetNumberOfDimensions(fibers) != 2))
       {   mexErrMsgTxt("Parameter 'fibers' has to be a double vector.");
       }
    }
      
    /* Verify validity of input argument 'values' */
    if (values != NULL)
    {  if (!mxIsDouble(values)  || ((mxGetM(values) > 1) &&
          (mxGetN(values) > 1)) || (mxGetNumberOfDimensions(values) != 2))
       {   mexErrMsgTxt("Parameter 'values' has to be a double vector.");
       }
    }   
   
    /* Verify validity of input argument 'D' */
    if (D != NULL)
    {  if (!mxIsDouble(D)  || ((mxGetM(D) > 1) &&
          (mxGetN(D) < 2)) || (mxGetNumberOfDimensions(D) != 2))
       {   mexErrMsgTxt("Parameter 'D' has to be a double matrix.");
       }
    }
   
   /* Verify all vectors have the same lenght */
   if  ((mxGetNumberOfElements(atoms) != mxGetNumberOfElements(fibers)) || (mxGetNumberOfElements(atoms) != mxGetNumberOfElements(values)))
       {   mexErrMsgTxt("Parameters 'atoms' ,'fibers' and 'values' have to be of equal length.");
       }


   /* Call the core subroutine */
    nTheta = mxGetM(D);
    nCoeffs = mxGetM(atoms);
   compute_diag_sub( mxGetPr(d), mxGetPr(atoms), mxGetPr(fibers), mxGetPr(values), mxGetPr(D), nFibers, nTheta, nCoeffs );
   plhs[0] = d;

    return ;
}
