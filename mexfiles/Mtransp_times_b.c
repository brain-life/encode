/* Mtransp_times_b.c

   ----------------------------------------------------------------------
   This file is part of LiFE toolbox

   Copyright (C) 2015 Cesar Caiafa & Franco Pestilli
   ----------------------------------------------------------------------
*/
#include "Mtransp_times_b.h"
#include "mex.h"


/* ----------------------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* ----------------------------------------------------------------------- */
{  const mxArray *atoms;
   const mxArray *voxels;
   const mxArray *fibers;
   const mxArray *values;
   const mxArray *D;
   const mxArray *Y;
   
   mxArray *w;
   
   
      
   int
       nFibers, nTheta, nCoeffs, k, i;
   
   /*unsigned int   dims[2]; */
     
   /* Free memory and exit if no parameters are given */
   if (nrhs == 0)
   {  if (nlhs != 0) { mexErrMsgTxt("No output arguments expected."); }
      return ;
   }

   /* Check for proper number of arguments */
   if (nrhs != 7)  { mexErrMsgTxt("Seven input arguments required."); }
   if ((nlhs  > 1)               ) { mexErrMsgTxt("Too many output arguments.");             }

    /* Extract the arguments */
   atoms = (mxArray *)prhs[0];
   voxels = (mxArray *)prhs[1];
   fibers = (mxArray *)prhs[2];
   values = (mxArray *)prhs[3];
   D = (mxArray *)prhs[4];
   Y = (mxArray *)prhs[5];
   
   nFibers = (int) mxGetScalar(prhs[6]);
   
   w = mxCreateNumericMatrix(nFibers, 1, mxDOUBLE_CLASS, mxREAL);
   
    /* Verify validity of input arguments */
    if (mxIsEmpty(atoms) || mxIsEmpty(voxels) || mxIsEmpty(fibers) || mxIsEmpty(values) || mxIsEmpty(D) || mxIsEmpty(Y))
    {   /* If arguments are empty, simply return an empty weights */
        plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        mexWarnMsgTxt("Returning empty weights.");
        return ;       
    }

    /* Verify validity of input argument 'atoms' */
    if (atoms != NULL)
    {  if (!mxIsDouble(atoms)  || ((mxGetM(atoms) > 1) &&
          (mxGetN(atoms) > 1)) || (mxGetNumberOfDimensions(atoms) != 2))
       {   mexErrMsgTxt("Parameter 'atoms' has to be a double vector.");
       }
    }
   
    /* Verify validity of input argument 'voxels' */
    if (voxels != NULL)
    {  if (!mxIsDouble(voxels)  || ((mxGetM(voxels) > 1) &&
          (mxGetN(voxels) > 1)) || (mxGetNumberOfDimensions(voxels) != 2))
       {   mexErrMsgTxt("Parameter 'voxels' has to be a double vector.");
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

    /* Verify validity of input argument 'Y' */
    if (Y != NULL)
    {  if (!mxIsDouble(Y)  || ((mxGetM(Y) > 1) &&
          (mxGetN(Y) < 2)) || (mxGetNumberOfDimensions(Y) != 2))
       {   mexErrMsgTxt("Parameter 'Y' has to be a double matrix.");
       }
    }   
   
   /* Verify all vectors have the same lenght */
   if ((mxGetNumberOfElements(atoms) != mxGetNumberOfElements(voxels)) || (mxGetNumberOfElements(atoms) != mxGetNumberOfElements(fibers)) || (mxGetNumberOfElements(atoms) != mxGetNumberOfElements(values)))
       {   mexErrMsgTxt("Parameters 'atoms', 'voxels', 'fibers' and 'values' have to be of equal length.");
       }
   
   /* Verify matrices D and Y have the same number of rows */
   if (mxGetM(D) != mxGetM(Y))
       {   mexErrMsgTxt("Matrices 'D' and 'Y' have to have equal number of rows.");
       }


   /* Call the core subroutine */  
   nTheta = mxGetM(Y);
   nCoeffs = mxGetM(atoms);
   Mtransp_times_b_sub( mxGetPr(w), mxGetPr(atoms), mxGetPr(voxels), mxGetPr(fibers), mxGetPr(values), mxGetPr(D), mxGetPr(Y), nFibers, nTheta, nCoeffs );
   plhs[0] = w;

    return ;
}
