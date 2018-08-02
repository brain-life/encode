/* compute_diag_sub.c

   ----------------------------------------------------------------------
   This file is part of LiFE toolbox

   Copyright (C) 2016 Cesar Caiafa & Franco Pestilli
   ----------------------------------------------------------------------
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>     /* provides DBL_EPSILON */
#include <sys/types.h>

#include "compute_diag.h"

void compute_diag_sub( double dPtr[], double atomsPtr[], double fibersPtr[], double valuesPtr[], double DPtr[] , int nFibers, int nTheta, int nCoeffs )
{   int
        k, i, atom_index;
    
    double
        val;
    
    for (k = 0; k < nCoeffs; k++)
    {
        val = 0;
        atom_index = (int)(atomsPtr[k]-1)*nTheta;
        
        for (i = 0; i < nTheta; i++)
        {
            val = val + DPtr[atom_index]*DPtr[atom_index];
            atom_index++;
        }
        val = val*valuesPtr[k]*valuesPtr[k];
        
        dPtr[(int)fibersPtr[k]-1] = dPtr[(int)fibersPtr[k]-1] + val;
            
    }
    
    return;


}