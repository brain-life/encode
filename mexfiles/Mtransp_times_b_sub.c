/* Mtransp_times_b_sub.c

   ----------------------------------------------------------------------
   This file is part of LiFE toolbox

   Copyright (C) 2015 Cesar Caiafa & Franco Pestilli
   ----------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>     /* provides DBL_EPSILON */
#include <sys/types.h>
#include <omp.h>

#include "Mtransp_times_b.h"

void Mtransp_times_b_sub( double wPtr[], double atomsPtr[], double voxelsPtr[], double fibersPtr[], double valuesPtr[], double DPtr[], double YPtr[], int nFibers, int nTheta, int nCoeffs )
{   int
        k, i, atom_index, voxel_index;
    
    double
        val;
#pragma omp parallel for private(i,k,atom_index,voxel_index) firstprivate(wPtr)    
    for (k = 0; k < nCoeffs; k++) 
    {
        val = 0;
        atom_index = (int)(atomsPtr[k]-1)*nTheta;
        voxel_index = (int)(voxelsPtr[k]-1)*nTheta;
        
        for (i = 0; i < nTheta; i++)
        {
            val = val + DPtr[i+atom_index]*YPtr[i+voxel_index];
        }
        val = val*valuesPtr[k];
        wPtr[(int)fibersPtr[k]-1] = wPtr[(int)fibersPtr[k]-1] + val;
            
    }
    
    return;


}
