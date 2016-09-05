/* M_times_w_sub.c

   ----------------------------------------------------------------------
   This file is part of LiFE toolbox

   Copyright (C) 2015 Cesar Caiafa & Franco Pestilli
   ----------------------------------------------------------------------
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>     /* provides DBL_EPSILON */
#include <sys/types.h>

#include "M_times_w.h"

void M_times_w_sub( double YPtr[], double atomsPtr[], double voxelsPtr[], double fibersPtr[], double valuesPtr[], double DPtr[], double wPtr[], int nTheta, int nVoxels, int nCoeffs )
{   int
        k, i, atom_index, voxel_index;
    
    double
        val;
    
    for (k = 0; k < nCoeffs; k++) 
    {
        
        atom_index = (int)(atomsPtr[k]-1)*nTheta;
        voxel_index = (int)(voxelsPtr[k]-1)*nTheta;

        for (i = 0; i < nTheta; i++)
        {
            YPtr[voxel_index] = YPtr[voxel_index] + DPtr[atom_index]*wPtr[(int)fibersPtr[k]-1]*valuesPtr[k];
            atom_index++;
            voxel_index++;
        }
            
    }
    
    return;

}



