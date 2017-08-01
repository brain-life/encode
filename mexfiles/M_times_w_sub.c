/* M_times_w_sub.c

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

#include "M_times_w.h"

void M_times_w_sub( double YPtr[], double atomsPtr[], double voxelsPtr[], double fibersPtr[], double valuesPtr[], double DPtr[], double wPtr[], int nTheta, int nVoxels, int nCoeffs )
{   int
        k, i, atom_index, voxel_index, thnum, maxnum;
    
    double
        val;
/*    printf("Number of iterations: %i\n", nCoeffs); */
#pragma omp parallel for private(i,k,atom_index,voxel_index) firstprivate(YPtr)
    for (k = 0; k < nCoeffs; k++) 
    {
        atom_index = (int)(atomsPtr[k]-1)*nTheta;
        voxel_index = (int)(voxelsPtr[k]-1)*nTheta;

/*	if(k % 100005 == 0){
	        thnum = omp_get_thread_num();
		maxnum = omp_get_max_threads();
		printf("Iteration: %i Thread Num: %i Thread max: %i\n", k, thnum, maxnum);
	}*/

        for (i = 0; i < nTheta; i++)
        {
            YPtr[i+voxel_index] = YPtr[i+voxel_index] + DPtr[i+atom_index]*wPtr[(int)fibersPtr[k]-1]*valuesPtr[k];
        }
            
    }
    
    return;

}

