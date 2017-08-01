/* M_times_w.h

   ----------------------------------------------------------------------
   This file is part of LiFE toolbox

   Copyright (C) 2015 Cesar Caiafa & Franco Pestilli
   ----------------------------------------------------------------------
*/
#ifndef __M_TIMES_W_H__
#define __M_TIMES_W_H__

/* The entries in b are non-negative, those in d strictly positive */
void M_times_w_sub( double YPtr[] , double atomsPtr[], double voxelsPtr[], double fibersPtr[], double valuesPtr[], double DPtr[], double wPtr[], int nTheta, int nVoxels, int nCoeffs );


#endif