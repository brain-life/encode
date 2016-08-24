/* Mtransp_times_b.h

   ----------------------------------------------------------------------
   This file is part of LiFE toolbox

   Copyright (C) 2015 Cesar Caiafa & Franco Pestilli
   ----------------------------------------------------------------------
*/
#ifndef __MTRANSP_TIMES_B_H__
#define __MTRANSP_TIMES_B_H__

/* The entries in b are non-negative, those in d strictly positive */
void Mtransp_times_b_sub( double wPtr[], double atomsPtr[], double voxelsPtr[], double fibersPtr[], double valuesPtr[], double DPtr[], double YPtr[], int nFibers, int nTheta, int nCoeffs );


#endif