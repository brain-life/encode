/*
 * Mtransp_times_b.h
 *
 * Code generation for function 'Mtransp_times_b'
 *
 */

#ifndef __MTRANSP_TIMES_B_H__
#define __MTRANSP_TIMES_B_H__

/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "Mtransp_times_b_types.h"

/* Function Declarations */
extern void Mtransp_times_b(const emlrtStack *sp, const emxArray_real_T *atoms,
  const emxArray_real_T *voxels, const emxArray_real_T *fibers, const
  emxArray_real_T *values, const emxArray_real_T *D, const emxArray_real_T *Y,
  real_T nFibers, emxArray_real_T *w);

#endif

/* End of code generation (Mtransp_times_b.h) */
