/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * M_times_w_LOOP.h
 *
 * Code generation for function 'M_times_w_LOOP'
 *
 */

#ifndef __M_TIMES_W_LOOP_H__
#define __M_TIMES_W_LOOP_H__

/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "M_times_w_LOOP_types.h"

/* Function Declarations */
extern void M_times_w_LOOP(const emlrtStack *sp, const emxArray_real_T *atoms,
  const emxArray_real_T *voxels, const emxArray_real_T *fibers, const
  emxArray_real_T *values, const emxArray_real_T *D, const emxArray_real_T *w,
  real_T nTheta, real_T nVoxels, emxArray_real_T *Y);

#endif

/* End of code generation (M_times_w_LOOP.h) */
