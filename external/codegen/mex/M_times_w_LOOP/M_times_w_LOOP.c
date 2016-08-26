/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * M_times_w_LOOP.c
 *
 * Code generation for function 'M_times_w_LOOP'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "M_times_w_LOOP.h"
#include "M_times_w_LOOP_emxutil.h"
#include "M_times_w_LOOP_data.h"

/* Variable Definitions */
static emlrtRTEInfo emlrtRTEI = { 1, 18, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m" };

static emlrtRTEInfo b_emlrtRTEI = { 1, 12, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m" };

static emlrtDCInfo emlrtDCI = { 5, 26, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 1 };

static emlrtBCInfo emlrtBCI = { -1, -1, 5, 26, "Y", "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 0 };

static emlrtDCInfo b_emlrtDCI = { 5, 43, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 1 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 5, 43, "D", "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 0 };

static emlrtDCInfo c_emlrtDCI = { 5, 53, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 1 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 5, 53, "w", "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 0 };

static emlrtBCInfo d_emlrtBCI = { -1, -1, 5, 66, "values", "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 0 };

static emlrtECInfo emlrtECI = { -1, 5, 22, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m" };

static emlrtDCInfo d_emlrtDCI = { 5, 9, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 1 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 5, 9, "Y", "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 0 };

static emlrtECInfo b_emlrtECI = { -1, 5, 5, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m" };

static emlrtDCInfo e_emlrtDCI = { 2, 11, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 1 };

static emlrtDCInfo f_emlrtDCI = { 2, 11, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 4 };

static emlrtDCInfo g_emlrtDCI = { 2, 18, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 1 };

static emlrtDCInfo h_emlrtDCI = { 2, 18, "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 4 };

static emlrtBCInfo f_emlrtBCI = { -1, -1, 5, 26, "voxels", "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 5, 43, "atoms", "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 0 };

static emlrtBCInfo h_emlrtBCI = { -1, -1, 5, 55, "fibers", "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 0 };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 5, 9, "voxels", "M_times_w_LOOP",
  "/N/dc2/projects/lifebid/code/lifebid/external/M_times_w_LOOP.m", 0 };

/* Function Definitions */
void M_times_w_LOOP(const emlrtStack *sp, const emxArray_real_T *atoms, const
                    emxArray_real_T *voxels, const emxArray_real_T *fibers,
                    const emxArray_real_T *values, const emxArray_real_T *D,
                    const emxArray_real_T *w, real_T nTheta, real_T nVoxels,
                    emxArray_real_T *Y)
{
  emxArray_real_T *b_Y;
  int32_T i0;
  real_T b_w;
  real_T b_values;
  real_T c_w;
  real_T c_values;
  int32_T loop_ub;
  int32_T k;
  emxArray_real_T *r0;
  emxArray_int32_T *r1;
  emxArray_real_T *c_Y;
  emxArray_real_T *d_Y;
  int32_T b_atoms;
  int32_T c_atoms;
  int32_T d_atoms;
  int32_T e_atoms;
  int32_T f_atoms;
  int32_T iv0[1];
  int32_T b_voxels;
  int32_T e_Y[1];
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &b_Y, 2, &b_emlrtRTEI, true);
  i0 = b_Y->size[0] * b_Y->size[1];
  if (nTheta >= 0.0) {
    b_w = nTheta;
  } else {
    b_w = emlrtNonNegativeCheckR2012b(nTheta, &f_emlrtDCI, sp);
  }

  if (b_w == (int32_T)muDoubleScalarFloor(b_w)) {
    b_Y->size[0] = (int32_T)b_w;
  } else {
    b_Y->size[0] = (int32_T)emlrtIntegerCheckR2012b(b_w, &e_emlrtDCI, sp);
  }

  if (nVoxels >= 0.0) {
    b_w = nVoxels;
  } else {
    b_w = emlrtNonNegativeCheckR2012b(nVoxels, &h_emlrtDCI, sp);
  }

  if (b_w == (int32_T)muDoubleScalarFloor(b_w)) {
    b_Y->size[1] = (int32_T)b_w;
  } else {
    b_Y->size[1] = (int32_T)emlrtIntegerCheckR2012b(b_w, &g_emlrtDCI, sp);
  }

  emxEnsureCapacity(sp, (emxArray__common *)b_Y, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  if (nTheta >= 0.0) {
    b_w = nTheta;
  } else {
    b_w = emlrtNonNegativeCheckR2012b(nTheta, &f_emlrtDCI, sp);
  }

  if (nVoxels >= 0.0) {
    b_values = nVoxels;
  } else {
    b_values = emlrtNonNegativeCheckR2012b(nVoxels, &h_emlrtDCI, sp);
  }

  if (b_w == (int32_T)muDoubleScalarFloor(b_w)) {
    c_w = b_w;
  } else {
    c_w = emlrtIntegerCheckR2012b(b_w, &e_emlrtDCI, sp);
  }

  if (b_values == (int32_T)muDoubleScalarFloor(b_values)) {
    c_values = b_values;
  } else {
    c_values = emlrtIntegerCheckR2012b(b_values, &g_emlrtDCI, sp);
  }

  loop_ub = (int32_T)c_w * (int32_T)c_values;
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_Y->data[i0] = 0.0;
  }

  k = 0;
  b_emxInit_real_T(sp, &r0, 1, &emlrtRTEI, true);
  emxInit_int32_T(sp, &r1, 1, &emlrtRTEI, true);
  b_emxInit_real_T(sp, &c_Y, 1, &emlrtRTEI, true);
  b_emxInit_real_T(sp, &d_Y, 1, &emlrtRTEI, true);
  while (k <= values->size[0] - 1) {
    /* k */
    i0 = voxels->size[0];
    b_atoms = 1 + k;
    if ((b_atoms >= 1) && (b_atoms < i0)) {
      c_atoms = b_atoms;
    } else {
      c_atoms = emlrtDynamicBoundsCheckR2012b(b_atoms, 1, i0, &f_emlrtBCI, sp);
    }

    b_w = voxels->data[c_atoms - 1];
    i0 = b_Y->size[1];
    if (b_w == (int32_T)muDoubleScalarFloor(b_w)) {
      b_atoms = (int32_T)b_w;
    } else {
      b_atoms = (int32_T)emlrtIntegerCheckR2012b(b_w, &emlrtDCI, sp);
    }

    emlrtDynamicBoundsCheckR2012b(b_atoms, 1, i0, &emlrtBCI, sp);
    i0 = atoms->size[0];
    b_atoms = 1 + k;
    if ((b_atoms >= 1) && (b_atoms < i0)) {
      d_atoms = b_atoms;
    } else {
      d_atoms = emlrtDynamicBoundsCheckR2012b(b_atoms, 1, i0, &g_emlrtBCI, sp);
    }

    b_w = atoms->data[d_atoms - 1];
    i0 = D->size[1];
    if (b_w == (int32_T)muDoubleScalarFloor(b_w)) {
      b_atoms = (int32_T)b_w;
    } else {
      b_atoms = (int32_T)emlrtIntegerCheckR2012b(b_w, &b_emlrtDCI, sp);
    }

    emlrtDynamicBoundsCheckR2012b(b_atoms, 1, i0, &b_emlrtBCI, sp);
    i0 = fibers->size[0];
    b_atoms = 1 + k;
    if ((b_atoms >= 1) && (b_atoms < i0)) {
      e_atoms = b_atoms;
    } else {
      e_atoms = emlrtDynamicBoundsCheckR2012b(b_atoms, 1, i0, &h_emlrtBCI, sp);
    }

    b_w = fibers->data[e_atoms - 1];
    i0 = w->size[0];
    if (b_w == (int32_T)muDoubleScalarFloor(b_w)) {
      b_atoms = (int32_T)b_w;
    } else {
      b_atoms = (int32_T)emlrtIntegerCheckR2012b(b_w, &c_emlrtDCI, sp);
    }

    emlrtDynamicBoundsCheckR2012b(b_atoms, 1, i0, &c_emlrtBCI, sp);
    i0 = values->size[0];
    b_atoms = k + 1;
    emlrtDynamicBoundsCheckR2012b(b_atoms, 1, i0, &d_emlrtBCI, sp);
    loop_ub = D->size[0];
    b_atoms = (int32_T)atoms->data[k];
    b_w = w->data[(int32_T)fibers->data[k] - 1];
    b_values = values->data[k];
    i0 = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r0, i0, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i0 = 0; i0 < loop_ub; i0++) {
      r0->data[i0] = D->data[i0 + D->size[0] * (b_atoms - 1)] * b_w * b_values;
    }

    i0 = b_Y->size[0];
    b_atoms = r0->size[0];
    if (i0 != b_atoms) {
      emlrtSizeEqCheck1DR2012b(i0, b_atoms, &emlrtECI, sp);
    }

    loop_ub = b_Y->size[0];
    i0 = r1->size[0];
    r1->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r1, i0, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i0 = 0; i0 < loop_ub; i0++) {
      r1->data[i0] = i0;
    }

    i0 = voxels->size[0];
    b_atoms = 1 + k;
    if ((b_atoms >= 1) && (b_atoms < i0)) {
      f_atoms = b_atoms;
    } else {
      f_atoms = emlrtDynamicBoundsCheckR2012b(b_atoms, 1, i0, &i_emlrtBCI, sp);
    }

    b_w = voxels->data[f_atoms - 1];
    i0 = b_Y->size[1];
    if (b_w == (int32_T)muDoubleScalarFloor(b_w)) {
      b_atoms = (int32_T)b_w;
    } else {
      b_atoms = (int32_T)emlrtIntegerCheckR2012b(b_w, &d_emlrtDCI, sp);
    }

    emlrtDynamicBoundsCheckR2012b(b_atoms, 1, i0, &e_emlrtBCI, sp);
    iv0[0] = r1->size[0];
    loop_ub = b_Y->size[0];
    b_voxels = (int32_T)voxels->data[k];
    i0 = c_Y->size[0];
    c_Y->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)c_Y, i0, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_Y->data[i0] = b_Y->data[i0 + b_Y->size[0] * (b_voxels - 1)];
    }

    e_Y[0] = c_Y->size[0];
    emlrtSubAssignSizeCheckR2012b(iv0, 1, e_Y, 1, &b_emlrtECI, sp);
    b_voxels = (int32_T)voxels->data[k];
    b_atoms = b_Y->size[0];
    loop_ub = (int32_T)voxels->data[k];
    i0 = d_Y->size[0];
    d_Y->size[0] = b_atoms;
    emxEnsureCapacity(sp, (emxArray__common *)d_Y, i0, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i0 = 0; i0 < b_atoms; i0++) {
      d_Y->data[i0] = b_Y->data[i0 + b_Y->size[0] * (loop_ub - 1)] + r0->data[i0];
    }

    loop_ub = d_Y->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_Y->data[r1->data[i0] + b_Y->size[0] * (b_voxels - 1)] = d_Y->data[i0];
    }

    k++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_real_T(&d_Y);
  emxFree_real_T(&c_Y);
  emxFree_int32_T(&r1);
  emxFree_real_T(&r0);
  i0 = Y->size[0];
  Y->size[0] = b_Y->size[0] * b_Y->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)Y, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = b_Y->size[0] * b_Y->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    Y->data[i0] = b_Y->data[i0];
  }

  emxFree_real_T(&b_Y);

  /*  vectorization */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (M_times_w_LOOP.c) */
