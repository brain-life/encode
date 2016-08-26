/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_M_times_w_LOOP_api.c
 *
 * Code generation for function '_coder_M_times_w_LOOP_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "M_times_w_LOOP.h"
#include "_coder_M_times_w_LOOP_api.h"
#include "M_times_w_LOOP_emxutil.h"
#include "M_times_w_LOOP_data.h"

/* Variable Definitions */
static emlrtRTEInfo c_emlrtRTEI = { 1, 1, "_coder_M_times_w_LOOP_api", "" };

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *D, const
  char_T *identifier, emxArray_real_T *y);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nTheta,
  const char_T *identifier);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *atoms, const
  char_T *identifier, emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  g_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *D, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  d_emlrt_marshallIn(sp, emlrtAlias(D), &thisId, y);
  emlrtDestroyArray(&D);
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  h_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nTheta,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nTheta), &thisId);
  emlrtDestroyArray(&nTheta);
  return y;
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *atoms, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(sp, emlrtAlias(atoms), &thisId, y);
  emlrtDestroyArray(&atoms);
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv1[1] = { 0 };

  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateNumericArray(1, iv1, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)u->data);
  emlrtSetDimensions((mxArray *)m0, u->size, 1);
  emlrtAssign(&y, m0);
  return y;
}

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv2[1];
  boolean_T bv0[1] = { true };

  static const int32_T iv3[1] = { -1 };

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 1U, iv3, &bv0[0],
    iv2);
  ret->size[0] = iv2[0];
  ret->allocatedSize = ret->size[0];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv4[2];
  int32_T i;
  int32_T iv5[2];
  boolean_T bv1[2] = { true, true };

  for (i = 0; i < 2; i++) {
    iv4[i] = -1;
  }

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv4, &bv1[0],
    iv5);
  ret->size[0] = iv5[0];
  ret->size[1] = iv5[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void M_times_w_LOOP_api(const mxArray * const prhs[8], const mxArray *plhs[1])
{
  emxArray_real_T *atoms;
  emxArray_real_T *voxels;
  emxArray_real_T *fibers;
  emxArray_real_T *values;
  emxArray_real_T *D;
  emxArray_real_T *w;
  emxArray_real_T *Y;
  real_T nTheta;
  real_T nVoxels;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  b_emxInit_real_T(&st, &atoms, 1, &c_emlrtRTEI, true);
  b_emxInit_real_T(&st, &voxels, 1, &c_emlrtRTEI, true);
  b_emxInit_real_T(&st, &fibers, 1, &c_emlrtRTEI, true);
  b_emxInit_real_T(&st, &values, 1, &c_emlrtRTEI, true);
  emxInit_real_T(&st, &D, 2, &c_emlrtRTEI, true);
  b_emxInit_real_T(&st, &w, 1, &c_emlrtRTEI, true);
  b_emxInit_real_T(&st, &Y, 1, &c_emlrtRTEI, true);

  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "atoms", atoms);
  emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "voxels", voxels);
  emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "fibers", fibers);
  emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "values", values);
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "D", D);
  emlrt_marshallIn(&st, emlrtAlias(prhs[5]), "w", w);
  nTheta = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "nTheta");
  nVoxels = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "nVoxels");

  /* Invoke the target function */
  M_times_w_LOOP(&st, atoms, voxels, fibers, values, D, w, nTheta, nVoxels, Y);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(Y);
  Y->canFreeData = false;
  emxFree_real_T(&Y);
  w->canFreeData = false;
  emxFree_real_T(&w);
  D->canFreeData = false;
  emxFree_real_T(&D);
  values->canFreeData = false;
  emxFree_real_T(&values);
  fibers->canFreeData = false;
  emxFree_real_T(&fibers);
  voxels->canFreeData = false;
  emxFree_real_T(&voxels);
  atoms->canFreeData = false;
  emxFree_real_T(&atoms);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_M_times_w_LOOP_api.c) */
