/*
 * _coder_Mtransp_times_b_api.c
 *
 * Code generation for function '_coder_Mtransp_times_b_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Mtransp_times_b.h"
#include "_coder_Mtransp_times_b_api.h"
#include "Mtransp_times_b_emxutil.h"

/* Variable Definitions */
static emlrtRTEInfo b_emlrtRTEI = { 1, 1, "_coder_Mtransp_times_b_api", "" };

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *D, const
  char_T *identifier, emxArray_real_T *y);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nFibers,
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

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nFibers,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nFibers), &thisId);
  emlrtDestroyArray(&nFibers);
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
  static const int32_T iv2[1] = { 0 };

  const mxArray *m1;
  y = NULL;
  m1 = emlrtCreateNumericArray(1, iv2, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m1, (void *)u->data);
  emlrtSetDimensions((mxArray *)m1, u->size, 1);
  emlrtAssign(&y, m1);
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
  int32_T iv3[1];
  boolean_T bv0[1];
  int32_T iv4[1];
  iv3[0] = -1;
  bv0[0] = true;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 1U, iv3, bv0, iv4);
  ret->size[0] = iv4[0];
  ret->allocatedSize = ret->size[0];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv5[2];
  boolean_T bv1[2];
  int32_T i;
  int32_T iv6[2];
  for (i = 0; i < 2; i++) {
    iv5[i] = -1;
    bv1[i] = true;
  }

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv5, bv1, iv6);
  ret->size[0] = iv6[0];
  ret->size[1] = iv6[1];
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

void Mtransp_times_b_api(const mxArray * const prhs[7], const mxArray *plhs[1])
{
  emxArray_real_T *atoms;
  emxArray_real_T *voxels;
  emxArray_real_T *fibers;
  emxArray_real_T *values;
  emxArray_real_T *D;
  emxArray_real_T *Y;
  emxArray_real_T *w;
  real_T nFibers;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  b_emxInit_real_T(&st, &atoms, 1, &b_emlrtRTEI, true);
  b_emxInit_real_T(&st, &voxels, 1, &b_emlrtRTEI, true);
  b_emxInit_real_T(&st, &fibers, 1, &b_emlrtRTEI, true);
  b_emxInit_real_T(&st, &values, 1, &b_emlrtRTEI, true);
  emxInit_real_T(&st, &D, 2, &b_emlrtRTEI, true);
  emxInit_real_T(&st, &Y, 2, &b_emlrtRTEI, true);
  b_emxInit_real_T(&st, &w, 1, &b_emlrtRTEI, true);

  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "atoms", atoms);
  emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "voxels", voxels);
  emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "fibers", fibers);
  emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "values", values);
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "D", D);
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[5]), "Y", Y);
  nFibers = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "nFibers");

  /* Invoke the target function */
  Mtransp_times_b(&st, atoms, voxels, fibers, values, D, Y, nFibers, w);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(w);
  w->canFreeData = false;
  emxFree_real_T(&w);
  Y->canFreeData = false;
  emxFree_real_T(&Y);
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

/* End of code generation (_coder_Mtransp_times_b_api.c) */
