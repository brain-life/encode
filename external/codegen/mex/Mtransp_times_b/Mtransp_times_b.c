/*
 * Mtransp_times_b.c
 *
 * Code generation for function 'Mtransp_times_b'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Mtransp_times_b.h"
#include "Mtransp_times_b_emxutil.h"
#include "Mtransp_times_b_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 23, "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m"
};

static emlrtRSInfo b_emlrtRSI = { 61, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

static emlrtRSInfo c_emlrtRSI = { 21, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

static emlrtRSInfo d_emlrtRSI = { 30, "eml_xdotu",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"
};

static emlrtMCInfo emlrtMCI = { 99, 13, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

static emlrtMCInfo b_emlrtMCI = { 98, 23, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

static emlrtMCInfo c_emlrtMCI = { 104, 13, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

static emlrtMCInfo d_emlrtMCI = { 103, 23, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

static emlrtRTEInfo emlrtRTEI = { 19, 18, "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m"
};

static emlrtDCInfo emlrtDCI = { 23, 40, "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  1 };

static emlrtBCInfo emlrtBCI = { -1, -1, 23, 40, "D", "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  0 };

static emlrtDCInfo b_emlrtDCI = { 23, 56, "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  1 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 23, 56, "Y", "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  0 };

static emlrtDCInfo c_emlrtDCI = { 20, 11, "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  1 };

static emlrtDCInfo d_emlrtDCI = { 20, 11, "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  4 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 23, 40, "atoms", "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  0 };

static emlrtBCInfo d_emlrtBCI = { -1, -1, 23, 56, "voxels", "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  0 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 23, 5, "w", "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  0 };

static emlrtDCInfo e_emlrtDCI = { 23, 5, "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  1 };

static emlrtBCInfo f_emlrtBCI = { -1, -1, 23, 7, "fibers", "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 23, 20, "w", "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  0 };

static emlrtDCInfo f_emlrtDCI = { 23, 20, "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  1 };

static emlrtBCInfo h_emlrtBCI = { -1, -1, 23, 22, "fibers", "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  0 };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 23, 67, "values", "Mtransp_times_b",
  "/Users/ccaiafa/Dropbox/DTI/life_BD/Beta_Version_2_0/life_BD-master/external/Mtransp_times_b.m",
  0 };

static emlrtRSInfo g_emlrtRSI = { 98, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

static emlrtRSInfo h_emlrtRSI = { 103, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

static emlrtRSInfo i_emlrtRSI = { 99, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

static emlrtRSInfo j_emlrtRSI = { 104, "eml_mtimes_helper",
  "/Applications/MATLAB_R2014b.app/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"
};

/* Function Declarations */
static void error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location);
static const mxArray *message(const emlrtStack *sp, const mxArray *b,
  emlrtMCInfo *location);

/* Function Definitions */
static void error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "error", true, location);
}

static const mxArray *message(const emlrtStack *sp, const mxArray *b,
  emlrtMCInfo *location)
{
  const mxArray *pArray;
  const mxArray *m2;
  pArray = b;
  return emlrtCallMATLABR2012b(sp, 1, &m2, 1, &pArray, "message", true, location);
}

void Mtransp_times_b(const emlrtStack *sp, const emxArray_real_T *atoms, const
                     emxArray_real_T *voxels, const emxArray_real_T *fibers,
                     const emxArray_real_T *values, const emxArray_real_T *D,
                     const emxArray_real_T *Y, real_T nFibers, emxArray_real_T
                     *w)
{
  int32_T i0;
  real_T d0;
  int32_T loop_ub;
  int32_T k;
  emxArray_real_T *a;
  emxArray_real_T *b;
  int32_T i1;
  int32_T i;
  boolean_T guard2 = false;
  const mxArray *y;
  static const int32_T iv0[2] = { 1, 21 };

  const mxArray *m0;
  char_T cv0[21];
  static const char_T cv1[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  const mxArray *b_y;
  static const int32_T iv1[2] = { 1, 45 };

  char_T cv2[45];
  static const char_T cv3[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D',
    'y', 'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p',
    'a', 'n', 's', 'i', 'o', 'n' };

  boolean_T guard1 = false;
  real_T c_y;
  ptrdiff_t n_t;
  ptrdiff_t incx_t;
  ptrdiff_t incy_t;
  double * xix0_t;
  double * yiy0_t;
  int32_T i2;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &b_st;
  d_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /*  function [ w ] = Mtransp_times_b( M, b ) */
  /*  [nFibers] = size(M.Phi,1); */
  /*  [nTheta]  = size(M.DictSig,1); */
  /*  [nAtoms] = size(M.Phi,2); */
  /*  [Nvoxels] = size(M.Phi,3); */
  /*   */
  /*  w = M.DictSig'*reshape(b,[nTheta,Nvoxels]); % This is still a little bit  */
  /*  % memory expensive because w is a (Nd x Nvoxels) matrix (Nd = # atoms) */
  /*  % See comments below for a very memory efficient implementation but very */
  /*  % slow without using a compiled version (MEX). */
  /*  M.Phi = reshape(M.Phi,[nFibers,nAtoms*Nvoxels]); */
  /*  w = double(ttv(M.Phi,w(:),2)); */
  /*  end */
  /* % */
  /*  The following is a memory efficient version that should be compiled in */
  /*  order to provide fast results */
  i0 = w->size[0];
  d0 = emlrtNonNegativeCheckFastR2012b(nFibers, &d_emlrtDCI, sp);
  w->size[0] = (int32_T)emlrtIntegerCheckFastR2012b(d0, &c_emlrtDCI, sp);
  emxEnsureCapacity(sp, (emxArray__common *)w, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  d0 = emlrtNonNegativeCheckFastR2012b(nFibers, &d_emlrtDCI, sp);
  loop_ub = (int32_T)emlrtIntegerCheckFastR2012b(d0, &c_emlrtDCI, sp);
  for (i0 = 0; i0 < loop_ub; i0++) {
    w->data[i0] = 0.0;
  }

  k = 1;
  emxInit_real_T(sp, &a, 2, &emlrtRTEI, true);
  b_emxInit_real_T(sp, &b, 1, &emlrtRTEI, true);
  while (k - 1 <= values->size[0] - 1) {
    /* k */
    st.site = &emlrtRSI;
    loop_ub = D->size[0];
    i0 = atoms->size[0];
    d0 = atoms->data[emlrtDynamicBoundsCheckFastR2012b(k, 1, i0, &c_emlrtBCI,
      &st) - 1];
    i0 = D->size[1];
    i1 = (int32_T)emlrtIntegerCheckFastR2012b(d0, &emlrtDCI, &st);
    i = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &emlrtBCI, &st);
    i0 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)a, i0, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i0 = 0; i0 < loop_ub; i0++) {
      a->data[a->size[0] * i0] = D->data[i0 + D->size[0] * (i - 1)];
    }

    loop_ub = Y->size[0];
    i0 = voxels->size[0];
    d0 = voxels->data[emlrtDynamicBoundsCheckFastR2012b(k, 1, i0, &d_emlrtBCI,
      &st) - 1];
    i0 = Y->size[1];
    i1 = (int32_T)emlrtIntegerCheckFastR2012b(d0, &b_emlrtDCI, &st);
    i = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &b_emlrtBCI, &st);
    i0 = b->size[0];
    b->size[0] = loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)b, i0, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i0 = 0; i0 < loop_ub; i0++) {
      b->data[i0] = Y->data[i0 + Y->size[0] * (i - 1)];
    }

    b_st.site = &c_emlrtRSI;
    i0 = Y->size[0];
    if (!(a->size[1] == i0)) {
      guard2 = false;
      if (a->size[1] == 1) {
        guard2 = true;
      } else {
        i0 = Y->size[0];
        if (i0 == 1) {
          guard2 = true;
        } else {
          y = NULL;
          m0 = emlrtCreateCharArray(2, iv0);
          for (i = 0; i < 21; i++) {
            cv0[i] = cv1[i];
          }

          emlrtInitCharArrayR2013a(&b_st, 21, m0, cv0);
          emlrtAssign(&y, m0);
          c_st.site = &h_emlrtRSI;
          d_st.site = &j_emlrtRSI;
          error(&c_st, message(&d_st, y, &c_emlrtMCI), &d_emlrtMCI);
        }
      }

      if (guard2) {
        b_y = NULL;
        m0 = emlrtCreateCharArray(2, iv1);
        for (i = 0; i < 45; i++) {
          cv2[i] = cv3[i];
        }

        emlrtInitCharArrayR2013a(&b_st, 45, m0, cv2);
        emlrtAssign(&b_y, m0);
        c_st.site = &g_emlrtRSI;
        d_st.site = &i_emlrtRSI;
        error(&c_st, message(&d_st, b_y, &emlrtMCI), &b_emlrtMCI);
      }
    }

    guard1 = false;
    if (a->size[1] == 1) {
      guard1 = true;
    } else {
      i0 = Y->size[0];
      if (i0 == 1) {
        guard1 = true;
      } else {
        b_st.site = &b_emlrtRSI;
        c_st.site = &d_emlrtRSI;
        if (a->size[1] < 1) {
          c_y = 0.0;
        } else {
          n_t = (ptrdiff_t)(a->size[1]);
          incx_t = (ptrdiff_t)(1);
          incy_t = (ptrdiff_t)(1);
          xix0_t = (double *)(&a->data[0]);
          yiy0_t = (double *)(&b->data[0]);
          c_y = ddot(&n_t, xix0_t, &incx_t, yiy0_t, &incy_t);
        }
      }
    }

    if (guard1) {
      c_y = 0.0;
      for (i0 = 0; i0 < a->size[1]; i0++) {
        c_y += a->data[a->size[0] * i0] * b->data[i0];
      }
    }

    i0 = w->size[0];
    i1 = fibers->size[0];
    d0 = fibers->data[emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &f_emlrtBCI,
      sp) - 1];
    i1 = (int32_T)emlrtIntegerCheckFastR2012b(d0, &e_emlrtDCI, sp);
    i = w->size[0];
    loop_ub = fibers->size[0];
    d0 = fibers->data[emlrtDynamicBoundsCheckFastR2012b(k, 1, loop_ub,
      &h_emlrtBCI, sp) - 1];
    loop_ub = (int32_T)emlrtIntegerCheckFastR2012b(d0, &f_emlrtDCI, sp);
    i2 = values->size[0];
    w->data[emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &e_emlrtBCI, sp) - 1] =
      w->data[emlrtDynamicBoundsCheckFastR2012b(loop_ub, 1, i, &g_emlrtBCI, sp)
      - 1] + c_y * values->data[emlrtDynamicBoundsCheckFastR2012b(k, 1, i2,
      &i_emlrtBCI, sp) - 1];
    k++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, sp);
  }

  emxFree_real_T(&b);
  emxFree_real_T(&a);

  /*   */
  /*  for f=1:Nf */
  /*      for l=1:Nelem(f) */
  /*          w(f) = w(f) + (D(:,atoms(f,l)))'*Y(:,voxels(f,l))*S0(voxels(f,l)); */
  /*      end */
  /*  end */
  /*  [nFibers] = size(M.Phi,1); */
  /*  [nTheta]  = size(M.DictSig,1); */
  /*  [Nvoxels] = size(M.Phi,3); */
  /*   */
  /*  b = reshape(b,[nTheta,Nvoxels]); % matrix reshape */
  /*   */
  /*  w = zeros(nFibers,1); % output */
  /*   */
  /*  for f=1:nFibers */
  /*      f */
  /*      A = sptenmat(M.Phi(f,:,:),1); % ->sptenmat */
  /*      aux = sum(b(:,A.subs(:,2)).*M.DictSig(:,A.subs(:,1)),1); */
  /*      w(f) = sum(aux'.*A.vals(:)); */
  /*  end */
  /*   */
  /*   */
  /*  end */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (Mtransp_times_b.c) */
