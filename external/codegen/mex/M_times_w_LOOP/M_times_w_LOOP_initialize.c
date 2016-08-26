/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * M_times_w_LOOP_initialize.c
 *
 * Code generation for function 'M_times_w_LOOP_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "M_times_w_LOOP.h"
#include "M_times_w_LOOP_initialize.h"
#include "M_times_w_LOOP_data.h"

/* Function Definitions */
void M_times_w_LOOP_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (M_times_w_LOOP_initialize.c) */
