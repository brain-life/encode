/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * M_times_w_LOOP_terminate.c
 *
 * Code generation for function 'M_times_w_LOOP_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "M_times_w_LOOP.h"
#include "M_times_w_LOOP_terminate.h"
#include "M_times_w_LOOP_data.h"

/* Function Definitions */
void M_times_w_LOOP_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void M_times_w_LOOP_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (M_times_w_LOOP_terminate.c) */
