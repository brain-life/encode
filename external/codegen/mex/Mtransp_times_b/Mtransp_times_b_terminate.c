/*
 * Mtransp_times_b_terminate.c
 *
 * Code generation for function 'Mtransp_times_b_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Mtransp_times_b.h"
#include "Mtransp_times_b_terminate.h"

/* Function Definitions */
void Mtransp_times_b_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void Mtransp_times_b_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (Mtransp_times_b_terminate.c) */
