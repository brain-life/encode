/*
 * Mtransp_times_b_initialize.c
 *
 * Code generation for function 'Mtransp_times_b_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Mtransp_times_b.h"
#include "Mtransp_times_b_initialize.h"
#include "Mtransp_times_b_data.h"

/* Function Definitions */
void Mtransp_times_b_initialize(emlrtContext *aContext)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (Mtransp_times_b_initialize.c) */
