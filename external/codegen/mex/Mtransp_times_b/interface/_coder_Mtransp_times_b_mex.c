/*
 * _coder_Mtransp_times_b_mex.c
 *
 * Code generation for function 'Mtransp_times_b'
 *
 */

/* Include files */
#include "mex.h"
#include "_coder_Mtransp_times_b_api.h"
#include "Mtransp_times_b_initialize.h"
#include "Mtransp_times_b_terminate.h"

/* Function Declarations */
static void Mtransp_times_b_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "Mtransp_times_b", NULL, false, {2045744189U,2170104910U,2743257031U,4284093946U}, NULL };
void *emlrtRootTLSGlobal = NULL;

/* Function Definitions */
static void Mtransp_times_b_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *outputs[1];
  const mxArray *inputs[7];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  int nInputs = nrhs;
  emlrtStack st = { NULL, NULL, NULL };
  /* Module initialization. */
  Mtransp_times_b_initialize(&emlrtContextGlobal);
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 7) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, mxINT32_CLASS, 7, mxCHAR_CLASS, 15, "Mtransp_times_b");
  } else if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, mxCHAR_CLASS, 15, "Mtransp_times_b");
  }
  /* Temporary copy for mex inputs. */
  for (n = 0; n < nInputs; ++n) {
    inputs[n] = prhs[n];
  }
  /* Call the function. */
  Mtransp_times_b_api(inputs, outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  Mtransp_times_b_terminate();
}

void Mtransp_times_b_atexit_wrapper(void)
{
   Mtransp_times_b_atexit();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(Mtransp_times_b_atexit_wrapper);
  /* Dispatch the entry-point. */
  Mtransp_times_b_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (_coder_Mtransp_times_b_mex.c) */
