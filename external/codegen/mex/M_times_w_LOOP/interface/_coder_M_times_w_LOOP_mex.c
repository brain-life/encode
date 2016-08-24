/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_M_times_w_LOOP_mex.c
 *
 * Code generation for function '_coder_M_times_w_LOOP_mex'
 *
 */

/* Include files */
#include "M_times_w_LOOP.h"
#include "_coder_M_times_w_LOOP_mex.h"
#include "M_times_w_LOOP_terminate.h"
#include "_coder_M_times_w_LOOP_api.h"
#include "M_times_w_LOOP_initialize.h"
#include "M_times_w_LOOP_data.h"

/* Function Declarations */
static void M_times_w_LOOP_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T
  nrhs, const mxArray *prhs[8]);

/* Function Definitions */
static void M_times_w_LOOP_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T
  nrhs, const mxArray *prhs[8])
{
  int32_T n;
  const mxArray *inputs[8];
  const mxArray *outputs[1];
  int32_T b_nlhs;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 8) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 8, 4,
                        14, "M_times_w_LOOP");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 14,
                        "M_times_w_LOOP");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  /* Call the function. */
  M_times_w_LOOP_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  M_times_w_LOOP_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(M_times_w_LOOP_atexit);

  /* Module initialization. */
  M_times_w_LOOP_initialize();

  /* Dispatch the entry-point. */
  M_times_w_LOOP_mexFunction(nlhs, plhs, nrhs, prhs);
}

/* End of code generation (_coder_M_times_w_LOOP_mex.c) */
