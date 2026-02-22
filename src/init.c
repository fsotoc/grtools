#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP grtools_matrixloglikC(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"grtools_matrixloglikC", (DL_FUNC) &grtools_matrixloglikC, 5},
  {NULL, NULL, 0}
};

void R_init_grtools(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
