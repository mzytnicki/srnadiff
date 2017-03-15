#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP srnadiff_rcpp_buildHmm(SEXP, SEXP, SEXP);
extern SEXP srnadiff_rcpp_viterbi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"srnadiff_rcpp_buildHmm", (DL_FUNC) &srnadiff_rcpp_buildHmm, 3},
    {"srnadiff_rcpp_viterbi",  (DL_FUNC) &srnadiff_rcpp_viterbi,  8},
    {NULL, NULL, 0}
};

void R_init_srnadiff(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
