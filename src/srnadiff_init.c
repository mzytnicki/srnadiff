#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _srnadiff_rcpp_buildHmm(SEXP, SEXP);
extern SEXP _srnadiff_rcpp_viterbi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _srnadiff_rcpp_normalization(SEXP, SEXP);
extern SEXP _srnadiff_rcpp_slice(SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP _srnadiff_rcpp_naive(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_srnadiff_rcpp_buildHmm",
        (DL_FUNC) &_srnadiff_rcpp_buildHmm,      2},
    {"_srnadiff_rcpp_viterbi",
        (DL_FUNC) &_srnadiff_rcpp_viterbi,       10},
    {"_srnadiff_rcpp_normalization",
        (DL_FUNC) &_srnadiff_rcpp_normalization, 2},
    {"_srnadiff_rcpp_slice",
        (DL_FUNC) &_srnadiff_rcpp_slice,         5},
/*
    {"_srnadiff_rcpp_naive",
        (DL_FUNC) &_srnadiff_rcpp_naive,         7},
*/
    {NULL, NULL, 0}
};

void R_init_srnadiff(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
