// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// kendall_weight_quickR
double kendall_weight_quickR(Rcpp::IntegerVector idxR, Rcpp::IntegerVector zR, int nR, int keyR, int kR, Rcpp::NumericVector uR);
RcppExport SEXP _kernrank_kendall_weight_quickR(SEXP idxRSEXP, SEXP zRSEXP, SEXP nRSEXP, SEXP keyRSEXP, SEXP kRSEXP, SEXP uRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idxR(idxRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type zR(zRSEXP);
    Rcpp::traits::input_parameter< int >::type nR(nRSEXP);
    Rcpp::traits::input_parameter< int >::type keyR(keyRSEXP);
    Rcpp::traits::input_parameter< int >::type kR(kRSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type uR(uRSEXP);
    rcpp_result_gen = Rcpp::wrap(kendall_weight_quickR(idxR, zR, nR, keyR, kR, uR));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kernrank_kendall_weight_quickR", (DL_FUNC) &_kernrank_kendall_weight_quickR, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_kernrank(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
