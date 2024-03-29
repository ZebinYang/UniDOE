// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CritEval
double CritEval(NumericMatrix x0, int nlevel, StringVector crit);
RcppExport SEXP _UniDOE_CritEval(SEXP x0SEXP, SEXP nlevelSEXP, SEXP critSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type nlevel(nlevelSEXP);
    Rcpp::traits::input_parameter< StringVector >::type crit(critSEXP);
    rcpp_result_gen = Rcpp::wrap(CritEval(x0, nlevel, crit));
    return rcpp_result_gen;
END_RCPP
}
// SATA_UD
List SATA_UD(int nsamp, int nv, int nlevel, StringVector init_method, NumericMatrix initX, StringVector crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed);
RcppExport SEXP _UniDOE_SATA_UD(SEXP nsampSEXP, SEXP nvSEXP, SEXP nlevelSEXP, SEXP init_methodSEXP, SEXP initXSEXP, SEXP critSEXP, SEXP maxiterSEXP, SEXP hits_ratioSEXP, SEXP levelpermtSEXP, SEXP rand_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< int >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< int >::type nlevel(nlevelSEXP);
    Rcpp::traits::input_parameter< StringVector >::type init_method(init_methodSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initX(initXSEXP);
    Rcpp::traits::input_parameter< StringVector >::type crit(critSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type hits_ratio(hits_ratioSEXP);
    Rcpp::traits::input_parameter< bool >::type levelpermt(levelpermtSEXP);
    Rcpp::traits::input_parameter< int >::type rand_seed(rand_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(SATA_UD(nsamp, nv, nlevel, init_method, initX, crit, maxiter, hits_ratio, levelpermt, rand_seed));
    return rcpp_result_gen;
END_RCPP
}
// SATA_AUD
List SATA_AUD(NumericMatrix xp, int nnew, int nv, int nlevel, StringVector init_method, NumericMatrix initX, StringVector crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed);
RcppExport SEXP _UniDOE_SATA_AUD(SEXP xpSEXP, SEXP nnewSEXP, SEXP nvSEXP, SEXP nlevelSEXP, SEXP init_methodSEXP, SEXP initXSEXP, SEXP critSEXP, SEXP maxiterSEXP, SEXP hits_ratioSEXP, SEXP levelpermtSEXP, SEXP rand_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< int >::type nnew(nnewSEXP);
    Rcpp::traits::input_parameter< int >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< int >::type nlevel(nlevelSEXP);
    Rcpp::traits::input_parameter< StringVector >::type init_method(init_methodSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initX(initXSEXP);
    Rcpp::traits::input_parameter< StringVector >::type crit(critSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type hits_ratio(hits_ratioSEXP);
    Rcpp::traits::input_parameter< bool >::type levelpermt(levelpermtSEXP);
    Rcpp::traits::input_parameter< int >::type rand_seed(rand_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(SATA_AUD(xp, nnew, nv, nlevel, init_method, initX, crit, maxiter, hits_ratio, levelpermt, rand_seed));
    return rcpp_result_gen;
END_RCPP
}
// SATA_AUD_COL
List SATA_AUD_COL(NumericMatrix xp, int nvnew, int nlevel, StringVector init_method, NumericMatrix initX, StringVector crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed);
RcppExport SEXP _UniDOE_SATA_AUD_COL(SEXP xpSEXP, SEXP nvnewSEXP, SEXP nlevelSEXP, SEXP init_methodSEXP, SEXP initXSEXP, SEXP critSEXP, SEXP maxiterSEXP, SEXP hits_ratioSEXP, SEXP levelpermtSEXP, SEXP rand_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< int >::type nvnew(nvnewSEXP);
    Rcpp::traits::input_parameter< int >::type nlevel(nlevelSEXP);
    Rcpp::traits::input_parameter< StringVector >::type init_method(init_methodSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initX(initXSEXP);
    Rcpp::traits::input_parameter< StringVector >::type crit(critSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type hits_ratio(hits_ratioSEXP);
    Rcpp::traits::input_parameter< bool >::type levelpermt(levelpermtSEXP);
    Rcpp::traits::input_parameter< int >::type rand_seed(rand_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(SATA_AUD_COL(xp, nvnew, nlevel, init_method, initX, crit, maxiter, hits_ratio, levelpermt, rand_seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_UniDOE_CritEval", (DL_FUNC) &_UniDOE_CritEval, 3},
    {"_UniDOE_SATA_UD", (DL_FUNC) &_UniDOE_SATA_UD, 10},
    {"_UniDOE_SATA_AUD", (DL_FUNC) &_UniDOE_SATA_AUD, 11},
    {"_UniDOE_SATA_AUD_COL", (DL_FUNC) &_UniDOE_SATA_AUD_COL, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_UniDOE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
