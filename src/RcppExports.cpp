// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rsatsne_cpp
Rcpp::List Rsatsne_cpp(NumericMatrix X1, NumericMatrix X2, NumericMatrix mat12, NumericMatrix mat21, int nk1, int nk2, int no_dims, double perplexity, double theta, bool verbose, int max_iter, bool distance_precomputed, NumericMatrix Y_in, bool init, int stop_lying_iter, int mom_switch_iter, double momentum, double final_momentum, double binding_force, int coupled_period, int uncoupled_period, int mnn_weight, double eta, double exaggeration_factor, unsigned int num_threads, bool return_logs);
RcppExport SEXP _Rsatsne_Rsatsne_cpp(SEXP X1SEXP, SEXP X2SEXP, SEXP mat12SEXP, SEXP mat21SEXP, SEXP nk1SEXP, SEXP nk2SEXP, SEXP no_dimsSEXP, SEXP perplexitySEXP, SEXP thetaSEXP, SEXP verboseSEXP, SEXP max_iterSEXP, SEXP distance_precomputedSEXP, SEXP Y_inSEXP, SEXP initSEXP, SEXP stop_lying_iterSEXP, SEXP mom_switch_iterSEXP, SEXP momentumSEXP, SEXP final_momentumSEXP, SEXP binding_forceSEXP, SEXP coupled_periodSEXP, SEXP uncoupled_periodSEXP, SEXP mnn_weightSEXP, SEXP etaSEXP, SEXP exaggeration_factorSEXP, SEXP num_threadsSEXP, SEXP return_logsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat12(mat12SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat21(mat21SEXP);
    Rcpp::traits::input_parameter< int >::type nk1(nk1SEXP);
    Rcpp::traits::input_parameter< int >::type nk2(nk2SEXP);
    Rcpp::traits::input_parameter< int >::type no_dims(no_dimsSEXP);
    Rcpp::traits::input_parameter< double >::type perplexity(perplexitySEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type distance_precomputed(distance_precomputedSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_in(Y_inSEXP);
    Rcpp::traits::input_parameter< bool >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type stop_lying_iter(stop_lying_iterSEXP);
    Rcpp::traits::input_parameter< int >::type mom_switch_iter(mom_switch_iterSEXP);
    Rcpp::traits::input_parameter< double >::type momentum(momentumSEXP);
    Rcpp::traits::input_parameter< double >::type final_momentum(final_momentumSEXP);
    Rcpp::traits::input_parameter< double >::type binding_force(binding_forceSEXP);
    Rcpp::traits::input_parameter< int >::type coupled_period(coupled_periodSEXP);
    Rcpp::traits::input_parameter< int >::type uncoupled_period(uncoupled_periodSEXP);
    Rcpp::traits::input_parameter< int >::type mnn_weight(mnn_weightSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type exaggeration_factor(exaggeration_factorSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type return_logs(return_logsSEXP);
    rcpp_result_gen = Rcpp::wrap(Rsatsne_cpp(X1, X2, mat12, mat21, nk1, nk2, no_dims, perplexity, theta, verbose, max_iter, distance_precomputed, Y_in, init, stop_lying_iter, mom_switch_iter, momentum, final_momentum, binding_force, coupled_period, uncoupled_period, mnn_weight, eta, exaggeration_factor, num_threads, return_logs));
    return rcpp_result_gen;
END_RCPP
}
// normalize_input_cpp
Rcpp::NumericMatrix normalize_input_cpp(Rcpp::NumericMatrix input);
RcppExport SEXP _Rsatsne_normalize_input_cpp(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(normalize_input_cpp(input));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _Rsatsne_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rsatsne_Rsatsne_cpp", (DL_FUNC) &_Rsatsne_Rsatsne_cpp, 26},
    {"_Rsatsne_normalize_input_cpp", (DL_FUNC) &_Rsatsne_normalize_input_cpp, 1},
    {"_Rsatsne_rcpp_hello", (DL_FUNC) &_Rsatsne_rcpp_hello, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rsatsne(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
