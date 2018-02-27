// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// safeSoftMax
arma::mat safeSoftMax(arma::mat logs, double expTrc);
RcppExport SEXP _sdetorus_safeSoftMax(SEXP logsSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type logs(logsSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(safeSoftMax(logs, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// solveTridiag
arma::vec solveTridiag(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int LU);
RcppExport SEXP _sdetorus_solveTridiag(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP LUSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type LU(LUSEXP);
    rcpp_result_gen = Rcpp::wrap(solveTridiag(a, b, c, d, LU));
    return rcpp_result_gen;
END_RCPP
}
// solveTridiagMatConsts
arma::mat solveTridiagMatConsts(arma::vec a, arma::vec b, arma::vec c, arma::mat d, int LU);
RcppExport SEXP _sdetorus_solveTridiagMatConsts(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP LUSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type LU(LUSEXP);
    rcpp_result_gen = Rcpp::wrap(solveTridiagMatConsts(a, b, c, d, LU));
    return rcpp_result_gen;
END_RCPP
}
// solvePeriodicTridiag
arma::vec solvePeriodicTridiag(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int LU);
RcppExport SEXP _sdetorus_solvePeriodicTridiag(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP LUSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type LU(LUSEXP);
    rcpp_result_gen = Rcpp::wrap(solvePeriodicTridiag(a, b, c, d, LU));
    return rcpp_result_gen;
END_RCPP
}
// forwardSweepTridiag
arma::mat forwardSweepTridiag(arma::vec a, arma::vec b, arma::vec c);
RcppExport SEXP _sdetorus_forwardSweepTridiag(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(forwardSweepTridiag(a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// forwardSweepPeriodicTridiag
arma::mat forwardSweepPeriodicTridiag(arma::vec a, arma::vec b, arma::vec c);
RcppExport SEXP _sdetorus_forwardSweepPeriodicTridiag(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(forwardSweepPeriodicTridiag(a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// crankNicolson1D
arma::mat crankNicolson1D(arma::mat u0, arma::vec b, arma::vec sigma2, arma::uvec N, double deltat, arma::uword Mx, double deltax, int imposePositive);
RcppExport SEXP _sdetorus_crankNicolson1D(SEXP u0SEXP, SEXP bSEXP, SEXP sigma2SEXP, SEXP NSEXP, SEXP deltatSEXP, SEXP MxSEXP, SEXP deltaxSEXP, SEXP imposePositiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type u0(u0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type deltat(deltatSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type Mx(MxSEXP);
    Rcpp::traits::input_parameter< double >::type deltax(deltaxSEXP);
    Rcpp::traits::input_parameter< int >::type imposePositive(imposePositiveSEXP);
    rcpp_result_gen = Rcpp::wrap(crankNicolson1D(u0, b, sigma2, N, deltat, Mx, deltax, imposePositive));
    return rcpp_result_gen;
END_RCPP
}
// crankNicolson2D
arma::mat crankNicolson2D(arma::mat u0, arma::mat bx, arma::mat by, arma::mat sigma2x, arma::mat sigma2y, arma::mat sigmaxy, arma::uvec N, double deltat, arma::uword Mx, double deltax, arma::uword My, double deltay, int imposePositive);
RcppExport SEXP _sdetorus_crankNicolson2D(SEXP u0SEXP, SEXP bxSEXP, SEXP bySEXP, SEXP sigma2xSEXP, SEXP sigma2ySEXP, SEXP sigmaxySEXP, SEXP NSEXP, SEXP deltatSEXP, SEXP MxSEXP, SEXP deltaxSEXP, SEXP MySEXP, SEXP deltaySEXP, SEXP imposePositiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type u0(u0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type bx(bxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type by(bySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma2x(sigma2xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma2y(sigma2ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigmaxy(sigmaxySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type deltat(deltatSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type Mx(MxSEXP);
    Rcpp::traits::input_parameter< double >::type deltax(deltaxSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type My(MySEXP);
    Rcpp::traits::input_parameter< double >::type deltay(deltaySEXP);
    Rcpp::traits::input_parameter< int >::type imposePositive(imposePositiveSEXP);
    rcpp_result_gen = Rcpp::wrap(crankNicolson2D(u0, bx, by, sigma2x, sigma2y, sigmaxy, N, deltat, Mx, deltax, My, deltay, imposePositive));
    return rcpp_result_gen;
END_RCPP
}
// driftWn1D
arma::vec driftWn1D(arma::vec x, double alpha, double mu, double sigma, int maxK, double expTrc);
RcppExport SEXP _sdetorus_driftWn1D(SEXP xSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(driftWn1D(x, alpha, mu, sigma, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// driftWn2D
arma::mat driftWn2D(arma::mat x, arma::mat A, arma::vec mu, arma::vec sigma, double rho, int maxK, double expTrc);
RcppExport SEXP _sdetorus_driftWn2D(SEXP xSEXP, SEXP ASEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(driftWn2D(x, A, mu, sigma, rho, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// euler1D
arma::mat euler1D(arma::vec x0, double alpha, double mu, double sigma, arma::uword N, double delta, int type, int maxK, double expTrc);
RcppExport SEXP _sdetorus_euler1D(SEXP x0SEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP NSEXP, SEXP deltaSEXP, SEXP typeSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(euler1D(x0, alpha, mu, sigma, N, delta, type, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// euler2D
arma::cube euler2D(arma::mat x0, arma::mat A, arma::vec mu, arma::vec sigma, double rho, arma::uword N, double delta, int type, int maxK, double expTrc);
RcppExport SEXP _sdetorus_euler2D(SEXP x0SEXP, SEXP ASEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP NSEXP, SEXP deltaSEXP, SEXP typeSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(euler2D(x0, A, mu, sigma, rho, N, delta, type, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// stepAheadWn1D
arma::mat stepAheadWn1D(arma::vec x0, double alpha, double mu, double sigma, arma::uword M, arma::uword N, double delta, int type, int maxK, double expTrc);
RcppExport SEXP _sdetorus_stepAheadWn1D(SEXP x0SEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP MSEXP, SEXP NSEXP, SEXP deltaSEXP, SEXP typeSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(stepAheadWn1D(x0, alpha, mu, sigma, M, N, delta, type, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// stepAheadWn2D
arma::cube stepAheadWn2D(arma::mat x0, arma::vec mu, arma::mat A, arma::vec sigma, double rho, arma::uword M, arma::uword N, double delta, int type, int maxK, double expTrc);
RcppExport SEXP _sdetorus_stepAheadWn2D(SEXP x0SEXP, SEXP muSEXP, SEXP ASEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP MSEXP, SEXP NSEXP, SEXP deltaSEXP, SEXP typeSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(stepAheadWn2D(x0, mu, A, sigma, rho, M, N, delta, type, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// linInterp
arma::vec linInterp(arma::vec x, arma::vec xGrid, arma::vec yGrid, bool equalSpaces);
RcppExport SEXP _sdetorus_linInterp(SEXP xSEXP, SEXP xGridSEXP, SEXP yGridSEXP, SEXP equalSpacesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yGrid(yGridSEXP);
    Rcpp::traits::input_parameter< bool >::type equalSpaces(equalSpacesSEXP);
    rcpp_result_gen = Rcpp::wrap(linInterp(x, xGrid, yGrid, equalSpaces));
    return rcpp_result_gen;
END_RCPP
}
// besselIExponScaled
arma::vec besselIExponScaled(arma::vec x, int nu, int maxK, bool equalSpaces);
RcppExport SEXP _sdetorus_besselIExponScaled(SEXP xSEXP, SEXP nuSEXP, SEXP maxKSEXP, SEXP equalSpacesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< bool >::type equalSpaces(equalSpacesSEXP);
    rcpp_result_gen = Rcpp::wrap(besselIExponScaled(x, nu, maxK, equalSpaces));
    return rcpp_result_gen;
END_RCPP
}
// dVmfCpp
arma::vec dVmfCpp(arma::mat x, arma::mat K, arma::mat M, arma::vec alpha, bool besselInterp, double l2pi);
RcppExport SEXP _sdetorus_dVmfCpp(SEXP xSEXP, SEXP KSEXP, SEXP MSEXP, SEXP alphaSEXP, SEXP besselInterpSEXP, SEXP l2piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type besselInterp(besselInterpSEXP);
    Rcpp::traits::input_parameter< double >::type l2pi(l2piSEXP);
    rcpp_result_gen = Rcpp::wrap(dVmfCpp(x, K, M, alpha, besselInterp, l2pi));
    return rcpp_result_gen;
END_RCPP
}
// clusterProbsVmf
arma::mat clusterProbsVmf(arma::mat cosData, arma::mat sinData, const arma::mat M, const arma::mat K, arma::rowvec alpha, double l2pi, bool besselInterp);
RcppExport SEXP _sdetorus_clusterProbsVmf(SEXP cosDataSEXP, SEXP sinDataSEXP, SEXP MSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP l2piSEXP, SEXP besselInterpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type cosData(cosDataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sinData(sinDataSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type l2pi(l2piSEXP);
    Rcpp::traits::input_parameter< bool >::type besselInterp(besselInterpSEXP);
    rcpp_result_gen = Rcpp::wrap(clusterProbsVmf(cosData, sinData, M, K, alpha, l2pi, besselInterp));
    return rcpp_result_gen;
END_RCPP
}
// weightedMuKappa
Rcpp::List weightedMuKappa(arma::mat cosData, arma::mat sinData, arma::mat weights, double kappaMax, bool isotropic);
RcppExport SEXP _sdetorus_weightedMuKappa(SEXP cosDataSEXP, SEXP sinDataSEXP, SEXP weightsSEXP, SEXP kappaMaxSEXP, SEXP isotropicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type cosData(cosDataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sinData(sinDataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type kappaMax(kappaMaxSEXP);
    Rcpp::traits::input_parameter< bool >::type isotropic(isotropicSEXP);
    rcpp_result_gen = Rcpp::wrap(weightedMuKappa(cosData, sinData, weights, kappaMax, isotropic));
    return rcpp_result_gen;
END_RCPP
}
// logLikWouPairs
double logLikWouPairs(arma::mat x, arma::vec t, arma::vec alpha, arma::vec mu, arma::vec sigma, double rho, int maxK, double expTrc);
RcppExport SEXP _sdetorus_logLikWouPairs(SEXP xSEXP, SEXP tSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikWouPairs(x, t, alpha, mu, sigma, rho, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// dWn1D
arma::vec dWn1D(arma::vec x, double mu, double sigma, int maxK, double expTrc, int vmApprox, double kt, double logConstKt);
RcppExport SEXP _sdetorus_dWn1D(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP maxKSEXP, SEXP expTrcSEXP, SEXP vmApproxSEXP, SEXP ktSEXP, SEXP logConstKtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    Rcpp::traits::input_parameter< int >::type vmApprox(vmApproxSEXP);
    Rcpp::traits::input_parameter< double >::type kt(ktSEXP);
    Rcpp::traits::input_parameter< double >::type logConstKt(logConstKtSEXP);
    rcpp_result_gen = Rcpp::wrap(dWn1D(x, mu, sigma, maxK, expTrc, vmApprox, kt, logConstKt));
    return rcpp_result_gen;
END_RCPP
}
// dTpdWou1D
arma::vec dTpdWou1D(arma::vec x, arma::vec x0, double t, double alpha, double mu, double sigma, int maxK, double expTrc, int vmApprox, double kt, double logConstKt);
RcppExport SEXP _sdetorus_dTpdWou1D(SEXP xSEXP, SEXP x0SEXP, SEXP tSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP maxKSEXP, SEXP expTrcSEXP, SEXP vmApproxSEXP, SEXP ktSEXP, SEXP logConstKtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    Rcpp::traits::input_parameter< int >::type vmApprox(vmApproxSEXP);
    Rcpp::traits::input_parameter< double >::type kt(ktSEXP);
    Rcpp::traits::input_parameter< double >::type logConstKt(logConstKtSEXP);
    rcpp_result_gen = Rcpp::wrap(dTpdWou1D(x, x0, t, alpha, mu, sigma, maxK, expTrc, vmApprox, kt, logConstKt));
    return rcpp_result_gen;
END_RCPP
}
// dTpdWou2D
arma::vec dTpdWou2D(arma::mat x, arma::mat x0, arma::vec t, arma::vec alpha, arma::vec mu, arma::vec sigma, double rho, int maxK, double expTrc);
RcppExport SEXP _sdetorus_dTpdWou2D(SEXP xSEXP, SEXP x0SEXP, SEXP tSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(dTpdWou2D(x, x0, t, alpha, mu, sigma, rho, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// rTpdWn2D
arma::cube rTpdWn2D(arma::uword n, arma::mat x0, arma::vec t, arma::vec mu, arma::vec alpha, arma::vec sigma, double rho, int maxK, double expTrc);
RcppExport SEXP _sdetorus_rTpdWn2D(SEXP nSEXP, SEXP x0SEXP, SEXP tSEXP, SEXP muSEXP, SEXP alphaSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(rTpdWn2D(n, x0, t, mu, alpha, sigma, rho, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// dStatWn2D
arma::vec dStatWn2D(arma::mat x, arma::vec alpha, arma::vec mu, arma::vec sigma, double rho, int maxK, double expTrc);
RcppExport SEXP _sdetorus_dStatWn2D(SEXP xSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP maxKSEXP, SEXP expTrcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    Rcpp::traits::input_parameter< double >::type expTrc(expTrcSEXP);
    rcpp_result_gen = Rcpp::wrap(dStatWn2D(x, alpha, mu, sigma, rho, maxK, expTrc));
    return rcpp_result_gen;
END_RCPP
}
// rStatWn2D
arma::mat rStatWn2D(arma::uword n, arma::vec mu, arma::vec alpha, arma::vec sigma, double rho);
RcppExport SEXP _sdetorus_rStatWn2D(SEXP nSEXP, SEXP muSEXP, SEXP alphaSEXP, SEXP sigmaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(rStatWn2D(n, mu, alpha, sigma, rho));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sdetorus_safeSoftMax", (DL_FUNC) &_sdetorus_safeSoftMax, 2},
    {"_sdetorus_solveTridiag", (DL_FUNC) &_sdetorus_solveTridiag, 5},
    {"_sdetorus_solveTridiagMatConsts", (DL_FUNC) &_sdetorus_solveTridiagMatConsts, 5},
    {"_sdetorus_solvePeriodicTridiag", (DL_FUNC) &_sdetorus_solvePeriodicTridiag, 5},
    {"_sdetorus_forwardSweepTridiag", (DL_FUNC) &_sdetorus_forwardSweepTridiag, 3},
    {"_sdetorus_forwardSweepPeriodicTridiag", (DL_FUNC) &_sdetorus_forwardSweepPeriodicTridiag, 3},
    {"_sdetorus_crankNicolson1D", (DL_FUNC) &_sdetorus_crankNicolson1D, 8},
    {"_sdetorus_crankNicolson2D", (DL_FUNC) &_sdetorus_crankNicolson2D, 13},
    {"_sdetorus_driftWn1D", (DL_FUNC) &_sdetorus_driftWn1D, 6},
    {"_sdetorus_driftWn2D", (DL_FUNC) &_sdetorus_driftWn2D, 7},
    {"_sdetorus_euler1D", (DL_FUNC) &_sdetorus_euler1D, 9},
    {"_sdetorus_euler2D", (DL_FUNC) &_sdetorus_euler2D, 10},
    {"_sdetorus_stepAheadWn1D", (DL_FUNC) &_sdetorus_stepAheadWn1D, 10},
    {"_sdetorus_stepAheadWn2D", (DL_FUNC) &_sdetorus_stepAheadWn2D, 11},
    {"_sdetorus_linInterp", (DL_FUNC) &_sdetorus_linInterp, 4},
    {"_sdetorus_besselIExponScaled", (DL_FUNC) &_sdetorus_besselIExponScaled, 4},
    {"_sdetorus_dVmfCpp", (DL_FUNC) &_sdetorus_dVmfCpp, 6},
    {"_sdetorus_clusterProbsVmf", (DL_FUNC) &_sdetorus_clusterProbsVmf, 7},
    {"_sdetorus_weightedMuKappa", (DL_FUNC) &_sdetorus_weightedMuKappa, 5},
    {"_sdetorus_logLikWouPairs", (DL_FUNC) &_sdetorus_logLikWouPairs, 8},
    {"_sdetorus_dWn1D", (DL_FUNC) &_sdetorus_dWn1D, 8},
    {"_sdetorus_dTpdWou1D", (DL_FUNC) &_sdetorus_dTpdWou1D, 11},
    {"_sdetorus_dTpdWou2D", (DL_FUNC) &_sdetorus_dTpdWou2D, 9},
    {"_sdetorus_rTpdWn2D", (DL_FUNC) &_sdetorus_rTpdWn2D, 9},
    {"_sdetorus_dStatWn2D", (DL_FUNC) &_sdetorus_dStatWn2D, 7},
    {"_sdetorus_rStatWn2D", (DL_FUNC) &_sdetorus_rStatWn2D, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_sdetorus(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
