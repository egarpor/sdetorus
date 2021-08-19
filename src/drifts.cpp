#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

// Declaration for safeSoftMax
arma::mat safeSoftMax(arma::mat logs, double expTrc = 30);

//' @title Drift of the WN diffusion in 1D
//'
//' @description Computes the drift of the WN diffusion in 1D in a vectorized way.
//'
//' @inheritParams dWn1D
//' @inheritParams dTpdWou1D
//' @inheritParams safeSoftMax
//' @return A vector of length \code{n} containing the drift evaluated at \code{x}.
//' @examples
//' driftWn1D(x = seq(0, pi, l = 10), alpha = 1, mu = 0, sigma = 1, maxK = 2,
//'           expTrc = 30)
//' @export
// [[Rcpp::export]]
arma::vec driftWn1D(arma::vec x, double alpha, double mu, double sigma, int maxK = 2, double expTrc = 30) {

  // Winding numbers
  int lk = 2 * maxK + 1;
  arma::mat twokpi = arma::linspace<arma::rowvec>(-maxK * 2 * M_PI, maxK * 2 * M_PI, lk);

  // mux
  arma::mat xtwokpi = repmat(x - mu, 1, lk);
  xtwokpi.each_row() += twokpi;

  // logs
  arma::mat w = -square(xtwokpi) / (sigma * sigma / alpha);

  // Weights
  w = safeSoftMax(w, expTrc);

  // Drift
  return alpha * sum(-xtwokpi % w, 1);

}


//' @title Drift of the WN diffusion in 2D
//'
//' @description Computes the drift of the WN diffusion in 2D in a vectorized way.
//'
//' @param A drift matrix of size \code{c(2, 2)}.
//' @inheritParams dTpdWou2D
//' @inheritParams dWn1D
//' @inheritParams safeSoftMax
//' @return A matrix of size \code{c(n, 2)} containing the drift evaluated at \code{x}.
//' @examples
//' alpha <- 3:1
//' mu <- c(0, 0)
//' sigma <- 1:2
//' rho <- 0.5
//' Sigma <- diag(sigma^2)
//' Sigma[1, 2] <- Sigma[2, 1] <- rho * prod(sigma)
//' A <- alphaToA(alpha = alpha, sigma = sigma, rho = rho)
//' x <- rbind(c(0, 1), c(1, 0.1), c(pi, pi), c(-pi, -pi), c(pi / 2, 0))
//' driftWn2D(x = x, A = A, mu = mu, sigma = sigma, rho = rho)
//' driftWn(x = x, A = A, mu = c(0, 0), Sigma = Sigma)
//' @export
// [[Rcpp::export]]
arma::mat driftWn2D(arma::mat x, arma::mat A, arma::vec mu, arma::vec sigma, double rho = 0, int maxK = 2, double expTrc = 30) {

  // Number of pairs
  arma::uword N = x.n_rows;

  // Sequence of winding numbers
  int lk = 2 * maxK + 1;
  arma::vec twokpi = arma::linspace<arma::vec>(-maxK * 2 * M_PI, maxK * 2 * M_PI, lk);

  // Bivariate vector (2 * K1 * M_PI, 2 * K2 * M_PI) for weighting
  arma::vec twokepivec(2);

  // Bivariate vector (2 * K1 * M_PI, 2 * K2 * M_PI) for wrapping
  arma::vec twokapivec(2);

  // Matrix of winding numbers combinations
  arma::mat winds(lk * lk, 2);

  // Create and initialize Sigma
  arma::mat Sigma = diagmat(square(sigma));
  Sigma(1, 0) = Sigma(0, 1) = rho * prod(sigma);

  // Inverse of 1/2 * A^(-1) * Sigma: 2 * Sigma^(-1) * A
  arma::mat invSigmaA = 2 * inv_sympd(Sigma) * A;

  // Log-normalizing constant for the Gaussian with covariance SigmaA
  double lognormconstSigmaA = -log(2 * M_PI) + 0.5 * log(det(invSigmaA));

  /*
   * Weights of the winding numbers for each data point
   */

  // We store the weights in a matrix
  arma::mat weightswinds(N, lk * lk);
  weightswinds.fill(lognormconstSigmaA);

  // Loop in the data
  for (arma::uword i = 0; i < N; i++) {

    // Compute the factors in the exponent that do not depend on the windings
    arma::vec xmu = x.submat(i, 0, i, 1).t() - mu;
    arma::vec xmuinvSigmaA = invSigmaA * xmu;
    double xmuinvSigmaAxmudivtwo = -0.5 * dot(xmuinvSigmaA, xmu);

    // Loop in the winding weight K1
    for (int wek1 = 0; wek1 < lk; wek1++) {

      // 2 * K1 * M_PI
      twokepivec(0) = twokpi(wek1);

      // Compute once the index
      int wekl1 = wek1 * lk;

      // Loop in the winding weight K2
      for (int wek2 = 0; wek2 < lk; wek2++) {

        // 2 * K2 * M_PI
        twokepivec(1) = twokpi(wek2);

        // Decomposition of the exponent
        weightswinds(i, wekl1 + wek2) += xmuinvSigmaAxmudivtwo - dot(xmuinvSigmaA, twokepivec) - dot(invSigmaA * twokepivec, twokepivec) / 2;

        // Matrix of winding number combinations
        winds(wekl1 + wek2, 0) = twokpi(wek1);
        winds(wekl1 + wek2, 1) = twokpi(wek2);

      }

    }

  }

  // Standardize weights
  weightswinds = safeSoftMax(weightswinds, expTrc);

  /*
   * Computation of the drift
   */

  // Create drift
  arma::mat drift = x;

  // Substract mu
  drift.each_row() -= mu.t();

  // Loop in the data
  for (arma::uword i = 0; i < N; i++) {

    // Add the winding effects
    drift.row(i) += weightswinds.row(i) * winds;

  }

  // Multiply drift by matrix and change sign
  drift = - drift * A.t();

  return drift;

}
