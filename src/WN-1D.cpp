#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

// Declaration for safeSoftMax
arma::mat safeSoftMax(arma::mat logs, double expTrc = 30);

//' @title WN density in 1D
//'
//' @description Computation of the WN density in 1D.
//'
//' @param x a vector of length \code{n} containing angles. They all must be in \eqn{[\pi,\pi)} so that the truncated wrapping by \code{maxK} windings is able to capture periodicity.
//' @param mu mean parameter. Must be in \eqn{[\pi,\pi)}.
//' @param sigma diffusion coefficient.
//' @param maxK maximum absolute value of the windings considered in the computation of the WN.
//' @inheritParams safeSoftMax
//' @param vmApprox whether to use the von Mises approximation to a wrapped normal (\code{1}) or not (\code{0}, default).
//' @param kt concentration for the von Mises, a suitable output from \code{\link{momentMatchWnVm}} (see examples).
//' @param logConstKt the logarithm of the von Mises normalizing constant associated to the concentration \code{kt} (see examples)
//' @return A vector of size \code{n} containing the density evaluated at \code{x}.
//' @examples
//' mu <- 0
//' sigma <- 1
//' dWn1D(x = seq(-pi, pi, l = 10), mu = mu, sigma = sigma, vmApprox = 0)
//'
//' # von Mises approximation
//' kt <- scoreMatchWnVm(sigma2 = sigma^2)
//' dWn1D(x = seq(-pi, pi, l = 10), mu = mu, sigma = sigma, vmApprox = 1, kt = kt,
//'       logConstKt = -log(2 * pi * besselI(x = kt, nu = 0, expon.scaled = TRUE)))
//' @export
// [[Rcpp::export]]
arma::vec dWn1D(arma::vec x, double mu, double sigma, int maxK = 2, double expTrc = 30, int vmApprox = 0, double kt = 0, double logConstKt = 0) {

  // Number of x's
  arma::uword N = x.n_elem;

  // Create result
  arma::vec densfinal(N); densfinal.zeros();

  // Winding numbers
  int lk = 2 * maxK + 1;
  arma::rowvec twokpi = arma::linspace<arma::rowvec>(-maxK * 2 * M_PI, maxK * 2 * M_PI, lk);

  // 2 * variance and inverse
  double sd2 = 2 * sigma * sigma;
  double isd2 = 1 / sd2;

  // Center x
  x -= mu;

  // Log-normalizing constant
  double lognormconstsd2 = -0.5 * log(M_PI * sd2);

  /*
   * Computation of the density: wrapping
   */

  // von Mises approximation?
  if (vmApprox == 1) {

    // Density
    densfinal = exp(logConstKt + kt * (cos(x) - 1));

  } else {

    // Loop in the winding wrapping
    for (int wak = 0; wak < lk; wak++) {

      // Decomposition of the negative exponent
      arma::vec exponent = x + twokpi(wak);
      exponent = -exponent % exponent * isd2 + lognormconstsd2;

      // Density
      densfinal += exp(exponent);

    }

  }

  return densfinal;

}


//' @title Approximation of the transition probability density of the WN diffusion in 1D
//'
//' @description Computation of the transition probability density (tpd) for a WN diffusion.
//'
//' @inheritParams dWn1D
//' @param x0 a vector of length \code{n} containing the starting angles. They all must be in \eqn{[\pi,\pi)}.
//' @param t a scalar containing the times separating \code{x} and \code{x0}.
//' @param alpha drift parameter.
//' @param sigma diffusion coefficient.
//' @inheritParams safeSoftMax
//' @return A vector of size \code{n} containing the tpd evaluated at \code{x}.
//' @details See Section 3.3 in García-Portugués et al. (2019) for details. See \code{\link{dTpdWou}} for the general case (less efficient for 2D).
//' @references
//' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \doi{10.1007/s11222-017-9790-2}
//' @examples
//' t <- 0.5
//' alpha <- 1
//' mu <- 0
//' sigma <- 1
//' x0 <- 0.1
//' dTpdWou1D(x = seq(-pi, pi, l = 10), x0 = rep(x0, 10), t = t, alpha = alpha,
//'           mu = mu, sigma = sigma, vmApprox = 0)
//'
//' # von Mises approximation
//' kt <- scoreMatchWnVm(sigma2 = sigma^2 * (1 - exp(-2 * alpha * t)) / (2 * alpha))
//' dTpdWou1D(x = seq(-pi, pi, l = 10), x0 = rep(x0, 10), t = t, alpha = alpha,
//'           mu = mu, sigma = sigma, vmApprox = 1, kt = kt,
//'           logConstKt = -log(2 * pi * besselI(x = kt, nu = 0,
//'                                              expon.scaled = TRUE)))
//' @export
// [[Rcpp::export]]
arma::vec dTpdWou1D(arma::vec x, arma::vec x0, double t, double alpha, double mu, double sigma, int maxK = 2, double expTrc = 30, int vmApprox = 0, double kt = 0, double logConstKt = 0) {

  // Number of x's
  arma::uword N = x.n_elem;

  // Create result
  arma::vec tpdfinal(N); tpdfinal.zeros();

  // Winding numbers
  int lk = 2 * maxK + 1;
  arma::rowvec twokpi = arma::linspace<arma::rowvec>(-maxK * 2 * M_PI, maxK * 2 * M_PI, lk);

  // Stationary variance (x 2)
  double sd2 = sigma * sigma / alpha;

  // Center x0
  x0 -= mu;

  // Dependent variance (x 2)
  double sd2t = sigma * sigma * (1 - exp(-2 * alpha * t)) / alpha;

  // Log-normalizing constant
  double lognormconstsd2t = -0.5 * log(M_PI * sd2t);

  /*
   * Weights of the winding numbers for each data point
   */

  // We store the weights in a matrix to skip the null later in the computation of the tpd
  arma::mat weightswindsinitial(N, lk);

  // Loop in the data
  for (arma::uword i = 0; i < N; i++) {

    // Negative exponent
    weightswindsinitial.row(i) = -square(x0(i) + twokpi) / sd2;

  }

  // Standardize weights for the tpd
  weightswindsinitial = safeSoftMax(weightswindsinitial, expTrc);

  /*
   * Computation of the tpd: wrapping + weighting
   */

  // Loop in the data
  for (arma::uword i = 0; i < N; i++) {

     // Loop on the winding weight
    for (int wek = 0; wek < lk; wek++) {

      // Skip zero weights
      if (weightswindsinitial(i, wek) > 0) {

        // mut
        double mut = mu + exp(-t * alpha) * (x0(i) + twokpi(wek));

        // von Mises approximation?
        if (vmApprox == 1) {

          // Decomposition of the exponent
          double exponent = logConstKt + kt * (cos(x(i) - mut) - 1);

          // Tpd
          tpdfinal(i) += exp(exponent) * weightswindsinitial(i, wek);

        } else {

          // Loop in the winding wrapping
          for (int wak = 0; wak < lk; wak++) {

            // Decomposition of the negative exponent
            double dif = x(i) + twokpi(wak) - mut;
            double exponent = -dif * dif / sd2t + lognormconstsd2t;

            // Tpd
            tpdfinal(i) += exp(exponent) * weightswindsinitial(i, wek);

          }

        }

      }

    }

  }

  return tpdfinal;

}
