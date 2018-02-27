#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

// Declaration for safeSoftMax
arma::mat safeSoftMax(arma::mat logs, double expTrc = 30);


//' @title Approximation of the transition probability density of the WN diffusion in 2D
//'
//' @description Computation of the transition probability density (tpd) for a WN diffusion (with diagonal diffusion matrix)
//'
//' @param x a matrix of dimension \code{c(n, 2)} containing angles. They all must be in \eqn{[\pi,\pi)} so that the truncated wrapping by \code{maxK} windings is able to capture periodicity.
//' @param x0 a matrix of dimension \code{c(n, 2)} containing the starting angles. They all must be in \eqn{[\pi,\pi)}. If all \code{x0} are the same, a matrix of dimension \code{c(1, 2)} can be passed for better performance.
//' @param alpha vector of length \code{3} parametrizing the \code{A} matrix as in \code{\link{alphaToA}}.
//' @param mu a vector of length \code{2} giving the mean.
//' @param sigma vector of length \code{2} containing the \strong{square root} of the diagonal of \eqn{\Sigma}, the diffusion matrix.
//' @param rho correlation coefficient of \eqn{\Sigma}.
//' @inheritParams dTpdWou1D
//' @inheritParams safeSoftMax
//' @return A vector of size \code{n} containing the tpd evaluated at \code{x}.
//' @details The function checks for positive definiteness. If violated, it resets \code{A} such that \code{solve(A) \%*\% Sigma} is positive definite.
//' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
//' @details See Section 3.3 in García-Portugués et al. (2017) for details. See \code{\link{dTpdWou}} for the general case (less efficient for 1D).
//' @references
//' García-Portugués, E., Sorensen, M., Mardia, K. V. and Hamelryck, T. (2017) Langevin diffusions on the torus: estimation and applications. \emph{Stat. Comput.}, \url{https://doi.org/10.1007/s11222-017-9790-2}.
//' @examples
//' set.seed(3455267)
//' alpha <- c(2, 1, -1)
//' sigma <- c(1.5, 2)
//' rho <- 0.9
//' Sigma <- diag(sigma^2)
//' Sigma[1, 2] <- Sigma[2, 1] <- rho * prod(sigma)
//' A <- alphaToA(alpha = alpha, sigma = sigma, rho = rho)
//' solve(A) %*% Sigma
//' mu <- c(pi, 0)
//' x <- t(euler2D(x0 = matrix(c(0, 0), nrow = 1), A = A, mu = mu,
//'                sigma = sigma, N = 500, delta = 0.1)[1, , ])
//' \dontrun{
//' sum(sapply(1:49, function(i) log(dTpdWou(x = matrix(x[i + 1, ], ncol = 2),
//'                                          x0 = x[i, ], t = 1.5, A = A,
//'                                          Sigma = Sigma, mu = mu))))
//' }
//' sum(log(dTpdWou2D(x = matrix(x[2:50, ], ncol = 2),
//'                   x0 = matrix(x[1:49, ], ncol = 2), t = 1.5, alpha = alpha,
//'                   mu = mu, sigma = sigma, rho = rho)))
//' \dontrun{
//' lgrid <- 56
//' grid <- seq(-pi, pi, l = lgrid + 1)[-(lgrid + 1)]
//' image(grid, grid, matrix(dTpdWou(x = as.matrix(expand.grid(grid, grid)),
//'                                  x0 = c(0, 0), t = 0.5, A = A,
//'                                  Sigma = Sigma, mu = mu),
//'                          nrow = 56, ncol = 56), zlim = c(0, 0.25),
//'       main = "dTpdWou")
//' image(grid, grid, matrix(dTpdWou2D(x = as.matrix(expand.grid(grid, grid)),
//'                                    x0 = matrix(0, nrow = 56^2, ncol = 2),
//'                                    t = 0.5, alpha = alpha, sigma = sigma,
//'                                    mu = mu),
//'                          nrow = 56, ncol = 56), zlim = c(0, 0.25),
//'       main = "dTpdWou2D")
//'
//' x <- seq(-pi, pi, l = 100)
//' t <- 1
//' image(x, x, matrix(dTpdWou2D(x = as.matrix(expand.grid(x, x)),
//'                              x0 = matrix(rep(0, 100 * 2), nrow = 100 * 100,
//'                                          ncol = 2),
//'                              t = t, alpha = alpha, mu = mu, sigma = sigma,
//'                              maxK = 2, expTrc = 30),
//'                              nrow = 100, ncol = 100),
//'       zlim = c(0, 0.25))
//' points(stepAheadWn2D(x0 = rbind(c(0, 0)), delta = t / 500,
//'                      A = alphaToA(alpha = alpha, sigma = sigma), mu = mu,
//'                      sigma = sigma, N = 500, M = 1000, maxK = 2,
//'                      expTrc = 30))
//' }
//' @export
// [[Rcpp::export]]
arma::vec dTpdWou2D(arma::mat x, arma::mat x0, arma::vec t, arma::vec alpha, arma::vec mu, arma::vec sigma, double rho = 0, int maxK = 2, double expTrc = 30) {

  /*
   * Create basic objects
   */

  // Number of pairs
  arma::uword N = x.n_rows;

  // Number of initial points
  arma::uword Nx0 = x0.n_rows;

  // Create and initialize A
  double quo = sigma(0) / sigma(1);
  double add = 0.5 * rho * (alpha(1) - alpha(0));
  arma::mat A(2, 2);
  A(0, 0) = alpha(0);
  A(1, 1) = alpha(1);
  A(0, 1) = (alpha(2) + add) * quo;
  A(1, 0) = (alpha(2) - add) / quo;

  // Create and initialize Sigma
  arma::mat Sigma = diagmat(square(sigma));
  Sigma(1, 0) = Sigma(0, 1) = rho * prod(sigma);

  // Sequence of winding numbers
  int lk = 2 * maxK + 1;
  arma::vec twokpi = arma::linspace<arma::vec>(-maxK * 2 * PI, maxK * 2 * PI, lk);

  // Bivariate vector (2 * K1 * PI, 2 * K2 * PI) for weighting
  arma::vec twokepivec(2);

  // Bivariate vector (2 * K1 * PI, 2 * K2 * PI) for wrapping
  arma::vec twokapivec(2);

  /*
   * Check if the t is common
   */

  int commonTime;
  if (t.n_elem == 1) {

    commonTime = 1;

  }else if (t.n_elem == N) {

    commonTime = 0;

  } else {

    stop("Length of t is neither 1 nor N");

  }

  /*
   * Check for symmetry and positive definiteness of A^(-1) * Sigma
   */

  // Only positive definiteness can be violated with the parametrization of A
  double prodDiagonal = add * add + alpha(0) * alpha(1);
  double testPosDef = prodDiagonal - alpha(2) * alpha(2);

  // Check positive definiteness
  if (testPosDef <= 0) {

    // Update alpha(2) such that testPosDef > 0
    alpha(2) = std::signbit(alpha(2)) * sqrt(prodDiagonal) * 0.999;

    // Reset A to a matrix with positive determinant
    A(0, 1) = (alpha(2) + add) * quo;
    A(1, 0) = (alpha(2) - add) / quo;

    // Update testPosDef
    testPosDef = prodDiagonal * 1e-3;

  }

  // A^(-1) * Sigma
  arma::mat invASigma(2, 2);
  invASigma(0, 0) = Sigma(0, 0) * (alpha(1) - rho * (add + alpha(2)));
  invASigma(0, 1) = invASigma(1, 0) = prod(sigma) * (rho * alpha(1) - add - alpha(2));
  invASigma(1, 1) = Sigma(1, 1) * (alpha(0) + rho * (add - alpha(2)));
  invASigma /= testPosDef;

  // Inverse of (1/2 * A^(-1) * Sigma): 2 * Sigma^(-1) * A
  arma::mat invSigmaA(2, 2);
  invSigmaA(0, 0) = invASigma(1, 1);
  invSigmaA(0, 1) = invSigmaA(1, 0) = -invASigma(0, 1);
  invSigmaA(1, 1) = invASigma(0, 0);
  double detInvSigmaA = testPosDef / det(Sigma);
  invSigmaA *= 2 * detInvSigmaA;

  // For normalizing constants
  double l2pi = log(2 * PI);

  // Log-normalizing constant for the Gaussian with covariance SigmaA
  double lognormconstSigmaA = -l2pi + 0.5 * log(4 * detInvSigmaA);

  // For calling log_det
  double sign = 1;

  /*
   * Computation of Gammat and exp(-t * A) analytically
   */

  // Quantities for computing exp(-t * A)
  double s = 0.5 * sum(alpha.head(2));
  double q = sqrt(fabs((alpha(0) - s) * (alpha(1) - s) - A(0, 1) * A(1, 0)));

  // Avoid indetermination in sinh(q * t) / q when q == 0
  if (q == 0) {

    q = 1e-6;

  }

  // s1(-t) and s2(-t)
  arma::vec est = exp(-s * t);
  arma::vec eqt = exp(q * t);
  arma::vec eqtinv = 1 / eqt;
  arma::vec cqt = 0.5 * (eqt + eqtinv);
  arma::vec sqt = (eqt - eqtinv) / (2 * q);
  arma::vec s1t = est % (cqt + s * sqt);
  arma::vec s2t = -est % sqt;

  // s1(-2t) and s2(-2t)
  est = est % est;
  eqt = eqt % eqt;
  eqtinv = eqtinv % eqtinv;
  cqt = 0.5 * (eqt + eqtinv);
  sqt = (eqt - eqtinv) / (2 * q);
  arma::vec s12t = est % (cqt + s * sqt);
  arma::vec s22t = -est % sqt;

  /*
   * Weights of the winding numbers for each data point
   */

  // We store the weights in a matrix to skip the null later in the computation of the tpd
  arma::mat weightswindsinitial(Nx0, lk * lk);
  weightswindsinitial.fill(lognormconstSigmaA);

  // Loop in the data
  for (arma::uword i = 0; i < Nx0; i++) {

    // Compute the factors in the exponent that do not depend on the windings
    arma::vec xmu = x0.row(i).t() - mu;
    arma::vec xmuinvSigmaA = invSigmaA * xmu;
    double xmuinvSigmaAxmudivtwo = -0.5 * dot(xmuinvSigmaA, xmu);

    // Loop in the winding weight K1
    for (int wek1 = 0; wek1 < lk; wek1++) {

      // 2 * K1 * PI
      twokepivec(0) = twokpi(wek1);

      // Compute once the index
      int wekl1 = wek1 * lk;

      // Loop in the winding weight K2
      for (int wek2 = 0; wek2 < lk; wek2++) {

        // 2 * K2 * PI
        twokepivec(1) = twokpi(wek2);

        // Negative exponent
        weightswindsinitial(i, wekl1 + wek2) += xmuinvSigmaAxmudivtwo - dot(xmuinvSigmaA, twokepivec) - 0.5 * dot(invSigmaA * twokepivec, twokepivec);

      }

    }

  }

  // Standardize weights for the tpd
  weightswindsinitial = safeSoftMax(weightswindsinitial, expTrc);

  /*
   * Computation of the tpd: wrapping + weighting
   */

  // The evaluations of the tpd are stored in a vector, no need to keep track of wrappings
  arma::vec tpdfinal(N); tpdfinal.zeros();

  // Variables inside the commonTime if-block
  arma::mat ExptiA(2, 2), invGammati(2, 2);
  double logDetInvGammati = 0, lognormconstGammati = 0;

  // If t is common, compute once
  if (commonTime == 1) {

    // Exp(-ti * A)
    ExptiA = s2t(0) * A;
    ExptiA.diag() += s1t(0);

    // Inverse and log-normalizing constant for the Gammat
    invGammati = 2 * inv_sympd((1 - s12t(0)) * invASigma - s22t(0) * Sigma);

    // Log-determinant of invGammati (assumed to be positive)
    arma::log_det(logDetInvGammati, sign, invGammati);

    // Log-normalizing constant for the Gaussian with covariance Gammati
    lognormconstGammati = -l2pi + 0.5 * logDetInvGammati;

  }

  // Variables inside the Nx0 if-block
  arma::vec x00(Nx0), muti(2);

  // If there is only one x0
  if (Nx0 == 1) {

    // Initial point x0 varying with i
    x00 = x0.row(0).t();

    // Common muti
    muti = mu + ExptiA * (x00 - mu);

    // Fill the rest of weightswindsinitial
    weightswindsinitial = arma::repmat(weightswindsinitial.row(0), N, 1);

  }

  // Loop in the data
  for (arma::uword i = 0; i < N; i++) {

    // Evaluation point x varying with i
    arma::vec xx = x.row(i).t();

    // If t is not common
    if (commonTime == 0) {

      // Exp(-ti * A)
      ExptiA = s2t(i) * A;
      ExptiA.diag() += s1t(i);

      // Inverse and log-normalizing constant for the Gammati
      invGammati = 2 * inv_sympd((1 - s12t(i)) * invASigma - s22t(i) * Sigma);

      // Log-determinant of invGammati (assumed to be positive)
      arma::log_det(logDetInvGammati, sign, invGammati);

      // Log-normalizing constant for the Gaussian with covariance Gammati
      lognormconstGammati = -l2pi + 0.5 * logDetInvGammati;

    }

    // If there are several x0's
    if (Nx0 > 1) {

      // Initial point x0 varying with i
      x00 = x0.row(i).t();

      // Common muti
      muti = mu + ExptiA * (x00 - mu);

    }

    // Loop on the winding weight K1
    for (int wek1 = 0; wek1 < lk; wek1++) {

      // 2 * K1 * PI
      twokepivec(0) = twokpi(wek1);

      // Compute once the index
      int wekl1 = wek1 * lk;

      // Loop on the winding weight K2
      for (int wek2 = 0; wek2 < lk; wek2++) {

        // Skip zero weights
        if (weightswindsinitial(i, wekl1 + wek2) > 0) {

          // 2 * K1 * PI
          twokepivec(1) = twokpi(wek2);

          // Compute the factors in the exponent that do not depend on the windings
          arma::vec xmuti = xx - muti - ExptiA * twokepivec;
          arma::vec xmutiInvGammati = invGammati * xmuti;
          double xmutiInvGammatixmutidiv2 = -0.5 * dot(xmutiInvGammati, xmuti);

          // Loop in the winding wrapping K1
          for (int wak1 = 0; wak1 < lk; wak1++) {

            // 2 * K1 * PI
            twokapivec(0) = twokpi(wak1);

            // Loop in the winding wrapping K2
            for (int wak2 = 0; wak2 < lk; wak2++) {

              // 2 * K2 * PI
              twokapivec(1) = twokpi(wak2);

              // Decomposition of the negative exponent
              double exponent = xmutiInvGammatixmutidiv2 - dot(xmutiInvGammati, twokapivec) - 0.5 * dot(invGammati * twokapivec, twokapivec) + lognormconstGammati;

              // Tpd
              tpdfinal(i) += exp(exponent) * weightswindsinitial(i, wekl1 + wek2);

            }

          }

        }

      }

    }

  }

  return tpdfinal;

}


//' @title Simulation from the approximated transition distribution of a WN diffusion in 2D
//'
//' @description Simulates from the approximate transition density of the WN diffusion in 2D.
//'
//' @param n sample size.
//' @param x0 a matrix of dimension \code{c(nx0, 2)} giving the starting values.
//' @param t vector of length \code{nx0} containing the times between observations.
//' @inheritParams dTpdWou2D
//' @inheritParams safeSoftMax
//' @return An array of dimension \code{c(n, 2, nx0)} containing the \code{n} samples of the transition distribution,
//' conditioned on that the process was at \code{x0} at \code{t} instants ago. The samples are all in \eqn{[\pi,\pi)}.
//' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
//' @examples
//' alpha <- c(3, 2, -1)
//' sigma <- c(0.5, 1)
//' mu <- c(pi, pi)
//' x <- seq(-pi, pi, l = 100)
//' t <- 0.5
//' image(x, x, matrix(dTpdWou2D(x = as.matrix(expand.grid(x, x)),
//'                             x0 = matrix(rep(0, 100 * 2),
//'                                         nrow = 100 * 100, ncol = 2),
//'                             t = t, mu = mu, alpha = alpha, sigma = sigma,
//'                             maxK = 2, expTrc = 30), nrow = 100, ncol = 100),
//'       zlim = c(0, 0.5))
//' points(rTpdWn2D(n = 500, x0 = rbind(c(0, 0)), t = t, mu = mu, alpha = alpha,
//'                 sigma = sigma)[, , 1], col = 3)
//' points(stepAheadWn2D(x0 = rbind(c(0, 0)), delta = t / 500,
//'                      A = alphaToA(alpha = alpha, sigma = sigma),
//'                      mu = mu, sigma = sigma, N = 500, M = 500, maxK = 2,
//'                      expTrc = 30), col = 4)
//' @export
// [[Rcpp::export]]
arma::cube rTpdWn2D(arma::uword n, arma::mat x0, arma::vec t, arma::vec mu, arma::vec alpha, arma::vec sigma, double rho = 0, int maxK = 2, double expTrc = 30) {

  /*
   * Create basic objects
   */

  // Number of different starting angles
  arma::uword nx0 = x0.n_rows;

  // Create and initialize A
  double quo = sigma(0) / sigma(1);
  double add = 0.5 * rho * (alpha(1) - alpha(0));
  arma::mat A(2, 2);
  A(0, 0) = alpha(0);
  A(1, 1) = alpha(1);
  A(0, 1) = (alpha(2) + add) * quo;
  A(1, 0) = (alpha(2) - add) / quo;

  // Create and initialize Sigma
  arma::mat Sigma = diagmat(square(sigma));
  Sigma(1, 0) = Sigma(0, 1) = rho * prod(sigma);

  // Sequence of winding numbers
  int lk = 2 * maxK + 1;
  arma::vec twokpi = arma::linspace<arma::vec>(-maxK * 2 * PI, maxK * 2 * PI, lk);

  // Bivariate vector (2 * K1 * PI, 2 * K2 * PI) for weighting
  arma::vec twokepivec(2);

  // Bivariate vector (2 * K1 * PI, 2 * K2 * PI) for wrapping
  arma::vec twokapivec(2);

  /*
   * Check if the t is common
   */

  int commonTime;
  if (t.n_elem == 1) {

    commonTime = 1;

  }else if (t.n_elem == nx0) {

    commonTime = 0;

  } else {

    stop("Length of t is neither 1 nor nx0");

  }

  /*
   * Check for symmetry and positive definiteness of A^(-1) * Sigma
   */

  // Only positive definiteness can be violated with the parametrization of A
  double prodDiagonal = add * add + alpha(0) * alpha(1);
  double testPosDef = prodDiagonal - alpha(2) * alpha(2);

  // Check positive definiteness
  if (testPosDef <= 0) {

    // Update alpha(2) such that testPosDef > 0
    alpha(2) = std::signbit(alpha(2)) * sqrt(prodDiagonal) * 0.999;

    // Reset A to a matrix with positive determinant
    A(0, 1) = (alpha(2) + add) * quo;
    A(1, 0) = (alpha(2) - add) / quo;

    // Update testPosDef
    testPosDef = prodDiagonal * 1e-3;

  }

  // A^(-1) * Sigma
  arma::mat invASigma(2, 2);
  invASigma(0, 0) = Sigma(0, 0) * (alpha(1) - rho * (add + alpha(2)));
  invASigma(0, 1) = invASigma(1, 0) = prod(sigma) * (rho * alpha(1) - add - alpha(2));
  invASigma(1, 1) = Sigma(1, 1) * (alpha(0) + rho * (add - alpha(2)));
  invASigma /= testPosDef;

  // Inverse of (1/2 * A^(-1) * Sigma): 2 * Sigma^(-1) * A
  arma::mat invSigmaA(2, 2);
  invSigmaA(0, 0) = invASigma(1, 1);
  invSigmaA(0, 1) = invSigmaA(1, 0) = -invASigma(0, 1);
  invSigmaA(1, 1) = invASigma(0, 0);
  double detInvSigmaA = testPosDef / det(Sigma);
  invSigmaA *= 2 * detInvSigmaA;

  // For normalizing constants
  double l2pi = log(2 * PI);

  // Log-normalizing constant for the Gaussian with covariance SigmaA
  double lognormconstSigmaA = -l2pi + 0.5 * log(4 * detInvSigmaA);

  /*
   * Computation of Gammat and exp(-t * A) analytically
   */

  // Quantities for computing exp(-t * A)
  double s = 0.5 * sum(alpha.head(2));
  double q = sqrt(fabs((alpha(0) - s) * (alpha(1) - s) - A(0, 1) * A(1, 0)));

  // Avoid indetermination in sinh(q * t) / q when q == 0
  if (q == 0) {

    q = 1e-6;

  }

  // s1(-t) and s2(-t)
  arma::vec est = exp(-s * t);
  arma::vec eqt = exp(q * t);
  arma::vec eqtinv = 1 / eqt;
  arma::vec cqt = 0.5 * (eqt + eqtinv);
  arma::vec sqt = (eqt - eqtinv) / (2 * q);
  arma::vec s1t = est % (cqt + s * sqt);
  arma::vec s2t = -est % sqt;

  // s1(-2t) and s2(-2t)
  est = est % est;
  eqt = eqt % eqt;
  eqtinv = eqtinv % eqtinv;
  cqt = 0.5 * (eqt + eqtinv);
  sqt = (eqt - eqtinv) / (2 * q);
  arma::vec s12t = est % (cqt + s * sqt);
  arma::vec s22t = -est % sqt;

  /*
   * Weights of the winding numbers for each data point
   */

  // We store the weights in a matrix to skip the null later in the computation of the tpd
  arma::mat weightswindsinitial(nx0, lk * lk);
  weightswindsinitial.fill(lognormconstSigmaA);

  // Loop in the data
  for (arma::uword i = 0; i < nx0; i++) {

    // Compute the factors in the exponent that do not depend on the windings
    arma::vec xmu = x0.row(i).t() - mu;
    arma::vec xmuinvSigmaA = invSigmaA * xmu;
    double xmuinvSigmaAxmudivtwo = -0.5 * dot(xmuinvSigmaA, xmu);

    // Loop in the winding weight K1
    for (int wek1 = 0; wek1 < lk; wek1++) {

      // 2 * K1 * PI
      twokepivec(0) = twokpi(wek1);

      // Compute once the index
      int wekl1 = wek1 * lk;

      // Loop in the winding weight K2
      for (int wek2 = 0; wek2 < lk; wek2++) {

        // 2 * K2 * PI
        twokepivec(1) = twokpi(wek2);

        // Negative exponent
        weightswindsinitial(i, wekl1 + wek2) += xmuinvSigmaAxmudivtwo - dot(xmuinvSigmaA, twokepivec) - 0.5 * dot(invSigmaA * twokepivec, twokepivec);

      }

    }

  }

  // Standardize weights for the tpd
  weightswindsinitial = safeSoftMax(weightswindsinitial, expTrc);

  /*
   * Sampling
   */

  // Probabilities of windings
  arma::mat probs = arma::cumsum(weightswindsinitial, 1);

  // Sample uniforms in [0, 1]. There is no weigthed sampling in Armadillo!
  arma::vec runif = arma::randu(n * nx0);

  // Matrix of muti for each x
  arma::mat mutix(n, 2);

  // Sample of independent N(0, 1)
  arma::cube x = arma::randn(n, 2, nx0);

  // Variables inside the commonTime if-block
  arma::mat ExptiA(2, 2), Gammati(2, 2), ch(2, 2);

  // If t is common, compute once
  if (commonTime == 1) {

    // Exp(-ti * A)
    ExptiA = s2t(0) * A;
    ExptiA.diag() += s1t(0);

    // Gammati
    Gammati = 0.5 * ((1 - s12t(0)) * invASigma - s22t(0) * Sigma);

    // Cholesky decomposition for correlate random sample
    ch = chol(Gammati);

  }

  // Loop throught the x0s
  for (arma::uword i = 0; i < nx0; i++) {

    // If t is not common
    if (commonTime == 0) {

      // Exp(-ti * A)
      ExptiA = s2t(i) * A;
      ExptiA.diag() += s1t(i);

      // Gammati
      Gammati = 0.5 * ((1 - s12t(i)) * invASigma - s22t(i) * Sigma);

      // Cholesky decomposition for correlate random sample
      ch = chol(Gammati);

    }

    // Common muti
    mutix.each_row() = (mu + ExptiA * (x0.row(i).t() - mu)).t();

    // Compute once the index
    arma::uword li = i * n;

    // Loop in the number of replicates
    for (arma::uword m = 0; m < n; m++) {

      // Choose windings with probabilities weightswindsinitial
      // Loop in the winding weight K1
      for (int wek1 = 0; wek1 < lk; wek1++) {

        // Compute once the index
        int wekl1 = wek1 * lk;

        // Loop in the winding weight K2
        for (int wek2 = 0; wek2 < lk; wek2++) {

          // Weighted sampling
          if (runif(li + m) <= probs(i, wekl1 + wek2)) {

            // 2 * K1 * PI
            twokepivec(0) = twokpi(wek1);

            // 2 * K2 * PI
            twokepivec(1) = twokpi(wek2);

            // Set mut for x
            mutix.row(m) += (ExptiA * twokepivec).t();

            // Skip the wek1 and wek2 loops
            goto sample;

          }

        }

      }

      // Goto identifier
      sample: ;

    }

    // Correlate random variables
    x.slice(i) *= ch;

    // Recenter depending on mut
    x.slice(i) += mutix;

  }

  // Wrap (convert to [-PI,PI) x [-PI,PI))
  x -= floor((x + PI) / (2 * PI)) * (2 * PI);

  return x;

}


//' @title Stationary density of a WN diffusion (with diagonal diffusion matrix) in 2D
//'
//' @description Stationary density of the WN diffusion.
//'
//' @inheritParams dTpdWou2D
//' @inheritParams safeSoftMax
//' @return A vector of size \code{n} containing the stationary density evaluated at \code{x}.
//' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
//' @examples
//' set.seed(345567)
//' alpha <- c(2, 1, -1)
//' sigma <- c(1.5, 2)
//' Sigma <- diag(sigma^2)
//' A <- alphaToA(alpha = alpha, sigma = sigma)
//' mu <- c(pi, pi)
//' dStatWn2D(x = toPiInt(matrix(1:20, nrow = 10, ncol = 2)), mu = mu,
//'           alpha = alpha, sigma = sigma)
//' dTpdWou(t = 10, x = toPiInt(matrix(1:20, nrow = 10, ncol = 2)), A = A,
//'          mu = mu, Sigma = Sigma, x0 = mu)
//' xth <- seq(-pi, pi, l = 100)
//' contour(xth, xth, matrix(dStatWn2D(x = as.matrix(expand.grid(xth, xth)),
//'                                    alpha = alpha, sigma = sigma, mu = mu),
//'                          nrow = length(xth), ncol = length(xth)), nlevels = 50)
//' points(rStatWn2D(n = 1000, mu = mu, alpha = alpha, sigma = sigma), col = 2)
//' @export
// [[Rcpp::export]]
arma::vec dStatWn2D(arma::mat x, arma::vec alpha, arma::vec mu, arma::vec sigma, double rho = 0, int maxK = 2, double expTrc = 30) {

  /*
   * Create basic objects
   */

  // Number of pairs
  arma::uword N = x.n_rows;

  // Create and initialize A
  double quo = sigma(0) / sigma(1);
  double add = 0.5 * rho * (alpha(1) - alpha(0));
  arma::mat A(2, 2);
  A(0, 0) = alpha(0);
  A(1, 1) = alpha(1);
  A(0, 1) = (alpha(2) + add) * quo;
  A(1, 0) = (alpha(2) - add) / quo;

  // Create and initialize Sigma
  arma::mat Sigma = diagmat(square(sigma));
  Sigma(1, 0) = Sigma(0, 1) = rho * prod(sigma);

  // Sequence of winding numbers
  int lk = 2 * maxK + 1;
  arma::vec twokpi = arma::linspace<arma::vec>(-maxK * 2 * PI, maxK * 2 * PI, lk);

  // Bivariate vector (2 * K1 * PI, 2 * K2 * PI) for weighting
  arma::vec twokepivec(2);

  /*
   * Check for symmetry and positive definiteness of A^(-1) * Sigma
   */

  // Only positive definiteness can be violated with the parametrization of A
  double prodDiagonal = add * add + alpha(0) * alpha(1);
  double testPosDef = prodDiagonal - alpha(2) * alpha(2);

  // Check positive definiteness
  if (testPosDef <= 0) {

    // Update alpha(2) such that testPosDef > 0
    alpha(2) = std::signbit(alpha(2)) * sqrt(prodDiagonal) * 0.999;

    // Reset A to a matrix with positive determinant
    A(0, 1) = (alpha(2) + add) * quo;
    A(1, 0) = (alpha(2) - add) / quo;

    // Update testPosDef
    testPosDef = prodDiagonal * 1e-3;

  }

  // Inverse of 1/2 * A^(-1) * Sigma: 2 * Sigma^(-1) * A
  arma::mat invSigmaA = 2 * inv_sympd(Sigma) * A;

  // For normalizing constants
  double l2pi = log(2 * PI);

  // Log-determinant of invSigmaA (assumed to be positive)
  double logDetInvSigmaA, sign;
  arma::log_det(logDetInvSigmaA, sign, invSigmaA);

  // Log-normalizing constant for the Gaussian with covariance SigmaA
  double lognormconstSigmaA = -l2pi + 0.5 * logDetInvSigmaA;

  /*
   * Evaluation of the density reusing the code from the weights of the winding numbers
   * in dTpdWou2D for each data point. Here we sum all the unstandarized weights
   * for each data point.
   */

  // We store the weights in a matrix to skip the null later in the computation of the tpd
  arma::mat weightswindsinitial(N, lk * lk);

  // Loop in the data
  for (arma::uword i = 0; i < N; i++) {

    // Compute the factors in the exponent that do not depend on the windings
    arma::vec xmu = x.row(i).t() - mu;
    arma::vec xmuinvSigmaA = invSigmaA * xmu;
    double xmuinvSigmaAxmudivtwo = 0.5 * dot(xmuinvSigmaA, xmu);

    // Loop in the winding weight K1
    for (int wek1 = 0; wek1 < lk; wek1++) {

      // 2 * K1 * PI
      twokepivec(0) = twokpi(wek1);

      // Compute once the index
      int wekl1 = wek1 * lk;

      // Loop in the winding weight K2
      for (int wek2 = 0; wek2 < lk; wek2++) {

        // 2 * K2 * PI
        twokepivec(1) = twokpi(wek2);

        // Decomposition of the exponent
        double exponent = xmuinvSigmaAxmudivtwo + dot(xmuinvSigmaA, twokepivec) + 0.5 * dot(invSigmaA * twokepivec, twokepivec) - lognormconstSigmaA;

        // Truncate the negative exponential
        if (exponent > expTrc) {

          weightswindsinitial(i, wekl1 + wek2) = 0;

        } else {

          weightswindsinitial(i, wekl1 + wek2) = exp(-exponent);

        }

      }

    }

  }

  // The density is the sum of the weights
  arma::vec dens = sum(weightswindsinitial, 1);

  return dens;

}


//' @title Simulation from the stationary density of a WN diffusion in 2D
//'
//' @description Simulates from the stationary density of the WN diffusion in 2D.
//'
//' @param n sample size.
//' @inheritParams dTpdWou2D
//' @return A matrix of dimension \code{c(n, 2)} containing the samples from the stationary distribution.
//' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
//' @examples
//' set.seed(345567)
//' alpha <- c(2, 1, -1)
//' sigma <- c(1.5, 2)
//' Sigma <- diag(sigma^2)
//' A <- alphaToA(alpha = alpha, sigma = sigma)
//' mu <- c(pi, pi)
//' plot(rStatWn2D(n = 1000, mu = mu, alpha = alpha, sigma = sigma))
//' points(toPiInt(mvtnorm::rmvnorm(n = 1000, mean = mu,
//'                                 sigma = solve(A) %*% Sigma / 2,
//'                                 method = "chol")), col = 2)
//' @export
// [[Rcpp::export]]
arma::mat rStatWn2D(arma::uword n, arma::vec mu, arma::vec alpha, arma::vec sigma, double rho = 0) {

  /*
   * Compute A^(-1) * Sigma and check for positive definiteness
   */

  // Create and initialize A
  double quo = sigma(0) / sigma(1);
  double add = 0.5 * rho * (alpha(1) - alpha(0));
  arma::mat A(2, 2);
  A(0, 0) = alpha(0);
  A(1, 1) = alpha(1);
  A(0, 1) = (alpha(2) + add) * quo;
  A(1, 0) = (alpha(2) - add) / quo;

  // Create and initialize Sigma
  arma::mat Sigma = diagmat(square(sigma));
  Sigma(1, 0) = Sigma(0, 1) = rho * prod(sigma);

  // Only positive definiteness can be violated with the parametrization of A
  double prodDiagonal = add * add + alpha(0) * alpha(1);
  double testPosDef = prodDiagonal - alpha(2) * alpha(2);

  // Check positive definiteness
  if (testPosDef <= 0) {

    // Update alpha(2) such that testPosDef > 0
    alpha(2) = std::signbit(alpha(2)) * sqrt(prodDiagonal) * 0.999;

    // Reset A to a matrix with positive determinant
    A(0, 1) = (alpha(2) + add) * quo;
    A(1, 0) = (alpha(2) - add) / quo;

    // Update testPosDef
    testPosDef = prodDiagonal * 1e-3;

  }

  // A^(-1) * Sigma / 2
  arma::mat invASigma(2, 2);
  invASigma(0, 0) = Sigma(0, 0) * (alpha(1) - rho * (add + alpha(2)));
  invASigma(0, 1) = invASigma(1, 0) = prod(sigma) * (rho * alpha(1) - add - alpha(2));
  invASigma(1, 1) = Sigma(1, 1) * (alpha(0) + rho * (add - alpha(2)));
  invASigma /= 2 * testPosDef;

  /*
   * Sample correlated normal variables
   */

  // Sample of independent N(0, 1)
  arma::mat x = arma::randn(n, 2);

  // Cholesky decomposition for correlate random sample
  arma::mat ch = chol(invASigma);

  // Correlate random variables
  x = x * ch;

  // Recenter
  x.each_row() += mu.t();

  // Wrap (convert to [-PI,PI) x [-PI,PI))
  x -= floor((x + PI) / (2 * PI)) * (2 * PI);

  return x;

}
