#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

// Declaration for safeSoftMax
arma::mat safeSoftMax(arma::mat logs, double expTrc = 30);

//' @title Loglikelihood of WN in 2D when only the initial and final points are observed
//'
//' @description Computation of the loglikelihood for a WN diffusion (with diagonal diffusion matrix) from a sample of initial and final pairs of angles.
//'
//' @param x a matrix of dimension \code{c(n, 4)} of initial and final pairs of angles. Each row is an observation containing \eqn{(\phi_0, \psi_0, \phi_t, \psi_t)}.
//' They all must be in \eqn{[\pi,\pi)} so that the truncated wrapping by \code{maxK} windings is able to capture periodicity.
//' @param t either a scalar or a vector of length \code{n} containing the times the initial and final dihedrals. If \code{t} is a scalar, a common time is assumed.
//' @inheritParams dTpdWou2D
//' @inheritParams safeSoftMax
//' @inheritParams dWn1D
//' @return A scalar giving the final loglikelihood, defined as the sum of the loglikelihood of the initial angles according to the stationary density
//' and the loglikelihood of the transitions from initial to final angles.
//' @details A negative penalty is added if positive definiteness is violated. If the output value is Inf, -100 * N is returned instead.
//' @examples
//' set.seed(345567)
//' x <- toPiInt(matrix(rnorm(200, mean = pi), ncol = 4, nrow = 50))
//' alpha <- c(2, 1, -0.5)
//' mu <- c(0, pi)
//' sigma <- sqrt(c(2, 1))
//'
//' # The same
//' logLikWouPairs(x = x, t = 0.5, alpha = alpha, mu = mu, sigma = sigma)
//' sum(
//'   log(dStatWn2D(x = x[, 1:2], alpha = alpha, mu = mu, sigma = sigma)) +
//'   log(dTpdWou2D(x = x[, 3:4], x0 = x[, 1:2], t = 0.5, alpha = alpha, mu = mu,
//'                  sigma = sigma))
//' )
//'
//' # Different times
//' logLikWouPairs(x = x, t = (1:50) / 50, alpha = alpha, mu = mu, sigma = sigma)
//' @export
// [[Rcpp::export]]
double logLikWouPairs(arma::mat x, arma::vec t, arma::vec alpha, arma::vec mu, arma::vec sigma, double rho = 0, int maxK = 2, double expTrc = 30) {

  /*
   * Create basic objects
   */

  // Number of pairs
  arma::uword N = x.n_rows;

  // Create log-likelihoods
  double loglikinitial = 0;
  double logliktpd = 0;
  double loglik = 0;

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
  arma::vec twokpi = arma::linspace<arma::vec>(-maxK * 2 * M_PI, maxK * 2 * M_PI, lk);

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

  // Add a penalty to the loglikelihood in case any assumption is violated
  double penalty = 0;

  // Only positive definiteness can be violated with the parametrization of A
  double prodDiagonal = add * add + alpha(0) * alpha(1);
  double testPosDef = prodDiagonal - alpha(2) * alpha(2);

  // Check positive definiteness
  if (testPosDef <= 0) {

    // Add a penalty
    penalty = -testPosDef * 10000 + 10;

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
  double l2pi = log(2 * M_PI);

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
  arma::mat weightswindsinitial(N, lk * lk);
  weightswindsinitial.fill(lognormconstSigmaA);

  // Loop in the data
  for (arma::uword i = 0; i < N; i++) {

    // Compute the factors in the exponent that do not depend on the windings
    arma::vec xmu = x.submat(i, 0, i, 1).t() - mu;
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

  // The unstandardized weights of the tpd give the required wrappings for the initial loglikelihood
  loglikinitial = accu(log(sum(exp(weightswindsinitial), 1)));

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

  // Loop in the data
  for (arma::uword i = 0; i < N; i++) {

    // Initial point x0 varying with i
    arma::vec x00 = x.submat(i, 0, i, 1).t();;

    // Evaluation point x varying with i
    arma::vec xx = x.submat(i, 2, i, 3).t();;

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

    // Common muti
    arma::vec muti = mu + ExptiA * (x00 - mu);

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

  // Logarithm of tpd
  tpdfinal = log(tpdfinal);

  // Set log(0) to -trunc, as this is the truncation of the negative exponentials
  tpdfinal.elem(find_nonfinite(tpdfinal)).fill(-expTrc);

  // Log-likelihood from tpd
  logliktpd = sum(tpdfinal);

  // Final loglikelihood
  loglik = loglikinitial + logliktpd;

  // Check if it is finite
  if (!std::isfinite(loglik)) loglik = -100 * N;

  return loglik - penalty;

}
