#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

//' @title Safe softmax function for computing weights
//'
//' @description Computes the weights \eqn{w_i = \frac{e^{p_i}}{\sum_{j=1}^k e^{p_j}}} from \eqn{p_i}, \eqn{i=1,\ldots,k}
//' in a safe way to avoid overflows and to truncate automatically to zero low values of \eqn{w_i}.
//'
//' @param logs matrix of logarithms where each row contains a set of \eqn{p_1,\ldots,p_k} to compute the weights from.
//' @param expTrc truncation for exponential: \code{exp(x)} with \code{x <= -expTrc} is set to zero. Defaults to \code{30}.
//' @return A matrix of the size as \code{logs} containing the weights for each row.
//' @details The \code{logs} argument must be always a matrix.
//' @examples
//' # A matrix
//' safeSoftMax(rbind(1:10, 20:11))
//' rbind(exp(1:10) / sum(exp(1:10)), exp(20:11) / sum(exp(20:11)))
//'
//' # A row-matrix
//' safeSoftMax(rbind(-100:100), expTrc = 30)
//' exp(-100:100) / sum(exp(-100:100))
//' @export
// [[Rcpp::export]]
arma::mat safeSoftMax(arma::mat logs, double expTrc = 30) {

  // Maximum of logs by rows to avoid overflows
  arma::vec m = max(logs, 1);

  // Recenter by columns
  logs.each_col() -= m;

  // Ratios by columns
  logs.each_col() -= log(sum(exp(logs), 1));

  // Truncate exponential by using a lambda function - requires C++ 11
  logs.transform([expTrc](double val) { return (val < -expTrc) ? double(0) : double(exp(val)); });

  return logs;

}


//' @title Thomas algorithm for solving tridiagonal matrix systems, with optional presaving of LU decomposition
//'
//' @description Implementation of the Thomas algorithm to solve efficiently the tridiagonal matrix system
//' \deqn{b_1 x_1 + c_1 x_2 + a_1x_n = d_1}{b[1] x[1] + c[1] x[2] + a[1]x[n] = d[1]}
//' \deqn{a_2 x_1 + b_2 x_2 + c_2x_3 = d_2}{a[2] x[1] + b[2] x[2] + c[2]x[3] = d[2]}
//' \deqn{\vdots \vdots \vdots}{...}
//' \deqn{a_{n-1} x_{n-2} + b_{n-1} x_{n-1} + c_{n-1}x_{n} = d_{n-1}}{a[n-1] x[n-2] + b[n-1] x[n-1] + c[n-1]x[n] = d[n-1]}
//' \deqn{c_n x_1 + a_{n} x_{n-1} + b_nx_n = d_n}{c[n] x[1] + a[n] x[n-1] + b[n]x[n] = d[n]}
//' with \eqn{a_1=c_n=0}{a[1]=c[n]=0} (usual tridiagonal matrix). If \eqn{a_1\neq0}{a[1]/=0} or \eqn{c_n\neq0}{c[n]/=0} (circulant tridiagonal matrix), then the Sherman--Morrison formula is employed.
//'
//' @param a,b,c subdiagonal (below main diagonal), diagonal and superdiagonal (above main diagonal), respectively. They all are vectors of length \code{n}.
//' @param d vector of constant terms, of length \code{n}. For \code{solveTridiagMatConsts}, it can be a matrix with \code{n} rows.
//' @param LU flag denoting if the forward sweep encoding the LU decomposition is supplied in vectors \code{b} and \code{c}. See details and examples.
//' @return
//' \itemize{
//' \item \code{solve*} functions: the solution, a vector of length \code{n} and a matrix with \code{n} rows for\cr \code{solveTridiagMatConsts}.
//' \item \code{forward*} functions: the matrix \code{cbind(b, c)} creating the suitable \code{b} and \code{c} arguments for calling \code{solve*} when \code{LU} is \code{TRUE}.
//' }
//' @details The Thomas algorithm is stable if the matrix is diagonally dominant.
//'
//' For the periodic case, two non-periodic tridiagonal systems with different constant terms (but same coefficients) are solved using \code{solveTridiagMatConsts}. These two solutions are combined by the Sherman--Morrison formula to obtain the solution to the periodic system.
//'
//' Note that the output of \code{solveTridiag} and \code{solveTridiagMatConsts} are independent from the values of \code{a[1]} and \code{c[n]}, but \code{solvePeriodicTridiag} is not.
//'
//' If \code{LU} is \code{TRUE}, then \code{b} and \code{c} must be precomputed with \code{forwardSweepTridiag} or\cr \code{forwardSweepPeriodicTridiag} for its use in the call of the appropriate solver, which will be slightly faster.
//' @references
//' Thomas, J. W. (1995). \emph{Numerical Partial Differential Equations: Finite Difference Methods}. Springer, New York. \doi{10.1007/978-1-4899-7278-1}
//'
//' Conte, S. D. and de Boor, C. (1980). \emph{Elementary Numerical Analysis: An Algorithmic Approach}. Third edition. McGraw-Hill, New York. \doi{10.1137/1.9781611975208}
//' @examples
//' # Tridiagonal matrix
//' n <- 10
//' a <- rnorm(n, 3, 1)
//' b <- rnorm(n, 10, 1)
//' c <- rnorm(n, 0, 1)
//' d <- rnorm(n, 0, 1)
//' A <- matrix(0, nrow = n, ncol = n)
//' diag(A) <- b
//' for (i in 1:(n - 1)) {
//'   A[i + 1, i] <- a[i + 1]
//'   A[i, i + 1] <- c[i]
//' }
//' A
//'
//' # Same solutions
//' drop(solveTridiag(a = a, b = b, c = c, d = d))
//' solve(a = A, b = d)
//'
//' # Presaving the forward sweep (encodes the LU factorization)
//' LU <- forwardSweepTridiag(a = a, b = b, c = c)
//' drop(solveTridiag(a = a, b = LU[, 1], c = LU[, 2], d = d, LU = 1))
//'
//' # With equal coefficient matrix
//' solveTridiagMatConsts(a = a, b = b, c = c, d = cbind(d, d + 1))
//' cbind(solve(a = A, b = d), solve(a = A, b = d + 1))
//' LU <- forwardSweepTridiag(a = a, b = b, c = c)
//' solveTridiagMatConsts(a = a, b = LU[, 1], c = LU[, 2], d = cbind(d, d + 1), LU = 1)
//'
//' # Periodic matrix
//' A[1, n] <- a[1]
//' A[n, 1] <- c[n]
//' A
//'
//' # Same solutions
//' drop(solvePeriodicTridiag(a = a, b = b, c = c, d = d))
//' solve(a = A, b = d)
//'
//' # Presaving the forward sweep (encodes the LU factorization)
//' LU <- forwardSweepPeriodicTridiag(a = a, b = b, c = c)
//' drop(solvePeriodicTridiag(a = a, b = LU[, 1], c = LU[, 2], d = d, LU = 1))
//' @export
// [[Rcpp::export]]
arma::vec solveTridiag(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int LU = 0) {

  // Check sizes of vectors
  arma::uword n = b.n_elem;
  if ((n != a.n_elem) | (n != c.n_elem) | (n != d.n_elem)) {

    stop("Incompatible lengths of a, b, c and d");

  }

  /*
   * Forward sweep
   */

  // n - 1
  arma::uword n1 = n - 1;

  if (LU == 0) {

    // First elements
    c(0) /= b(0);
    d(0) /= b(0);

    // Loop skipping first elements of c and d and last element of d
    for (arma::uword i = 1; i < n1; i++) {

      // Common factor
      double m = 1.0 / (b(i) - a(i) * c(i - 1));

      // Update c and d
      c(i) *= m;
      d(i) = (d(i) - a(i) * d(i - 1)) * m;

    }

    // Final update for d
    d(n1) = (d(n1) - a(n1) * d(n1 - 1)) / (b(n1) - a(n1) * c(n1 - 1));

  } else {

    // First element
    d(0) *= b(0);

    // Loop skipping first element of d
    for (arma::uword i = 1; i < n; i++) {

      // Update only d
      d(i) = (d(i) - a(i) * d(i - 1)) * b(i);

    }

  }

  /*
   * Back substitution
   */

  for (arma::uword i = n1; i-- > 0;) {

    d(i) -= c(i) * d(i + 1);

  }

  return d;

}


//' @rdname solveTridiag
//' @export
// [[Rcpp::export]]
arma::mat solveTridiagMatConsts(arma::vec a, arma::vec b, arma::vec c, arma::mat d, int LU = 0) {

  // Check sizes of vectors
  arma::uword n = b.n_elem;
  if ((n != a.n_elem) | (n != c.n_elem) | (n != d.n_rows)) {

    stop("Incompatible sizes of a, b, c and d");

  }

  /*
   * Forward sweep
   */

  // n - 1
  arma::uword n1 = n - 1;

  if (LU == 0) {

    // First elements
    c(0) /= b(0);
    d.row(0) /= b(0);

    // Loop skipping first elements of c and d and last element of d
    for (arma::uword i = 1; i < n1; i++) {

      // Common factor
      double m = 1.0 / (b(i) - a(i) * c(i - 1));

      // Update c and d
      c(i) *= m;
      d.row(i) = (d.row(i) - a(i) * d.row(i - 1)) * m;

    }

    // Final update for d
    d.row(n1) = (d.row(n1) - a(n1) * d.row(n1 - 1)) / (b(n1) - a(n1) * c(n1 - 1));

  } else {

    // First element
    d.row(0) *= b(0);

    // Loop skipping first element of d
    for (arma::uword i = 1; i < n; i++) {

      // Update only d
      d.row(i) = (d.row(i) - a(i) * d.row(i - 1)) * b(i);

    }

  }

  /*
   * Back substitution
   */

  for (arma::uword i = n1; i-- > 0;) {

    d.row(i) -= c(i) * d.row(i + 1);

  }

  return d;

}


//' @rdname solveTridiag
//' @export
// [[Rcpp::export]]
arma::vec solvePeriodicTridiag(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int LU = 0) {

  // Check sizes of vectors
  arma::uword n = b.n_elem;
  if ((n != a.n_elem) | (n != c.n_elem) | (n != d.n_elem)) {

    stop("Incompatible lengths of a, b, c and d");

  }

  // Create w
  arma::vec w = arma::zeros(n);

  // n - 1
  arma::uword n1 = n - 1;

  // beta
  double beta = -a(0);

  if (LU == 0) {

    // Fill w
    w(0) = b(0);
    w(n1) = -c(n1);

    // Prepare coefficients
    beta /= w(0);
    b(n1) -= beta * c(n1);
    b(0) *= 2;

  } else {

    // Fill w
    w(0) = 0.5 / b(0);
    w(n1) = -c(n1) / b(n1);

    // Prepare coefficients
    beta /= w(0);

  }

  // Solve two tridiagonal systems
  arma::mat y = solveTridiagMatConsts(a, b, c, join_horiz(d, w), LU);

  // Sherman--Morrison formula
  beta = (y(0, 0) + beta * y(n1, 0)) / (1.0 - y(0, 1) - beta * y(n1, 1));
  d = y.col(0) + beta * y.col(1);

  return d;

}


//' @rdname solveTridiag
//' @export
// [[Rcpp::export]]
arma::mat forwardSweepTridiag(arma::vec a, arma::vec b, arma::vec c) {

  // Check sizes of vectors
  arma::uword n = b.n_elem;
  if ((n != a.n_elem) | (n != c.n_elem)) {

    stop("Incompatible lengths of a, b and c");

  }

  // First element
  b(0) = 1 / b(0);
  c(0) *= b(0);

  // Loop skipping first elements of c and d and last element of d
  for (arma::uword i = 1; i < n; i++) {

    // Common factor
    b(i) = 1.0 / (b(i) - a(i) * c(i - 1));

    // Update c and d
    c(i) *= b(i);

  }

  return join_horiz(b, c);

}


//' @rdname solveTridiag
//' @export
// [[Rcpp::export]]
arma::mat forwardSweepPeriodicTridiag(arma::vec a, arma::vec b, arma::vec c) {

  // Check sizes of vectors
  arma::uword n = b.n_elem;
  if ((n != a.n_elem) | (n != c.n_elem)) {

    stop("Incompatible lengths of a, b and c");

  }

  // Common
  n--;

  // Prepare coefficients
  b(n) += a(0) / b(0) * c(n);
  b(0) *= 2;

  return forwardSweepTridiag(a, b, c);

}
