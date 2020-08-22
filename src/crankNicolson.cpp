# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
# include <Rmath.h>
# include <R.h>
using namespace Rcpp;

// Declarations for solvePeriodicTridiag and forwardSweepPeriodicTridiag
arma::vec solvePeriodicTridiag(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int LU);
arma::vec forwardSweepPeriodicTridiag(arma::vec a, arma::vec b, arma::vec c);

//' @title Crank--Nicolson finite difference scheme for the 1D Fokker--Planck equation with periodic boundaries
//'
//' @description Implementation of the Crank--Nicolson scheme for solving the Fokker--Planck equation
//' \deqn{p(x, t)_t = -(p(x, t) b(x))_x + \frac{1}{2}(\sigma^2(x) p(x, t))_{xx},}{p(x, t)_t = -(p(x, t) * b(x))_x + 1/2 * (\sigma^2(x) p(x, t))_{xx},}
//' where \eqn{p(x, t)} is the transition probability density of the circular diffusion
//' \deqn{dX_t=b(X_t)dt+\sigma(X_t)dW_t}{dX_t=b(X_t)dt+\sigma(X_t)dW_t}.
//'
//' @param u0 matrix of size \code{c(Mx, 1)} giving the initial condition. Typically, the evaluation of a density highly concentrated at a given point. If \code{nt == 1}, then \code{u0} can be a matrix \code{c(Mx, nu0)} containing different starting values in the columns.
//' @param b vector of length \code{Mx} containing the evaluation of the drift.
//' @param sigma2 vector of length \code{Mx} containing the evaluation of the squared diffusion coefficient.
//' @param N increasing integer vector of length \code{nt} giving the indexes of the times at which the solution is desired. The times of the solution are \code{delta * c(0:max(N))[N + 1]}.
//' @param deltat time step.
//' @param Mx size of the equispaced spatial grid in \eqn{[-\pi,\pi)}.
//' @param deltax space grid discretization.
//' @param imposePositive flag to indicate whether the solution should be transformed in order to be always larger than a given tolerance. This prevents spurious negative values. The tolerance will be taken as \code{imposePositiveTol} if this is different from \code{FALSE} or \code{0}.
//' @return
//' \itemize{
//' \item If \code{nt > 1}, a matrix of size \code{c(Mx, nt)} containing the discretized solution at the required times.
//' \item If \code{nt == 1}, a matrix of size \code{c(Mx, nu0)} containing the discretized solution at a fixed time for different starting values.
//' }
//' @details The function makes use of \code{\link{solvePeriodicTridiag}} for obtaining implicitly the next step in time of the solution.
//'
//' If \code{imposePositive = TRUE}, the code implicitly assumes that the solution integrates to one at any step. This might b unrealistic if the initial condition is not properly represented in the grid (for example, highly concentrated density in a sparse grid).
//' @references
//' Thomas, J. W. (1995). \emph{Numerical Partial Differential Equations: Finite Difference Methods}. Springer, New York. \url{https://doi.org/10.1007/978-1-4899-7278-1}
//' @examples
//' # Parameters
//' Mx <- 200
//' N <- 200
//' x <- seq(-pi, pi, l = Mx + 1)[-c(Mx + 1)]
//' times <- seq(0, 1, l = N + 1)
//' u0 <- dWn1D(x, pi/2, 0.05)
//' b <- driftWn1D(x, alpha = 1, mu = pi, sigma = 1)
//' sigma2 <- rep(1, Mx)
//'
//' # Full trajectory of the solution (including initial condition)
//' u <- crankNicolson1D(u0 = cbind(u0), b = b, sigma2 = sigma2, N = 0:N,
//'                      deltat = 1 / N, Mx = Mx, deltax = 2 * pi / Mx)
//'
//' # Mass conservation
//' colMeans(u) * 2 * pi
//'
//' # Visualization of tpd
//' plotSurface2D(times, x, z = t(u), levels = seq(0, 3, l = 50))
//'
//' # Only final time
//' v <- crankNicolson1D(u0 = cbind(u0), b = b, sigma2 = sigma2, N = N,
//'                      deltat = 1 / N, Mx = Mx, deltax = 2 * pi / Mx)
//' sum(abs(u[, N + 1] - v))
//' @export
// [[Rcpp::export]]
arma::mat crankNicolson1D(arma::mat u0, arma::vec b, arma::vec sigma2, arma::uvec N, double deltat, arma::uword Mx, double deltax, int imposePositive = 0) {

  // Take elements or rows depending on N
  arma::uword nt = N.n_elem;
  arma::uword nu0;
  if (nt > 1) {

    nu0 = u0.n_elem;

  } else {

    nu0 = u0.n_rows;

  }

  // Check if N is sorted
  if (!N.is_sorted()) {

    stop("N is not sorted increasingly");

  }

  // Check dimensions
  if ((Mx != nu0) | (Mx != b.n_elem) | (Mx != sigma2.n_elem)) {

    stop("Incompatible lengths of u0, b and sigma2");

  }

  // Common factors
  double r = deltat / (4 * deltax * deltax);
  b = b * deltax;

  // Coefficients
  arma::vec gamma = shift(-b + sigma2, -1) * r;
  arma::vec beta = 2 * sigma2 * r - 1;
  arma::vec alpha = shift(b + sigma2, +1) * r;

  // Save LU decomposition as two vectors
  arma::mat LU = forwardSweepPeriodicTridiag(alpha, -beta - 2, gamma);
  arma::vec betaLU = LU.col(0);
  arma::vec gammaLU = LU.col(1);

  // Save the solution vector at a selected time points or just at the final time (with possible different starting values)?
  if (nt > 1) {

    // Solution
    arma::mat u = arma::zeros(Mx, nt);

    // Counter for saving solution at N(m)
    arma::uword m = 0;

    // Loop in time
    for (arma::uword n = 0; n < ((N(nt - 1)) + 1); n++) {

      // Save solution?
      if (n == N(m)) {

        u.col(m) = u0;
        m++;

      }

      // Intercept
      arma::vec d = -gamma % shift(u0, -1) + beta % u0 - alpha % shift(u0, +1);

      // Solve tridiagonal system with presaved LU decomposition
      u0 = solvePeriodicTridiag(alpha, betaLU, gammaLU, d, 1);

    }

    // Ensure positiveness of the solution
    if (imposePositive != 0) {

      // Compute minimums and check which are negative
      arma::rowvec minU = min(u, 0);
      arma::uvec ind = find(minU < 0);

      // Standarize solution: (f - m) / (If + 2 * pi * m)
      for (arma::uword i = 0; i < ind.n_elem; i++) {

        // Make column positive and reweight
        arma::uword j = ind(i);
        double m = minU(j);
        u.col(j) -= m;
        u.col(j) /= (1 - m * 2 * PI);

      }

    }

    return u;

  } else {

    // Loop in different starting values
    for (arma::uword j = 0; j < u0.n_cols; j++) {

      // Loop in time
      for (arma::uword n = 0; n < N(0); n++) {

        // Previous column
        arma::vec u0j = u0.col(j);

        // Intercept
        arma::vec d = -gamma % shift(u0j, -1) + beta % u0j - alpha % shift(u0j, +1);

        // Solve tridiagonal system with presaved LU decomposition
        u0.col(j) = solvePeriodicTridiag(alpha, betaLU, gammaLU, d, 1);

      }

    }

    // Ensure positiveness of the solution
    if (imposePositive != 0) {

      // Compute minimums and check which are negative
      arma::rowvec minU = min(u0, 0);
      arma::uvec ind = find(minU < 0);

      // Standarize solution: (f - m) / (If + 2 * pi * m)
      for (arma::uword i = 0; i < ind.n_elem; i++) {

        // Make column positive and reweight
        arma::uword j = ind(i);
        double m = minU(j);
        u0.col(j) -= m;
        u0.col(j) /= (1 - m * 2 * PI);

      }

    }

    return u0;

  }

}


//' @title Crank--Nicolson finite difference scheme for the 2D Fokker--Planck equation with periodic boundaries
//'
//' @description Implementation of the Crank--Nicolson scheme for solving the Fokker--Planck equation
//' \deqn{p(x, y, t)_t = -(p(x, y, t) b_1(x, y))_x -(p(x, y, t) b_2(x, y))_y+}
//' \deqn{+ \frac{1}{2}(\sigma_1^2(x, y) p(x, y, t))_{xx} + \frac{1}{2}(\sigma_2^2(x, y) p(x, y, t))_{yy} + (\sigma_{12}(x, y) p(x, y, t))_{xy},}{p(x, y, t)_t = -(p(x, y, t) * b_1(x, y))_x -(p(x, y, t) * b_2(x, y))_y + 1/2 * (\sigma_1^2(x, y) *p(x, y, t))_{xx} + 1/2 * (\sigma_2^2(x, y) p(x, y, t))_{yy} + (\sigma_{12}(x, y) p(x, y, t))_{xy},}
//' where \eqn{p(x, y, t)} is the transition probability density of the toroidal diffusion
//' \deqn{dX_t=b_1(X_t,Y_t)dt+\sigma_1(X_t,Y_t)dW^1_t+\sigma_{12}(X_t,Y_t)dW^2_t,}{dX_t=b_1(X_t,Y_t)dt+\sigma_1(X_t,Y_t)dW^1_t+\sigma_{12}(X_t,Y_t)dW^2_t,}
//' \deqn{dY_t=b_2(X_t,Y_t)dt+\sigma_{12}(X_t,Y_t)dW^1_t+\sigma_2(X_t,Y_t)dW^2_t.}{dY_t=b_2(X_t,Y_t)dt+\sigma_{12}(X_t,Y_t)dW^1_t+\sigma_2(X_t,Y_t)dW^2_t.}
//'
//' @param u0 matrix of size \code{c(Mx * My, 1)} giving the initial condition matrix column-wise stored. Typically, the evaluation of a density highly concentrated at a given point. If \code{nt == 1}, then \code{u0} can be a matrix \code{c(Mx * My, nu0)} containing different starting values in the columns.
//' @param bx,by matrices of size \code{c(Mx, My)} containing the evaluation of the drift in the first and second space coordinates, respectively.
//' @param sigma2x,sigma2y,sigmaxy matrices of size \code{c(Mx, My)} containing the evaluation of the entries of the diffusion matrix (it has to be positive definite)\cr
//' \code{rbind(c(sigma2x, sigmaxy),
//'             c(sigmaxy, sigma2y))}.
//' @inheritParams crankNicolson1D
//' @param Mx,My sizes of the equispaced spatial grids in \eqn{[-\pi,\pi)} for each component.
//' @param deltax,deltay space grid discretizations for each component.
//' @param imposePositive flag to indicate whether the solution should be transformed in order to be always larger than a given tolerance. This prevents spurious negative values. The tolerance will be taken as \code{imposePositiveTol} if this is different from \code{FALSE} or \code{0}.
//' @return
//' \itemize{
//' \item If \code{nt > 1}, a matrix of size \code{c(Mx * My, nt)} containing the discretized solution at the required times with the \code{c(Mx, My)} matrix stored column-wise.
//' \item If \code{nt == 1}, a matrix of size \code{c(Mx * My, nu0)} containing the discretized solution at a fixed time for different starting values.
//' }
//' @details The function makes use of \code{\link{solvePeriodicTridiag}} for obtaining implicitly the next step in time of the solution.
//'
//' If \code{imposePositive = TRUE}, the code implicitly assumes that the solution integrates to one at any step. This might b unrealistic if the initial condition is not properly represented in the grid (for example, highly concentrated density in a sparse grid).
//' @references
//' Thomas, J. W. (1995). \emph{Numerical Partial Differential Equations: Finite Difference Methods}. Springer, New York. \url{https://doi.org/10.1007/978-1-4899-7278-1}
//' @examples
//' # Parameters
//' Mx <- 100
//' My <- 100
//' N <- 200
//' x <- seq(-pi, pi, l = Mx + 1)[-c(Mx + 1)]
//' y <- seq(-pi, pi, l = My + 1)[-c(My + 1)]
//' m <- c(pi / 2, pi)
//' p <- c(0, 1)
//' u0 <- c(outer(dWn1D(x, p[1], 0.5), dWn1D(y, p[2], 0.5)))
//' bx <- outer(x, y, function(x, y) 5 * sin(m[1] - x))
//' by <- outer(x, y, function(x, y) 5 * sin(m[2] - y))
//' sigma2 <- matrix(1, nrow = Mx, ncol = My)
//' sigmaxy <- matrix(0.5, nrow = Mx, ncol = My)
//'
//' # Full trajectory of the solution (including initial condition)
//' u <- crankNicolson2D(u0 = cbind(u0), bx = bx, by = by, sigma2x = sigma2,
//'                      sigma2y = sigma2, sigmaxy = sigmaxy,
//'                      N = 0:N, deltat = 1 / N, Mx = Mx, deltax = 2 * pi / Mx,
//'                      My = My, deltay = 2 * pi / My)
//'
//' # Mass conservation
//' colMeans(u) * 4 * pi^2
//'
//' # Only final time
//' v <- crankNicolson2D(u0 = cbind(u0), bx = bx, by = by, sigma2x = sigma2,
//'                      sigma2y = sigma2, sigmaxy = sigmaxy,
//'                      N = N, deltat = 1 / N, Mx = Mx, deltax = 2 * pi / Mx,
//'                      My = My, deltay = 2 * pi / My)
//' sum(abs(u[, N + 1] - v))
//'
//' \dontrun{
//' # Visualization of tpd
//' library(manipulate)
//' manipulate({
//'   plotSurface2D(x, y, z = matrix(u[, j + 1], Mx, My),
//'                 main = round(mean(u[, j + 1]) * 4 * pi^2, 4),
//'                 levels = seq(0, 2, l = 21))
//'   points(p[1], p[2], pch = 16)
//'   points(m[1], m[2], pch = 16)
//' }, j = slider(0, N))
//' }
//' @export
// [[Rcpp::export]]
arma::mat crankNicolson2D(arma::mat u0, arma::mat bx, arma::mat by, arma::mat sigma2x, arma::mat sigma2y, arma::mat sigmaxy, arma::uvec N, double deltat, arma::uword Mx, double deltax, arma::uword My, double deltay, int imposePositive = 0) {

  // Take elements or rows depending on N
  arma::uword nt = N.n_elem;
  arma::uword nu0;
  if (nt > 1) {

    nu0 = u0.n_elem;

  } else {

    nu0 = u0.n_rows;

  }

  // Check if N is sorted
  if (!N.is_sorted()) {

    stop("N is not sorted increasingly");

  }

  // Mx * My
  arma::uword M2 = Mx * My;

  // Check dimensions
  if ((Mx != bx.n_rows) | (Mx != by.n_rows) | (Mx != sigma2x.n_rows) | (Mx != sigma2y.n_rows) | (Mx != sigmaxy.n_rows) |
      (My != bx.n_cols) | (My != by.n_cols) | (My != sigma2x.n_cols) | (My != sigma2y.n_cols) | (My != sigmaxy.n_cols) |
      (nu0 != M2)) {

    stop("Incompatible sizes of u0, bx, by, sigma2x, sigma2y, sigmaxy");

  }

  // Zero shift
  arma::uvec idxy = arma::linspace<arma::uvec>(0, M2 - 1, M2);

  // Column-ordering to row-ordering (and viceversa)
  arma::uvec col2row = floor(idxy / My);
  col2row = (idxy - col2row * My) * Mx + col2row;
  arma::uvec row2col = floor(idxy / Mx);
  row2col = (idxy - row2col * Mx) * My + row2col;

  // Order diagonals
  arma::umat idmat = reshape(idxy, Mx, My);
  arma::uvec p1p1 = vectorise(shift(shift(idmat, -1, 1), -1, 0), 0);
  arma::uvec p1m1 = vectorise(shift(shift(idmat, -1, 1), +1, 0), 0);
  arma::uvec m1p1 = vectorise(shift(shift(idmat, +1, 1), -1, 0), 0);
  arma::uvec m1m1 = vectorise(shift(shift(idmat, +1, 1), +1, 0), 0);
  // arma::uvec p1p1 = shift(shift(idxy, -Mx), -1);
  // arma::uvec p1m1 = shift(shift(idxy, +Mx), -1);
  // arma::uvec m1p1 = shift(shift(idxy, -Mx), +1);
  // arma::uvec m1m1 = shift(shift(idxy, +Mx), +1);

  // Common factors
  double rx = deltat / (4 * deltax * deltax);
  double ry = deltat / (4 * deltay * deltay);

  // Coefficients x
  bx = bx * deltax;
  arma::vec gammax = vectorise(shift(-bx + sigma2x, -1, 0) * rx, 0);
  arma::vec betax = vectorise(2 * sigma2x * rx, 0);
  arma::vec alphax = vectorise(shift(bx + sigma2x, 1, 0) * rx, 0);

  // Coefficients y
  by = by * deltay;
  arma::vec gammay = vectorise(shift(-by + sigma2y, -1, 1) * ry, 1).t();
  arma::vec betay = vectorise(2 * sigma2y * ry, 1).t();
  arma::vec alphay = vectorise(shift(by + sigma2y, 1, 1) * ry, 1).t();

  // Coefficients xy
  sigmaxy = sigmaxy * deltat / (8 * deltax * deltay);
  arma::vec Cp1p1 = vectorise(shift(shift(sigmaxy, -1, 1), -1, 0), 0);
  arma::vec Cp1m1 = -vectorise(shift(shift(sigmaxy, -1, 1), +1, 0), 0);
  arma::vec Cm1p1 = -vectorise(shift(shift(sigmaxy, +1, 1), -1, 0), 0);
  arma::vec Cm1m1 = vectorise(shift(shift(sigmaxy, +1, 1), +1, 0), 0);
  // arma::vec Cp1p1 = shift(shift(vectorise(sigmaxy, 0), -Mx), -1);
  // arma::vec Cp1m1 = -shift(shift(vectorise(sigmaxy, 0), +Mx), -1);
  // arma::vec Cm1p1 = -shift(shift(vectorise(sigmaxy, 0), -Mx), +1);
  // arma::vec Cm1m1 = shift(shift(vectorise(sigmaxy, 0), +Mx), +1);

  // Save LU decomposition as two vectors for x
  arma::mat LU = forwardSweepPeriodicTridiag(-alphax, betax + 1, -gammax);
  arma::vec betaxLU = LU.col(0);
  arma::vec gammaxLU = LU.col(1);

  // Save LU decomposition as two vectors for y
  LU = forwardSweepPeriodicTridiag(-alphay, betay + 1, -gammay);
  arma::vec betayLU = LU.col(0);
  arma::vec gammayLU = LU.col(1);

  // Save the solution vector at a selected time points or just at the final time (with possible different starting values)?
  if (nt > 1) {

    // Solution
    arma::mat u = arma::zeros(M2, nt);

    // Counter for saving solution at N(m)
    arma::uword m = 0;

    // Loop in time
    for (arma::uword n = 0; n < ((N(nt - 1)) + 1); n++) {

      // Save solution?
      if (n == N(m)) {

        u.col(m) = u0;
        m++;

      }

      // Y = u0 by column order
      arma::vec Y = u0;

      // Y = u0 by row order
      arma::vec Yr = Y(col2row);

      // Intercept x
      arma::vec dx = gammax % shift(Y, -1) - betax % Y + alphax % shift(Y, +1);

      // Intercept y
      arma::vec dy = gammay % shift(Yr, -1) - betay % Yr + alphay % shift(Yr, +1);

      // Intercept xy
      arma::vec dxy = Y(p1p1) % Cp1p1 + Y(p1m1) % Cp1m1 + Y(m1p1) % Cm1p1 + Y(m1m1) % Cm1m1;

      /*
       * Y0: Explicit step
       */

      // Y0 = U(n) + deltat * F(U(n))
      Y += 2 * (dx + dy(row2col) + dxy);

      /*
       * Y1: Implicit step for x, tridiagonal *column-wise*
       */

      // Y1 = Y0 + 0.5 * deltat * (F(Y0) - F1(U(n)))
      Y = solvePeriodicTridiag(-alphax, betaxLU, gammaxLU, Y - dx, 1);

      /*
       * Y2: Implicit step for y, tridiagonal *row-wise*
       */

      // Y2 = Y1 + 0.5 * deltat * (F(Y1) - F2(U(n)))
      Y = solvePeriodicTridiag(-alphay, betayLU, gammayLU, Y(col2row) - dy, 1);

      /*
       * Store solution
       */

      u0 = Y(row2col);

    }

    // Ensure positiveness of the solution
    if (imposePositive != 0) {

      // Compute minimums and check which are negative
      arma::rowvec minU = min(u, 0);
      arma::uvec ind = find(minU < 0);

      // Standarize solution: (f - m) / (If + 4 * pi^2 * m)
      for (arma::uword i = 0; i < ind.n_elem; i++) {

        // Make column positive and reweight
        arma::uword j = ind(i);
        double m = minU(j);
        u.col(j) -= m;
        u.col(j) /= (1 - m * 4 * PI * PI);

      }

    }

    return u;

  } else {

    // Loop in different starting values
    for (arma::uword j = 0; j < u0.n_cols; j++) {

      // Loop in time
      for (arma::uword n = 0; n < N(0); n++) {

        // Previous column Y = u0 by column order
        arma::vec Y = u0.col(j);

        // Y = u0 by row order
        arma::vec Yr = Y(col2row);

        // Intercept x
        arma::vec dx = gammax % shift(Y, -1) - betax % Y + alphax % shift(Y, +1);

        // Intercept y
        arma::vec dy = gammay % shift(Yr, -1) - betay % Yr + alphay % shift(Yr, +1);

        // Intercept xy
        arma::vec dxy = Y(p1p1) % Cp1p1 + Y(p1m1) % Cp1m1 + Y(m1p1) % Cm1p1 + Y(m1m1) % Cm1m1;

        /*
         * Y0: Explicit step
         */

        // Y0 = U(n) + deltat * F(U(n))
        Y += 2 * (dx + dy(row2col) + dxy);

        /*
         * Y1: Implicit step for x, tridiagonal *column-wise*
         */

        // Y1 = Y0 + 0.5 * deltat * (F(Y0) - F1(U(n)))
        Y = solvePeriodicTridiag(-alphax, betaxLU, gammaxLU, Y - dx, 1);

        /*
         * Y2: Implicit step for y, tridiagonal *row-wise*
         */

        // Y2 = Y1 + 0.5 * deltat * (F(Y1) - F2(U(n)))
        Y = solvePeriodicTridiag(-alphay, betayLU, gammayLU, Y(col2row) - dy, 1);

        /*
         * Store solution
         */

        u0.col(j) = Y(row2col);

      }

    }

    // Ensure positiveness of the solution
    if (imposePositive != 0) {

      // Compute minimums and check which are negative
      arma::rowvec minU = min(u0, 0);
      arma::uvec ind = find(minU < 0);

      // Standarize solution: (f - m) / (If + 4 * pi^2 * m)
      for (arma::uword i = 0; i < ind.n_elem; i++) {

        // Make column positive and reweight
        arma::uword j = ind(i);
        double m = minU(j);
        u0.col(j) -= m;
        u0.col(j) /= (1 - m * 4 * PI * PI);

      }

    }

    return u0;

  }

}
