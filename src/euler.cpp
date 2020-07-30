#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

// Declaration for safeSoftMax and drifts
arma::mat safeSoftMax(arma::mat logs, double expTrc = 30);
arma::vec driftWn1D(arma::vec x, double alpha, double mu, double sigma, int maxK = 2, double expTrc = 30);
arma::mat driftWn2D(arma::mat x, arma::mat A, arma::vec mu, arma::vec sigma, double rho = 0, int maxK = 2, double expTrc = 30);

//' @title Simulation of trajectories of the WN or vM diffusion in 1D
//'
//' @description Simulation of the Wrapped Normal (WN) diffusion or von Mises (vM) diffusion by the Euler method in 1D, for several starting values.
//'
//' @param x0 vector of length \code{nx0} giving the initial points.
//' @inheritParams driftWn1D
//' @inheritParams dTpdWou1D
//' @param N number of discretization steps.
//' @param delta discretization step.
//' @param type integer giving the type of diffusion. Currently, only \code{1} for WN and \code{2} for vM are supported.
//' @return A matrix of size \code{c(nx0, N + 1)} containing the \code{nx0} discretized trajectories. The first column corresponds to the vector \code{x0}.
//' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
//' @examples
//' N <- 100
//' nx0 <- 20
//' x0 <- seq(-pi, pi, l = nx0 + 1)[-(nx0 + 1)]
//' set.seed(12345678)
//' samp <- euler1D(x0 = x0, mu = 0, alpha = 3, sigma = 1, N = N,
//'                 delta = 0.01, type = 2)
//' tt <- seq(0, 1, l = N + 1)
//' plot(rep(0, nx0), x0, pch = 16, col = rainbow(nx0), xlim = c(0, 1))
//' for (i in 1:nx0) linesCirc(tt, samp[i, ], col = rainbow(nx0)[i])
//' @export
// [[Rcpp::export]]
arma::mat euler1D(arma::vec x0, double alpha, double mu, double sigma, arma::uword N = 100, double delta = 0.01, int type = 1, int maxK = 2, double expTrc = 30) {

  // Number of x0's
  arma::uword nx0 = x0.n_elem;

  // Create result
  arma::mat x(nx0, N + 1);

  // Initialize by rows
  x.col(0) = x0;

  // Create drift
  arma::vec drift(nx0);

  // Create diffusion
  double diffusion = sigma;

  // Noise
  arma::vec Z = rnorm(N * nx0) * sqrt(delta) * diffusion;

  // Loop for simulation
  for (arma::uword i = 1; i < N + 1; i++) {

    // Drift
    switch(type) {

      // WN
      case 1 :
      drift = driftWn1D(x.col(i - 1), alpha, mu, sigma, maxK, expTrc);
      break;

      // vM
      case 2 :
      drift = alpha * sin(mu - x.col(i - 1));
      break;

    }

    // Index
    arma::uword inx0 = nx0 * (i - 1);

    // Update
    x.col(i) = x.col(i - 1) + drift * delta + Z.subvec(inx0, inx0 + nx0 - 1);

    // Convert result to [-PI,PI)
    x.col(i) -= floor((x.col(i) + PI) / (2 * PI)) * (2 * PI);

  }

  return x;

}


//' @title Simulation of trajectories of the WN or MvM diffusion in 2D
//'
//' @description Simulation of the Wrapped Normal (WN) diffusion or Multivariate von Mises (MvM) diffusion by the Euler method in 2D, for several starting values.
//'
//' @param x0 matrix of size \code{c(nx0, 2)} giving the initial points.
//' @inheritParams driftWn2D
//' @inheritParams euler1D
//' @inheritParams safeSoftMax
//' @return An array of size \code{c(nx0, 2, N + 1)} containing the \code{nx0} discretized trajectories. The first slice corresponds to the matrix \code{x0}.
//' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
//' @examples
//' N <- 100
//' nx0 <- 5
//' x0 <- seq(-pi, pi, l = nx0 + 1)[-(nx0 + 1)]
//' x0 <- as.matrix(expand.grid(x0, x0))
//' nx0 <- nx0^2
//' set.seed(12345678)
//' samp <- euler2D(x0 = x0, mu = c(0, 0), A = rbind(c(3, 1), 1:2),
//'                 sigma = c(1, 1), N = N, delta = 0.01, type = 2)
//' plot(x0[, 1], x0[, 2], xlim = c(-pi, pi), ylim = c(-pi, pi), pch = 16,
//'      col = rainbow(nx0))
//' for (i in 1:nx0) linesTorus(samp[i, 1, ], samp[i, 2, ],
//'                            col = rainbow(nx0, alpha = 0.5)[i])
//' @export
// [[Rcpp::export]]
arma::cube euler2D(arma::mat x0, arma::mat A, arma::vec mu, arma::vec sigma, double rho = 0, arma::uword N = 100, double delta = 0.01, int type = 1, int maxK = 2, double expTrc = 30) {

  // Number of x0's
  arma::uword nx0 = x0.n_rows;

  // Create and initialize result
  arma::cube x = arma::zeros(nx0, 2, N + 1);
  x.each_slice() = x0;

  // Create drifts
  arma::mat drift(nx0, 2);

  // Uncorrelated noise
  arma::mat Z = arma::zeros(N * nx0, 2);
  Z.col(0) = as<arma::vec>(rnorm(N * nx0));
  Z.col(1) = as<arma::vec>(rnorm(N * nx0));

  // Correlate noise
  Z.col(1) = (rho * Z.col(0) + sqrt(1 - rho * rho) * Z.col(1));
  Z.col(1) *= sqrt(delta) * sigma(1);
  Z.col(0) *= sqrt(delta) * sigma(0);

  /*
   * vM quantities
   */

  // mut, alphat
  arma::rowvec mut = mu.t();
  arma::rowvec alphat = diagvec(A).t();

  if (type == 2) {

    // Set diagonal to zero
    A.diag().zeros();

  }

  // Loop on the discretization
  for (arma::uword i = 1; i < N + 1; i++) {

    // Drift
    if (type == 1) { // WN

      drift = driftWn2D(x.slice(i - 1), A, mu, sigma, rho, maxK, expTrc);

    } else if (type == 2) { // vM

      arma::mat difmu = x.slice(i - 1);
      difmu.each_row() -= mut;
      arma::mat si = sin(difmu), asi = si;
      asi.each_row() %= alphat;
      arma::mat ci = cos(difmu);
      arma::mat aci = si * A;
      aci %= ci;
      drift = -asi + aci;

    }

    // Index
    arma::uword inx0 = nx0 * (i - 1);

    // Update
    x.slice(i) = x.slice(i - 1) + drift * delta + Z.rows(inx0, inx0 + nx0 - 1);

    // Convert result to [-PI,PI)
    x.slice(i) -= floor((x.slice(i) + PI) / (2 * PI)) * (2 * PI);

  }

  return x;

}


//' @title Multiple simulation of trajectory ends of the WN or vM diffusion in 1D
//'
//' @description Simulates \code{M} trajectories starting from different initial values \code{x0} of the WN or vM diffusion in 1D, by the Euler method, and returns their ends.
//'
//' @inheritParams euler1D
//' @inheritParams dTpdWou1D
//' @param M number of Monte Carlo replicates.
//' @return A matrix of size \code{c(nx0, M)} containing the \code{M} trajectory ends for each starting value \code{x0}.
//' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
//' @examples
//' N <- 100
//' nx0 <- 20
//' x0 <- seq(-pi, pi, l = nx0 + 1)[-(nx0 + 1)]
//' set.seed(12345678)
//' samp1 <- euler1D(x0 = x0, mu = 0, alpha = 3, sigma = 1, N = N,
//'                  delta = 0.01, type = 2)
//' tt <- seq(0, 1, l = N + 1)
//' plot(rep(0, nx0), x0, pch = 16, col = rainbow(nx0), xlim = c(0, 1))
//' for (i in 1:nx0) linesCirc(tt, samp1[i, ], col = rainbow(nx0)[i])
//' set.seed(12345678)
//' samp2 <- stepAheadWn1D(x0 = x0, mu = 0, alpha = 3, sigma = 1, M = 1,
//'                        N = N, delta = 0.01, type = 2)
//' points(rep(1, nx0), samp2[, 1], pch = 16, col = rainbow(nx0))
//' samp1[, N + 1]
//' samp2[, 1]
//' @export
// [[Rcpp::export]]
arma::mat stepAheadWn1D(arma::vec x0, double alpha, double mu, double sigma, arma::uword M, arma::uword N = 100, double delta = 0.01, int type = 1, int maxK = 2, double expTrc = 30) {

  // Number of x0's
  arma::uword nx0 = x0.n_elem;
  arma::uword Mnx0 = M * nx0;

  // Create and initialize result
  arma::vec x = repmat(x0, M, 1);

  // Create drift
  arma::vec drift(Mnx0);

  // Diffusion coefficient
  double diffusion = sigma;

  // Noise
  arma::vec Z = rnorm(N * Mnx0) * sqrt(delta) * diffusion;

  // Loop
  for (arma::uword i = 0; i < N; i++) {

    // Drift
    switch(type) {

      // WN
      case 1 :
      drift = driftWn1D(x, alpha, mu, sigma, maxK, expTrc);
      break;

      // vM
      case 2 :
      drift = alpha * sin(mu - x);
      break;

    }

    // Index
    arma::uword iM = Mnx0 * i;

    // Update
    x += drift * delta + Z.subvec(iM, iM + Mnx0 - 1);

    // Convert result to [-PI, PI)
    x -= floor((x + PI) / (2 * PI)) * (2 * PI);

  }

  return reshape(x, nx0, M);

}


//' @title Multiple simulation of trajectory ends of the WN or MvM diffusion in 2D
//'
//' @description Simulates \code{M} trajectories starting from different initial values \code{x0} of the WN or MvM diffusion in 2D, by the Euler method, and returns their ends.
//'
//' @inheritParams euler2D
//' @inheritParams dTpdWou2D
//' @inheritParams stepAheadWn1D
//' @return An array of size \code{c(nx0, 2, M)} containing the \code{M} trajectory ends for each starting value \code{x0}.
//' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
//' @examples
//' N <- 100
//' nx0 <- 3
//' x0 <- seq(-pi, pi, l = nx0 + 1)[-(nx0 + 1)]
//' x0 <- as.matrix(expand.grid(x0, x0))
//' nx0 <- nx0^2
//' set.seed(12345678)
//' samp1 <- euler2D(x0 = x0, mu = c(0, 0), A = rbind(c(3, 1), 1:2),
//'                  sigma = c(1, 1), N = N, delta = 0.01, type = 2)
//' plot(x0[, 1], x0[, 2], xlim = c(-pi, pi), ylim = c(-pi, pi), pch = 16,
//'      col = rainbow(nx0))
//' for (i in 1:nx0) linesTorus(samp1[i, 1, ], samp1[i, 2, ],
//'                            col = rainbow(nx0, alpha = 0.75)[i])
//' set.seed(12345678)
//' samp2 <- stepAheadWn2D(x0 = x0, mu = c(0, 0), A = rbind(c(3, 1), 1:2),
//'                        sigma = c(1, 1), M = 2, N = N, delta = 0.01,
//'                        type = 2)
//' points(samp2[, 1, 1], samp2[, 2, 1], pch = 16, col = rainbow(nx0))
//' samp1[, , N + 1]
//' samp2[, , 1]
//' @export
// [[Rcpp::export]]
arma::cube stepAheadWn2D(arma::mat x0, arma::vec mu, arma::mat A, arma::vec sigma, double rho = 0, arma::uword M = 100, arma::uword N = 100, double delta = 0.01, int type = 1, int maxK = 2, double expTrc = 30) {

  // Number of x0's
  arma::uword nx0 = x0.n_rows;
  arma::uword Mnx0 = M * nx0;

  // Create and initialize result
  arma::mat x = repmat(x0, M, 1);

  // Create drift
  arma::mat drift(Mnx0, 2);

  // Uncorrelated noise
  arma::mat Z = arma::zeros(N * Mnx0, 2);
  Z.col(0) = as<arma::vec>(rnorm(N * Mnx0));
  Z.col(1) = as<arma::vec>(rnorm(N * Mnx0));

  // Correlate noise
  Z.col(1) = (rho * Z.col(0) + sqrt(1 - rho * rho) * Z.col(1));
  Z.col(1) *= sqrt(delta) * sigma(1);
  Z.col(0) *= sqrt(delta) * sigma(0);

  /*
   * vM quantities
   */

  // mut, alphat
  arma::rowvec mut = mu.t();
  arma::rowvec alphat = diagvec(A).t();

  if (type == 2) {

    // Set diagonal to zero
    A.diag().zeros();

  }

  // Loop
  for (arma::uword i = 0; i < N; i++) {

    // Drift
    if (type == 1) { // WN

      drift = driftWn2D(x, A, mu, sigma, rho, maxK, expTrc);

    } else if (type == 2) { // vM

      arma::mat difmu = x;
      difmu.each_row() -= mut;
      arma::mat si = sin(difmu), asi = si;
      asi.each_row() %= alphat;
      arma::mat ci = cos(difmu);
      arma::mat aci = si * A;
      aci %= ci;
      drift = -asi + aci;

    }

    // Index
    arma::uword iM = Mnx0 * i;

    // Update
    x += drift * delta + Z.rows(iM, iM + Mnx0 - 1);

    // Convert result to [-PI,PI)
    x -= floor((x + PI) / (2 * PI)) * (2 * PI) ;

  }

  // Return a cube
  arma::cube cx = arma::zeros(nx0, 2, M);
  cx.tube(0, 0, nx0 - 1, 0) = reshape(x.col(0), nx0, M);
  cx.tube(0, 1, nx0 - 1, 1) = reshape(x.col(1), nx0, M);
  return cx;

}
