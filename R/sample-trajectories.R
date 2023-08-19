

#' @title Simulation of trajectories for the WN diffusion in 1D
#'
#' @description Simulation of the Wrapped Normal (WN) diffusion in 1D by subsampling a fine trajectory obtained by the Euler discretization.
#'
#' @param x0 initial point.
#' @inheritParams dTpdWou1D
#' @param N number of discretization steps in the resulting trajectory.
#' @param delta discretization step.
#' @param NFine number of discretization steps for the fine trajectory. Must be larger than \code{N}.
#' @param deltaFine discretization step for the fine trajectory. Must be smaller than \code{delta}.
#' @return A vector of length \code{N + 1} containing \code{x0} in the first entry and the discretized trajectory.
#' @details The fine trajectory is subsampled using the indexes \code{seq(1, NFine + 1, by = NFine / N)}.
#' @examples
#' isRStudio <- identical(.Platform$GUI, "RStudio")
#' if (isRStudio) {
#'   manipulate::manipulate({
#'     x <- seq(0, N * delta, by = delta)
#'     plot(x, x, ylim = c(-pi, pi), type = "n", ylab = expression(X[t]), xlab = "t")
#'     linesCirc(x, rTrajWn1D(x0 = 0, alpha = alpha, mu = 0, sigma = sigma, N = N,
#'                              delta = 0.01))
#'     }, delta = slider(0.01, 5.01, step = 0.1),
#'     N = manipulate::slider(10, 500, step = 10, initial = 200),
#'     alpha = manipulate::slider(0.01, 5, step = 0.1, initial = 1),
#'     sigma = manipulate::slider(0.01, 5, step = 0.1, initial = 1))
#' }
#' @export
rTrajWn1D <- function(x0, alpha, mu, sigma, N = 100, delta = 0.01,
                      NFine = ceiling(N * delta / deltaFine),
                      deltaFine = min(delta / 100, 1e-3)) {

  # Sample by Euler
  samp <- drop(euler1D(x0 = x0, alpha = alpha, mu = mu, sigma = sigma,
                       N = NFine, delta = deltaFine))

  # Subsample
  samp <- samp[seq(1, NFine + 1, by = NFine / N)]

  return(samp)

}


#' @title Simulation of trajectories for the WN diffusion in 2D
#'
#' @description Simulation of the Wrapped Normal (WN) diffusion in 2D by subsampling a fine trajectory obtained by the Euler discretization.
#'
#' @param x0 vector of length \code{2} giving the initial point.
#' @inheritParams dTpdWou2D
#' @inheritParams rTrajWn1D
#' @return A matrix of size \code{c(N + 1, 2)} containing \code{x0} in the first entry and the discretized trajectory.
#' @details The fine trajectory is subsampled using the indexes \code{seq(1, NFine + 1, by = NFine / N)}.
#' @examples
#' samp <- rTrajWn2D(x0 = c(0, 0), alpha = c(1, 1, -0.5), mu = c(pi, pi),
#'                     sigma = c(1, 1), N = 1000, delta = 0.01)
#' plot(samp, xlim = c(-pi, pi), ylim = c(-pi, pi), pch = 19, cex = 0.25,
#'      xlab = expression(X[t]), ylab = expression(Y[t]), col = rainbow(1000))
#' linesTorus(samp[, 1], samp[, 2], col = rainbow(1000))
#' @export
rTrajWn2D <- function(x0, alpha, mu, sigma, rho = 0, N = 100, delta = 0.01,
                      NFine = ceiling(N * delta / deltaFine),
                      deltaFine = min(delta / 100, 1e-3)) {

  # Sample by Euler
  samp <- t(euler2D(x0 = rbind(x0), A = alphaToA(alpha = alpha, sigma = sigma,
                                                 rho = rho),
                    mu = mu, sigma = sigma, rho = rho, N = NFine,
                    delta = deltaFine)[1, , ])

  # Subsample
  samp <- samp[seq(1, NFine + 1, by = NFine / N), ]

  return(samp)

}


#' @title Simulation of trajectories of a Langevin diffusion
#'
#' @description Simulation of an arbitrary Langevin diffusion in dimension \code{p} by subsampling a fine trajectory obtained by the Euler discretization.
#'
#' @param x0 vector of length \code{p} giving the initial point.
#' @param drift drift for the diffusion.
#' @param SigDif matrix of size \code{c(p, p)} giving the infinitesimal (constant) covariance matrix of the diffusion.
#' @inheritParams rTrajWn2D
#' @param ... parameters to be passed to \code{drift}.
#' @param circular whether to wrap the resulting trajectory to \eqn{[-\pi,\pi)^p}.
#' @return A vector of length \code{N + 1} containing \code{x0} in the first entry and the discretized trajectory.
#' @details The fine trajectory is subsampled using the indexes \code{seq(1, NFine + 1, by = NFine / N)}.
#' @examples
#' isRStudio <- identical(.Platform$GUI, "RStudio")
#' if (isRStudio) {
#'   # 1D
#'   manipulate::manipulate({
#'     x <- seq(0, N * delta, by = delta)
#'     plot(x, x, ylim = c(-pi, pi), type = "n", ylab = expression(X[t]), xlab = "t")
#'     linesCirc(x, rTrajLangevin(x0 = 0, drift = driftJp, SigDif = sigma,
#'                                alpha = alpha, mu = 0, psi = psi, N = N,
#'                                delta = 0.01))
#'     }, delta = manipulate::slider(0.01, 5.01, step = 0.1),
#'     N = manipulate::slider(10, 500, step = 10, initial = 200),
#'     alpha = manipulate::slider(0.01, 5, step = 0.1, initial = 1),
#'     psi = manipulate::slider(-2, 2, step = 0.1, initial = 1),
#'     sigma = manipulate::slider(0.01, 5, step = 0.1, initial = 1))
#'
#'   # 2D
#'   samp <- rTrajLangevin(x0 = c(0, 0), drift = driftMvm, alpha = c(1, 1),
#'                         mu = c(2, -1), A = diag(rep(0, 2)),
#'                         SigDif = diag(rep(1, 2)), N = 1000, delta = 0.1)
#'   plot(samp, xlim = c(-pi, pi), ylim = c(-pi, pi), pch = 19, cex = 0.25,
#'        xlab = expression(X[t]), ylab = expression(Y[t]), col = rainbow(1000))
#'   linesTorus(samp[, 1], samp[, 2], col = rainbow(1000))
#' }
#' @export
rTrajLangevin <- function(x0, drift, SigDif, N = 100, delta = 0.01,
                          NFine = ceiling(N * delta / deltaFine),
                          deltaFine = min(delta / 100, 1e-3),
                          circular = TRUE, ...) {

  # Dimension
  p <- length(x0)

  # Correlate noise
  if (p == 1) {

    Z <- matrix(rnorm(NFine, mean = 0, sd = sqrt(deltaFine * SigDif)), ncol = 1)

  } else {

    Z <- mvtnorm::rmvnorm(NFine, mean = rep(0, p), sigma = deltaFine * SigDif)

  }

  # Sample by Euler
  samp <- matrix(nrow = NFine + 1, ncol = p)
  samp[1, ] <- x0
  for (i in 1:NFine) {

    samp[i + 1, ] <- samp[i, ] +
      drift(samp[i, , drop = FALSE], ...) * deltaFine + Z[i, ]

  }

  # Subsample
  samp <- samp[seq(1, NFine + 1, by = NFine / N), ]

  # Wrap it
  if (circular) {

    samp <- toPiInt(samp)

  }

  return(samp)

}
