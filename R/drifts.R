

#' @title Drift for the JP diffusion
#'
#' @description Drift for the Langevin diffusion associated to the Jones and Pewsey (JP) family of circular distributions.
#'
#' @param x vector with the evaluation points for the drift.
#' @param alpha strength of the drift.
#' @param mu unconditional mean of the diffusion.
#' @param psi shape parameter, see details.
#' @return A vector of the same length as \code{x} containing the drift.
#' @details Particular interesting choices for the shape parameter are:
#' \itemize{
#' \item \code{psi = -1}: gives the Wrapped Cauchy as stationary density.
#' \item \code{psi = 0}: is the sinusoidal drift of the vM diffusion.
#' \item \code{psi = 1}: gives the Cardioid as stationary density.
#' }
#' See Section 2.2.3 in García-Portugués et al. (2019) for details.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \url{https://doi.org/10.1007/s11222-017-9790-2}
#'
#' Jammalamadaka, S. R. and SenGupta, A. (2001) \emph{Topics in Circular Statistics}. World Scientific, Singapore. \url{https://doi.org/10.1142/4031}
#'
#' Jones, M. C. and Pewsey, A. (2005). A family of symmetric distributions on the circle. \emph{Journal of the American Statistical Association}, 100(472):1422--1428. \url{https://doi.org/10.1198/016214505000000286}
#' @examples
#' x <- seq(-pi, pi, l = 200)
#' plot(x, x, type = "n", ylab = "drift")
#' for (i in 0:20) {
#'   lines(x, driftJp(x = x, alpha = 1, mu = 0, psi = -1 + 2 * i / 20),
#'         col = rainbow(21)[i + 1])
#' }
#' @export
driftJp <- function(x, alpha, mu, psi) {

  if (psi != 0) {

    r <- 2 * alpha * psi
    er <- exp(r)
    erm <- 1 / er
    (er - erm) / (psi * ((er + erm) + (er - erm) * cos(mu - x))) * sin(mu - x)

  } else {

    alpha * sin(mu - x)

  }

}


#' @title Drift for the MvM diffusion
#'
#' @description Drift for the Langevin diffusion associated to the Multivariate von Mises (MvM) in dimension \code{p}.
#'
#' @param x matrix of size \code{c(n, p)} with the evaluation points for the drift.
#' @param mu vector of length \code{p} with the unconditional mean of the diffusion.
#' @param alpha vector of length \code{p} with the strength of the drift in the diagonal (\eqn{\sin}{sin} terms).
#' @param A matrix of size \code{c(p, p)} with the strength of the drift in cross terms (\eqn{\cos}{cos}-\eqn{\sin}{sin} terms). The diagonal has to be zero.
#' @return A matrix of the same size as \code{x} containing the drift.
#' @details See Section 2.2.1 in García-Portugués et al. (2019) for details.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \url{https://doi.org/10.1007/s11222-017-9790-2}
#' @examples
#' # 1D
#' x <- seq(-pi, pi, l = 200)
#' plot(x, x, type = "n", ylab = "drift")
#' for (i in 0:20) {
#'   lines(x, driftMvm(x = x, alpha = 3 * i / 20, mu = 0, A = 0),
#'         col = rainbow(21)[i + 1])
#' }
#'
#' # 2D
#' x <- seq(-pi, pi, l = 100)
#' plotSurface2D(x, x, f = function(x) sqrt(rowSums(driftMvm(x = x,
#'               alpha = c(2, 2), mu = c(-1, -1), A = rbind(c(0, 0), c(0, 0)))^2)),
#'               fVect = TRUE)
#' @export
driftMvm <- function(x, alpha, mu, A = 0) {

  # Dimension
  p <- length(mu)

  # Choice between vM and MvM
  if (p == 1) {

    drift <- alpha * sin(mu - x)

  } else {

    # Difference of x and mu
    x <- sweep(x, 2, mu, "-")

    # Cosine and sine
    cosx <- cos(x)
    sinx <- sin(x)

    # Drift: d(t(u) %*% A %*% v) = t(u) %*% A %*% d(v) + t(v) %*% t(A) %*% d(u)
    drift <- -sweep(sinx, 2, alpha, "*") + (sinx %*% A) * cosx

  }

  return(drift)

}


#' @title Drift for the mivM diffusion
#'
#' @description Drift for the Langevin diffusion associated to a mixture of \code{m} independent (multivariate) von Mises (mivM) of dimension \code{p}.
#'
#' @inheritParams driftMvm
#' @param A matrix of size \code{c(m, p)} giving the strengths of the drifts.
#' @param M matrix of size \code{c(m, p)} giving the means.
#' @param sigma diffusion coefficient.
#' @param p vector of length \code{m} giving the proportions. Must add to one.
#' @inheritParams safeSoftMax
#' @return A matrix of the same size as \code{x} containing the drift.
#' @details \code{\link{driftMixVm}} is more efficient for the circular case. The diffusion matrix is \eqn{\sigma\bold{I}}{sigma}. See Section 2.2.4 in García-Portugués et al. (2019) for details.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \url{https://doi.org/10.1007/s11222-017-9790-2}
#' @examples
#' # 1D
#' x <- seq(-pi, pi, l = 200)
#' plot(x, x, type = "n", ylab = "drift")
#' for (i in 1:10) {
#'   lines(x, driftMixIndVm(x = cbind(x), A = cbind(c(2, 2)),
#'         M = cbind(c(0, -pi + 2 * pi * i / 10)), sigma = 1, p = c(0.5, 0.5)),
#'         col = rainbow(10)[i])
#' }
#'
#' # 2D
#' x <- seq(-pi, pi, l = 100)
#' plotSurface2D(x, x, f = function(x) sqrt(rowSums(driftMixIndVm(x = x,
#'               A = rbind(c(1, 1), c(1, 1)), M = rbind(c(1, 1), c(-1, -1)),
#'               sigma = 1, p = c(0.25, 0.75))^2)), fVect = TRUE)
#' @export
driftMixIndVm <- function(x, A, M, sigma, p, expTrc = 30) {

  # Components
  m <- length(p)

  # Number of evaluation points
  n <- nrow(x)

  # Transpose x for using recyclying by columns
  x <- t(x)

  # Concentrations for mixtures
  K <- 2 * A / sigma^2

  # Compute weights
  logs <- matrix(0, nrow = m, ncol = n)
  for (j in 1:m) {

    # Logs of components (mu - x using recyclying by columns)
    logs[j, ] <- log(p[j]) + colSums(K[j, ] * (cos(M[j, ] - x) - 1)) - sum(logBesselI0Scaled(K[j, ]))

  }
  w <- safeSoftMax(logs = t(logs), expTrc = expTrc)

  # Compute weighted drift
  drift <- matrix(0, nrow = n, ncol = nrow(x))
  for (j in 1:m) {

    drift <- drift + t(A[j, ] * sin(M[j, ] - x)) * w[, j]

  }

  return(drift)

}


#' @title Drift for the mivM diffusion (circular case)
#'
#' @description Drift for the Langevin diffusion associated to a mixture of \code{m} independent von Mises (mivM) of dimension one.
#'
#' @inheritParams driftJp
#' @param alpha vector of length \code{m} giving the strengths of the drifts.
#' @param mu vector of length \code{m} giving the means.
#' @inheritParams driftMixIndVm
#' @inheritParams safeSoftMax
#' @details \code{\link{driftMixIndVm}} is more general, but less efficient for the circular case. See Section 2.2.4 in García-Portugués et al. (2019) for details.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \url{https://doi.org/10.1007/s11222-017-9790-2}
#' @return A vector of the same length as \code{x} containing the drift.
#' @examples
#' x <- seq(-pi, pi, l = 200)
#' plot(x, x, type = "n", ylab = "drift")
#' for (i in 1:10) {
#'   lines(x, driftMixVm(x = x, alpha = c(2, 2), mu = c(0, -pi + 2 * pi * i / 10),
#'         sigma = 1, p = c(0.5, 0.5)), col = rainbow(10)[i])
#' }
#' @export
driftMixVm <- function(x, alpha, mu, sigma, p, expTrc = 30) {

  # x as a matrix c(m, length(x)) for recyclying columns
  x <- repRow(x, length(p))

  # mu - x using recyclying by columns
  x <- mu - x

  # Concentrations for mixtures
  kappa <- 2 * alpha / sigma^2

  # Logs of components
  logs <- log(p) + kappa * (cos(x) - 1) - logBesselI0Scaled(x = kappa)

  # Weights
  w <- t(safeSoftMax(logs = t(logs), expTrc = 30))

  # Drift
  colSums(alpha * sin(x) * w)

}


#' @title Drift for the WN diffusion
#'
#' @description Drift for the Langevin diffusion associated to the (multivariate) Wrapped Normal (WN) in dimension \code{p}.
#'
#' @inheritParams driftMixIndVm
#' @inheritParams driftMvm
#' @inheritParams dTpdWou2D
#' @param A matrix of size \code{c(p, p)} giving the drift strength.
#' @param Sigma diffusion matrix, of size \code{c(p, p)}.
#' @param invSigmaA the matrix \code{solve(Sigma) \%*\% A} (optional).
#' @details See Section 2.2.2 in García-Portugués et al. (2019) for details.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \url{https://doi.org/10.1007/s11222-017-9790-2}
#' @return A matrix of the same size as \code{x} containing the drift.
#' @details \code{\link{driftWn1D}} and \code{\link{driftWn2D}} are more efficient for the 1D and 2D cases.
#' @examples
#' # 1D
#' x <- seq(-pi, pi, l = 200)
#' plot(x, x, type = "n", ylab = "drift")
#' for (i in 1:20) {
#'   lines(x, driftWn(x = cbind(x), A = 1 * i / 20, mu = 0, Sigma = 1),
#'         col = rainbow(20)[i])
#' }
#'
#' # 2D
#' x <- seq(-pi, pi, l = 100)
#' plotSurface2D(x, x, f = function(x) sqrt(rowSums(
#'               driftWn(x = x, A = alphaToA(alpha = c(1, 1, 0.5),
#'                                           sigma = c(1.5, 1.5)), mu = c(0, 0),
#'               Sigma = diag(c(1.5^2, 1.5^2)))^2)), fVect = TRUE)
#' @export
driftWn <- function(x, A, mu, Sigma, invSigmaA = NULL, maxK = 2, expTrc = 30) {

  # Winding numbers
  sK <- seq(-maxK, maxK, by = 1)

  # Inverse of the Gaussian exponent
  if (is.null(invSigmaA)) {

    invSigmaA <- solve(Sigma) %*% A

  }

  # Dimension
  p <- ncol(x)

  # Number of evaluation points
  n <- nrow(x)

  # Transpose x for using recyclying by columns
  x <- t(x)

  # mu - x using recyclying by columns
  x <- mu - x

  # Grid of windings
  twokpi <- as.matrix(do.call(what = expand.grid,
                              args = rep(list(2 * pi * sK), p)))
  lk <- nrow(twokpi)

  # Compute weights
  logs <- matrix(0, nrow = lk, ncol = n)
  for (j in 1:lk) {

    # x - twokpi using recyclying by columns
    v <- t(x - twokpi[j, ])

    # Logs of components
    logs[j, ] <- -rowSums((v %*% invSigmaA) * v)

  }
  w <- safeSoftMax(logs = t(logs), expTrc = expTrc)

  # Compute weighted drift
  drift <- matrix(0, nrow = n, ncol = p)
  for (j in 1:lk) {

    drift <- drift + t(x - twokpi[j, ]) * w[, j]

  }

  # Multiply by A
  drift <- drift %*% t(A)

  return(drift)

}
