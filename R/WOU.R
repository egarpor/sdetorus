

#' @title Conditional probability density of the WOU process
#'
#' @description Conditional probability density of the Wrapped
#' Ornstein--Uhlenbeck (WOU) process.
#'
#' @param x matrix of size \code{c(n, p)} with the evaluation points in
#' \eqn{[-\pi,\pi)^p}.
#' @inheritParams dTpdWou1D
#' @inheritParams driftWn
#' @param x0 vector of length \code{p} with the initial point in
#' \eqn{[-\pi,\pi)^p}.
#' @inheritParams dTpdMou
#' @param invASigma the matrix \code{solve(Sigma) \%*\% A} (optional).
#' @return A vector of length \code{n} with the density evaluated at \code{x}.
#' @details See Section 3.3 in García-Portugués et al. (2019) for details.
#' \code{\link{dTpdWou1D}} and \code{\link{dTpdWou2D}} are more efficient
#' implementations for the 1D and 2D cases, respectively.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019)
#' Langevin diffusions on the torus: estimation and applications.
#' \emph{Statistics and Computing}, 29(2):1--22. \doi{10.1007/s11222-017-9790-2}
#' @examples
#' # 1D
#' t <- 0.5
#' alpha <- 1
#' mu <- 0
#' sigma <- 1
#' x0 <- pi
#' x <- seq(-pi, pi, l = 10)
#' dTpdWou(x = cbind(x), x0 = x0, t = t, A = alpha, mu = 0, Sigma = sigma^2) -
#' dTpdWou1D(x = cbind(x), x0 = rep(x0, 10), t = t, alpha = alpha, mu = 0,
#'           sigma = sigma)
#'
#' # 2D
#' t <- 0.5
#' alpha <- c(2, 1, -1)
#' sigma <- c(1.5, 2)
#' rho <- 0.9
#' Sigma <- diag(sigma^2)
#' Sigma[1, 2] <- Sigma[2, 1] <- rho * prod(sigma)
#' A <- alphaToA(alpha = alpha, sigma = sigma, rho = rho)
#' mu <- c(pi, 0)
#' x0 <- c(0, 0)
#' x <- seq(-pi, pi, l = 5)
#' x <- as.matrix(expand.grid(x, x))
#' dTpdWou(x = x, x0 = x0, t = t, A = A, mu = mu, Sigma = Sigma) -
#' dTpdWou2D(x = x, x0 = rbind(x0), t = t, alpha = alpha, mu = mu,
#'           sigma = sigma, rho = rho)
#' @export
dTpdWou <- function(x, t, A, mu, Sigma, x0, maxK = 2, eigA = NULL,
                    invASigma = NULL) {

  # Read dimension and number of evaluation points
  p <- ncol(x)
  nx <- nrow(x)
  if (is.null(p)) {
    p <- length(x)
    nx <- 1
  }

  # Winding numbers
  sk <- seq(-maxK, maxK, by = 1)

  # Grid of windings
  largs <- vector(mode = "list", length = p)
  for (i in 1:p) largs[[i]] <- sk
  grid.K <- as.matrix(do.call(expand.grid, args = largs))
  length.grid.K <- nrow(grid.K)

  # Eigen decomposition of A
  if (is.null(eigA)) {

    eigA <- eigen(A)

  }

  # Inverse of the product of matrices
  if (is.null(invASigma)) {

    invASigma <- solve(Sigma) %*% A

  }

  # Matrix exponential of -t*A
  ExptA <- eigA$vectors %*% diag(exp(-eigA$values * t), nrow = p, ncol = p) %*%
    solve(eigA$vectors)

  # Difference between x and mu by rows
  x0.minus.mu <- x0 - mu

  # Computation of conditional covariance
  Gamma <- covtMou(t = t, A = A, Sigma = Sigma, eigA = eigA)

  # Inverse of conditional covariance
  inv.Gamma <- solve(Gamma)

  # Initialize weights
  w.m <- numeric(length.grid.K)

  # Initialize shifted conditional means
  mutVmm <- matrix(mu + ExptA %*% (x0 - mu), nrow = length.grid.K,
                   ncol = p, byrow = TRUE)

  # Compute the weights and the shifted conditional means
  for (ki in 1:length.grid.K) {

    # The k
    k <- grid.K[ki, ]

    # Vector of shifted th-mu
    v <- x0.minus.mu + 2 * k * pi

    # Exponential of the weights (try to save computations later)
    w.m[ki] <- exp(-(t(v) %*% invASigma) %*% v)

    # Shifted conditional mean
    mutVmm[ki, ] <- mutVmm[ki, ] + 2 * pi * ExptA %*% k

  }

  # Standardize
  w.m <- w.m / sum(w.m)

  # Final result
  result <- numeric(nx)

  # Loop for the winding numbers weights
  for (mi in 1:length.grid.K) {

    # Difference between x and mutVmm[mi,] by rows
    x.minus.mutVmm <- sweep(x, 2, mutVmm[mi, ], "-")

    # Cumulant in the inner loop
    result.k <- numeric(nx)

    # Loop for the wrapping
    for (ki in 1:length.grid.K) {

      # The k
      k <- grid.K[ki, ]

      # Vector of shifted th-mu
      v <- sweep(x.minus.mutVmm, 2, 2 * k * pi, "+")

      # Exponential of the conditional Gaussian
      result.k <- result.k + exp(-rowSums((v %*% inv.Gamma) * v) / 2)

    }

    # Outer cumulant
    result <- result + result.k * w.m[mi]

  }

  # Standardize
  result <- result / sqrt((2 * pi)^p * det(Gamma))

  return(result)

}
