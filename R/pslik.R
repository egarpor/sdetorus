

#' @title Wrapped Euler and Shoji--Ozaki pseudo-transition probability densities
#'
#' @description Wrapped pseudo-transition probability densities.
#'
#' @param x a matrix of dimension \code{c(n, p)}. If a vector is provided, is assumed that \code{p = 1}.
#' @param x0 a matrix of dimension \code{c(n, p)}. If all \code{x0} are the same, a matrix of dimension \code{c(1, p)} can be passed for better performance. If a vector is provided, is assumed that \code{p = 1}.
#' @param t time step between \code{x} and \code{x0}.
#' @param method a string for choosing \code{"E"} (Euler), \code{"SO"} (Shoji--Ozaki) or \code{"SO2"} (Shoji--Ozaki with Ito's expansion in the drift) method.
#' @param b drift function. Must return a matrix of the same size as \code{x}.
#' @param jac.b jacobian of the drift function.
#' @param b1 first derivative of the drift function (univariate). Must return a vector of the same length as \code{x}.
#' @param b2 second derivative of the drift function (univariate). Must return a vector of the same length as \code{x}.
#' @param sigma2 diagonal of the diffusion matrix (if univariate, this is the square of the diffusion coefficient). Must return an object of the same size as \code{x}.
#' @param circular flag to indicate circular data.
#' @param maxK maximum absolute winding number used if \code{circular = TRUE}.
#' @param vmApprox flag to indicate von Mises approximation to wrapped normal. See\cr \code{\link{momentMatchWnVm}} and \code{\link{scoreMatchWnBvm}}.
#' @param twokpi optional matrix of winding numbers to avoid its recomputation. See details.
#' @param ... additional parameters passed to \code{b}, \code{b1}, \code{b2}, \code{jac.b} and \code{sigma2}.
#' @return Output from \code{\link{mleOptimWrapper}}.
#' @details See Section 3.2 in García-Portugués et al. (2019) for details. \code{"SO2"} implements Shoji and Ozai (1998)'s expansion with for \code{p = 1}. \code{"SO"} is the same expansion, for arbitrary \code{p}, but considering null second derivatives.
#'
#' \code{twokpi} is \code{repRow(2 * pi * c(-maxK:maxK), n = n)} if \code{p = 1} and\cr \code{as.matrix(do.call(what = expand.grid, args = rep(list(2 * pi * c(-maxK:maxK)), p)))} otherwise.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \url{https://doi.org/10.1007/s11222-017-9790-2}
#'
#' Shoji, I. and Ozaki, T. (1998) A statistical method of estimation and simulation for systems of stochastic differential equations. \emph{Biometrika}, 85(1):240--243. \url{https://doi.org/10.1093/biomet/85.1.240}
#' @examples
#' # 1D
#' grid <- seq(-pi, pi, l = 501)[-501]
#' alpha <- 1
#' sigma <- 1
#' t <- 0.5
#' x0 <- pi/2
#' # manipulate::manipulate({
#'
#'   # Drifts
#'   b <- function(x) driftWn1D(x = x, alpha = alpha, mu = 0, sigma = sigma)
#'   b1 <- function(x, h = 1e-4) {
#'     l <- length(x)
#'     res <- driftWn1D(x = c(x + h, x - h), alpha = alpha, mu = 0, sigma = sigma)
#'     drop(res[1:l] - res[(l + 1):(2 * l)])/(2 * h)
#'   }
#'   b2 <- function(x, h = 1e-4) {
#'     l <- length(x)
#'     res <- driftWn1D(x = c(x + h, x, x - h), alpha = alpha, mu = 0,
#'                      sigma = sigma)
#'     drop(res[1:l] - 2 * res[(l + 1):(2 * l)] + res[(2 * l + 1):(3 * l)])/(h^2)
#'   }
#'
#'   # Squared diffusion
#'   sigma2 <- function(x) rep(sigma^2, length(x))
#'
#'   # Plot
#'   plot(grid, dTpdPde1D(Mx = length(grid), x0 = x0, t = t, alpha = alpha,
#'                        mu = 0, sigma = sigma), type = "l",
#'        ylab = "Density", xlab = "", ylim = c(0, 0.75), lwd = 2)
#'   lines(grid, dTpdWou1D(x = grid, x0 = rep(x0, length(grid)), t = t,
#'                        alpha = alpha, mu = 0, sigma = sigma), col = 2)
#'   lines(grid, dPsTpd(x = grid, x0 = x0, t = t, method = "E", b = b,
#'                      b1 = b1, b2 = b2, sigma2 = sigma2), col = 3)
#'   lines(grid, dPsTpd(x = grid, x0 = x0, t = t, method = "SO", b = b,
#'                      b1 = b1, b2 = b2, sigma2 = sigma2), col = 4)
#'   lines(grid, dPsTpd(x = grid, x0 = x0, t = t, method = "SO2", b = b,
#'                      b1 = b1, b2 = b2, sigma2 = sigma2),
#'         col = 5)
#'   lines(grid, dPsTpd(x = grid, x0 = x0, t = t, method = "E", b = b,
#'                      b1 = b1, b2 = b2, sigma2 = sigma2, vmApprox = TRUE),
#'         col = 6)
#'   lines(grid, dPsTpd(x = grid, x0 = x0, t = t, method = "SO", b = b,
#'                      b1 = b1, b2 = b2, sigma2 = sigma2, vmApprox = TRUE),
#'         col = 7)
#'   lines(grid, dPsTpd(x = grid, x0 = x0, t = t, method = "SO2", b = b,
#'                      b1 = b1, b2 = b2, sigma2 = sigma2, vmApprox = TRUE),
#'         col = 8)
#'   legend("topright", legend = c("PDE", "WOU", "E", "SO1", "SO2", "EvM",
#'                                 "SO1vM", "SO2vM"), lwd = 2, col = 1:8)
#'
#' # }, x0 = manipulate::slider(-pi, pi, step = 0.1, initial = -pi),
#' # alpha = manipulate::slider(0.1, 5, step = 0.1, initial = 1),
#' # sigma = manipulate::slider(0.1, 5, step = 0.1, initial = 1),
#' # t = manipulate::slider(0.1, 5, step = 0.1, initial = 1))
#'
#' # 2D
#' grid <- seq(-pi, pi, l = 76)[-76]
#' alpha1 <- 2
#' alpha2 <- 1
#' alpha3 <- 0.5
#' sig1 <- 1
#' sig2 <- 2
#' t <- 0.5
#' x01 <- pi/2
#' x02 <- -pi/2
#' # manipulate::manipulate({
#'
#'   alpha <- c(alpha1, alpha2, alpha3)
#'   sigma <- c(sig1, sig2)
#'   x0 <- c(x01, x02)
#'
#'   # Drifts
#'   b <- function(x) driftWn2D(x = x, A = alphaToA(alpha = alpha, sigma = sigma),
#'                              mu = rep(0, 2), sigma = sigma)
#'   jac.b <- function(x, h = 1e-4) {
#'     l <- nrow(x)
#'     res <- driftWn2D(x = rbind(cbind(x[, 1] + h, x[, 2]),
#'                                cbind(x[, 1] - h, x[, 2]),
#'                                cbind(x[, 1], x[, 2] + h),
#'                                cbind(x[, 1], x[, 2] - h)),
#'                      A = alphaToA(alpha = alpha, sigma = sigma), mu = rep(0, 2),
#'                      sigma = sigma)
#'     cbind(res[1:l, ] - res[(l + 1):(2 * l), ],
#'           res[2 * l + 1:l, ] - res[2 * l + (l + 1):(2 * l), ]) / (2 * h)
#'   }
#'
#'   # Squared diffusion
#'   sigma2 <- function(x) matrix(sigma^2, nrow = length(x) / 2L, ncol = 2)
#'
#'   # Plot
#'   old_par <- par(mfrow = c(3, 2))
#'   plotSurface2D(grid, grid, z = dTpdPde2D(Mx = length(grid), My = length(grid),
#'                                           x0 = x0, t = t, alpha = alpha,
#'                                           mu = rep(0, 2), sigma = sigma),
#'                 levels = seq(0, 1, l = 20), main = "Exact")
#'   plotSurface2D(grid, grid,
#'                 f = function(x) drop(dTpdWou2D(x = x, x0 = repRow(x0, nrow(x)),
#'                                                 t = t, alpha = alpha,
#'                                                 mu = rep(0, 2), sigma = sigma)),
#'                 levels = seq(0, 1, l = 20), fVect = TRUE, main = "WOU")
#'   plotSurface2D(grid, grid,
#'                 f = function(x) dPsTpd(x = x, x0 = rbind(x0), t = t,
#'                                        method = "E", b = b, jac.b = jac.b,
#'                                        sigma2 = sigma2),
#'                 levels = seq(0, 1, l = 20), fVect = TRUE, main = "E")
#'   plotSurface2D(grid, grid,
#'                 f = function(x) dPsTpd(x = x, x0 = rbind(x0), t = t,
#'                                        method = "SO", b = b, jac.b = jac.b,
#'                                        sigma2 = sigma2),
#'                 levels = seq(0, 1, l = 20), fVect = TRUE, main = "SO")
#'   plotSurface2D(grid, grid,
#'                 f = function(x) dPsTpd(x = x, x0 = rbind(x0), t = t,
#'                                        method = "E", b = b, jac.b = jac.b,
#'                                        sigma2 = sigma2, vmApprox = TRUE),
#'                 levels = seq(0, 1, l = 20), fVect = TRUE, main = "EvM")
#'   plotSurface2D(grid, grid,
#'                 f = function(x) dPsTpd(x = x, x0 = rbind(x0), t = t,
#'                                        method = "SO", b = b, jac.b = jac.b,
#'                                        sigma2 = sigma2, vmApprox = TRUE),
#'                 levels = seq(0, 1, l = 20), fVect = TRUE, main = "SOvM")
#'   par(old_par)
#'
#' # }, x01 = manipulate::slider(-pi, pi, step = 0.1, initial = -pi),
#' # x02 = manipulate::slider(-pi, pi, step = 0.1, initial = -pi),
#' # alpha1 = manipulate::slider(0.1, 5, step = 0.1, initial = 1),
#' # alpha2 = manipulate::slider(0.1, 5, step = 0.1, initial = 1),
#' # alpha3 = manipulate::slider(-5, 5, step = 0.1, initial = 0),
#' # sig1 = manipulate::slider(0.1, 5, step = 0.1, initial = 1),
#' # sig2 = manipulate::slider(0.1, 5, step = 0.1, initial = 1),
#' # t = manipulate::slider(0.01, 5, step = 0.01, initial = 1))
#' @export
dPsTpd <- function(x, x0, t, method = c("E", "SO", "SO2"), b, jac.b, sigma2, b1,
                   b2, circular = TRUE, maxK = 2, vmApprox = FALSE, twokpi = NULL,
                   ...) {

  # Get n and p
  if (is.matrix(x) & is.matrix(x0)) {

    n <- dim(x)[1]
    p <- dim(x)[2]

    # Has x0 the same size as x?
    if (ncol(x0) != p) {

      stop("Incompatible dimensions for x and x0")

    } else {

      nx0 <- nrow(x0)
      nx0 <- ifelse(nx0 == n, 1, n)

    }

  } else if (!is.matrix(x) & !is.matrix(x0)) {

    n <- length(x)
    p <- 1
    x <- matrix(x, nrow = n, ncol = 1)
    nx0 <- length(x0)
    nx0 <- ifelse(nx0 == n, 1, n)

  } else {

    stop("Incompatible dimensions for x and x0")

  }

  # Set maxK = 0 if the process is linear
  if (!circular) {

    maxK <- 0

  }

  # Series truncation
  sK <- -maxK:maxK

  # Translate data
  y <- x
  x <- x0

  if (method[1] == "E") {

    # Diffusion and drift, matrices of n x p, as we only encode the diagonal
    # of the diffusion coefficient
    bx <- b(x = x, ...)
    s2x <- sigma2(x = x, ...)

    if (vmApprox & circular) {

      # Get the closest von Mises kappa to the sigma^2 of the WN
      kx <- repRow(matrix(momentMatchWnVm(sigma2 = t * s2x), ncol = p), nx0)

      # Difference on the exponent
      cyx <- cos(y - repRow(x + bx * t, nx0))

      # Sum of independent-wise densities (because sigma2 is diagonal)
      dens <- log(2 * pi) + matrix(logBesselI0Scaled(x = kx), ncol = p) -
        kx * (cyx - 1)
      dens <- exp(-rowSums(dens))

    } else {

      # Windings
      if (is.null(twokpi)) {

        twokpi <- repRow(2 * pi * sK, n = n)

      }

      # Ensure they are matrices
      bx <- matrix(bx, ncol = p)
      s2x <- repRow(matrix(t * s2x, ncol = p), nx0)

      # Difference in the exponent
      y <- y - repRow(x + bx * t, nx0)

      # Wrap it for better accuracy of the truncated series if circular
      if (circular) {

        y <- toPiInt(y)

      }

      # Product of independent-wise densities (because sigma2 is diagonal)
      dens <- sapply(1:p, function(j) {

        # Wrapping
        dist <- sweep(twokpi, 1, y[, j], "+")

        # Final loglikelihood
        rowSums(exp(-(dist^2 / s2x[, j] + log(2 * pi * s2x[, j])) / 2))

      })
      dens <- apply(rbind(dens), 1, prod)

    }

  } else if (method[1] == "SO") {

    if (p == 1) {

      # Drift, first derivative and second derivative and diffusion
      bx <- b(x = x, ...)
      b1x <- b1(x = x, ...)
      s2x <- sigma2(x = x, ...)

      # Exponential
      ebx <- exp(b1x * t)

      # Expectation and variance
      Ex <- rep(x + bx / b1x * (ebx - 1), nx0)
      Vx <- rep(s2x / (2 * b1x) * (ebx^2 - 1), nx0)

      if (vmApprox & circular) {

        # Get the closest von Mises kappa to the sigma^2 of the WN
        kx <- momentMatchWnVm(sigma2 = Vx)

        # Log-density
        dens <- log(2 * pi) + logBesselI0Scaled(x = kx) - kx * (cos(y - Ex) - 1)
        dens <- drop(exp(-dens))

      } else {

        # Windings
        if (is.null(twokpi)) {

          twokpi <- repRow(2 * pi * sK, n = n)

        }

        # Difference in the exponent
        y <- y - Ex

        # Wrap it for better accuracy of the truncated series if circular
        if (circular) {

          y <- toPiInt(y)

        }

        # Wrappings
        dist <- sweep(twokpi, 1, y, "+")

        # Final density
        dens <- rowSums(exp(-(dist^2 / Vx + log(2 * pi * Vx)) / 2))

      }

    } else {

      # Windings
      if (is.null(twokpi)) {

        twokpi <- as.matrix(do.call(what = expand.grid,
                                    args = rep(list(2 * pi * sK), p)))

      }

      # Drift
      bx <- repRow(b(x = x, ...), nx0)

      # Identity matrix
      I <- diag(rep(1, p), ncol = p)

      if (nx0 == 1) {

        dens <- sapply(1:n, function(i) {

          # Diffusion
          s2xInv2 <- diag(2 / c(sigma2(x = x[i, ], ...)), nrow = p)

          # Jacobian and eigendecomposition
          jac.bx <- jac.b(x = x[i, , drop = FALSE], ...)
          eig.jac.bx <- eigen(x = jac.bx, symmetric = FALSE)
          eig.jac.bx$svectors <- solve(eig.jac.bx$vectors)

          # Expectation (x + P * D * diag(e^(lambda * t) - 1) * P^(-1) * b)
          Ex <- x[i, ] +
            eig.jac.bx$vectors %*% diag((exp(eig.jac.bx$values * t) - 1) /
                                          eig.jac.bx$values, nrow = p) %*%
            eig.jac.bx$svectors %*% bx[i, ]

          # Inverse covariance matrix (2 * V^(-1) * P *
          # diag(1 / (e^(2 * lambda * t) - 1)) * D * P^(-1))
          SVx <- s2xInv2 %*% eig.jac.bx$vectors %*%
            diag(eig.jac.bx$values / (exp(2 * eig.jac.bx$values * t) - 1),
                 nrow = p) %*% eig.jac.bx$svectors
          det.SVx <- det(SVx)

          if (p == 2 & vmApprox & circular) {

            # Density
            kx <- scoreMatchWnBvm(invSigma = SVx)
            kx[3] <- 2 * kx[3]
            dens <- dBvm(x = y[i, ], mu = Ex, kappa = kx)

          } else {

            # Recenter
            y[i, ] <- drop(y[i, ] - Ex)

            # Wrap it for better accuracy of the truncated series if circular
            if (circular) {

              y[i, ] <- toPiInt(y[i, ])

            }

            # Density
            dens <- sum(apply(twokpi, 1, function(wind)
              exp(-((y[i, ] + wind) %*% SVx %*% (y[i, ] + wind) -
                      log(det.SVx)) / 2)
              ))

          }

          dens

        })

      } else {

        # Diffusion
        s2xInv2 <- diag(2 / c(sigma2(x = x[1, ], ...)), nrow = p)

        # Jacobian and eigendecomposition
        jac.bx <- jac.b(x = x[1, , drop = FALSE], ...)
        eig.jac.bx <- eigen(x = jac.bx, symmetric = FALSE)
        eig.jac.bx$svectors <- solve(eig.jac.bx$vectors)

        # Expectation (x + P * D * diag(e^(lambda * t) - 1) * P^(-1) * b)
        Ex <- x[1, ] + eig.jac.bx$vectors %*%
          diag((exp(eig.jac.bx$values * t) - 1) /
                 eig.jac.bx$values, nrow = p) %*%
          eig.jac.bx$svectors %*% bx[1, ]

        # Inverse covariance matrix (2 * V^(-1) * P *
        # diag(1 / (e^(2 * lambda * t) - 1)) * D * P^(-1))
        SVx <- s2xInv2 %*% eig.jac.bx$vectors %*%
          diag(eig.jac.bx$values / (exp(2 * eig.jac.bx$values * t) - 1),
               nrow = p) %*% eig.jac.bx$svectors
        det.SVx <- det(SVx)

        if (p == 2 & vmApprox & circular) {

          # Density
          kx <- scoreMatchWnBvm(invSigma = SVx)
          kx[3] <- 2 * kx[3]
          dens <- dBvm(x = y, mu = Ex, kappa = kx)

        } else {

          # Recenter and wrap
          y <- sweep(y, 2, drop(Ex), "-")

          # Wrap it for better accuracy of the truncated series if circular
          if (circular) {

            y <- toPiInt(y)

          }

          # Wrapping of final density
          dens <- rowSums(apply(twokpi, 1, function(wind) {

            z <- sweep(y, 2, wind, "+")
            exp(-(rowSums((z %*% SVx) * z) - log(det.SVx) +
                    p * log(2 * pi)) / 2)

          }))

        }

      }

    }

  } else if (method[1] == "SO2" & p == 1) {

      # Drift, first derivative and second derivative and diffusion
      bx <- b(x = x, ...)
      b1x <- b1(x = x, ...)
      b2x <- b2(x = x, ...)
      s2x <- sigma2(x = x, ...)

      # Exponential
      ebx <- exp(b1x * t)

      # Expectation and variance
      Ex <- rep(x + bx / b1x * (ebx - 1) + s2x / 2 * b2x /
                  (b1x)^2 * (ebx - 1 - b1x * t), nx0)
      Vx <- rep(s2x / (2 * b1x) * (ebx^2 - 1), nx0)

      if (circular & vmApprox) {

        # Get the closest von Mises kappa to the sigma^2 of the WN
        kx <- momentMatchWnVm(sigma2 = Vx)

        # Log-density
        dens <- log(2 * pi) + logBesselI0Scaled(x = kx) - kx * (cos(y - Ex) - 1)
        dens <- drop(exp(-dens))

      } else {

        # Windings
        if (is.null(twokpi)) {

          twokpi <- repRow(2 * pi * sK, n = n)

        }

        # Difference in the exponent
        y <- y - Ex

        # Wrap it for better accuracy of the truncated series if circular
        if (circular) {

          y <- toPiInt(y)

        }

        # Wrappings
        dist <- sweep(twokpi, 1, y, "+")

        # Final density
        dens <- rowSums(exp(-(dist^2 / Vx + log(2 * pi * Vx)) / 2))

      }

  } else {

    stop("Only \"E\", \"SO\" and \"SO2\" (with p = 1) are supported")

  }

  return(dens)

}
