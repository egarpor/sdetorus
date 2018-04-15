

#' @title Efficient computation of Bessel related functions
#'
#' @description Computation of \eqn{\log(I_0(x))-x}{log(I0(x))-x} and the inverse of \eqn{A_1(k)=\frac{I_0(k)}{I_1(k)}}{A1(k)=I0(k)/I1(k)}.
#'
#' @param x evaluation vector. For \code{logBesselI0Scaled}, \code{x} must contain non-negative values. For \code{a1Inv}, \code{x} must be in \eqn{[0,1]}.
#' @param splineApprox whether to use a pre-computed spline approximation (faster) or not.
#' @return A vector of the same length as \code{x}.
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
#' @details Both functions may rely on pre-computed spline interpolations (\code{logBesselI0ScaledSpline} and \code{a1InvSpline}) that are created in \code{\link{.GlobalEnv}} the first time the functions are used. Otherwise, a call to \code{besselI} is done for \eqn{\log(I_0(x))-x}{log(I0(x))-x} and \eqn{A_1(k)=x}{A1(k)=x} is solved numerically. The data in which the interpolation is based is given in the examples.
#'
#' For \code{x} larger than \code{5e4}, the asymptotic expansion of \code{\link[Bessel]{besselIasym}} is employed.
#' @examples
#' \dontrun{
#' # Save grid for log besselI0 scaled
#' x1 <- c(seq(0, 1, by = 1e-4), seq(1 + 1e-2, 10, by = 1e-3),
#'         seq(10 + 1e-1, 100, by = 1e-2), seq(100 + 1e0, 1e3, by = 1e0),
#'         seq(1000 + 1e1, 5e4, by = 2e1))
#' logBesselI0ScaledEvalGrid <- log(besselI(x = x1, nu = 0, expon.scaled = TRUE))
#' save(list = "logBesselI0ScaledEvalGrid", file = "logBesselI0ScaledEvalGrid.rda",
#'      compress = TRUE)
#' }
#' \dontrun{
#' # Save grid for A1 inverse
#' x2 <- rev(c(seq(1e-04, 0.9 - 1e-4, by = 1e-4), seq(0.9, 1 - 1e-05, by = 1e-5)))
#' a1InvEvalGrid <- sapply(x2, function(k) {
#'   uniroot(f = function(x) k - besselI(x, nu = 1, expon.scaled = TRUE) /
#'           besselI(x, nu = 0, expon.scaled = TRUE),
#'           lower = 1e-06, upper = 1e+05, tol = 1e-15)$root
#' })
#' save(list = "a1InvEvalGrid", file = "a1InvEvalGrid.rda", compress = TRUE)
#' }
#' # Accuracy logBesselI0Scaled
#' x <- seq(0, 1e3, l = 1e3)
#' summary(logBesselI0Scaled(x = x, splineApprox = TRUE) -
#'         logBesselI0Scaled(x = x, splineApprox = FALSE))
#'
#' # Accuracy a1Inv
#' y <- seq(0, 1 - 1e-4, l = 1e3)
#' summary(a1Inv(x = y, splineApprox = TRUE) -
#'         a1Inv(x = y, splineApprox = FALSE))
#' @export
logBesselI0Scaled <- function(x, splineApprox = TRUE) {

  if (splineApprox) {

    if (is.null(logBesselI0ScaledSpline)) {

      assign(x = "logBesselI0ScaledSpline",
             value = splinefun(x = x1, y = logBesselI0ScaledEvalGrid),
             pos = environment(logBesselI0Scaled))

    }

    res <- numeric(length(x))
    indAsymp <- x >= 5e4
    indNoAsymp <- !indAsymp
    res[indNoAsymp] <- logBesselI0ScaledSpline(x[indNoAsymp])

    if (any(indAsymp)) {

      res[indAsymp] <- Bessel::besselIasym(x = x[indAsymp], nu = 0, 
                                           expon.scaled = TRUE, log = TRUE)

    }

  } else {

    res <- log(besselI(x = x, nu = 0, expon.scaled = TRUE))

  }

  return(res)

}


#' @rdname logBesselI0Scaled
#' @export
a1Inv <- function(x, splineApprox = TRUE) {

  if (splineApprox) {

    if (is.null(a1InvSpline)) {

      assign(x = "a1InvSpline",
             value = splinefun(x = x2, y = a1InvEvalGrid),
             envir = environment(a1Inv))

    }

    res <- pmax(a1InvSpline(x), 0)
    indOne <- x >= 1

    if (any(indOne)) {

      message("a1Inv: x larger or equal to 1, thus not on the image of A1. Setting a1Inv(x) = Inf.")
      res[indOne] <- Inf

    }

  } else {

    res <- sapply(x, function(y) uniroot(f = function(k) besselI(x = k, nu = 1, expon.scaled = TRUE) / besselI(x = k, nu = 0, expon.scaled = TRUE) - y, interval = c(0, 1e4), tol = 1e-15)$root)

  }

  return(res)

}


logBesselI0ScaledSpline <- NULL
a1InvSpline <- NULL


#' @title Score and moment matching of a univariate or bivariate wrapped normal by a von Mises
#'
#' @description Given a wrapped normal density, find the parameters of a von Mises that matches it according to two characteristics: moments and scores. Score matching estimators are available for univariate and bivariate cases and moment matching only for the former.
#'
#' @param sigma,sigma2 standard deviation or variance of the wrapped normal.
#' @param Sigma,invSigma covariance or precision matrix of the bivariate wrapped normal.
#' @return Vector of parameters \eqn{(\kappa_1,\kappa_2,\lambda)}, where \eqn{(\kappa_1,\kappa_2,2\lambda)} is a suitable input for \code{kappa} in \code{dBvm}.
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
#' @details If the precision matrix is singular or if there are no solutions for the score matching estimator, \code{c(0, 0, 0)} is returned.
#' @references
#' Mardia, K. V., Kent, J. T., Laha, A. K. (2016). Score matching estimators for directional distributions. \emph{arXiv:1604.0847}.
#' @examples
#' # Univariate WN approximation
#' sigma <- 0.5
#' curve(dWn1D(x = x, mu = 0, sigma = sigma), from = -pi, to = pi, ylab = "Density",
#'       ylim = c(0, 1))
#' curve(dVm(x = x, mu = 0, kappa = momentMatchWnVm(sigma = sigma)), from = -pi,
#'       to = pi, col = "red", add = TRUE)
#' curve(dVm(x = x, mu = 0, kappa = scoreMatchWnVm(sigma = sigma)), from = -pi,
#'       to = pi, col = "green", add = TRUE)
#'
#' # Bivariate WN approximation
#'
#' # WN
#' alpha <- c(2, 1, 1)
#' sigma <- c(1, 1)
#' mu <- c(pi / 2, pi / 2)
#' x <- seq(-pi, pi, l = 101)[-101]
#' plotSurface2D(x, x, f = function(x) dStatWn2D(x = x, alpha = alpha, mu = mu,
#'                                               sigma = sigma), fVect = TRUE)
#' A <- alphaToA(alpha = alpha, sigma = sigma)
#' S <- 0.5 * solve(A) %*% diag(sigma)
#'
#' # Score matching
#' kappa <- scoreMatchWnBvm(Sigma = S)
#'
#' # dBvm uses lambda / 2 in the exponent
#' plotSurface2D(x, x, f = function(x)
#'               sdetorus:::dBvm(x = x, mu = mu,
#'                               kappa = c(kappa[1:2], 2 * kappa[3])),
#'               fVect = TRUE)
#'
#' # With singular Sigma
#' invSigma <- matrix(c(1, sqrt(0.999), sqrt(0.999), 1), nrow = 2, ncol = 2)
#' scoreMatchWnBvm(invSigma = invSigma)
#' invSigma <- matrix(1, nrow = 2, ncol = 2)
#' scoreMatchWnBvm(invSigma = invSigma)
#' @export
scoreMatchWnBvm <- function(Sigma, invSigma) {

  # Compute Sigma
  if (missing(Sigma)) {

    Sigma <- tryCatch(solve(invSigma), error = function(e) NA)

  }

  # Check for regularity
  if (anyNA(Sigma)) {

    message("scoreMatchWnBvm: non-invertible invSigma, setting kappa = c(0, 0, 0)")
    kappa <- c(0, 0, 0)

  } else {

    # Vector d - extra exp's to avoid overflows
    d <- c(exp(-0.5 * Sigma[1, 1]), exp(-0.5 * Sigma[2, 2]))
    d <- c(d, 0.5 * (exp(-0.5 * (Sigma[1, 1] + Sigma[2, 2]) + Sigma[1, 2]) -
                     exp(-0.5 * (Sigma[1, 1] + Sigma[2, 2]) - Sigma[1, 2])))

    # Matrix W computed robustly
    W <- matrix(0, nrow = 3, ncol = 3)
    W[1, 1] <- 0.5 * (1 - exp(-2 * Sigma[1, 1]))
    W[2, 2] <- 0.5 * (1 - exp(-2 * Sigma[2, 2]))
    W[3, 3] <- 0.5 * (1 - 0.5 * (
      exp(2 * (-(Sigma[1, 1] + Sigma[2, 2]) + 2 * Sigma[1, 2])) +
      exp(2 * (-(Sigma[1, 1] + Sigma[2, 2]) - 2 * Sigma[1, 2]))
      ))
    W[1, 3] <- W[3, 1] <- -0.25 * (
      exp(-(2 * Sigma[1, 1] + 0.5 * Sigma[2, 2]) + 2 * Sigma[1, 2]) -
      exp(-(2 * Sigma[1, 1] + 0.5 * Sigma[2, 2]) - 2 * Sigma[1, 2])
    )
    W[2, 3] <- W[3, 2] <- -0.25 * (
      exp(-(2 * Sigma[2, 2] + 0.5 * Sigma[1, 1]) + 2 * Sigma[1, 2]) -
      exp(-(2 * Sigma[2, 2] + 0.5 * Sigma[1, 1]) - 2 * Sigma[1, 2])
    )

    # Get (kappa1, kappa2, lambda)
    if (any(!is.finite(W) | !is.finite(d))) {

      message("scoreMatchWnBvm: non-finite values in W and d, setting kappa = c(0, 0, 0)")
      kappa <- c(0, 0, 0)

    } else {

      kappa <- tryCatch(solve(a = W, b = d), error = function(e) {
        message("scoreMatchWnBvm: W singular, setting kappa = c(0, 0, 0)")
        NA})

    }

  }

  return(kappa)

}


#' @rdname scoreMatchWnBvm
#' @export
scoreMatchWnVm <- function(sigma, sigma2) {

  if (missing(sigma2)) {

    sigma2 <- sigma^2

  }
  c <- exp(-0.5 * sigma2)
  2 * c / (1 - c^4)

}


#' @rdname scoreMatchWnBvm
#' @export
momentMatchWnVm <- function(sigma, sigma2) {

  if (missing(sigma2)) {

    sigma2 <- sigma^2

  }

  a1Inv(x = exp(-0.5 * sigma2))

}
