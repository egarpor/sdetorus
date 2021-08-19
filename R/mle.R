

#' @title Maximum pseudo-likelihood estimation by wrapped pseudo-likelihoods
#'
#' @description Maximum pseudo-likelihood using the Euler and Shoji--Ozaki pseudo-likelihoods.
#'
#' @param data a matrix of dimension \code{c(n, p)}.
#' @param delta discretization step.
#' @inheritParams dPsTpd
#' @inheritParams mleOptimWrapper
#' @param ... further parameters passed to \code{\link{mleOptimWrapper}}.
#' @return Output from \code{\link{mleOptimWrapper}}.
#' @details See Section 3.2 in García-Portugués et al. (2019) for details. \code{"SO2"} implements Shoji and Ozai (1998)'s expansion with for \code{p = 1}. \code{"SO"} is the same expansion, for arbitrary \code{p}, but considering null second derivatives.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \doi{10.1007/s11222-017-9790-2}
#'
#' Shoji, I. and Ozaki, T. (1998) A statistical method of estimation and simulation for systems of stochastic differential equations. \emph{Biometrika}, 85(1):240--243. \doi{10.1093/biomet/85.1.240}
#' @examples
#' \donttest{
#' # Example in 1D
#'
#' delta <- 0.5
#' pars <- c(0.25, 0, 2)
#' set.seed(12345678)
#' samp <- rTrajWn1D(x0 = 0, alpha = pars[1], mu = pars[2], sigma = pars[3],
#'                   N = 100, delta = delta)
#' b <- function(x, pars) driftWn1D(x = x, alpha = pars[1], mu = pars[2],
#'                                  sigma = pars[3], maxK = 2, expTrc = 30)
#' b1 <- function(x, pars, h = 1e-4) {
#'   l <- length(x)
#'   res <- b(x = c(x + h, x - h), pars = pars)
#'   drop(res[1:l] - res[(l + 1):(2 * l)])/(2 * h)
#' }
#' b2 <- function(x, pars, h = 1e-4) {
#'   l <- length(x)
#'   res <- b(x = c(x + h, x, x - h), pars = pars)
#'   drop(res[1:l] - 2 * res[(l + 1):(2 * l)] + res[(2 * l + 1):(3 * l)])/(h^2)
#' }
#' sigma2 <- function(x, pars) rep(pars[3]^2, length(x))
#' lower <- c(0.1, -pi, 0.1)
#' upper <- c(10, pi, 10)
#' psMle(data = samp, delta = delta, method = "E", b = b, sigma2 = sigma2,
#'       start = pars, lower = lower, upper = upper)
#' psMle(data = samp, delta = delta, method = "E", b = b, sigma2 = sigma2,
#'       start = pars, lower = lower, upper = upper, vmApprox = TRUE)
#' psMle(data = samp, delta = delta, method = "SO2", b = b, b1 = b1,
#'       b2 = b2, sigma2 = sigma2, start = pars, lower = lower, upper = upper)
#' psMle(data = samp, delta = delta, method = "SO2", b = b, b1 = b1,
#'       b2 = b2, sigma2 = sigma2, start = pars, lower = lower,
#'       upper = upper, vmApprox = TRUE)
#' psMle(data = samp, delta = delta, method = "SO", b = b, b1 = b1,
#'       lower = lower, upper = upper, sigma2 = sigma2, start = pars)
#' approxMleWn1D(data = samp, delta = delta, start = pars)
#' mlePde1D(data = samp, delta = delta, b = b, sigma2 = sigma2,
#'          start = pars, lower = lower, upper = upper)
#'
#' # Example in 2D
#'
#' delta <- 0.5
#' pars <- c(1, 0.5, 0, 0, 0, 1, 2)
#' set.seed(12345678)
#' samp <- rTrajWn2D(x0 = c(0, 0), alpha = pars[1:3], mu = pars[4:5],
#'                   sigma = pars[6:7], N = 100, delta = delta)
#' b <- function(x, pars) driftWn2D(x = x, A = alphaToA(alpha = pars[1:3],
#'                                                      sigma = pars[6:7]),
#'                                  mu = pars[4:5], sigma = pars[6:7], maxK = 2,
#'                                  expTrc = 30)
#' jac.b <- function(x, pars, h = 1e-4) {
#'   l <- nrow(x)
#'   res <- b(x = rbind(cbind(x[, 1] + h, x[, 2]),
#'                      cbind(x[, 1] - h, x[, 2]),
#'                      cbind(x[, 1], x[, 2] + h),
#'                      cbind(x[, 1], x[, 2] - h)), pars = pars)
#'   cbind(res[1:l, ] - res[(l + 1):(2 * l), ],
#'         res[2 * l + 1:l, ] - res[2 * l + (l + 1):(2 * l), ]) / (2 * h)
#' }
#' sigma2 <- function(x, pars) matrix(pars[6:7]^2, nrow = length(x) / 2L, ncol = 2)
#' lower <- c(0.01, 0.01, -25, -pi, -pi, 0.01, 0.01)
#' upper <- c(25, 25, 25, pi, pi, 25, 25)
#' psMle(data = samp, delta = delta, method = "E", b = b, sigma2 = sigma2,
#'       start = pars, lower = lower, upper = upper)
#' psMle(data = samp, delta = delta, method = "E", b = b, sigma2 = sigma2,
#'       start = pars, lower = lower, upper = upper, vmApprox = TRUE)
#' psMle(data = samp, delta = delta, method = "SO", b = b, jac.b = jac.b,
#'       sigma2 = sigma2, start = pars, lower = lower, upper = upper)
#' approxMleWn2D(data = samp, delta = delta, start = c(pars, 0))
#' # Set maxit = 5 to test and avoid a very long evaluation
#' mlePde2D(data = samp, delta = delta, b = b, sigma2 = sigma2, start = pars,
#'          lower = lower, upper = upper, maxit = 5)
#' }
#' @export
psMle <- function(data, delta, method = c("E", "SO", "SO2"), b, jac.b, sigma2,
                  b1, b2, start, lower, upper, circular = TRUE, maxK = 2,
                  vmApprox = FALSE, ...) {

  # Get N
  if (is.matrix(data)) {

    N <- dim(data)[1]
    p <- dim(data)[2]

  } else {

    N <- length(data)
    p <- 1
    data <- matrix(data, ncol = 1, nrow = N)

  }

  # Set maxK = 0 if the process is linear
  if (!circular) {

    maxK <- 0

  }

  # Translate data
  y <- data[-1, , drop = FALSE]
  x <- data[-N, , drop = FALSE]

  # Series truncation
  sK <- -maxK:maxK

  # twokpi to avoid recomputing
  if (!vmApprox) {

    if (method == "SO" & p > 1) {

      twokpi <- as.matrix(do.call(what = expand.grid,
                                  args = rep(list(2 * pi * sK), p)))

    } else {

      twokpi <- sapply(sK, function(k) rep(2 * pi * k, N - 1))

    }

  }

  minusLogLik <- function(pars) {

    -sum(log(dPsTpd(x = y, x0 = x, t = delta, method = method, b = b,
                    jac.b = jac.b, b1 = b1, b2 = b2, sigma2 = sigma2,
                    circular = circular, maxK = maxK, vmApprox = vmApprox,
                    twokpi = twokpi, pars = pars)))

  }

  # Optimization
    mleOptimWrapper(minusLogLik = minusLogLik, start = start, lower = lower,
                    upper = upper, ...)

}


#' @title Approximate MLE of the WN diffusion in 1D
#'
#' @description Approximate Maximum Likelihood Estimation (MLE) for the Wrapped Normal (WN) in 1D using the wrapped Ornstein--Uhlenbeck diffusion.
#'
#' @inheritParams psMle
#' @param alpha,mu,sigma if their values are provided, the likelihood function is optimized with respect to the rest of unspecified parameters. The number of elements in \code{start}, \code{lower} and \code{upper} has to be modified accordingly (see examples).
#' @inheritParams safeSoftMax
#' @return Output from \code{\link{mleOptimWrapper}}.
#' @details See Section 3.3 in García-Portugués et al. (2019) for details.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \doi{10.1007/s11222-017-9790-2}
#' @examples
#' alpha <- 0.5
#' mu <- 0
#' sigma <- 2
#' samp <- rTrajWn1D(x0 = 0, alpha = alpha, mu = mu, sigma = sigma, N = 1000,
#'                     delta = 0.1)
#' approxMleWn1D(data = samp, delta = 0.1, start = c(alpha, mu, sigma))
#' approxMleWn1D(data = samp, delta = 0.1, sigma = sigma, start = c(alpha, mu),
#'                 lower = c(0.01, -pi), upper = c(25, pi))
#' approxMleWn1D(data = samp, delta = 0.1, mu = mu, start = c(alpha, sigma),
#'                 lower = c(0.01, 0.01), upper = c(25, 25))
#' @export
approxMleWn1D <- function(data, delta, start, alpha = NA, mu = NA, sigma = NA,
                          lower = c(0.01, -pi, 0.01), upper = c(25, pi, 25),
                          vmApprox = FALSE, maxK = 2, ...) {

  # Get N
  N <- length(data)

  # Translate data
  y <- data[-1]
  x <- data[-N]

  # Specfied parameters
  specPars <- c(alpha, mu, sigma)
  indUnSpecPars <- is.na(specPars)

  minusLogLik <- function(pars) {

    # Change only unspecified parameters
    specPars[indUnSpecPars] <- pars

    -sum(log(dTpdWou1D(x = y, t = delta, alpha = specPars[1], mu = specPars[2],
                      sigma = specPars[3], x0 = x, maxK = maxK)))

  }

  # Optimization
    mleOptimWrapper(minusLogLik = minusLogLik, start = start, lower = lower,
                    upper = upper, ...)

}


#' @title Approximate MLE of the WN diffusion in 2D
#'
#' @description Approximate Maximum Likelihood Estimation (MLE) for the Wrapped Normal (WN) in 2D using the wrapped Ornstein--Uhlenbeck diffusion.
#'
#' @inheritParams psMle
#' @inheritParams approxMleWn1D
#' @param alpha,mu,sigma,rho if their values are provided, the likelihood function is optimized with respect to the rest of unspecified parameters. The number of elements in \code{start}, \code{lower} and \code{upper} has to be modified accordingly (see examples).
#' @inheritParams safeSoftMax
#' @return Output from \code{\link{mleOptimWrapper}}.
#' @details See Section 3.3 in García-Portugués et al. (2019) for details.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \doi{10.1007/s11222-017-9790-2}
#' @examples
#' alpha <- c(2, 2, -0.5)
#' mu <- c(0, 0)
#' sigma <- c(1, 1)
#' rho <- 0.2
#' samp <- rTrajWn2D(x0 = c(0, 0), alpha = alpha, mu = mu, sigma = sigma,
#'                   rho = rho, N = 1000, delta = 0.1)
#' approxMleWn2D(data = samp, delta = 0.1, start = c(alpha, mu, sigma, rho))
#' approxMleWn2D(data = samp, delta = 0.1, alpha = alpha,
#'               start = c(mu, sigma), lower = c(-pi, -pi, 0.01, 0.01),
#'               upper = c(pi, pi, 25, 25))
#' mleMou(data = samp, delta = 0.1, start = c(alpha, mu, sigma),
#'        optMethod = "Nelder-Mead")
#' @export
approxMleWn2D <- function(data, delta, start, alpha = rep(NA, 3),
                          mu = rep(NA, 2), sigma = rep(NA, 2), rho = NA,
                          lower = c(0.01, 0.01, -25, -pi, -pi, 0.01, 0.01, -0.99),
                          upper = c(rep(25, 3), pi, pi, 25, 25, 0.99),
                          maxK = 2, ...) {

  # Get N
  if (!is.matrix(data)) {

    stop("data must be a matrix")

  }
  N <- nrow(data)
  p <- ncol(data)

  # Translate data
  y <- data[-1, ]
  x <- data[-N, ]

  # Specfied parameters
  specPars <- c(alpha, mu, sigma, rho)
  indUnSpecPars <- is.na(specPars)

  minusLogLik <- function(pars) {

    # Change only unspecified parameters
    specPars[indUnSpecPars] <- pars

    -sum(log(dTpdWou2D(x = y, x0 = x, t = delta, alpha = specPars[1:3],
                      mu = specPars[4:5], sigma = specPars[6:7],
                      rho = specPars[8], maxK = maxK)))

  }

  # Check positvedefiniteness (if alpha is not provided)
  if (any(is.na(alpha))) {

    region <- function(pars) {

      # Test
      prodDiagonal = 0.25 * (pars[8] * (pars[2] - pars[1]))^2 +
        alpha[1] * alpha[2] + pars[1] * pars[2]
      testPosDef <- prodDiagonal - pars[3]^2
      if (testPosDef < 0) {

        pars[3] <- sign(pars[3]) *  sqrt(prodDiagonal) * 0.999
        penalty <- 1e3 * testPosDef

      } else {

        penalty <- 0

      }

      return(list("pars" = pars, "penalty" = -penalty))

    }

  } else {

    region <- function(pars) {

      list("pars" = pars, "penalty" = 0)

    }

  }

  # Optimization
  mleOptimWrapper(minusLogLik = minusLogLik, start = start, lower = lower,
                  upper = upper, ...)

}


#' @title High-frequency estimate of the diffusion matrix
#'
#' @description Estimation of the \eqn{\Sigma} in the multivariate diffusion
#' \deqn{dX_t=b(X_t)dt+\Sigma dW_t}{dX_t=b(X_t)dt+\Sigma dW_t}
#' by the high-frequency estimate
#' \deqn{\hat\Sigma = \frac{1}{N\Delta}\sum_{i=1}^N(X_i-X_{i-1})(X_i-X_{i-1})^T}{\hat\Sigma = \frac{1}{N\Delta}\sum_{i=1}^N(X_i-X_{i-1})(X_i-X_{i-1})^T.}
#'
#' @param data vector or matrix of size \code{c(N, p)} containing the discretized process.
#' @param delta discretization step.
#' @param circular whether the process is circular or not.
#' @param diagonal,isotropic enforce different constraints for the diffusion matrix.
#' @return The estimated diffusion matrix of size \code{c(p, p)}.
#' @details See Section 3.1 in García-Portugués et al. (2019) for details.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019) Langevin diffusions on the torus: estimation and applications. \emph{Statistics and Computing}, 29(2):1--22. \doi{10.1007/s11222-017-9790-2}
#' @examples
#' # 1D
#' x <- drop(euler1D(x0 = 0, alpha = 1, mu = 0, sigma = 1, N = 1000,
#'                   delta = 0.01))
#' sigmaDiff(x, delta = 0.01)
#'
#' # 2D
#' x <- t(euler2D(x0 = rbind(c(pi, pi)), A = rbind(c(2, 1), c(1, 2)),
#'                mu = c(pi, pi), sigma = c(1, 1), N = 1000,
#'                delta = 0.01)[1, , ])
#' sigmaDiff(x, delta = 0.01)
#' sigmaDiff(x, delta = 0.01, circular = FALSE)
#' sigmaDiff(x, delta = 0.01, diagonal = TRUE)
#' sigmaDiff(x, delta = 0.01, isotropic = TRUE)
#' @export
sigmaDiff <- function(data, delta, circular = TRUE, diagonal = FALSE,
                      isotropic = FALSE) {

  # Dimension and number of points
  p <- ncol(data)
  if (is.null(p)) {

    N <- length(data) - 1
    p <- 1

  } else {

    N <- nrow(data) - 1

  }

  # Covariance matrix
  dif <- diffCirc(data, circular = circular)
  sigmaDiff <- crossprod(dif) / (delta * N)

  # Constrain the matrix if requested
  if (isotropic) {

    sigmaDiff <- diag(mean(diag(sigmaDiff)), nrow = p, ncol = p)

  } else if (diagonal) {

    sigmaDiff <- diag(diag(sigmaDiff), nrow = p, ncol = p)

  }

  return(sigmaDiff)

}


#' @title Approximate MLE of the WN diffusion in 2D from a sample of initial and final pairs of angles.
#'
#' @description Approximate Maximum Likelihood Estimation (MLE) for the Wrapped Normal (WN) diffusion, using the wrapped Ornstein--Uhlenbeck diffusion and assuming initial stationarity.
#'
#' @inheritParams logLikWouPairs
#' @inheritParams approxMleWn2D
#' @inheritParams safeSoftMax
#' @inheritParams psMle
#' @return Output from \code{\link{mleOptimWrapper}}.
#' @examples
#' mu <- c(0, 0)
#' alpha <- c(1, 2, 0.5)
#' sigma <- c(1, 1)
#' rho <- 0.5
#' set.seed(4567345)
#' begin <- rStatWn2D(n = 200, mu = mu, alpha = alpha, sigma = sigma)
#' end <- t(apply(begin, 1, function(x) rTrajWn2D(x0 = x, alpha = alpha,
#'                                                mu = mu, sigma = sigma,
#'                                                rho = rho, N = 1,
#'                                                delta = 0.1)[2, ]))
#' data <- cbind(begin, end)
#' approxMleWnPairs(data = data, delta = 0.1,
#'                  start = c(2, pi/2, 2, 0.5, 0, 2, 1, 0.5))
#' @export
approxMleWnPairs <- function(data, delta, start = c(0, 0, 1, 1, 0, 1, 1),
                             alpha = rep(NA, 3), mu = rep(NA, 2),
                             sigma = rep(NA, 2), rho = NA,
                             lower = c(-pi, -pi, 0.01, 0.01, -25, 0.01, 0.01, -0.99),
                             upper = c(pi, pi, 25, 25, 25, 25, 25, 0.99),
                             maxK = 2, expTrc = 30, ...) {

  # Specfied parameters
  specPars <- c(alpha, mu, sigma, rho)
  indUnSpecPars <- is.na(specPars)

  minusLogLik <- function(pars) {

    # Change only unspecified parameters
    specPars[indUnSpecPars] <- pars
    -logLikWouPairs(x = data, t = delta, alpha = specPars[1:3],
                   mu = specPars[4:5], sigma = specPars[6:7], rho = specPars[8],
                   maxK = maxK, expTrc = expTrc)

  }

  # Optimization
  mleOptimWrapper(minusLogLik = minusLogLik, start = start, lower = lower,
                  upper = upper, ...)

}
