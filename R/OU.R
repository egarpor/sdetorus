

#' @title Simulation of trajectories for the univariate OU diffusion
#'
#' @description Simulation of trajectories of the \emph{univariate} Ornstein-Uhlenbeck (OU) diffusion
#' \deqn{dX_t=\alpha(\mu - X_t)dt+\sigma dW_t, X_0=x_0} using the exact transition probability density.
#'
#' @param x0 initial point.
#' @param alpha strength of the drift.
#' @param mu unconditional mean of the diffusion.
#' @param sigma diffusion coefficient.
#' @param N number of discretization steps in the resulting trajectory.
#' @param delta time discretization step.
#' @return A vector of length \code{N + 1} containing \code{x0} in the first entry and the exact discretized trajectory on the remaining elements.
#' @details The law of the discretized trajectory is a multivariate normal with mean \code{\link{meantOu}} and covariance matrix \code{\link{covstOu}}. See \code{\link{rTrajMou}} for the multivariate case (less efficient for dimension one).
#' @examples
#' \dontrun{
#' require(manipulate)
#' manipulate({
#'  set.seed(345678);
#'  plot(seq(0, N * delta, by = delta), rTrajOu(x0 = 0, alpha = alpha, mu = 0,
#'       sigma = sigma, N = N, delta = delta), ylim = c(-4, 4), type = "l",
#'       ylab = expression(X[t]), xlab = "t")
#'  }, delta = slider(0.01, 5.01, step = 0.1),
#'  N = slider(10, 500, step = 10, initial = 200),
#'  alpha = slider(0.01, 5, step = 0.1, initial = 1),
#'  sigma = slider(0.01, 5, step = 0.1, initial = 1))
#'  }
#' @export
rTrajOu <- function(x0, alpha, mu, sigma, N = 100, delta = 1e-3) {

  # Vector of sampling times
  times <- seq(delta, N * delta, length.out = N)

  # Sample using the exact transition density
  samp <- mvtnorm::rmvnorm(n = 1, mean = meantOu(t = times, alpha = alpha,
                                                 mu = mu, x0 = x0),
                           sigma = covstOu(s = times, t = times, alpha = alpha,
                                           sigma = sigma))

  return(c(x0, drop(samp)))

}


#' @title Transition probability density of the univariate OU diffusion
#'
#' @description Transition probability density of the \emph{univariate} Ornstein-Uhlenbeck (OU) diffusion
#' \deqn{dX_t=\alpha(\mu - X_t)dt+\sigma dW_t, X_0=x_0.}{dX_t=alpha(mu - X_t)dt+sigma dW_t, X0=x0.}
#'
#' @param x vector with the evaluation points.
#' @param t,s time between observations.
#' @param log flag to indicate whether to compute the logarithm of the density.
#' @inheritParams rTrajOu
#' @return A vector of the same length as \code{x} containing the evaluation of the density.
#' @details The transition probability density is a normal density with mean \code{\link{meantOu}} and variance \code{\link{vartOu}}. See \code{\link{dTpdMou}} for the multivariate case (less efficient for dimension one).
#' @examples
#' x <- seq(-4, 4, by = 0.01)
#' plot(x, dTpdOu(x = x, x0 = 3, t = 0.1, alpha = 1, mu = -1, sigma = 1),
#'      type = "l", ylim = c(0, 1.5), xlab = "x", ylab = "Density",
#'      col = rainbow(20)[1])
#' for (i in 2:20) {
#'   lines(x, dTpdOu(x = x, x0 = 3, t = i / 10, alpha = 1, mu = -1, sigma = 1),
#'         col = rainbow(20)[i])
#' }
#' @export
dTpdOu <- function(x, x0, t, alpha, mu, sigma, log = FALSE) {

  dnorm(x, mean = meantOu(t = t, alpha = alpha, mu = mu, x0 = x0),
        sd = sqrt(vartOu(t = t, alpha = alpha, sigma = sigma)), log = log)

}


#' @rdname dTpdOu
#' @export
meantOu <- function(x0, t, alpha, mu) {

  mu + (x0 - mu) * exp(-alpha * t)

}


#' @rdname dTpdOu
#' @export
vartOu <- function(t, alpha, sigma) {

  sigma^2/(2 * alpha) * (1 - exp(-2 * alpha * t))

}


#' @rdname dTpdOu
#' @export
covstOu <- function(s, t, alpha, sigma) {

  sigma^2/(2 * alpha) * outer(s, t, function(ss, tt)
    exp(alpha * (2 * pmin(ss, tt) - (ss + tt))) - exp(-alpha * (ss + tt)))

}


#' @title Maximum likelihood estimation of the OU diffusion
#'
#' @description Computation of the maximum likelihood estimator of the parameters of the \emph{univariate} Ornstein-Uhlenbeck (OU) diffusion from a discretized trajectory \eqn{\{X_{\Delta i}\}_{i=1}^N}{\{X_{Delta * i}\}_{i=1}^N}. The objective function to minimize is
#' \deqn{\sum_{i=2}^n\log p_{\Delta}(X_{\Delta i} | X_{\Delta (i - 1)}).}{\sum_{i=2}^N log p_{Delta}(X_{Delta * i} | X_{Delta * (i - 1)}).}
#'
#' @param data a vector of size \code{N} with the discretized trajectory of the diffusion.
#' @param alpha,mu,sigma arguments to fix a parameter to a given value and perform the estimation on the rest. Defaults to \code{NA}, meaning that the parameter is estimated. Note that \code{start}, \code{lower} and \code{upper} must be changed accordingly if parameters are fixed, see examples.
#' @inheritParams rTrajOu
#' @inheritParams mleOptimWrapper
#' @param ... further arguments to be passed to \code{\link{mleOptimWrapper}}.
#' @return Output from \code{\link{mleOptimWrapper}}.
#' @details The first element in \code{data} is not taken into account for estimation. See \code{\link{mleMou}} for the multivariate case (less efficient for dimension one).
#' @examples
#' set.seed(345678)
#' data <- rTrajOu(x0 = 0, alpha = 1, mu = 0, sigma = 1, N = 100, delta = 0.1)
#' mleOu(data = data, delta = 0.1, start = c(2, 1, 2), lower = c(0.1, -10, 0.1),
#'       upper = c(25, 10, 25))
#'
#' # Fixed sigma and mu
#' mleOu(data = data, delta = 0.1, mu = 0, sigma = 1, start = 2, lower = 0.1,
#'       upper = 25, optMethod = "nlm")
#' @export
mleOu <- function(data, delta, alpha = NA, mu = NA, sigma = NA, start,
                  lower = c(0.01, -5, 0.01), upper = c(25, 5, 25), ...) {

  # Get N
  N <- length(data)

  # Translate data
  y <- data[-1]
  x <- data[-N]

  # Specified parameters
  specPars <- c(alpha, mu, sigma)
  indUnSpecPars <- is.na(specPars)

  minusLogLik <- function(pars) {

    # Change only unspecified parameters
    specPars[indUnSpecPars] <- pars

    -sum(dTpdOu(x = y, x0 = x, t = delta, alpha = specPars[1], mu = specPars[2],
                sigma = specPars[3], log = TRUE))

  }

  # Optimization
  mleOptimWrapper(minusLogLik = minusLogLik, start = start, lower = lower,
                  upper = upper, ...)

}


#' @title Simulation of trajectories for the multivariate OU diffusion
#'
#' @description Simulation of trajectories of the \emph{multivariate} Ornstein-Uhlenbeck (OU) diffusion
#' \deqn{dX_t=A(\mu - X_t)dt+\Sigma^\frac{1}{2}dW_t, X_0=x_0}{dX_t=A(mu - X_t)dt+Sigma^(1/2) dW_t, X0=x0} using the exact transition probability density.
#'
#' @param x0 a vector of length \code{p} containing initial point.
#' @param A the drift matrix, of size \code{c(p, p)}.
#' @param mu unconditional mean of the diffusion, a vector of length \code{p}.
#' @param Sigma square of the diffusion matrix, a matrix of size \code{c(p, p)}.
#' @inheritParams rTrajOu
#' @return A matrix of size \code{c(N + 1, p)} containing \code{x0} in the first row and  the exact discretized trajectory on the remaining rows.
#' @details The law of the discretized trajectory at \emph{each} time step is a multivariate normal with mean \code{\link{meantMou}} and covariance matrix \code{\link{covtMou}}. See \code{\link{rTrajOu}} for the univariate case (more efficient).
#'
#' \code{solve(A) \%*\% Sigma} has to be a covariance matrix (symmetric and positive definite) in order to have a proper transition density. For the bivariate case, this can be ensured with the \code{\link{alphaToA}} function. In the multivariate case, it is ensured if \code{Sigma} is isotropic and \code{A} is a covariance matrix.
#' @examples
#' set.seed(987658)
#' data <- rTrajMou(x0 = c(0, 0), A = alphaToA(alpha = c(1, 2, 0.5),
#'                  sigma = 1:2), mu = c(1, 1), Sigma = diag((1:2)^2),
#'                  N = 200, delta = 0.1)
#' plot(data, pch = 19, col = rainbow(201), cex = 0.25)
#' arrows(x0 = data[-201, 1], y0 = data[-201, 2], x1 = data[-1, 1],
#'        y1 = data[-1, 2], col = rainbow(201), angle = 10, length = 0.1)
#' @export
rTrajMou <- function(x0, A, mu, Sigma, N = 100, delta = 1e-3) {

  # Get p
  p <- ncol(A)

  # Eigendecomposition of A
  eigA <- eigen(A)

  # Covariance matrix
  covt <- covtMou(t = delta, Sigma = Sigma, eigA = eigA)

  # Sample using the exact transition density
  samp <- matrix(x0, nrow = N + 1, ncol = p, byrow = TRUE)
  sapply(2:(N + 1), function(i) {
    samp[i, ] <<-
      mvtnorm::rmvnorm(n = 1, mean = meantMou(t = delta, x0 = samp[i - 1, ],
                                              mu = mu, eigA = eigA),
                       sigma = covt)
  })

  return(samp)

}


#' @title Transition probability density of the multivariate OU diffusion
#'
#' @description Transition probability density of the \emph{multivariate} Ornstein-Uhlenbeck (OU) diffusion
#' \deqn{dX_t=A(\mu - X_t)dt+\Sigma^\frac{1}{2}dW_t, X_0=x_0.}{dX_t=A(mu - X_t)dt+Sigma^(1/2) dW_t, X0=x0.}
#'
#' @param x matrix of with \code{p} columns containing the evaluation points.
#' @inheritParams dTpdOu
#' @inheritParams rTrajMou
#' @param eigA optional argument containing \code{eigen(A)} for reuse.
#' @return A matrix of the same size as \code{x} containing the evaluation of the density.
#' @details The transition probability density is a multivariate normal with mean \code{\link{meantMou}} and covariance \code{\link{covtMou}}. See \code{\link{dTpdOu}} for the univariate case (more efficient).
#'
#' \code{solve(A) \%*\% Sigma} has to be a covariance matrix (symmetric and positive definite) in order to have a proper transition density. For the bivariate case, this can be ensured with the \code{\link{alphaToA}} function. In the multivariate case, it is ensured if \code{Sigma} is isotropic and \code{A} is a covariance matrix.
#' @examples
#' x <- seq(-4, 4, by = 0.1)
#' xx <- as.matrix(expand.grid(x, x))
#' \dontrun{
#' require(manipulate)
#' manipulate(image(x, x, matrix(dTpdMou(x = xx, x0 = c(1, 2), t = t,
#'                                       A = alphaToA(alpha = c(1, 2, 0.5),
#'                                                    sigma = 1:2),
#'                                       mu = c(0, 0), Sigma = diag((1:2)^2)),
#'                               nrow = length(x), ncol = length(x)),
#'                  zlim = c(0, 0.25)), t = slider(0.1, 5, step = 0.1))
#' }
#' @export
dTpdMou <- function(x, x0, t, A, mu, Sigma, eigA, log = FALSE) {

  # Eigendecomposition
  if (missing(eigA)) {

    eigA <- eigen(A)

  }

  # Density
  mvtnorm::dmvnorm(x = x, mean = meantMou(t = t, x0 = x0, mu = mu, eigA = eigA),
                   sigma = covtMou(t = t, Sigma = Sigma, eigA = eigA),
                   log = log)

}


#' @rdname dTpdMou
#' @export
meantMou <- function(t, x0, A, mu, eigA) {

  # Dimension
  p <- length(mu)

  # Eigendecompositon
  if (missing(eigA)) {

    eigA <- eigen(A)

  }

  # Mean
  mu + (eigA$vectors %*% diag(exp(-eigA$values * t), nrow = p, ncol = p) %*%
          solve(eigA$vectors)) %*% (x0 - mu)

}


#' @rdname dTpdMou
#' @export
covtMou <- function(t, A, Sigma, eigA) {

  # Eigendecomposition
  if (missing(eigA)) {

    eigA <- eigen(A)

  }

  # Eigenvalues and eigenvectors
  lambda <- eigA$values
  v <- eigA$vectors
  sv <- solve(v)

  # Dimension
  p <- ncol(v)

  # Covariance matrix
  0.5 * v %*% diag((1- exp(-2 * t * lambda)) / lambda, nrow = p, ncol = p) %*%
    sv %*% Sigma

}


#' @title Maximum likelihood estimation of the multivariate OU diffusion
#'
#' @description Computation of the maximum likelihood estimator of the parameters of the \emph{multivariate} Ornstein-Uhlenbeck (OU) diffusion from a discretized trajectory \eqn{\{X_{\Delta i}\}_{i=1}^N}{\{X_{Delta * i}\}_{i=1}^N}. The objective function to minimize is
#' \deqn{\sum_{i=2}^n\log p_{\Delta}(X_{\Delta i} | X_{\Delta (i - 1)}).}{\sum_{i=2}^N log p_{Delta}(X_{Delta * i} | X_{Delta * (i - 1)}).}
#'
#' @param data a matrix of size \code{c(N, p)} with the discretized trajectory of the diffusion.
#' @param alpha,mu,sigma arguments to fix a parameter to a given value and perform the estimation on the rest. Defaults to \code{NA}, meaning that the parameter is estimated. Note that \code{start}, \code{lower} and \code{upper} must be changed accordingly if parameters are fixed, see examples.
#' @inheritParams rTrajOu
#' @inheritParams mleOptimWrapper
#' @param ... further arguments to be passed to \code{\link{mleOptimWrapper}}.
#' @return Output from \code{\link{mleOptimWrapper}}.
#' @details The first row in \code{data} is not taken into account for estimation. See \code{\link{mleOu}} for the univariate case (more efficient).
#'
#' \code{mleMou} only handles \code{p = 2} currently. It imposes that \code{Sigma} is diagonal and handles the parametrization of \code{A} by \code{\link{alphaToA}}.
#' @examples
#' set.seed(345678)
#' data <- rTrajMou(x0 = c(0, 0), A = alphaToA(alpha = c(1, 1, 0.5),
#'                                             sigma = 1:2), mu = c(1, 1),
#'                  Sigma = diag((1:2)^2), N = 200, delta = 0.5)
#' mleMou(data = data, delta = 0.5, start = c(1, 1, 0, 1, 1, 1, 2),
#'        lower = c(0.1, 0.1, -25, -10, -10, 0.1, 0.1),
#'        upper = c(25, 25, 25, 10, 10, 25, 25), maxit = 500)
#'
#' # Fixed sigma and mu
#' mleMou(data = data, delta = 0.5, mu = c(1, 1), sigma = 1:2, start = c(1, 1, 0),
#'        lower = c(0.1, 0.1, -25), upper = c(25, 25, 25))
#' @export
mleMou <- function(data, delta, alpha = rep(NA, 3), mu = rep(NA, 2),
                   sigma = rep(NA, 2), start,
                   lower = c(0.01, 0.01, -25, -pi, -pi, 0.01, 0.01),
                   upper = c(25, 25, 25, pi, pi, 25, 25), ...) {

  # Get N and p
  N <- nrow(data)
  p <- ncol(data)
  if (p > 2) {

    stop("mleMou only handles p = 2 currently\n")

  }

  # Specified parameters
  specPars <- c(alpha, mu, sigma)
  indUnSpecPars <- is.na(specPars)

  minusLogLik <- function(pars) {

    # Change only unspecified parameters
    specPars[indUnSpecPars] <- pars

    # Eigendecomposition
    A <- alphaToA(alpha = specPars[1:3], sigma = specPars[6:7])
    eigA <- eigen(A)

    # Covariance
    covt <- covtMou(t = delta, Sigma = diag(specPars[6:7]^2), eigA = eigA)

    # Exponential
    exptA <- eigA$vectors %*% diag(exp(-eigA$values * delta),
                                   nrow = p, ncol = p) %*% solve(eigA$vectors)

    # x - mut
    xMut <- sweep(data[-1, ] - sweep(data[-N, ], 2, specPars[4:5], "-") %*%
                    t(exptA), 2, specPars[4:5], "-")

    # Loglikelihood
    -sum(mvtnorm::dmvnorm(x = xMut, mean = rep(0, p), sigma = covt, log = TRUE))

  }

  region <- function(pars) {

    # Check positvedefiniteness
    testPosDef <- pars[1] * pars[2] - pars[3]^2

    if (testPosDef < 0) {

      pars[3] <- sign(pars[3]) * sqrt(pars[1] * pars[2]) * 0.99
      penalty <- 1e3 * testPosDef - 10

    } else {

      penalty <- 0

    }

    return(list("pars" = pars, "penalty" = -penalty))

  }

  # Optimization
  mleOptimWrapper(minusLogLik = minusLogLik, start = start, region = region,
                  lower = lower, upper = upper, ...)

}


#' @title Valid drift matrices for the Ornstein-Uhlenbeck diffusion in 2D
#'
#' @description Constructs drift matrices \eqn{A} such that \code{solve(A) \%*\% Sigma} is symmetric.
#'
#' @param A matrix of size \code{c(2, 2)}.
#' @param alpha vector of length \code{3} containing the \code{A} matrix. The first two elements are the diagonal.
#' @param sigma vector of length \code{2} containing the \strong{square root} of the diagonal of \code{Sigma}.
#' @param rho correlation of \code{Sigma}.
#' @param Sigma the diffusion matrix of size \code{c(2, 2)}.
#' @return The drift matrix \code{A} or the \code{alpha} vector.
#' @details The parametrization enforces that \code{solve(A) \%*\% Sigma} is symmetric. Positive definiteness is guaranteed if \code{alpha[3]^2 < rho^2 * (alpha[1] - alpha[2])^2 / 4 + alpha[1] * alpha[2]}.
#' @examples
#' # Parameters
#' alpha <- 3:1
#' Sigma <- rbind(c(1, 0.5), c(0.5, 4))
#'
#' # Covariance matrix
#' A <- alphaToA(alpha = alpha, Sigma = Sigma)
#' S <- 0.5 * solve(A) %*% Sigma
#' det(S)
#'
#' # Check
#' aToAlpha(A = alphaToA(alpha = alpha, Sigma = Sigma), Sigma = Sigma)
#' alphaToA(alpha = aToAlpha(A = A, Sigma = Sigma), Sigma = Sigma)
#' @export
alphaToA <- function(alpha, sigma, rho = 0, Sigma) {

  # Compute sigma and rho if not provided
  if (missing(sigma)) {

    if (missing(Sigma)) {

      stop("Must provide either sigma and rho, or Sigma")

    } else {

      sigma <- diag(sqrt(Sigma))
      rho <- Sigma[1, 2] / prod(sigma)

    }

  }

  # Compute A
  quo <- sigma[1] / sigma[2]
  add <- 0.5 * rho * (alpha[2] - alpha[1])
  A <- matrix(c(alpha[1], (alpha[3] + add) * quo,
                (alpha[3] - add)/ quo, alpha[2]),
              nrow = 2, ncol = 2, byrow = TRUE)

  return(A)

}


#' @rdname alphaToA
#' @export
aToAlpha <- function(A, sigma, rho = 0, Sigma) {

  # Compute sigma and rho if not provided
  if (missing(sigma)) {

    if (missing(Sigma)) {

      stop("Must provide either sigma and rho, or Sigma")

    } else {

      sigma <- sqrt(diag(Sigma))
      rho <- Sigma[1, 2] / prod(sigma)

    }

  }

  # Compute alpha
  quo <- sigma[2] / sigma[1]
  add <- 0.5 * rho * (A[2, 2] - A[1, 1])
  alpha <- c(diag(A), A[1, 2] * quo - add)

  return(alpha)

}
