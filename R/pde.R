

#' @title Transition probability density in 1D by PDE solving
#'
#' @description Computation of the transition probability density (tpd) of the Wrapped Normal (WN) or von Mises (vM) diffusion, by solving its associated Fokker-Planck Partial Differential Equation (PDE) in 1D.
#'
#' @inheritParams crankNicolson1D
#' @param x0 point giving the mean of the initial circular density, a WN with standard deviation equal to \code{sdInitial}.
#' @param t time separating \code{x0} and the evaluation of the tpd.
#' @param type either \code{"WN"} or \code{"vM"}.
#' @inheritParams dTpdWou1D
#' @param Mt size of the time grid in \eqn{[0, t]}.
#' @param sdInitial the standard deviation of the concentrated WN giving the initial condition.
#' @param ... Further parameters passed to \code{\link{crankNicolson1D}}.
#' @return A vector of length \code{Mx} with the tpd evaluated at \code{seq(-pi, pi, l = Mx + 1)[-(Mx + 1)]}.
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
#' @details A combination of small \code{sdInitial} and coarse space-time discretization (small \code{Mx} and \code{Mt}) is prone to create numerical instabilities. See Sections 3.4.1, 2.2.1 and 2.2.2 in García-Portugués et al. (2017) for details.
#' @references 
#' García-Portugués, E., Sorensen, M., Mardia, K. V. and Hamelryck, T. (2017) Langevin diffusions on the torus: estimation and applications. \emph{Stat. Comput.}, \url{https://doi.org/10.1007/s11222-017-9790-2}.
#' @examples
#' Mx <- 100
#' x <- seq(-pi, pi, l = Mx + 1)[-c(Mx + 1)]
#' x0 <- pi
#' t <- 0.5
#' alpha <- 1
#' mu <- 0
#' sigma <- 1
#' \dontrun{
#' require(manipulate)
#' manipulate({
#' plot(x, dTpdPde1D(Mx = Mx, x0 = x0, t = t, alpha = alpha, mu = 0,
#'                   sigma = sigma), type = "l", ylab = "Density",
#'      xlab = "", ylim = c(0, 0.75))
#' lines(x, dTpdWou1D(x = x, x0 = rep(x0, Mx), t = t, alpha = alpha, mu = 0,
#'                     sigma = sigma), col = 2)
#' }, x0 = slider(-pi, pi, step = 0.01, initial = 0),
#' alpha = slider(0.01, 5, step = 0.01, initial = 1),
#' sigma = slider(0.01, 5, step = 0.01, initial = 1),
#' t = slider(0.01, 5, step = 0.01, initial = 1))
#' }
#' @export
dTpdPde1D <- function(Mx = 500, x0, t, alpha, mu, sigma, type = "WN",
                      Mt = ceiling(1e2 * t), sdInitial = 0.1, ...) {

  # Initialize list of parameters
  grid <- seq(-pi, pi, l = Mx + 1)[-c(Mx + 1)]

  # Initial conditions
  u0 <- matrix(dWn1D(x = grid, mu = x0, sigma = sdInitial), ncol = 1)

  # Drift and diffusion
  if (type == "WN") {

    bGrid <- driftWn1D(x = grid, alpha = alpha, mu = mu, sigma = sigma)
    sigma2Grid <- rep(sigma^2, Mx)

  } else if (type == "vM") {

    bGrid <- alpha * sin(mu - grid)
    sigma2Grid <- rep(sigma^2, Mx)

  }

  # Solution Tpd
  crankNicolson1D(u0 = u0, b = bGrid, sigma2 = sigma2Grid, Mx = Mx,
                  deltax = grid[2] - grid[1], N = Mt, deltat = t / Mt, ...)

}


#' @title Transition probability density in 2D by PDE solving
#'
#' @description Computation of the transition probability density (tpd) of the Wrapped Normal (WN) or Multivariate von Mises (MvM) diffusion, by solving its associated Fokker-Planck Partial Differential Equation (PDE) in 2D.
#'
#' @inheritParams crankNicolson2D
#' @inheritParams dTpdPde1D
#' @param x0 point giving the mean of the initial circular density, an isotropic WN with standard deviations equal to \code{sdInitial}.
#' @param alpha for \code{"WN"}, a vector of length \code{3} parametrizing the \code{A} matrix as in \code{\link{alphaToA}}. For \code{"vM"}, a vector of length \code{3} containing \code{c(alpha[1:2], A[1, 2])}, from the arguments \code{alpha} and \code{A} in \code{\link{driftMvm}}.
#' @param mu vector of length \code{2} giving the mean.
#' @param sigma for \code{"WN"}, a vector of length \code{2} containing the \strong{square root} of the diagonal of the diffusion matrix. For \code{"vM"}, the standard deviation giving the isotropic diffusion matrix.
#' @param rho for \code{"WN"}, the correlation of the diffusion matrix.
#' @param sdInitial standard deviations of the concentrated WN giving the initial condition.
#' @param ... Further parameters passed to \code{\link{crankNicolson2D}}.
#' @return A matrix of size \code{c(Mx, My)} with the tpd evaluated at the combinations of \code{seq(-pi, pi, l = Mx + 1)[-(Mx + 1)]} and \code{seq(-pi, pi, l = My + 1)[-(My + 1)]}.
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
#' @details A combination of small \code{sdInitial} and coarse space-time discretization (small \code{Mx} and \code{Mt}) is prone to create numerical instabilities. See Sections 3.4.2, 2.2.1 and 2.2.2 in García-Portugués et al. (2017) for details.
#' @references 
#' García-Portugués, E., Sorensen, M., Mardia, K. V. and Hamelryck, T. (2017) Langevin diffusions on the torus: estimation and applications. \emph{Stat. Comput.}, \url{https://doi.org/10.1007/s11222-017-9790-2}.
#' @examples
#' M <- 100
#' x <- seq(-pi, pi, l = M + 1)[-c(M + 1)]
#' image(x, x, dTpdPde2D(Mx = M, My = M, x0 = c(0, pi), t = 1,
#'                       alpha = c(1, 1, 0.5), mu = c(pi / 2, 0), sigma = 1:2),
#'       zlim = c(0, 0.25), col = colorRamps::matlab.like(20),
#'       xlab = "x", ylab = "y")
#' @export
dTpdPde2D <- function(Mx = 50, My = 50, x0, t, alpha, mu, sigma, rho = 0, 
                      type = "WN", Mt = ceiling(1e2 * t), sdInitial = 0.1, 
                      ...) {

  # Drifts
  gridX <- seq(-pi, pi, l = Mx + 1)[-c(Mx + 1)]
  gridY <- seq(-pi, pi, l = My + 1)[-c(My + 1)]
  grid <- as.matrix(expand.grid(gridX, gridY))

  # Drift and diffusion
  if (type == "WN") {

    bGrid <- driftWn2D(x = grid, A = alphaToA(alpha = alpha, sigma = sigma),
                       mu = mu, sigma = sigma, rho = rho)
    bxGrid <- matrix(bGrid[, 1], nrow = Mx, ncol = My)
    byGrid <- matrix(bGrid[, 2], nrow = Mx, ncol = My)
    sigma2xGrid <- matrix(sigma[1]^2, nrow = Mx, ncol = My)
    sigma2yGrid <- matrix(sigma[2]^2, nrow = Mx, ncol = My)
    sigmaxyGrid <- matrix(0, nrow = Mx, ncol = My)

  } else if (type == "vM") {

    bxGrid <- matrix(alpha[1] * sin(mu[1] - grid[, 1]) -
                       alpha[3] * cos(mu[1] - grid[, 1]) *
                       sin(mu[2] - grid[, 2]), nrow = Mx, ncol = My)
    byGrid <- matrix(alpha[2] * sin(mu[2] - grid[, 2]) -
                       alpha[3] * sin(mu[1] - grid[, 1]) *
                       cos(mu[2] - grid[, 2]), nrow = Mx, ncol = My)
    sigma2xGrid <- matrix(sigma[1]^2, nrow = Mx, ncol = My)
    sigma2yGrid <- sigma2xGrid
    sigmaxyGrid <- matrix(0, nrow = Mx, ncol = My)

  }

  # Initial condition
  u0 <- matrix(dWn1D(x = gridX, mu = x0[1], sigma = sdInitial) %*%
                 t(dWn1D(x = gridY, mu = x0[2], sigma = sdInitial)), ncol = 1)

  # Go only throgout the unique indices (tpd[x1, x0])
  matrix(crankNicolson2D(u0 = u0, bx = bxGrid, by = byGrid,
                         sigma2x = sigma2xGrid, sigma2y = sigma2yGrid,
                         sigmaxy = sigmaxyGrid, N = Mt, deltat = t / Mt,
                         Mx = Mx, deltax = gridX[2] - gridX[1], My = My,
                         deltay = gridY[2] - gridY[1], ...), nrow = Mx,
         ncol = My)

}


#' @title MLE for toroidal process via PDE solving in 1D
#'
#' @description Maximum Likelihood Estimation (MLE) for arbitrary diffusions, based on the transition probability density (tpd) obtained as the numerical solution of the Fokker-Planck Partial Differential Equation (PDE) in 1D.
#'
#' @param b drift function. Must return a vector of the same size as its argument.
#' @param sigma2 function giving the squared diffusion coefficient. Must return a vector of the same size as its argument.
#' @inheritParams dTpdPde1D
#' @inheritParams mleOu
#' @param linearBinning flag to indicate whether linear binning should be applied for the initial values of the tpd, instead of usual simple binning (cheaper). Linear binning is always done in the evaluation of the tpd.
#' @return Output from \code{\link{mleOptimWrapper}}.
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
#' @details See Sections 3.4.1 and 3.4.4 in García-Portugués et al. (2017) for details.
#' @references 
#' García-Portugués, E., Sorensen, M., Mardia, K. V. and Hamelryck, T. (2017) Langevin diffusions on the torus: estimation and applications. \emph{Stat. Comput.}, \url{https://doi.org/10.1007/s11222-017-9790-2}.
#' @examples
#' \dontrun{
#' # Test in OU
#' alpha <- 2
#' mu <- 0
#' sigma <- 1
#' set.seed(234567)
#' traj <- rTrajOu(x0 = 0, alpha = alpha, mu = mu, sigma = sigma, N = 500,
#'                 delta = 0.5)
#' b <- function(x, pars) pars[1] * (pars[2] - x)
#' sigma2 <- function(x, pars) rep(pars[3]^2, length(x))
#'
#' exactOu <- mleOu(traj, delta = 0.5, start = c(1, 1, 2),
#'                  lower = c(0.1, -pi, 0.1), upper = c(10, pi, 10))
#' pdeOu <- mlePde1D(data = traj, delta = 0.5, Mx = 100, Mt = 100, b = b,
#'                   sigma2 = sigma2, start = c(1, 1, 2),
#'                   lower = c(0.1, -pi, -10), upper = c(10, pi, 10),
#'                   verbose = 2)
#' pdeOuLin <- mlePde1D(data = traj, delta = 0.5, Mx = 100, Mt = 100, b = b,
#'                      sigma2 = sigma2, start = c(1, 1, 2),
#'                      lower = c(0.1, -pi, -10), upper = c(10, pi, 10),
#'                      linearBinning = TRUE, verbose = 2)
#' head(exactOu)
#' head(pdeOu)
#' head(pdeOuLin)
#'
#' # Test in WN diffusion
#' alpha <- 2
#' mu <- 0
#' sigma <- 1
#' set.seed(234567)
#' traj <- rTrajWn1D(x0 = 0, alpha = alpha, mu = mu, sigma = sigma, N = 500,
#'                  delta = 0.5)
#'
#' exactOu <- mleOu(traj, delta = 0.5, start = c(1, 1, 2),
#'                  lower = c(0.1, -pi, 0.1), upper = c(10, pi, 10))
#' pdeWn <- mlePde1D(data = traj, delta = 0.5, Mx = 100, Mt = 100, b = function(x, pars)
#'                   driftWn1D(x = x, alpha = pars[1], mu = pars[2], sigma = pars[3]),
#'                   sigma2 = function(x, pars) rep(pars[3]^2, length(x)),
#'                   start = c(1, 1, 2), lower = c(0.1, -pi, -10),
#'                   upper = c(10, pi, 10), verbose = 2)
#' pdeWnLin <- mlePde1D(data = traj, delta = 0.5, Mx = 100, Mt = 100, b = function(x, pars)
#'                      driftWn1D(x = x, alpha = pars[1], mu = pars[2], sigma = pars[3]),
#'                      sigma2 = function(x, pars) rep(pars[3]^2, length(x)),
#'                      start = c(1, 1, 2), lower = c(0.1, -pi, -10),
#'                      upper = c(10, pi, 10), linearBinning = TRUE, verbose = 2)
#' head(exactOu)
#' head(pdeWn)
#' head(pdeWnLin)
#' }
#' @export
mlePde1D <- function(data, delta, b, sigma2, Mx = 500,
                     Mt = ceiling(1e2 * delta), sdInitial = 0.1,
                     linearBinning = FALSE, start, lower, upper, ...) {

  # Size
  N <- length(data)

  # Number of parameters
  npar <- ifelse(is.matrix(start), ncol(start), length(start))

  # Circular grid for solution
  grid <- seq(-pi, pi, l = Mx + 1)[-c(Mx + 1)]

  # State step
  delx <- grid[2] - grid[1]

  # Time step
  delt <- delta / Mt

  # Lower and upper indexes for linear binning of the evaluation points
  # grid[lowerIndBins] <= data <= grid[upperIndBins] (safer using findInterval
  # than ceiling((data - pi) / delx))
  lowerIndBins <- findInterval(x = data, vec = grid)
  upperIndBins <- toInt(lowerIndBins + 1, a = 1, b = Mx + 1)

  # Evaluation indexes (indEval for x[i])
  indEval <- cbind(lowerIndBins[-1], upperIndBins[-1])

  # Evaluation weights for linear binning
  wEval <- weightsLinearInterp1D(x = data[-1], g1 = grid[indEval[, 1]],
                                 g2 = grid[indEval[, 2]], circular = TRUE)

  if (linearBinning) {

    # Unique indexes for both lower and upper bins
    uniqueInd <- unique(c(lowerIndBins, upperIndBins))

    # Link from lower and upper bins to the aggregated (lower + upper)
    # unique values
    matchLowerIndBins <- match(x = lowerIndBins, table = uniqueInd)
    matchUpperIndBins <- match(x = upperIndBins, table = uniqueInd)

    # Start indexes
    indStart <- cbind(matchLowerIndBins[-N], matchUpperIndBins[-N])

    # Start weights
    wStart <- weightsLinearInterp1D(x = data[-N], g1 = grid[lowerIndBins[-N]],
                                    g2 = grid[upperIndBins[-N]],
                                    circular = TRUE)

    # Initial conditions (normalized to ensure they are a density in the grid)
    u0 <- sapply(uniqueInd, function(i) dWn1D(x = grid, mu = grid[i],
                                              sigma = sdInitial))
    u0 <- sweep(u0, 2, apply(u0, 2, periodicTrapRule1D), FUN = "/")

    # Loglikelihood (expensive part)
    minusLogLik <- function(x) {

      # Update drifts and squared diffusion
      bGrid <- b(x = grid, pars = x)
      sigma2Grid <- sigma2(x = grid, pars = x)

      # Go only throgout the unique indices (tpd[xEval, xStart])
      tpd <- crankNicolson1D(u0 = u0, b = bGrid, sigma2 = sigma2Grid, N = Mt,
                             deltat = delt, Mx = Mx, deltax = delx,
                             imposePositive = TRUE)

      # Linear binning in evaluation and start
      res <- -sum(sapply(1:(N - 1), function(i)
        log(max(crossprod(wEval[i, ], tpd[indEval[i, ], indStart[i, ]])
                %*% wStart[i, ], 1e-08))))

      return(res)

    }

  } else {

    # closestIndBins <- sapply(seq_along(data), function(i)
    #   which.min(abs(toPiInt(data[i] - grid))))
    closestIndBins <- toInt(round((data + pi) / delx) + 1, a = 1, b = Mx + 1)

    # Unique indexes of closests bins
    simpleUniqueInd <- unique(closestIndBins)

    # Link from closest bin to the unique values
    matchClosestIndBins <- match(x = closestIndBins, table = simpleUniqueInd)

    # Start indexes
    indStart <- matchClosestIndBins[-N]

    # Initial conditions (normalized to ensure they are a density in the grid)
    u0 <- sapply(simpleUniqueInd, function(i) dWn1D(x = grid, mu = grid[i],
                                                    sigma = sdInitial))
    u0 <- sweep(u0, 2, apply(u0, 2, periodicTrapRule1D), FUN = "/")

    # Loglikelihood (expensive part)
    minusLogLik <- function(x) {

      # Update drifts and squared diffusion
      bGrid <- b(x = grid, pars = x)
      sigma2Grid <- sigma2(x = grid, pars = x)

      # Go only throgout the unique indices (tpd[xEval, xStart])
      tpd <- crankNicolson1D(u0 = u0, b = bGrid, sigma2 = sigma2Grid, N = Mt,
                             deltat = delt, Mx = Mx, deltax = delx)

      # Linear binning in evaluation
      res <- -sum(sapply(1:(N - 1), function(i)
        log(max(crossprod(wEval[i, ], tpd[indEval[i, ], indStart[i]]), 1e-08))))

      return(res)

    }

  }

  # Optimization
  mleOptimWrapper(minusLogLik = minusLogLik, start = start, lower = lower,
                  upper = upper, ...)

}


#' @title MLE for toroidal process via PDE solving in 2D
#'
#' @description Maximum Likelihood Estimation (MLE) for arbitrary diffusions, based on the transition probability density (tpd) obtained as the numerical solution of the Fokker-Planck Partial Differential Equation (PDE) in 2D.
#'
#' @param b drift function. Must return a vector of the same size as its argument.
#' @param sigma2 function giving the diagonal of the diffusion matrix. Must return a vector of the same size as its argument.
#' @inheritParams psMle
#' @inheritParams dTpdPde2D
#' @inheritParams mlePde1D
#' @return Output from \code{\link{mleOptimWrapper}}.
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
#' @details See Sections 3.4.2 and 3.4.4 in García-Portugués et al. (2017) for details. The function currently includes the \code{region} function for imposing a feasibility region on the parameters of the bivariate WN diffusion.
#' @references 
#' García-Portugués, E., Sorensen, M., Mardia, K. V. and Hamelryck, T. (2017) Langevin diffusions on the torus: estimation and applications. \emph{Stat. Comput.}, \url{https://doi.org/10.1007/s11222-017-9790-2}.
#' @examples
#' \dontrun{
#' # Test in OU process
#' alpha <- c(1, 2, -0.5)
#' mu <- c(0, 0)
#' sigma <- c(0.5, 1)
#' set.seed(2334567)
#' data <- rTrajMou(x0 = c(0, 0), A = alphaToA(alpha = alpha, sigma = sigma),
#'                  mu = mu, Sigma = diag(sigma^2), N = 500, delta = 0.5)
#' b <- function(x, pars) sweep(-x, 2, pars[4:5], "+") %*%
#'                        t(alphaToA(alpha = pars[1:3], sigma = sigma))
#' sigma2 <- function(x, pars) repRow(sigma^2, nrow(x))
#'
#' exactOu <- mleMou(data = data, delta = 0.5, sigma = sigma,
#'                   start = c(1, 1, 0, 2, 2),
#'                   lower = c(0.1, 0.1, -25, -10, -10),
#'                   upper = c(25, 25, 25, 10, 10))
#' head(exactOu, 2)
#' pdeOu <- mlePde2D(data = data, delta = 0.5, b = b, sigma2 = sigma2,
#'                   Mx = 10, My = 10, Mt = 10, start = rbind(c(1, 1, 0, 2, 2)),
#'                   lower = c(0.1, 0.1, -25, -10, -10),
#'                   upper = c(25, 25, 25, 10, 10), verbose = 2)
#' head(pdeOu, 2)
#' pdeOuLin <- mlePde2D(data = data, delta = 0.5, b = b, sigma2 = sigma2,
#'                      Mx = 10, My = 10, Mt = 10, start = rbind(c(1, 1, 0, 2, 2)),
#'                      lower = c(0.1, 0.1, -25, -10, -10),
#'                      upper = c(25, 25, 25, 10, 10), verbose = 2,
#'                      linearBinning = TRUE)
#' head(pdeOuLin, 2)
#'
#' # Test in WN diffusion
#' alpha <- c(1, 0.5, 0.25)
#' mu <- c(0, 0)
#' sigma <- c(2, 1)
#' set.seed(234567)
#' data <- rTrajWn2D(x0 = c(0, 0), alpha = alpha, mu = mu, sigma = sigma,
#'                     N = 200, delta = 0.5)
#' b <- function(x, pars) driftWn2D(x = x, A = alphaToA(alpha = pars[1:3],
#'                                                      sigma = sigma),
#'                                  mu = pars[4:5], sigma = sigma)
#' sigma2 <- function(x, pars) repRow(sigma^2, nrow(x))
#'
#' exactOu <- mleMou(data = data, delta = 0.5, sigma = sigma,
#'                   start = c(1, 1, 0, 1, 1),
#'                   lower = c(0.1, 0.1, -25, -25, -25),
#'                   upper = c(25, 25, 25, 25, 25), optMethod = "nlm")
#' pdeWn <- mlePde2D(data = data, delta = 0.5, b = b, sigma2 = sigma2,
#'                   Mx = 20, My = 20, Mt = 10, start = rbind(c(1, 1, 0, 1, 1)),
#'                   lower = c(0.1, 0.1, -25, -25, -25),
#'                   upper = c(25, 25, 25, 25, 25), verbose = 2, optMethod = "nlm")
#' pdeWnLin <- mlePde2D(data = data, delta = 0.5, b = b, sigma2 = sigma2,
#'                      Mx = 20, My = 20, Mt = 10, start = rbind(c(1, 1, 0, 1, 1)),
#'                      lower = c(0.1, 0.1, -25, -25, -25),
#'                      upper = c(25, 25, 25, 25, 25), verbose = 2,
#'                      linearBinning = TRUE)
#'
#' head(exactOu)
#' head(pdeOu)
#' head(pdeOuLin)
#' }
#' @export
mlePde2D <- function(data, delta, b, sigma2, Mx = 50, My = 50,
                     Mt = ceiling(1e2 * delta), sdInitial = 0.1,
                     linearBinning = FALSE, start, lower, upper, ...) {

  # Get N
  if (!is.matrix(data)) {

    stop("data must be a matrix")

  }
  N <- nrow(data)

  # Number of parameters
  npar <- ifelse(is.matrix(start), ncol(start), length(start))

  # Circular grids
  gridX <- seq(-pi, pi, l = Mx + 1)[-c(Mx + 1)]
  gridY <- seq(-pi, pi, l = My + 1)[-c(My + 1)]
  gridXY <- as.matrix(expand.grid(gridX, gridY))

  # State steps
  delx <- gridX[2] - gridX[1]
  dely <- gridY[2] - gridY[1]

  # Time step
  delt <- delta / Mt

  # Lower and upper indexes for linear binning of the evaluation points
  # gridX[lowerIndBinsX] <= data[, 1] <= gridX[upperIndBinsX]
  # gridY[lowerIndBinsY] <= data[, 2] <= gridY[upperIndBinsY]
  lowerIndBinsX <- findInterval(x = data[, 1], vec = gridX)
  upperIndBinsX <- toInt(lowerIndBinsX + 1, a = 1, b = Mx + 1)
  lowerIndBinsY <- findInterval(x = data[, 2], vec = gridY)
  upperIndBinsY <- toInt(lowerIndBinsY + 1, a = 1, b = My + 1)

  # Evaluation indexes (indEval for x[i]), lower-lower, upper-lower,
  # lower-upper, upper-upper
  indEval <- matrix(kIndex(i = rep(c(lowerIndBinsX[-1], upperIndBinsX[-1]), 2),
                           j = c(rep(lowerIndBinsY[-1], 2),
                                 rep(upperIndBinsY[-1], 2)),
                           nr = Mx, nc = My, byRows = TRUE),
                    ncol = 4, nrow = N - 1)

  # Evaluation weights for linear binning, sorted as lower-lower, upper-lower,
  # lower-upper, upper-upper
  wEval <- weightsLinearInterp2D(x = data[-1, 1], y = data[-1, 2],
                                 gx1 = gridX[lowerIndBinsX[-1]],
                                 gx2 = gridX[upperIndBinsX[-1]],
                                 gy1 = gridY[lowerIndBinsY[-1]],
                                 gy2 = gridY[upperIndBinsY[-1]],
                                 circular = TRUE)

  if (linearBinning) {

    # Unique indexes for both (lower, lower), (lower, upper), (upper, lower)
    # and (upper, upper) bins
    uniqueInd <- unique(cbind(rep(c(lowerIndBinsX, upperIndBinsX), 2),
                              c(rep(lowerIndBinsY, 2), rep(upperIndBinsY, 2))))

    # Link from lower and upper bins to the aggregated (lower + upper)
    # unique values
    matchLowerLowerIndBins <- matMatch(x = cbind(lowerIndBinsX, lowerIndBinsY),
                                       mat = uniqueInd)
    matchUpperLowerIndBins <- matMatch(x = cbind(upperIndBinsX, lowerIndBinsY),
                                       mat = uniqueInd)
    matchLowerUpperIndBins <- matMatch(x = cbind(lowerIndBinsX, upperIndBinsY),
                                       mat = uniqueInd)
    matchUpperUpperIndBins <- matMatch(x = cbind(upperIndBinsX, upperIndBinsY),
                                       mat = uniqueInd)

    # Start indexes
    indStart <- cbind(matchLowerLowerIndBins[-N], matchUpperLowerIndBins[-N],
                      matchLowerUpperIndBins[-N], matchUpperUpperIndBins[-N])

    # Start weights
    wStart <- weightsLinearInterp2D(x = data[-N, 1], y = data[-N, 2],
                                    gx1 = gridX[lowerIndBinsX[-N]],
                                    gx2 = gridX[upperIndBinsX[-N]],
                                    gy1 = gridY[lowerIndBinsY[-N]],
                                    gy2 = gridY[upperIndBinsY[-N]],
                                    circular = TRUE)

    # Initial conditions (normalized to ensure they are a density in the grid)
    u0 <- apply(uniqueInd, 1, function(ij) c(tcrossprod(
      dWn1D(x = gridX, mu = gridX[ij[1]], sigma = sdInitial),
      dWn1D(x = gridY, mu = gridY[ij[2]], sigma = sdInitial))))
    u0 <- sweep(u0, 2, apply(u0, 2, periodicTrapRule2D), FUN = "/")

    # Loglikelihood (expensive part)
    minusLogLik <- function(x) {

      # Update drifts and parameters
      bGrid <- b(x = gridXY, pars = x)
      bxGrid <- matrix(bGrid[, 1], nrow = Mx, ncol = My)
      byGrid <- matrix(bGrid[, 2], nrow = Mx, ncol = My)
      sigmaGrid <- sigma2(x = gridXY, pars = x)
      sigma2xGrid <- matrix(sigmaGrid[, 1], nrow = Mx, ncol = My)
      sigma2yGrid <- matrix(sigmaGrid[, 2], nrow = Mx, ncol = My)
      sigmaxyGrid <- matrix(0, nrow = Mx, ncol = My)

      # Go only throgout the unique indices (tpd[x1, x0])
      tpd <- crankNicolson2D(u0 = u0, bx = bxGrid, by = byGrid,
                             sigma2x = sigma2xGrid, sigma2y = sigma2yGrid,
                             sigmaxy = sigmaxyGrid, N = Mt, deltat = delt,
                             Mx = Mx, deltax = delx, My = My, deltay = dely,
                             imposePositive = TRUE)

      # Linear binning in evaluation and start
      res <- -sum(sapply(1:(N - 1), function(i)
        log(max(crossprod(wEval[i, ], tpd[indEval[i, ], indStart[i, ]]) %*%
                  wStart[i, ], 1e-08))))

      return(res)

    }

  } else {

    # closestIndBins <- cbind(
    #   sapply(1:nrow(data), function(i)
    #     which.min(abs(toPiInt(data[i, 1] - gridX)))),
    #   sapply(1:nrow(data), function(i)
    #     which.min(abs(toPiInt(data[i, 2] - gridY)))))
    closestIndBins <- toInt(round(sweep(data + pi, 2, c(delx, dely), "/")) + 1,
                            a = 1, b = Mx + 1)

    # Unique indexes of closests bins
    simpleUniqueInd <- unique(closestIndBins)

    # Link from closest bin to the unique values
    matchClosestIndBins <- matMatch(x = closestIndBins, mat = simpleUniqueInd)

    # Start indexes
    indStart <- matchClosestIndBins[-N]

    # Initial conditions (normalized to ensure they are a density in the grid)
    u0 <- apply(simpleUniqueInd, 1, function(ij) c(tcrossprod(
      dWn1D(x = gridX, mu = gridX[ij[1]], sigma = sdInitial),
      dWn1D(x = gridY, mu = gridY[ij[2]], sigma = sdInitial))))
    u0 <- sweep(u0, 2, apply(u0, 2, periodicTrapRule2D), FUN = "/")

    # Loglikelihood (expensive part)
    minusLogLik <- function(x) {

      # Update drifts and parameters
      bGrid <- b(x = gridXY, pars = x)
      bxGrid <- matrix(bGrid[, 1], nrow = Mx, ncol = My)
      byGrid <- matrix(bGrid[, 2], nrow = Mx, ncol = My)
      sigmaGrid <- sigma2(x = gridXY, pars = x)
      sigma2xGrid <- matrix(sigmaGrid[, 1], nrow = Mx, ncol = My)
      sigma2yGrid <- matrix(sigmaGrid[, 2], nrow = Mx, ncol = My)
      sigmaxyGrid <- matrix(0, nrow = Mx, ncol = My)

      # Go only throgout the unique indices (tpd[x1, x0])
      tpd <- crankNicolson2D(u0 = u0, bx = bxGrid, by = byGrid,
                             sigma2x = sigma2xGrid, sigma2y = sigma2yGrid,
                             sigmaxy = sigmaxyGrid, N = Mt, deltat = delt,
                             Mx = Mx, deltax = delx, My = My, deltay = dely,
                             imposePositive = TRUE)

      # Linear binning in evaluation
      res <- -sum(sapply(1:(N - 1), function(i)
        log(max(crossprod(wEval[i, ], tpd[indEval[i, ], indStart[i]]), 1e-08))))

      return(res)

    }

  }

  # Check positvedefiniteness
  region <- function(pars) {

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
  mleOptimWrapper(minusLogLik = minusLogLik, region = region, start = start,
                  lower = lower, upper = upper, ...)

}
