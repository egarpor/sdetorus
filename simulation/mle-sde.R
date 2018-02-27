
# Load packages
library(circular)
library(DirStats) # Remove dependency
library(doSNOW)
library(sdetorus)


### WN, WC, mivM and OU simulation functions
{

# Generate WN process 1D
rWn1D <- function(alpha, mu, sigma, NFine = 5e5, deltaFine = 1e-3,
                  seed = as.integer(runif(1) * 2e9), maxK = 2) {

  # Save seed for posterior analysis
  set.seed(seed)

  # Initial point
  x0 <- toPiInt(rnorm(n = 1, mean = mu, sd = sigma / sqrt(2 * alpha)))

  # Fine trajectory
  data <- euler1D(x0 = x0, mu = mu, alpha = alpha, sigma = sigma, N = NFine,
                  type = 1L, delta = deltaFine, maxK = maxK)[1, ]
  # twokpi <- 2 * pi * seq(-3, 3, by = 1)
  # sda <- sigma / sqrt(2 * alpha)
  # sdt <- sqrt(vartOu(t = deltaFine, alpha = alpha, sigma = sigma))
  # data <- numeric(NFine + 1); data[1] <- x0
  # for (i in 1:NFine) {
  #   w <- dnorm(data[i] + twokpi, mean = mu, sd = sda)
  #   twompi <- sample(x = twokpi, size = 1, prob = w)
  #   data[i + 1] <- toPiInt(rnorm(n = 1, mean = meantOu(x0 = data[i] + twompi,
  #                                                      t = deltaFine,
  #                                                      alpha = alpha,
  #                                                       mu = mu), sd = sdt))
  # }

  return(list("data" = data, "seed" = seed,
              "pars" = c("alpha" = alpha, "mu" = mu, "sigma" = sigma,
                         "deltaFine" = deltaFine, "NFine" = NFine)))

}

# Generate WC process 1D
rWc1D <- function(alpha, mu, sigma, NFine = 5e5, deltaFine = 1e-3,
                  seed = as.integer(runif(1) * 2e9)) {

  # Save seed for posterior analysis
  set.seed(seed)

  # Initial point
  x0 <- toPiInt(rcauchy(n = 1, location = mu, scale = sigma / sqrt(2 * alpha)))

  # Fine trajectory
  data <- rTrajLangevin(x0 = x0, drift = driftJp, SigDif = sigma^2, N = NFine,
                        delta = deltaFine, NFine = NFine, deltaFine = deltaFine,
                        circular = TRUE, alpha = alpha, mu = mu, psi = -1)

  return(list("data" = data, "seed" = seed,
              "pars" = c("alpha" = alpha, "mu" = mu, "sigma" = sigma,
                         "deltaFine" = deltaFine, "NFine" = NFine)))

}

# Fit WN and WC 1D
fit1D <- function(listData, NSub, deltaSub, type = "WN", start = "stat",
                  vmApprox = FALSE, maxK = 2, pde = FALSE, SO = TRUE, lower,
                  upper, optMethod = "Nelder-Mead", alphaFix = NA, muFix = NA,
                  sigmaFix = NA, sdInitial = 0.1, ...) {

  ## Startup ##
  {

  lp <- length(listData$pars)
  data <- listData$data[seq.int(1, listData$pars[lp],
                                by = deltaSub / listData$pars[lp - 1])[1:NSub]]
  delta <- deltaSub

  }

  ## Fixed parameters ##
  {

  # Specfied parameters
  specPars <- listData$pars[1:(lp - 2)]
  indUnSpecPars <- is.na(c(alphaFix, muFix, sigmaFix))
  lp <- sum(indUnSpecPars)
  alpha <- ifelse(is.na(alphaFix), NA, specPars[1])
  mu <- ifelse(is.na(muFix), NA, specPars[2])
  sigma <- ifelse(is.na(sigmaFix), NA, specPars[3])

  # Lower and upper values
  if (missing(lower)) {

    lower <- c(0.01, -pi, 0.01)[indUnSpecPars]

  }
  if (missing(upper)) {

    upper <- c(25, pi, 25)[indUnSpecPars]

  }

  }

  ## Assign results ##

  assignRes <- function(res, suffix) {

    res <- c(res[1], toPiInt(res[2]), res[3], res[4], time)
    names(res) <- paste(c("alpha", "mu", "sigma", "loglik", "time"), suffix,
                        sep = "")

    return(res)

  }

  ## Diffusion-specific ##
  if (type == "WN") {

    # Drift and derivatives of the drift
    b <- function(x, pars){

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      driftWn1D(x = x, alpha = specPars[1], mu = specPars[2],
                sigma = specPars[3], maxK = maxK)

    }
    b1 <- function(x, pars, h = 1e-5){

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      l <- length(x)
      res <- driftWn1D(x = c(x + h, x - h), alpha = specPars[1],
                       mu = specPars[2], sigma = specPars[3], maxK = maxK)
      drop(res[1:l] - res[(l + 1):(2 * l)]) / (2 * h)

    }

    # Diffusion
    sigma2 <- function(x, pars){

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      rep(specPars[3]^2, length(x))

    }

  } else if (type == "WC") {

    # Drift and derivatives of the drift
    b <- function(x, pars){

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      driftJp(x = x, alpha = specPars[1], mu = specPars[2], psi = -1)

    }
    b1 <- function(x, pars, h = 1e-5){

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      l <- length(x)
      res <- driftJp(x = c(x + h, x - h), alpha = specPars[1], mu = specPars[2],
                     psi = -1)

      drop(res[1:l] - res[(l + 1):(2 * l)]) / (2 * h)

    }

    # Diffusion
    sigma2 <- function(x, pars){

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      rep(specPars[3]^2, length(x))

    }

  }

  ## Stationary distribution ##
  {

  # Start time
  time <- proc.time()[3]

  # Diffusion matrix
  sigmaStat <- sqrt(sigmaDiff(data = data, delta = delta))
  if (!is.na(sigmaFix)) {
    sigmaStat <- specPars[3]
  }

  if (type == "WN") {

    # MLE
    ml <- mle.wrappednormal(x = circular(data), max.iter = 1e3, tol = 1e-5)
    muStat <- toPiInt(as.numeric(ml$mu))
    alphaStat <- sigmaStat^2 / (2 * ml$sd^2)

    # Loglikelihood
    loglikStat <- -sum(log(dwrappednormal(circular(data), mu = circular(muStat),
                                          sd = ml$sd)))

  } else if (type == "WC") {

    # MLE
    ml <- mle.wrappedcauchy(x = circular(data), max.iter = 1e3, tol = 1e-5)
    muStat <- toPiInt(as.numeric(ml$mu))
    kappaStat <- 2 * tanh(ml$rho)
    alphaStat <- kappaStat * sigmaStat^2 / 2

    # Loglikelihood and ratio
    loglikStat <- -sum(log(dJp(data, mu = muStat, kappa = kappaStat, psi = -1)))

  }

  # End time
  time <- proc.time()[3] - time
  names(time) <- NULL

  # Object
  resStat <- assignRes(c(alphaStat, toPiInt(muStat), sigmaStat, loglikStat),
                       suffix = "Stat")

  }

  ## Starting values ##
  {

  # Choose between stationary estimates or true parameters
  if (start == "stat") {

    start <- rbind(c(alphaStat, muStat, sigmaStat)[indUnSpecPars])

  } else if (start == "true") {

    start <- rbind(specPars[indUnSpecPars])

  } else if (start == "rand") {

    start <- rbind(c(alphaStat, muStat, sigmaStat),
                   cbind(runif(9, 0.01, 5), runif(9, -pi, pi),
                         runif(9, 0.01, 5)))[, indUnSpecPars]

  }

  }

  ## Euler and Ozaki ##
  {

  # Function to avoid repeating code
  ps <- function(method, circular, suffix) {

    # Start time
    time <- proc.time()[3]

    # Procedure
    resAux <- tryCatch({

      fit <- psMle(data = data, delta = delta, method = method, b = b, b1 = b1,
                   sigma2 = sigma2, start = start, lower = lower, upper = upper,
                   circular = circular, vmApprox = vmApprox, maxK = maxK,
                   optMethod = optMethod, ...)

      if (fit$convergence & fit$localMinimumGuaranteed) {

        c(fit$par, fit$value)

      } else {

        rep(NA, lp + 1)

      }

    }, error = function(e) {show(e); rep(NaN, lp + 1)})

    # End time
    time <- proc.time()[3] - time

    # Object
    res <- c(specPars, NA)
    res[c(which(indUnSpecPars), length(res))] <- resAux
    res <- assignRes(res, suffix = suffix)

    return(res)

  }

  resE <- ps(method = "E", circular = TRUE, suffix = "E")
  resSO <- ps(method = "SO", circular = TRUE, suffix = "SO")

  }

  ## Approximate tpd WN ##
  if (type == "WN") {

    # Start time
    time <- proc.time()[3]

    # Procedure
    resAux <- tryCatch({

      fit <- approxMleWn1D(data = data, delta = delta, alpha = alpha, mu = mu,
                           sigma = sigma, start = start, lower = lower,
                           upper = upper, vmApprox = vmApprox, maxK = maxK,
                           optMethod = optMethod, ...)

      if (fit$convergence & fit$localMinimumGuaranteed) {

        c(fit$par, fit$value)

      } else {

        rep(NA, lp + 1)

      }

    }, error = function(e) {show(e); rep(NaN, lp + 1)})

    # End time
    time <- proc.time()[3] - time

    # Object
    resAp <- c(specPars, NA)
    resAp[c(which(indUnSpecPars), length(resAp))] <- resAux
    resAp <- assignRes(res = resAp, suffix = "Ap")

  }

  ## Exact MLE ##
  if (pde) {

    # Start time
    time <- proc.time()[3]

    # Procedure
    resAux <- tryCatch({

      fitInitial <- mlePde1D(data = data, delta = delta, b = b, sigma2 = sigma2,
                             Mx = 50, linearBinning = FALSE,
                             sdInitial = sdInitial, start = start,
                             lower = lower, upper = upper,
                             optMethod = optMethod, ...)

      fit <- mlePde1D(data = data, delta = delta, b = b, sigma2 = sigma2,
                      linearBinning = FALSE, sdInitial = sdInitial,
                      start = t(sapply(fitInitial$solutionsOutput,
                                       function(x) x$par)),
                      lower = lower, upper = upper, optMethod = optMethod, ...)

      if (fit$convergence & fit$localMinimumGuaranteed) {

        c(fit$par, fit$value)

      } else {

        rep(NA, lp + 1)

      }

    }, error = function(e) {show(e); rep(NaN, lp + 1)})

    # End time
    time <- proc.time()[3] - time

    # Object
    resPde <- c(specPars, NA)
    resPde[c(which(indUnSpecPars), length(resPde))] <- resAux
    resPde <- assignRes(res = resPde, suffix = "Pde")

  }

  ## Return ##
  if (type == "WN") {

    if (pde) {

      if (SO) {

        return(c(resStat, resE, resSO, resAp, resPde))

      } else {

        return(c(resStat, resE, resAp, resPde))

      }

    } else {

      if (SO) {

        return(c(resStat, resE, resSO, resAp))

      } else {

        return(c(resStat, resE, resAp))

      }

    }

  } else if (type == "WC") {

    if (pde) {

      if (SO) {

        return(c(resStat, resE, resSO, resPde))

      } else {

        return(c(resStat, resE, resPde))

      }

    } else {

      if (SO) {

        return(c(resStat, resE, resSO))

      } else {

        return(c(resStat, resE))

      }

    }

  }

}

# Generate WN process 2D
rWn2D <- function(alpha1, alpha2, alpha3, mu1, mu2, sigma1, sigma2, NFine = 5e5,
                  deltaFine = 1e-3, seed = as.integer(runif(1) * 2e9),
                  maxK = 2) {

  # Save seed for posterior analysis
  set.seed(seed)

  # Initial point
  x0 <- toPiInt(mvtnorm::rmvnorm(n = 1, mean = c(mu1, mu2),
                                 sigma = 1/2 * solve(alphaToA(
                                   alpha = c(alpha1, alpha2, alpha3),
                                   sigma = c(sigma1, sigma2)))
                                 %*% diag(c(sigma1^2, sigma2^2))))

  # Fine trajectory
  data <- t(euler2D(x0 = x0, mu = c(mu1, mu2),
                    A = alphaToA(alpha = c(alpha1, alpha2, alpha3),
                                 sigma = c(sigma1, sigma2)),
                    sigma = c(sigma1, sigma2), N = NFine, delta = deltaFine,
                    type = 1L, maxK = maxK)[1, , ])

  return(list("data" = data, "seed" = seed,
              pars = c("alpha" = c(alpha1, alpha2, alpha3), "mu" = c(mu1, mu2),
                       "sigma" = c(sigma1, sigma2), "deltaFine" = deltaFine,
                       "NFine" = NFine)))

}

# Generate MivM process 2D
rMivM2D <- function(alpha1, alpha2, mu1, mu2, sigma, p, NFine = 5e5,
                    deltaFine = 1e-3, seed = as.integer(runif(1) * 2e9)) {

  # Save seed for posterior analysis
  set.seed(seed)

  # Initial point
  if (runif(1) < p) {

    x0 <- circular:::RvonmisesRad(n = 2, mu = mu1, kappa = 2 * alpha1 / sigma^2)

  } else {

    x0 <- circular:::RvonmisesRad(n = 2, mu = mu2, kappa = 2 * alpha2 / sigma^2)

  }

  # Fine trajectory
  data <- rTrajLangevin(x0 = x0, drift = driftMixIndVm,
                        SigDif = diag(rep(sigma^2, 2)), N = NFine,
                        delta = deltaFine, NFine = NFine, deltaFine = deltaFine,
                        A = rbind(rep(alpha1, 2), rep(alpha2, 2)),
                        M = rbind(rep(mu1, 2), rep(mu2, 2)), sigma = sigma,
                        p = c(p, 1 - p), circular = TRUE)

  return(list("data" = data, "seed" = seed,
              pars = c("alpha" = c(alpha1, alpha2), "mu" = c(mu1, mu2), "p" = p,
                       "sigma" = sigma, "deltaFine" = deltaFine,
                       "NFine" = NFine)))

}

# Fit WN and MivM 2D
fit2D <- function(listData, NSub, deltaSub, type = "WN", start = "stat",
                  vmApprox = FALSE, maxK = 2, pde = FALSE, SO = TRUE, lower, upper,
                  optMethod = "Nelder-Mead", alphaFix = rep(NA, 3),
                  muFix = rep(NA, 2), pFix = NA, sigmaFix = rep(NA, 2),
                  sdInitial = 0.1, ...) {

  ## Startup ##
  {

  lp <- length(listData$pars)
  data <- listData$data[seq.int(1, listData$pars[lp],
                                by = deltaSub/listData$pars[lp - 1])[1:NSub], ]
  delta <- deltaSub
  pars <- listData$pars[1:(lp - 2)]

  }

  ## Fixed parameters ##
  {

  # Specfied parameters
  specPars <- listData$pars[1:(lp - 2)]

  if (type == "WN") {

    indUnSpecPars <- is.na(c(alphaFix, muFix, sigmaFix))
    lp <- sum(indUnSpecPars)
    if (any(is.na(alphaFix))) {

      alpha <- rep(NA, 3)

    } else {

      alpha <- specPars[1:3]

    }
    if (any(is.na(muFix))) {

      mu <- rep(NA, 2)

    } else {

      mu <- specPars[4:5]

    }
    if (any(is.na(sigmaFix))) {

      sigma <- rep(NA, 2)

    } else {

      sigma <- specPars[6:7]

    }

    # Lower and upper values
    if (missing(lower)) {

      lower <- c(0.01, 0.01, -25, -pi, -pi, 0.01, 0.01)[indUnSpecPars]

    }
    if (missing(upper)) {

      upper <- c(25, 25, 25, pi, pi, 25, 25)[indUnSpecPars]

    }

  } else if (type == "mivM") {

    indUnSpecPars <- is.na(c(alphaFix[1:2], muFix, pFix, sigmaFix[1]))
    lp <- sum(indUnSpecPars)
    if (any(is.na(alphaFix))) {

      alpha <- rep(NA, 2)

    } else {

      alpha <- specPars[1:2]

    }
    if (any(is.na(muFix))) {

      mu <- rep(NA, 2)

    } else {

      mu <- specPars[3:4]

    }
    if (any(is.na(pFix))) {

      p <- NA

    } else {

      p <- specPars[5]

    }
    if (any(is.na(sigmaFix))) {

      sigma <- NA

    } else {

      sigma <- specPars[6]

    }

    # Lower and upper values
    if (missing(lower)) {

      lower <- c(0.01, 0.01, -pi, -pi, -Inf, 0.01)[indUnSpecPars]

    }
    if (missing(upper)) {

      upper <- c(25, 25, pi, pi, Inf, 25)[indUnSpecPars]

    }

  }

  }

  ## Assign results ##

  assignRes <- function(res, suffix, time) {

    if (type == "WN") {

      res <- c(res[1:3], toPiInt(res[4:5]), res[6:7], res[8], time)
      names(res) <- paste(c("alpha1", "alpha2", "alpha3", "mu1", "mu2",
                            "sigma1", "sigma2", "loglik", "time"),
                          suffix, sep = "")

    } else if (type == "mivM") {

      res <- c(res[1:2], toPiInt(res[3:4]), res[5:6], res[7], time)
      names(res) <- paste(c("alpha1", "alpha2", "mu1", "mu2", "p",
                            "sigma", "loglik", "time"), suffix, sep = "")

    }

    return(res)

  }

  ## Diffusion-specific ##
  if (type == "WN") {

    # Drift and derivatives of the drift
    b <- function(x, pars) {

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      driftWn2D(x = x, A = alphaToA(alpha = specPars[1:3],
                                    sigma = specPars[6:7]), mu = specPars[4:5],
                sigma = specPars[6:7], maxK = maxK)

    }
    jac.b <- function(x, pars, h = 1e-5) {

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      l <- nrow(x)
      res <- driftWn2D(x = rbind(cbind(x[, 1] + h, x[, 2]),
                                 cbind(x[, 1] - h, x[, 2]),
                                 cbind(x[, 1], x[, 2] + h),
                                 cbind(x[, 1], x[, 2] - h)),
                       A = alphaToA(alpha = specPars[1:3],
                                    sigma = specPars[6:7]),
                       mu = specPars[4:5], sigma = specPars[6:7], maxK = maxK)
      cbind(res[1:l, ] - res[(l + 1):(2 * l), ],
            res[2 * l + 1:l, ] - res[2 * l + (l + 1):(2 * l), ]) / (2 * h)
    }
    sigma2 <- function(x, pars) {

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      matrix(specPars[6:7]^2, nrow = length(x) / 2L, ncol = 2)

    }

  } else if (type == "mivM") {

    # Drift and derivatives of the drift
    b <- function(x, pars) {

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      p <- (tanh(specPars[5]) + 1) / 2
      driftMixIndVm(x = x, A = matrix(specPars[1:2], nrow = 2, ncol = 2),
                    M = matrix(specPars[3:4], nrow = 2, ncol = 2),
                    p = c(p, 1 - p), sigma = specPars[6])

    }
    jac.b <- function(x, pars, h = 1e-5) {

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      l <- nrow(x)
      p <- (tanh(specPars[5]) + 1) / 2
      res <- driftMixIndVm(x = rbind(cbind(x[, 1] + h, x[, 2]),
                                     cbind(x[, 1] - h, x[, 2]),
                                     cbind(x[, 1], x[, 2] + h),
                                     cbind(x[, 1], x[, 2] - h)),
                       A = matrix(specPars[1:2], nrow = 2, ncol = 2),
                       M = matrix(specPars[3:4], nrow = 2, ncol = 2),
                       p = c(p, 1 - p), sigma = specPars[6])
      cbind(res[1:l, ] - res[(l + 1):(2 * l), ],
            res[2 * l + 1:l, ] - res[2 * l + (l + 1):(2 * l), ]) / (2 * h)
    }
    sigma2 <- function(x, pars) {

      # Change only unspecified parameters
      specPars[indUnSpecPars] <- pars

      matrix(specPars[6]^2, nrow = length(x) / 2L, ncol = 2)

    }

  }

  ## Stationary distribution ##

  # Start time
  time <- proc.time()[3]

  if (type == "WN") {

    # Diffusion matrix
    sigmaStat <- sigmaDiff(data = data, delta = delta, diagonal = TRUE)
    if (!any(is.na(sigma))) {

      sigmaStat <- diag(sigma^2, nrow = 2)

    }

    # MLE
    ml1 <- mle.wrappednormal(x = circular(data[, 1]))
    ml2 <- mle.wrappednormal(x = circular(data[, 2]))
    ml <- mleOptimWrapper(minusLogLik = function(x) {
      -sum(log(DirStats::dwntd(th1 = data[, 1], th2 = data[, 2], mu1 = x[1],
                               mu2 = x[2], sigma1 = x[3], sigma2 = x[4],
                               rho = x[5], K = 3)))
    }, start = c(toPiInt(as.numeric(ml1$mu)), toPiInt(as.numeric(ml2$mu)),
                 ml1$sd, ml2$sd, 0), lower = c(-pi, -pi, 0.01, 0.01, -1),
    upper = c(pi, pi, 20, 20, 1), optMethod = "Nelder-Mead")
    mu1Stat <- ml$par[1]
    mu2Stat <- ml$par[2]
    SWN <- matrix(c(ml$par[3]^2, ml$par[5] * ml$par[3] * ml$par[4],
                    ml$par[5] * ml$par[3] * ml$par[4], ml$par[4]^2),
                  nrow = 2, ncol = 2, byrow = TRUE)

    A <- 0.5 * sigmaStat %*% solve(SWN)
    alpha1Stat <- A[1, 1]
    alpha2Stat <- A[2, 2]
    alpha3Stat <- A[1, 2] * sqrt(sigmaStat[2, 2] / sigmaStat[1, 1])
    sigma1Stat <- sqrt(sigmaStat[1, 1])
    sigma2Stat <- sqrt(sigmaStat[2, 2])

    # End time
    time <- proc.time()[3] - time
    names(time) <- NULL

    # Object
    resStat <- c(alpha1Stat = alpha1Stat,
                 alpha2Stat = alpha2Stat,
                 alpha3Stat = alpha3Stat,
                 mu1Stat = toPiInt(mu1Stat),
                 mu2Stat = toPiInt(mu2Stat),
                 sigma1Stat = sigma1Stat,
                 sigma2Stat = sigma2Stat,
                 loglikStat = ml$value,
                 timeStat = time)

  } else if (type == "mivM") {

    # Diffusion matrix
    sigmaStat <- sqrt(sigmaDiff(data = data, delta = delta,
                                isotropic = TRUE)[1, 1])
    if (!any(is.na(sigma))) {

      sigmaStat <- sigma
      names(sigmaStat) <- NULL

    }

    # MLE
    ml1 <- sdetorus:::emVmf(data = data, k = 2, kappaMax = 100, maxIter = 2000)
    ml <- mleOptimWrapper(minusLogLik = function(x) {
      p <- (tanh(x[5]) + 1) / 2
      -sum(log(sdetorus:::dVmf(x = data, K = matrix(x[1:2], 2, 2),
                                M = matrix(x[3:4], 2, 2),
                                alpha = c(p, 1 - p))))
    }, start = c(rowMeans(ml1$K), circular:::MeanCircularRad(ml1$M[1, ]),
                 circular:::MeanCircularRad(ml1$M[2, ]),
                 atanh(2 * ml1$alpha[1] - 1)),
    lower = c(0.01, 0.01, -2 * pi, -2 * pi, -Inf),
    upper = c(20, 20, 2 * pi, 2 * pi, Inf), optMethod = "Nelder-Mead",
    selectSolution = "lowest")
    kappa1Stat <- ml$par[1]
    kappa2Stat <- ml$par[2]
    mu1Stat <- ml$par[3]
    mu2Stat <- ml$par[4]
    pStat <- ml$par[5]
    alpha1Stat <- kappa1Stat * sigmaStat^2 / 2
    alpha2Stat <- kappa2Stat * sigmaStat^2 / 2

    # End time
    time <- proc.time()[3] - time
    names(time) <- NULL

    # Object
    resStat <- c(alpha1Stat = alpha1Stat,
                 alpha2Stat = alpha2Stat,
                 mu1Stat = toPiInt(mu1Stat),
                 mu2Stat = toPiInt(mu2Stat),
                 pStat = pStat,
                 sigmaStat = sigmaStat,
                 loglikStat = ml$value,
                 timeStat = time)

  }

  ## Starting values ##

  # Choose between stationary estimates or true parameters
  if (start == "stat") {

    if (type == "WN") {

      start <- rbind(c(alpha1Stat, alpha2Stat, alpha3Stat, mu1Stat, mu2Stat,
                       sigma1Stat, sigma2Stat)[indUnSpecPars])

    } else if (type == "mivM") {

      start <- rbind(c(alpha1Stat, alpha2Stat, mu1Stat, mu2Stat,
                       pStat, sigmaStat)[indUnSpecPars])

    }

  } else if (start == "true") {

    start <- rbind(pars[indUnSpecPars])

  } else if (start == "rand") {

    a1 <- runif(9, 0.01, 5)
    a2 <- runif(9, 0.01, 5)
    m1 <- runif(9, -pi, pi)
    m2 <- runif(9, -pi, pi)
    s1 <- runif(9, 0.01, 5)
    s2 <- runif(9, 0.01, 5)
    p <- 2 * atanh(runif(9, 0, 1) - 1)

    if (type == "WN") {

      start <- rbind(c(alpha1Stat, alpha2Stat, alpha3Stat, mu1Stat, mu2Stat,
                       sigma1Stat, sigma2Stat),
                     cbind(a1, a2, runif(9, -sqrt(abs(a1 * a2)),
                                         sqrt(abs(a1 * a2))),
                           m1, m2, s1, s2))[, indUnSpecPars]

    } else if (type == "mivM") {

      start <- rbind(c(alpha1Stat, alpha2Stat, mu1Stat, mu2Stat, pStat,
                       sigmaStat),
                     cbind(a1, a2, m1, m2, p, s1))[, indUnSpecPars]

    }

  }

  ## Euler and Ozaki ##
  {

  # Function to avoid repeating code
  ps <- function(method, circular, suffix) {

    # Start time
    time <- proc.time()[3]

    # Procedure
    resAux <- tryCatch({

      fit <- psMle(data = data, delta = delta, method = method, b = b,
                   jac.b = jac.b, sigma2 = sigma2, start = start,
                   lower = lower, upper = upper, circular = circular,
                   vmApprox = vmApprox, maxK = maxK, optMethod = optMethod, ...)

      if (fit$convergence & fit$localMinimumGuaranteed) {

        c(fit$par, fit$value)

      } else {

        rep(NA, lp + 1)

      }

    }, error = function(e) {show(e); rep(NaN, lp + 1)})

    # End time
    time <- proc.time()[3] - time

    # Object
    res <- c(specPars, NA)
    res[c(which(indUnSpecPars), length(res))] <- resAux
    assignRes(res, suffix, time)

  }

  resE <- ps(method = "E", circular = TRUE, suffix = "E")
  if (SO) {

    resSO <- ps(method = "SO", circular = TRUE, suffix = "SO")

  }

  }

  ## Approximate tpd WN ##
  if (type == "WN") {

    # Start time
    time <- proc.time()[3]

    # Procedure
    resAux <- tryCatch({

      fit <- approxMleWn2D(data = data, delta = delta, alpha = alpha, mu = mu,
                           sigma = sigma, start = start, lower = lower,
                           upper = upper, maxK = maxK, optMethod = optMethod,
                           rho = 0, ...)

      if (fit$convergence & fit$localMinimumGuaranteed) {

        c(fit$par, fit$value)

      } else {

        rep(NA, lp + 1)

      }

    }, error = function(e) {show(e); rep(NaN, lp + 1)})

    # End time
    time <- proc.time()[3] - time

    # Object
    resAp <- c(specPars, NA)
    resAp[c(which(indUnSpecPars), length(resAp))] <- resAux
    resAp <- assignRes(resAp, "Ap", time)

  }

  ## Exact MLE ##
  if (pde) {

    # Start time
    time <- proc.time()[3]

    # Procedure
    resAux <- tryCatch({

      fitInitial <- mlePde2D(data = data, delta = delta, b = b, sigma2 = sigma2,
                             Mx = 15, My = 15, linearBinning = FALSE,
                             sdInitial = sdInitial, start = start,
                             lower = lower, upper = upper,
                             optMethod = optMethod, ...)

      fit <- mlePde2D(data = data, delta = delta, b = b, sigma2 = sigma2,
                      start = t(sapply(fitInitial$solutionsOutput,
                                       function(x) x$par)),
                      lower = lower, upper = upper, sdInitial = sdInitial,
                      linearBinning = FALSE, optMethod = optMethod, ...)

      if (fit$convergence & fit$localMinimumGuaranteed) {

        c(fit$par, fit$value)

      } else {

        rep(NA, lp + 1)

      }

    }, error = function(e) {show(e); rep(NaN, lp + 1)})

    # End time
    time <- proc.time()[3] - time

    # Object
    resPde <- c(specPars, NA)
    resPde[c(which(indUnSpecPars), length(resPde))] <- resAux
    resPde <- assignRes(resPde, "Pde", time)

  }

  if (type == "WN") {

    if (pde) {

      if (SO) {

        return(c(resStat, resE, resSO, resAp, resPde))

      } else {

        return(c(resStat, resE, resAp, resPde))

      }

    } else {

      if (SO) {

        return(c(resStat, resE, resSO, resAp))

      } else {

        return(c(resStat, resE, resAp))

      }

    }

  } else if (type == "mivM") {

    if (pde) {

      if (SO) {

        return(c(resStat, resE, resSO, resPde))

      } else {

        return(c(resStat, resE, resPde))

      }

    } else {

      if (SO) {

        return(c(resStat, resE, resSO))

      } else {

        return(c(resStat, resE))

      }

    }

  }

}

}


#' @title Wrapper for creating a LaTeX table with the significant best estimate
#'
#' @param errors data frame with column names containing \code{factors} and \code{estimators}.
#' @param factors names of the factors of the data frame to split the summary according to.
#' @param estimators names of the estimators whose errors are to be summarized.
#' @param alpha significance level for the one-sided \code{\link[stats]{t.test}} for paired samples.
#' @value A list with the entries:
#' \itemize{
#' \item \code{uniqueFactors}: the unique factors in the data frame.
#' \item \code{tab}: vector of character strings, suitable for a LaTeX table, containing
#' \code{meanErrors} with bold fonts for the significant lowest MSEs.
#' \item \code{meanErrors}: matrix of size \code{c(length(uniqueFactors), length(estimators))}
#' containing the average MSEs.
#' }
aggregateSignif <- function(errors, factors, estimators, alpha = 0.05) {

  # Columns of the factors and estimators
  indFactors <- na.omit(match(x = factors, table = colnames(errors)))
  indEst <- na.omit(match(x = estimators, table = colnames(errors)))

  # Go through unique factors
  uniqueFactors <- unique(errors[, indFactors])
  nu <- nrow(uniqueFactors)
  tab <- character(length = length(uniqueFactors))
  meanErrors <- matrix(nrow = nu, ncol = length(indEst))
  for (i in 1:nu) {

    # Index subset with factors == uniqueFactors
    ind <- apply(errors[, indFactors] == repRow(uniqueFactors[i, ],
                                                nrow(errors)), 1, all)

    # Subset
    subErrors <- errors[ind, indEst]

    # MSE
    meanErrors[i, ] <- colMeans(subErrors, na.rm = TRUE)

    # Index of minimum MSE
    indMin <- which.min(meanErrors[i, ])

    # Are the other MSEs equal?
    signifs <- sapply(seq_along(indEst)[-indMin], function(j) {
      t.test(x = subErrors[, indMin], y = subErrors[, j], alternative = "less",
             paired = TRUE, var.equal = FALSE)$p.value
    }) > alpha

    # Boldfaces
    bold <- rep("", length(indEst))
    bold[-indMin][signifs] <- bold[indMin] <- "\\bf "

    # Print table
    tab[i] <- paste(paste(bold, sprintf("%.4f", meanErrors[i, ]), sep = ""),
                    collapse = " & ")

    # Assign names
    names(tab)[i] <- paste(paste(factors, uniqueFactors[i, ], sep = ""),
                           collapse = "-")

  }

  return(list("uniqueFactors" = uniqueFactors, "tab" = tab,
              "meanErrors" = meanErrors))

}


#' @title Wrapper for creating a LaTeX table with the average of the relative efficiencies of a set of multivariate estimators
#'
#' @param errors data frame with column names containing \code{factors} and \code{estimators}.
#' @param factors names of the factors of the data frame to split the summary according to.
#' @param estimators list containing in each entry the names of the estimators
#' whose errors are to be summarized.
#' @value A list with the entries:
#' \itemize{
#' \item \code{uniqueFactors}: the unique factors in the data frame.
#' \item \code{tab}: vector of character strings, suitable for a LaTeX table, containing
#' \code{meanErrors} with bold fonts for the significant lowest MSEs.
#' \item \code{meanErrors}: matrix of size \code{c(length(uniqueFactors), length(estimators[[1]]))}
#' containing the average MSEs.
#' }
aggregateRes <- function(errors, factors, estimators) {

  # Columns of the factors and estimators
  indFactors <- na.omit(match(x = factors, table = colnames(errors)))

  # Check estimators list
  nEstimators <- length(estimators)
  nMethods <- length(estimators[[1]])
  stopifnot(all(nMethods == lapply(estimators, length)))

  # Go through unique factors
  uniqueFactors <- unique(errors[, indFactors])
  nu <- nrow(uniqueFactors)
  tab <- character(length = length(uniqueFactors))
  meanErrors <- matrix(0, nrow = nu, ncol = nMethods)
  for (i in 1:nu) {

    # Index subset with factors == uniqueFactors
    ind <- apply(errors[, indFactors] == repRow(uniqueFactors[i, ],
                                                nrow(errors)), 1, all)
    # Subset
    subErrors <- errors[ind, ]

    # Loop in the estimators
    for (l in 1:nEstimators) {

      # Match estimators
      indVars <- na.omit(match(x = estimators[[l]], table = colnames(errors)))

      # MSE
      meanErrorsAux <- colMeans(subErrors[, indVars], na.rm = TRUE)

      # Index of minimum MSE
      indMin <- which.min(meanErrorsAux)

      # RE with respect to the lowest MSE
      meanErrorsAux <- (meanErrorsAux[indMin] / nEstimators) / meanErrorsAux

      # Accumulate RE
      meanErrors[i, ] <- meanErrors[i, ] + meanErrorsAux

    }

    # Boldfaces
    bold <- rep("", nMethods)
    bold[which.max(meanErrors[i, ])] <- "\\bf "

    # Print table for average RE
    tab[i] <- paste(paste(bold, sprintf("%.4f", meanErrors[i, ])), collapse = " & ")

    # Assign names
    names(tab)[i] <- paste(paste(factors, uniqueFactors[i, ], sep = ""),
                           collapse = "-")

  }

  return(list("uniqueFactors" = uniqueFactors, "tab" = tab,
              "meanErrors" = meanErrors))

}


### Tests for implementations
dontrun <- TRUE
if (!dontrun) {

  ## 1D

  # Parameters
  NSub <- 250; deltaSub <- .5; maxK <- 2
  alpha <- 1; mu <- 0; sigma <- 2
  optMethod <- "Nelder-Mead"

  # WN
  f <- matrix(nrow = 100, ncol = 20)
  pb <- txtProgressBar(style = 3)
  for (i in 1:100) {

    listData <- rWn1D(alpha = alpha, mu = mu, sigma = sigma, NFine = 5e5,
                      deltaFine = 1e-3, seed = 34567 + i * 345, maxK = maxK)
    f[i, ] <- fit1D(listData = listData, NSub = NSub, deltaSub = deltaSub,
                    type = "WN", maxK = maxK, pde = FALSE, start = "rand",
                    optMethod = optMethod, sigmaFix = sigma, eigTol = 1e-3,
                    condTol = 1e3)

    setTxtProgressBar(pb = pb, i / 100)

  }
  colnames(f) <- names(fit1D(listData = listData, NSub = NSub,
                             deltaSub = deltaSub, type = "WN", maxK = maxK,
                             pde = FALSE, start = "true", optMethod = optMethod,
                             sigmaFix = sigma))
  colMeans(f, na.rm = TRUE)
  colMeans(t(t(f) - colMeans(f, na.rm = TRUE))^2, na.rm = TRUE)
  f <- as.data.frame(f)
  rbind(summary((f$alphaE - alpha)^2, na.rm = TRUE),
        summary((f$alphaSO - alpha)^2, na.rm = TRUE),
        summary((f$alphaAp - alpha)^2, na.rm = TRUE))
  rbind(summary((f$muE - mu)^2, na.rm = TRUE),
        summary((f$muSO - mu)^2, na.rm = TRUE),
        summary((f$muAp - mu)^2, na.rm = TRUE))

  hist((f$alphaE - alpha)^2, breaks = seq(0, 0.25, l = 20), ylim = c(0, 80))
  rug((f$alphaE - alpha)^2, col = 1)
  hist((f$alphaSO - alpha)^2, add = TRUE, border = 2,
       breaks = seq(0, 0.25, l = 20))
  rug((f$alphaSO - alpha)^2, col = 2)
  hist((f$alphaAp - alpha)^2, add = TRUE, border = 3,
       breaks = seq(0, 0.25, l = 20))
  rug((f$alphaAp - alpha)^2, col = 3)

  # WC
  listData <- rWc1D(alpha = alpha, mu = mu, sigma = sigma, NFine = 1e6,
                    deltaFine = 1e-3)

  fit1D(listData = listData, NSub = NSub, deltaSub = deltaSub, type = "WC",
        maxK = maxK, pde = TRUE, start = "true", optMethod = optMethod,
        sigmaFix = sigma)

  ## 2D

  # Parameters
  NSub <- 500; deltaSub <- 0.5; maxK <- 2
  alpha <- c(5, 5, 0); mu <- c(2, -2); sigma <- c(2, 2); p <- 0.5
  optMethod <- "Nelder-Mead"

  # WN
  listData <- rWn2D(alpha1 = alpha[1], alpha2 = alpha[2], alpha3 = alpha[3],
                    mu1 = mu[1], mu2 = mu[2], sigma1 = sigma[1],
                    sigma2 = sigma[2], NFine = 5e5, deltaFine = 1e-3)

  fit2D(listData = listData, NSub = NSub, deltaSub = deltaSub, type = "WN",
        maxK = maxK, pde = FALSE, start = "stat", sigmaFix = sigma,
        optMethod = optMethod)

  # mivM
  listData <- rMivM2D(alpha1 = alpha[1], alpha2 = alpha[2], mu1 = mu[1],
                      mu2 = mu[2], sigma = sigma[1], p = p, NFine = 5e5,
                      deltaFine = 1e-3, seed = 1434809263)

  fit2D(listData = listData, NSub = NSub, deltaSub = deltaSub, type = "mivM",
        maxK = maxK, pde = FALSE, start = "stat", sigmaFix = sigma[1],
        pFix = atanh(2 * p - 1), optMethod = optMethod)

}


### Simulation setting
dontrun <- TRUE
if (!dontrun) {

# Number of cores
nCores <- 30 # min(parallel::detectCores() - 1, 15)

# Kind of process
# p <- 1; type <- "WN"
# p <- 1; type <- "WC"
# p <- 2; type <- "WN"
p <- 2; type <- "mivM"

# Common parameters
M <- 1e3
delta <- c(0.05, 0.20, 0.5, 1)
NSub <- 250
deltaFine <- 1e-3
NFine <- NSub * max(delta) / deltaFine
pde <- FALSE
SO <- TRUE
start <- "stat"
K <- ifelse(p == 1, 5, 3)
eigTol <- 1e-3
condTol <- 1e3
sdInitial <- 0.1

}


### Simulations
dontrun <- TRUE
if (!dontrun) {

# Parameters process
if (type == "WN" & p == 1) {

  alpha <- c(0.5, 1)
  mu <- pi/2
  sigma <- 1:2
  cases <- expand.grid(alpha = alpha, mu = mu, sigma = sigma, delta = delta,
                       rep = 1:M)

} else if (type == "WC" & p == 1) {

  alpha <- c(0.5, 1)
  mu <- pi/2
  sigma <- 1:2
  cases <- expand.grid(alpha = alpha, mu = mu, sigma = sigma, delta = delta,
                       rep = 1:M)

} else if (type == "WN" & p == 2) {

  alpha1 <- 1:2
  mu1 <- pi/2
  sigma1 <- 1:2
  cases <- expand.grid(alpha1 = alpha1, mu1 = mu1, sigma1 = sigma1,
                       delta = delta, rep = 1:M)
  cases$alpha2 <- cases$alpha1
  cases$alpha3 <- cases$alpha1 / 2
  cases$mu2 <- -cases$mu1
  cases$sigma2 <- cases$sigma1
  cases <- cases[, c(1, 6:7, 2, 8, 3, 9, 4:5)]

} else if (type == "mivM" & p == 2) {

  alpha1 <- 0.75
  mu1 <- pi/2
  alpha2 <- 1.5
  mu2 <- -pi/2
  prop <- c(0.25, 0.50, 0.75)
  sigma <- 1
  cases <- expand.grid(alpha1 = alpha1, alpha2 = alpha2, mu1 = mu1, mu2 = mu2,
                       p = prop, sigma = sigma, delta = delta, rep = 1:M)

}

# Set seed
set.seed(2345678)
seed <- rep(round(runif(M, 1, 9876541)) + 91919, each = nrow(cases) / M)

# Create a cluster
cl <- makeSOCKcluster(nCores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nrow(cases), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
res <- foreach(i = 1:nrow(cases), .combine = rbind,
               .packages = c("Rcpp", "circular", "colorRamps","sdetorus"),
               .options.snow = opts) %dopar% {

  # Time
  time <- proc.time()[3]

  # Sample
  if (p == 1) {

    if (type == "WN") {

      samp <- rWn1D(alpha = cases[i, 1], mu = cases[i, 2], sigma = cases[i, 3],
                    NFine = NFine, deltaFine = deltaFine, seed = seed[i],
                    maxK = maxK)

    } else if (type == "WC") {

      samp <- rWc1D(alpha = cases[i, 1], mu = cases[i, 2], sigma = cases[i, 3],
                    NFine = NFine, deltaFine = deltaFine, seed = seed[i])

    }

  } else if (p == 2) {

    if (type == "WN") {

      samp <- rWn2D(alpha1 = cases[i, 1], alpha2 = cases[i, 2],
                    alpha3 = cases[i, 3], mu1 = cases[i, 4], mu2 = cases[i, 5],
                    sigma1 = cases[i, 6], sigma2 = cases[i, 7], NFine = NFine,
                    deltaFine = deltaFine, seed = seed[i], maxK = maxK)

    } else if (type == "mivM") {

      samp <- rMivM2D(alpha1 = cases[i, 1], alpha2 = cases[i, 2],
                      mu1 = cases[i, 3], mu2 = cases[i, 4], p = cases[i, 5],
                      sigma = cases[i, 6], NFine = NFine, deltaFine = deltaFine,
                      seed = seed[i])

    }

  }

  # Fit
  if (p == 1) {

    fit <- fit1D(listData = samp, NSub = NSub, deltaSub = cases[i, 4],
                 type = type, maxK = maxK, pde = pde, start = start,
                 optMethod = "Nelder-Mead", sigmaFix = cases[i, 3],
                 eigTol = eigTol, condTol = condTol)

  } else if (p == 2) {

    if (type == "WN") {

    fit <- fit2D(listData = samp, NSub = NSub, deltaSub = cases[i, 8],
                 type = type, maxK = maxK, pde = pde, start = start,
                 optMethod = "Nelder-Mead", sigmaFix = cases[i, 6:7],
                 eigTol = eigTol, condTol = condTol)

    } else if (type == "mivM") {

      fit <- fit2D(listData = samp, NSub = NSub, deltaSub = cases[i, 7],
                   type = type, maxK = maxK, pde = pde, SO = SO, start = start,
                   optMethod = "Nelder-Mead", sigmaFix = cases[i, 6],
                   pFix = atanh(2 * cases[i, 5] - 1), eigTol = eigTol,
                   condTol = condTol)

    }

  }
  cat(i, "\n")

  # Time
  if (i %% 500 == 0) {

    cat("Time", proc.time()[3] - time, "\n")

  }

  # Result
  c(seed = seed[i], fit)

}
res <- cbind(cases, res)
close(pb)

# Close cluster
stopCluster(cl)

# Save data
save(list = "res", file = paste(type, "-", p, "D-", ifelse(pde, "pde", "nopde"),
                                "025.RData", sep = ""))

}


### Analysis 1D
{

# Load data
# p <- 1; type <- "WN"
p <- 1; type <- "WC"
pde <- TRUE
load(file = paste(type, "-", p, "D-", ifelse(pde, "pde", "nopde"),
                  ".RData", sep = ""))

# Parameters of each scenario
indPars <- match(x = c("delta", "alpha", "mu", "sigma"), table = colnames(res))

# Columns of the estimates
ap <- c("E", "SO", "Ap", "Pde")
if (type == "WC") {

  ap <- ap[-3]

}

indAlpha <- na.omit(match(x = paste("alpha", ap, sep = ""),
                          table = colnames(res)))
indMu <- na.omit(match(x = paste("mu", ap, sep = ""), table = colnames(res)))

# Estimates
estimates <- res[, c(indPars, indAlpha, indMu)]

# Columns of the estimates in the new object
indAlpha <- na.omit(match(x = paste("alpha", ap, sep = ""),
                          table = colnames(estimates)))
indMu <- na.omit(match(x = paste("mu", ap, sep = ""),
                       table = colnames(estimates)))

# Squared errors
estimates[, indAlpha] <- sweep(estimates[, indAlpha], 1, estimates$alpha,
                               FUN = "-")^2
estimates[, indMu] <- toPiInt(sweep(estimates[, indMu], 1, estimates$mu,
                                    FUN = "-"))^2

# Percentage of NA's
nas <- aggregate(x = estimates, by = list(res$delta, res$alpha, res$mu,
                                          res$sigma),
                 FUN = function(x) mean(is.na(x)))[, -c(1:4)] * 100
colMeans(nas)

# Table for average REs
fac <- c("alpha", "sigma", "delta")
a <- aggregateRes(errors = estimates, factors = fac,
                  estimators = list(paste("alpha", ap, sep = ""),
                                    paste("mu", ap, sep = "")))
A <- matrix(a$tab, ncol = 4, byrow = TRUE)

# Total average of average REs
colMeans(a$meanErrors)

# Average REs by delta (averages above the scenarios)
t(sapply(1:4, function(k) colMeans(a$meanErrors[(4 * (k - 1) + 1):(4 * k), ])))

# LaTeX tables with average REs
cat(paste(apply(cbind(c("$0.05$", "$0.20$", "$0.50$", "$1.00$"), A[, 1:2]), 1,
                function(x) paste(x, collapse = " & ")), collapse = "\\\\\n"))
cat("\n\n")
cat(paste(apply(cbind(c("$0.05$", "$0.20$", "$0.50$", "$1.00$"), A[, 3:4]), 1,
                function(x) paste(x, collapse = " & ")), collapse = "\\\\\n"))

# # Table for MSEs
# fac <- c("alpha", "sigma", "delta")
# a1 <- aggregateSignif(errors = estimates, factors = fac,
#                       estimators = paste("alpha", ap, sep = ""))
# a2 <- aggregateSignif(errors = estimates, factors = fac,
#                       estimators = paste("mu", ap, sep = ""))
# A1 <- matrix(a1$tab, ncol = 4, byrow = TRUE)
# A2 <- matrix(a2$tab, ncol = 4, byrow = TRUE)
# A <- rbind(A1, A2)
# A <- as.matrix(A[c(rbind(1:nrow(A1), nrow(A1) + 1:nrow(A2))), ])
# cat(paste(apply(A[, 1:2], 1, function(x) paste(x, collapse = " & ")),
#           collapse = "\\\\\n"))
# cat("\n\n")
# cat(paste(apply(A[, 3:4], 1, function(x) paste(x, collapse = " & ")),
#           collapse = "\\\\\n"))

}


### Analysis 2D
{

# Load data
# p <- 2; type <- "WN"
p <- 2; type <- "mivM"
pde <- FALSE
load(file = paste(type, "-", p, "D-", ifelse(pde, "pde", "nopde"), ".RData",
                  sep = ""))

# Parameters of each scenario
if (type == "WN") {

  indPars <- match(x = c("delta", "alpha1", "alpha2", "alpha3", "mu1", "mu2",
                         "sigma1", "sigma2"), table = colnames(res))

} else if (type == "mivM") {

  indPars <- match(x = c("delta", "alpha1", "alpha2", "mu1", "mu2", "p",
                         "sigma"), table = colnames(res))

}

# Columns of the estimates
ap <- c("E", "SO", "Ap", "Pde")
if (!pde) {

  ap <- ap[-4]

}
if (type == "mivM") {

  ap <- ap[-3]

}
indAlpha1 <- na.omit(match(x = paste(c("alpha1"), ap, sep = ""),
                           table = colnames(res)))
indAlpha2 <- na.omit(match(x = paste(c("alpha2"), ap, sep = ""),
                           table = colnames(res)))
indMu1 <- na.omit(match(x = paste("mu1", ap, sep = ""),
                        table = colnames(res)))
indMu2 <- na.omit(match(x = paste("mu2", ap, sep = ""),
                        table = colnames(res)))
if (type == "WN") {

  indAlpha3 <- na.omit(match(x = paste("alpha3", ap, sep = ""),
                             table = colnames(res)))

} else if (type == "mivM") {

  indP <- na.omit(match(x = paste("p", ap, sep = ""), table = colnames(res)))

}

# Estimates
if (type == "WN") {

  estimates <- res[, c(indPars, indAlpha1, indAlpha2, indAlpha3,
                       indMu1, indMu2)]

} else if (type == "mivM") {

  estimates <- res[, c(indPars, indAlpha1, indAlpha2, indMu1, indMu2, indP)]

}

# Columns of the estimates in the new object
indAlpha1 <- na.omit(match(x = paste("alpha1", ap, sep = ""),
                           table = colnames(estimates)))
indAlpha2 <- na.omit(match(x = paste("alpha2", ap, sep = ""),
                           table = colnames(estimates)))
indMu1 <- na.omit(match(x = paste("mu1", ap, sep = ""),
                        table = colnames(estimates)))
indMu2 <- na.omit(match(x = paste("mu2", ap, sep = ""),
                        table = colnames(estimates)))
if (type == "WN") {

  indAlpha3 <- na.omit(match(x = paste("alpha3", ap, sep = ""),
                             table = colnames(estimates)))

} else if (type == "mivM") {

  indP <- na.omit(match(x = paste("p", ap, sep = ""),
                        table = colnames(estimates)))

}

# Squared errors
if (type == "WN") {

  estimates[, indAlpha1] <- sweep(estimates[, indAlpha1], 1, estimates$alpha1,
                                  FUN = "-")^2
  estimates[, indAlpha2] <- sweep(estimates[, indAlpha2], 1, estimates$alpha2,
                                  FUN = "-")^2
  estimates[, indMu1] <- toPiInt(sweep(estimates[, indMu1], 1, estimates$mu1,
                                       FUN = "-"))^2
  estimates[, indMu2] <- toPiInt(sweep(estimates[, indMu2], 1, estimates$mu2,
                                       FUN = "-"))^2
  estimates[, indAlpha3] <- toPiInt(sweep(estimates[, indAlpha3], 1,
                                          estimates$alpha3, FUN = "-"))^2

} else if (type == "mivM") {

  # Break unidentifiability
  indMin <- sapply(indMu1, function(i)
    apply(cbind(toPiInt(estimates[, i] - estimates$mu1)^2,
                       toPiInt(estimates[, i] - estimates$mu2)^2),
                 1, function(x) {w <- which.min(x); ifelse(length(w), w, NA)})
    )
  A1 <- cbind(estimates$alpha1, estimates$alpha2)
  A2 <- cbind(estimates$alpha2, estimates$alpha1)
  M1 <- cbind(estimates$mu1, estimates$mu2)
  M2 <- cbind(estimates$mu2, estimates$mu1)
  P <- cbind(estimates$p, 1 - estimates$p)
  for (i in 1:ncol(indMin)) {

    estimates[, indAlpha1[i]] <- (estimates[, indAlpha1[i]] -
                                    A1[cbind(1:nrow(indMin), indMin[, i])])^2
    estimates[, indAlpha2[i]] <- (estimates[, indAlpha2[i]] -
                                    A2[cbind(1:nrow(indMin), indMin[, i])])^2
    estimates[, indMu1[i]] <- toPiInt(estimates[, indMu1[i]] -
                                        M1[cbind(1:nrow(indMin), indMin[, i])])^2
    estimates[, indMu2[i]] <- toPiInt(estimates[, indMu2[i]] -
                                        M2[cbind(1:nrow(indMin), indMin[, i])])^2
    estimates[, indP[i]] <- ((tanh(estimates[, indP[i]]) + 1) / 2 -
                               P[cbind(1:nrow(indMin), indMin[, i])])^2

  }

}

# Percentage of NA's
if (type == "WN") {

  nas <- aggregate(x = estimates[, -c(1:8)],
                   by = list(res$delta, res$alpha1, res$alpha2, res$alpha3,
                             res$mu1, res$mu2, res$sigma1, res$sigma2),
                   FUN = function(x) mean(is.na(x)))[, -c(1:8)] * 100
  nas[, -c(1:8)] <- nas[, -c(1:8)] * 100

} else if (type == "mivM") {

  nas <- aggregate(x = estimates[, -c(1:7)],
                   by = list(res$delta, res$alpha1, res$alpha2, res$mu1,
                             res$mu2, res$p, res$sigma),
                   FUN = function(x) mean(is.na(x)))
  nas[, -c(1:7)] <- nas[, -c(1:7)] * 100

}
colMeans(nas)

# Table for average REs
if (type == "WN") {

  fac <- c("alpha1", "alpha2", "alpha3", "sigma1", "sigma2", "delta")
  listEstimators <- list(paste("alpha1", ap, sep = ""),
                         paste("alpha2", ap, sep = ""),
                         paste("alpha3", ap, sep = ""),
                         paste("mu1", ap, sep = ""),
                         paste("mu2", ap, sep = ""))

} else if (type == "mivM") {

  fac <- c("alpha1", "alpha2", "p", "sigma", "delta")
  listEstimators <- list(paste("alpha1", ap, sep = ""),
                         paste("alpha2", ap, sep = ""),
                         paste("mu1", ap, sep = ""),
                         paste("mu2", ap, sep = ""),
                         paste("p", ap, sep = ""))

}
a <- aggregateRes(errors = estimates, factors = fac,
                  estimators = listEstimators)
nScenarios <- ifelse(type == "WN", 4, 3)
A <- matrix(a$tab, ncol = nScenarios, byrow = TRUE)

# Total average of average REs
colMeans(a$meanErrors)

# Average REs by delta (averages above the scenarios)
t(sapply(1:4, function(k)
  colMeans(a$meanErrors[(nScenarios * (k - 1) + 1):(nScenarios * k), ])))

# LaTeX tables with average REs
if (type == "WN") {

  cat(paste(apply(cbind(c("$0.05$", "$0.20$", "$0.50$", "$1.00$"), A[, 1:2]), 1,
                  function(x) paste(x, collapse = " & ")), collapse = "\\\\\n"))
  cat("\n\n")
  cat(paste(apply(cbind(c("$0.05$", "$0.20$", "$0.50$", "$1.00$"), A[, 3:4]), 1,
                  function(x) paste(x, collapse = " & ")), collapse = "\\\\\n"))

} else if (type == "mivM") {

  cat(paste(apply(cbind(c("$0.05$", "$0.20$", "$0.50$", "$1.00$"), A), 1,
                  function(x) paste(x, collapse = " & ")), collapse = "\\\\\n"))

}

# # Table for MSEs
# if (type == "WN") {
#
#   fac <- c("alpha1", "alpha2", "alpha3", "sigma1", "sigma2", "delta")
#   A1 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("alpha1", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#   A2 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("alpha2", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#   A3 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("alpha3", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#   A4 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("mu1", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#   A5 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("mu2", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#
# } else if (type == "mivM") {
#
#   fac <- c("alpha1", "alpha2", "sigma", "delta")
#   A1 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("alpha1", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#   A2 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("alpha2", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#   A3 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("mu1", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#   A4 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("mu2", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#   A5 <- matrix(aggregateSignif(errors = estimates, factors = fac,
#                                estimators = paste("p", ap, sep = ""))$tab,
#                ncol = 1, byrow = TRUE)
#
# }
# A <- rbind(A1, A2, A3, A4, A5)
# A <- as.matrix(A[c(rbind(1:nrow(A1), nrow(A1) + 1:nrow(A2),
#                          2 * nrow(A1) + 1:nrow(A2),
#                          3 * nrow(A1) + 1:nrow(A2),
#                          4 * nrow(A1) + 1:nrow(A2))), ])
# if (type == "WN") {
#
#   cat(paste(apply(A[, 1:2], 1, function(x) paste(x, collapse = " & ")),
#             collapse = "\\\\\n"))
#   cat(paste(apply(A[, 3:4], 1, function(x) paste(x, collapse = " & ")),
#             collapse = "\\\\\n"))
#
# } else if (type == "mivM") {
#
#   cat(paste(apply(cbind(c("$0.05$", "$0.20$", "$0.50$", "$1.00$"), A), 1,
#                   function(x) paste(x, collapse = " & ")), collapse = "\\\\\n"))
#
# }

}


### Visualization tpds
{

library(manipulate)

## WC

# tpd
Mx <- 500
Mt <- 200
x <- seq(-pi, pi, l = M + 1)[-c(Mx + 1)]
times <- seq(0, 1, l = Mt + 1)
u0 <- dWn1D(x, mu = pi/2, sigma = 0.1)
b <- driftJp(x, alpha = 1, mu = 0, psi = -1)
sigma2 <- rep(2^2, M)
u <- crankNicolson1D(u0 = cbind(u0), b = b, sigma2 = sigma2, N = 0:Mt,
                     deltat = 1 / Mt, Mx = Mx, deltax = 2 * pi / Mx)

# Visualization of tpd
plotSurface2D(times, x, z = t(u), levels = seq(0, 3, l = 50))

# Time cuts
manipulate(plot(x, u[, i], type = "l", ylab = "tpd",
                main = sprintf("t = %.3f", (i - 1)/200),
                ylim = c(0, 1)), i = slider(1, 201))

## mivM

# tpd
Mx <- 200
My <- 200
Mt <- 200
x <- seq(-pi, pi, l = Mx + 1)[-c(Mx + 1)]
y <- seq(-pi, pi, l = My + 1)[-c(My + 1)]
u0 <- c(outer(dWn1D(x, mu = 0, sigma = 0.1), dWn1D(y, mu = pi, sigma = 0.1)))
p <- 0.25
b <- driftMixIndVm(x = expand.grid(x, y),
                   A = rbind(rep(1, 2), rep(2, 2)) * 3/4,
                   M = rbind(rep(1, 2), rep(-1, 2)) * pi/2, sigma = 1,
                   p = c(p, 1 - p))
bx <- matrix(b[, 1], nrow = Mx, ncol = My)
by <- matrix(b[, 2], nrow = Mx, ncol = My)
sigma2 <- matrix(1, nrow = Mx, ncol = My)
sigmaxy <- matrix(0, nrow = Mx, ncol = My)
u <- crankNicolson2D(u0 = cbind(u0), bx = bx, by = by, sigma2x = sigma2,
                     sigma2y = sigma2, sigmaxy = sigmaxy,
                     N = 0:Mt, deltat = 1 / Mt, Mx = Mx, deltax = 2 * pi / Mx,
                     My = My, deltay = 2 * pi / My)

# Time cuts
manipulate({
  plotSurface2D(x, y, z = matrix(u[, j + 1], Mx, My),
                main = sprintf("t = %.3f", (i - 1)/200),
                levels = seq(0, 1, l = 50))
}, j = slider(0, Mt))

}
