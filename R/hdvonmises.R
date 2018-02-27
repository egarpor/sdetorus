
dBvm <- function(x, mu, kappa, logConst) {
  
  # Convert to matrices if n == 1
  if (!is.matrix(x)) {
    
    x <- matrix(x, nrow = 1)
    
  }
  
  # Normalizing constant
  if (missing(logConst)) {
    
    logConst <- log(constBvm(kappa = kappa))
    
  }
  
  # Differences
  thMu1 <- x[, 1] - mu[1]
  thMu2 <- x[, 2] - mu[2]
  
  # Density
  exp(-logConst + kappa[1] * cos(thMu1) + kappa[2] * cos(thMu2) +
        kappa[3] * sin(thMu1) * sin(thMu2) / 2)
  
}

constBvm <- function(M = 25, kappa) {
  
  # Truncation series
  m <- 0:M
  
  if (all(kappa > 0)) {
    
    const <- 4 * pi^2 * sum(exp(
      lchoose(n = 2 * m, k = m) +
        m * (2 * (log(kappa[3]) - log(4)) - log(kappa[1]) - log(kappa[2])) +
        kappa[1] + log(besselI(x = kappa[1], nu = m, expon.scaled = TRUE)) +
        kappa[2] + log(besselI(x = kappa[2], nu = m, expon.scaled = TRUE))
    ))
    
  } else if (kappa[3] < 1e-15) {
    
    const <- 4 * pi^2 * besselI(x = kappa[1], nu = 0) * besselI(x = kappa[2], nu = 0)
    
  } else {
    
    warning("No series expansion for kappa1 == 0 or kappa2 == 0 and kappa3 != 0. This is a bimodal density!")
    const <- 1 / mcTorusIntegrate(f = function(x)
      dBvm(x = x, mu = c(0, 0), kappa = kappa, logConst = 0), p = 2, M = 1e3,
      fVect = TRUE, seed = 3456789)
    
  }
  
  return(const)
  
}

dVmf <- function(x, M, K, alpha, besselInterp = FALSE) {
  
  # Convert to matrices
  if (!is.matrix(M)) {
    
    M <- matrix(M, nrow = 1, byrow = TRUE)
    
  }
  if (!is.matrix(K)) {
    
    K <- matrix(K, nrow = 1, byrow = TRUE)
    
  }
  if (!is.matrix(x)) {
    
    x <- matrix(x, ncol = 1)
    
  }
  if (missing(alpha)) {
    
    alpha <- 1
    
  }
  
  
  # Check if x, M, K and alpha agree
  dM <- dim(M)
  dK <- dim(K)
  if (!(all(dM == dK) & (dM[2] == ncol(x)) & (dM[1] == length(alpha)))) {
    
    stop("Incompatible sizes of x, M, K and alpha")
    
  }
  
  # Density
  drop(.dVmfCpp(x = x, K = K, M = M, alpha = alpha, besselInterp = besselInterp))
  
}

emVmf <- function(data, k, M, K, alpha, tol = c(0.001, 0.001, 0.001 / k),
                  kappaMax = 500, maxIter = 100, isotropic = FALSE,
                  besselInterp = FALSE, verbose = 0) {
  
  # Sample size
  n <- nrow(data)
  
  # Dimension
  p <- ncol(data)
  
  # Number of parameters: M and K are k * p, alpha are k - 1
  npar <- 2 * k * p + k - 1
  if (npar > n) {
    
    warning(paste("Trying to fit more parameters than sample size.
                  EM will overfit and likely break down due to high values of kappa.
                  Maximum number of mixtures allowed is k =",
                  ceiling((n + 1) / (2 * p + 1)), sep = ""), immediate. = TRUE)
    
  }
  
  # Initialize means, concentrations and proportions
  if (missing(M)) {
    
    M <- matrix(runif(k * p, -pi, pi), nrow = k, ncol = p, byrow = TRUE)
    
  }
  if (missing(K)) {
    
    if (isotropic) {
      
      K <- matrix(runif(k * p, 0, 1), nrow = k, ncol = p, byrow = TRUE)
      
    } else {
      
      K <- matrix(runif(k, 0, 1), nrow = k, ncol = p, byrow = TRUE)
      
    }
    
  }
  if (missing(alpha)) {
    
    alpha <- runif(k)
    alpha <- alpha / sum(alpha)
    
  }
  MNew <- M
  KNew <- K
  alphaNew <- alpha
  
  # The evaluations of f in the data and in each component
  fih <- matrix(0, nrow = n, ncol = k)
  
  # The probability of the hidden states
  pih <- fih
  
  # Positions instead of angles
  cosData <- cos(data)
  sinData <- sin(data)
  
  # Factor vMF density
  l2pi <- -p * log(2 * pi)
  
  # EM algorithm
  j <- 0
  repeat {
    
    # Time
    t0 <- proc.time()[3]
    
    ## Expectation step
    
    pih <- .clusterProbsVmf(cosData = cosData, sinData = sinData, K = K, M = M,
                            alpha = alpha, l2pi = l2pi, besselInterp = besselInterp)
    
    ## Maximization step
    
    # Weights mixtures
    alphaNew <- colMeans(pih)
    
    # Circular means and concentrations
    MK <- .weightedMuKappa(cosData = cosData, sinData = sinData, weights = pih,
                           isotropic = isotropic)
    
    # MNew and KNew
    MNew <- MK$M
    KNew <- MK$K
    
    # Increments
    deltaM <- abs(MNew - M)
    deltaM <- max(pmin(deltaM, abs(2 * pi - deltaM)))
    deltaK <- max(abs(KNew - K))
    deltaAlpha <- max(abs(alphaNew - alpha))
    j <- j + 1
    
    # Display progress
    if (verbose > 0) {
      
      # Increments and time
      cat(paste("Iteration: ", j, ". tolM: ", round(deltaM, 4), ". tolK: ",
                round(deltaK, 4), ". tolAlpha: ", round(deltaAlpha, 4),
                ". Time: ", round(proc.time()[3] - t0, 4), ".\n", sep = ""))
      
      if (verbose > 1) {
        
        # Values of parameters
        cat("M:", round(MNew, 3), "\nK:", round(KNew, 3), "\na:",
            round(alphaNew, 3), "\n\n")
        
        if (verbose > 2 & p < 10 & j %% 5 == 0) {
          
          # Plot assignments
          pairs(data, col = apply(pih, 1, which.max),
                main = paste("Assignments at iteration", j), lower.panel = NULL)
          
        }
        
      }
      
    }
    
    # Check convergence
    if (deltaM < tol[1] & deltaK < tol[2] & deltaAlpha < tol[3]) {
      
      # Convergence
      if (verbose) {
        
        cat("Succesful convergence\n")
        
      }
      convergence <- TRUE
      break
      
    } else if (j > maxIter) {
      
      # No convergence
      if (verbose) {
        
        cat("EM algorithm with k =", k, "components did not converge\n")
        
      }
      convergence <- FALSE
      break
      
    } else {
      
      # Update parameters
      M <- MNew
      K <- KNew
      alpha <- alphaNew
      
    }
    
  }
  
  # Loglikelihood and result
  logLik <- sum(log(dVmf(x = data, M = M, K = K, alpha = alpha,
                         besselInterp = besselInterp)))
  return(list(M = M, K = K, alpha = alpha, pih = pih,
              BIC = -2 * logLik + log(n) * npar, convergence = convergence))
  
}

mcTorusIntegrate <- function (f, p, M = 1e5, seed, fVect = TRUE, ...) {
  
  # Set seed if present
  if (!missing(seed)) {
    
    old <- .Random.seed
    set.seed(seed)
    on.exit({.Random.seed <<- old})
    
  }
  
  # Sample uniformly on the torus
  sampleUnifTorus <- matrix(runif(n = M * p, min = -pi, max = pi), nrow = M, ncol = p)
  
  # Evaluations
  if (fVect) {
    
    fEvals <- f(sampleUnifTorus, ...)
    
  } else {
    
    fEvals <- apply(sampleUnifTorus, 1, f, ...)
    
  }
  
  # Torus area
  area <- (2 * pi)^p
  
  return(mean(fEvals) * area)
  
}
