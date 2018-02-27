
# Sort objects by memory size, expressed in Gb
sizeObjectsSession <- function(k = 10) {

  sort(x = sapply(ls(envir = .GlobalEnv), function(x) object.size(get(x))),
       decreasing = TRUE)[1:k] / 1e9

}


## Computations 1D
{

# Load packages
library(sdetorus)
library(colorRamps)
library(doSNOW)
library(abind)
library(ggplot2)
library(reshape2)
library(Matrix)
packages <- c("Rcpp", "mvtnorm", "foreach", "colorRamps", "sdetorus",
              "abind", "Matrix")

# Type of process
# type <- "WN"
type <- "vM"

# Parameters for accuray of PDE solution. Accuracy is needed because we work
# with logs. For example, with a coarser grid, SoVm seemed to outperform So in
# KL divergence because of the inaccuracy of the PDE solution. The accuracy is
# primarly affected by the state grid if the initial condition is not very
# sharp.
Mx <- 3000
Mt <- 1500
sdInitial <- 0.1
xFullGrid <- seq(-pi, pi, l = Mx + 1)[-(Mx + 1)]

# State subgrids for approximating the integrals
MxSub <- Mx / 3
indSubGrid <- seq.int(1, Mx, by = round(Mx / MxSub))
xSubGrid <- c(seq(-pi, pi, l = Mx + 1)[-(Mx + 1)])[indSubGrid]

# Times for measurement of divergences
tx <- c(seq(0.05, 0.30, by = 0.05), seq(0.4, 1.5, by = 0.1), 2, 2.5, 3)
indtx <- round(tx * Mt / max(tx))

# Tolerance and number of winding numbers used (2 * K + 1)
tol <- 1e-10
xMax <- sqrt(.Machine$double.xmax)
maxK <- 3

# Control over PDE
loadPde <- FALSE
onlyPde <- FALSE & !loadPde
savePde <- FALSE & !loadPde

# Smoothing matrix required to properly compare approximations and PDE solution
s <- dWn1D(x = xSubGrid, mu = 0, sigma = sdInitial, maxK = maxK)
s[s < tol] <- 0
s <- (s / periodicTrapRule1D(s)) / MxSub * 2 * pi
S <- matrix(0, MxSub, MxSub)
ind0 <- 1:MxSub + MxSub/2 + 1
for (i in 1:MxSub) {

  S[i, ] <- s[toInt(ind0 - i, 1, MxSub + 1)]

}

# Declare it as sparse
S <- Matrix(S, sparse = TRUE)

# Array of parameters
arrayPars <- rbind(c(0.5, pi/2, 0.5),
                   c(1, pi/2, 0.5),
                   c(2, pi/2, 0.5),
                   c(0.5, pi/2, 1),
                   c(1, pi/2, 1),
                   c(2, pi/2, 1),
                   c(0.5, pi/2, 2),
                   c(1, pi/2, 2),
                   c(2, pi/2, 2))

# Drifts
if (type == "WN") {

  b <- function(x, pars) driftWn1D(x = x, alpha = pars[1], mu = pars[2],
                                   sigma = pars[3], maxK = maxK)
  b1 <- function(x, pars, h = 1e-5) {
    l <- length(x)
    res <- b(x = c(x + h, x - h), pars = pars)
    drop(res[1:l] - res[(l + 1):(2 * l)])/(2 * h)
  }
  b2 <- function(x, pars, h = 1e-5) {
    l <- length(x)
    res <- b(x = c(x + h, x, x - h), pars = pars)
    drop(res[1:l] - 2 * res[(l + 1):(2 * l)] + res[(2 * l + 1):(3 * l)])/(h^2)
  }
  sigma2 <- function(x, pars) rep(pars[3]^2, length(x))

} else if (type == "vM") {

  b <- function(x, pars) pars[1] * sin(pars[2] - x)
  b1 <- function(x, pars) -pars[1] * cos(pars[2] - x)
  b2 <- function(x, pars, h = 1e-5) -pars[1] * sin(pars[2] - x)
  sigma2 <- function(x, pars) rep(pars[3]^2, length(x))

}

# Parallel execution. We use this strategy:
# - Automatic determination of number of chunks and cores is obtained from
#   a maximum number of cores. This maximizes the number of jobs per node.
# - Sequential loop over parameters, computation of drifts and stationary
#   densities.
# - Parallel loop for computing PDE solution for the full trajectory from t = 0
#   to t = max(tx) and for different starting values (this is the part
#   parallelized). The parallel job is broken into into chunks of tasks that are
#   equally assigned  to each core.
# - Store a subset of the solution, which is used for approximating the
#   integrals, as an array (subgrid, times, subinitialgrid).
# - Parallel loop (can be turned into sequential) in the times for computing the
#   smoothed approximations. We access there the solution to the PDE.
# - Sequential loops on the starting values for computing the approximations.
#   (can be turned into chunked parallel loops)
maxCores <- 30 # min(12, detectCores() - 1)
chunkSize <- 50 # MxSub %/% rev(which((MxSub %% 1:maxCores) == 0))[1]
nTasks <- MxSub %/% chunkSize
nCores <- maxCores
cl <- makeSOCKcluster(maxCores)
registerDoSNOW(cl)
res1D <- array(dim = c(nrow(arrayPars), length(tx), 8))

# Sequential loop over parameters
for (l in 1:nrow(arrayPars)) {

  # Selection of parameters
  pars <- arrayPars[l, ]

  # Drifts for PDE
  if (type == "WN") {

    bGrid <- driftWn1D(x = xFullGrid, alpha = pars[1], mu = pars[2],
                       sigma = pars[3], maxK = maxK)
    sigma2Grid <- rep(pars[3]^2, Mx)

  } else if (type == "vM") {

    bGrid <- pars[1] * sin(pars[2] - xFullGrid)
    sigma2Grid <- rep(pars[3]^2, Mx)

  }

  # Stationary density
  if (type == "WN") {

    dStat <- dWn1D(x = xSubGrid, mu = pars[2],
                   sigma = pars[3] / sqrt(2 * pars[1]), maxK = maxK)

  } else {

    dStat <- dVm(x = xSubGrid, mu = pars[2], kappa = 2 * pars[1] / (pars[3]^2))

  }
  dStat[dStat > xMax] <- xMax
  dStat[dStat < tol | !is.finite(dStat)] <- tol
  dStat <- dStat / periodicTrapRule1D(fx = dStat)

  # Progress
  cat("\n PDE for", l, "/", nrow(arrayPars), "\n")

  # Load PDE?
  if (loadPde) {

    if (exists("p")) {

      rm(p)
      gc()

    }
    load(file = paste("pde-", type, "-", l, "-1D.RData", sep = ""))

  } else {

    # Parallel loop for starting values with task chunking
    # Result stored as sparse matrix (initial + evaluation, ti)
    pb <- txtProgressBar(max = nTasks, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    p <- foreach(k = 1:nTasks, .combine = rBind, .packages = packages,
                 .noexport = "S", .options.snow = opts) %dopar% {

      # Chunk of tasks
      foreach(iChunk = 1:chunkSize, .combine = rBind) %do% {

        # Index in the complete task
        i <- (k - 1) * chunkSize + iChunk

        # Solve PDE from initial condition and get subsolution
        u0 <- matrix(dWn1D(x = xFullGrid, mu = xSubGrid[i], sigma = sdInitial,
                           maxK = maxK), ncol = 1)
        u0 <- u0 / periodicTrapRule1D(fx = u0)
        q <- crankNicolson1D(u0 = u0, b = bGrid, sigma2 = sigma2Grid, Mx = Mx,
                             deltax = xFullGrid[2] - xFullGrid[1], N = indtx,
                             deltat = max(tx) / Mt)[indSubGrid, ]

        # Result as a sparse matrix to save space - normalization comes in
        # later!
        q[q < tol] <- 0
        Matrix(q, sparse = TRUE)

      }

    }
    close(pb)

    # Save PDE?
    if (savePde) {

      save(list = "p", file = paste("pde-", type, "-", l, "-1D.RData",
                                    sep = ""))

    }

    # Free memory
    gc()

  }

  # Compute approximations?
  if (!onlyPde) {

    # Progress
    cat("\n Approxs for", l, "/", nrow(arrayPars), "\n")

    # Parallel/sequential loop over times
    pb <- txtProgressBar(max = length(tx), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    res1D[l, , ] <- foreach(ti = 1:length(tx), .combine = rbind,
                            .packages = packages,
                            .options.snow = opts) %dopar% {

      # Progress
      cat("\n Time", ti, "/", length(tx))

      # Time
      t <- tx[ti]

      # Only the i-th time and as a dense matrix (eval, intial)
      pti <- matrix(p[, ti], nrow = MxSub, ncol = MxSub)

      # Normalize PDE: rounding to tolerance, set zeros to tolerance and ensure
      # it is a proper density
      if (any(pti < -tol)) {

        cat(paste("PDE negative values in i =", i, "with val =",
                  signif(min(pti), 2), "\n"))

      }
      pti <- round(pti, digits = -log10(tol))
      pti[pti > xMax] <- xMax
      pti[pti < tol | !is.finite(pti)] <- tol
      pti <- sweep(pti, 2, apply(pti, 2,
                                 function(x) periodicTrapRule1D(fx = x)), "/")

      # Stationary
      klS <- repCol(dStat, MxSub)
      klS <- apply(log(pti / klS) * pti, 2, periodicTrapRule1D) * dStat

      # Parallel KL computation for pseudo-likelihoods and approximation (with
      # smooth initial condition)
      kl <- function(method, circular, vmApprox) {

        # Progress
        cat("\n", paste(ifelse(circular, "", "U"), method,
                        ifelse(vmApprox, "vM", ""), sep = ""), "\n")

        # Pseudo-likelihood with task chunking (export objects not defined
        # inside kl function)
        ap <- foreach(k = 1:nTasks, .combine = cBind,
                      .export = c("xSubGrid", "t", "b", "b1", "sigma2", "pars",
                                  "maxK", "tol", "xMax", "chunkSize", "MxSub",
                                  "nTasks"), .packages = packages) %do% {

          # Chunk of tasks
          foreach(iChunk = 1:chunkSize, .combine = cbind) %do% {

            # Index in the complete task
            i <- (k - 1) * chunkSize + iChunk

            # Result
            if (method == "WOU") {

              q <- tryCatch(dTpdWou1D(x = xSubGrid, x0 = rep(xSubGrid[i], MxSub),
                                     t = t, alpha = pars[1], mu = pars[2],
                                     sigma = pars[3], maxK = maxK),
                            error = function(e) {
                message(paste("Error in method = ", method, ", circular = ",
                              circular, ", vmApprox = ", vmApprox, sep = ""))
                rep(tol, MxSub)
              })

            } else {

              q <- tryCatch(dPsTpd(x = xSubGrid, x0 = xSubGrid[i], t = t, b = b,
                                   b1 = b1, sigma2 = sigma2, method = method,
                                   pars = pars, maxK = maxK, circular = circular,
                                   vmApprox = vmApprox), error = function(e) {
                message(paste("Error in method = ", method, ", circular = ",
                              circular, ", vmApprox = ", vmApprox, sep = ""))
                rep(tol, MxSub)
              })

            }

            # Normalize
            q[q > xMax] <- xMax
            q[q < tol | !is.finite(q)] <- 0
            Matrix(q, sparse = TRUE)

          }

        }

        # Smoothing
        ap <- ap %*% S

        # Normalize solution: rounding to tolerance, set zeros to tolerance and
        # ensure it is a proper density
        ap <- round(ap, digits = -log10(tol))
        ap[ap < tol] <- tol
        ap <- sweep(ap, 2, apply(ap, 2,
                                 function(x) periodicTrapRule1D(fx = x)), "/")

        # KL
        apply(log(pti / ap) * pti, 2, periodicTrapRule1D) * dStat

      }

      # Pseudo-likelihoods
      klEu <- kl(method = "E", circular = TRUE, vmApprox = FALSE)
      klSo <- kl(method = "SO", circular = TRUE, vmApprox = FALSE)
      klUEu <- kl(method = "E", circular = FALSE, vmApprox = FALSE)
      klUSo <- kl(method = "SO", circular = FALSE, vmApprox = FALSE)
      klEuVm <- kl(method = "E", circular = TRUE, vmApprox = TRUE)
      klSoVm <- kl(method = "SO", circular = TRUE, vmApprox = TRUE)
      klAp <- kl(method = "WOU", circular = TRUE, vmApprox = FALSE)

      # Expectation
      eKlS <- periodicTrapRule1D(fx = klS)
      eKlEu <- periodicTrapRule1D(fx = klEu)
      eKlSo <- periodicTrapRule1D(fx = klSo)
      eKlUEu <- periodicTrapRule1D(fx = klUEu)
      eKlUSo <- periodicTrapRule1D(fx = klUSo)
      eKlEuVm <- periodicTrapRule1D(fx = klEuVm)
      eKlSoVm <- periodicTrapRule1D(fx = klSoVm)
      eKlAp <- periodicTrapRule1D(fx = klAp)

      # Result
      c(eKlS, eKlEu, eKlSo, eKlUEu, eKlUSo, eKlEuVm, eKlSoVm, eKlAp)

    }
    close(pb)

    # Save data
    save(list = "res1D", file = paste("kl-tpd-", type, "-1D.RData", sep = ""))

  }

}

# Save data
save(list = "res1D", file = paste("kl-tpd-", type, "-1D.RData", sep = ""))

# Stop cluster
stopCluster(cl)

}


## Analysis 1D
{

# Load packages
library(ggplot2)
library(reshape2)

# Type of process
# type <- "WN"
type <- "vM"

# Load plot
load(file = paste("kl-tpd-", type, "-1D.RData", sep = ""))

# Format array into a data frame
dfRes1D <- melt(res1D)
colnames(dfRes1D) <- c("model", "t", "Approx", "KL")
dfRes1D$model <- factor(dfRes1D$model, levels = 1:9)
dfRes1D$t <- c(seq(0.05, 0.30, by = 0.05), seq(0.4, 1.5, by = 0.1),
               2, 2.5, 3)[dfRes1D$t]
dfRes1D$Approx <- factor(dfRes1D$Approx, levels = 1:8,
                         labels = c("S", "E", "SO", "UE", "USO", "EvM",
                                    "SOvM", "WOU"))

# Captions
labels <- c(
  "1" = "list(alpha == 0.5, sigma == 0.5)",
  "2" = "list(alpha == 1, sigma == 0.5)",
  "3" = "list(alpha == 2, sigma == 0.5)",
  "4" = "list(alpha == 0.5, sigma == 1)",
  "5" = "list(alpha == 1, sigma == 1)",
  "6" = "list(alpha == 2, sigma == 1)",
  "7" = "list(alpha == 0.5, sigma == 2)",
  "8" = "list(alpha == 1, sigma == 2)",
  "9" = "list(alpha == 2, sigma == 2)"
)

# Recreate ggplot2's colors manually
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Manual legends
breaks <- c("S", "E", "UE", "EvM", "SO", "USO", "SOvM", "WOU")
colours <- gg_color_hue(4)[c(1, rep(2, 3), rep(3, 3), 4)]
shapes <- c(16, rep(c(16, 1, 6), 2), 16)
names(colours) <- names(shapes) <- breaks

# Beautiful plot
s <- 12
m <- ggplot(data = dfRes1D, aes(x = t, y = KL, colour = Approx, shape = Approx))
m <- m + geom_line(alpha = 0.75) + geom_point(alpha = 0.75, size = 2) +
  scale_y_log10() +  facet_wrap(~ model, ncol = 3, nrow = 3,
                                labeller = as_labeller(labels, label_parsed)) +
  scale_color_manual(name = "Approx.", labels = breaks, breaks = breaks,
                     values = colours) +
  scale_shape_manual(name = "Approx.", labels = breaks, breaks = breaks,
                     values = shapes) +
  theme(aspect.ratio = 1, strip.text = element_text(size = s),
        legend.text = element_text(size = s),
        legend.title = element_text(size = s + 2),
        axis.text.x = element_text(size = s),
        axis.text.y = element_text(size = s),
        axis.title.x = element_text(size = s),
        axis.title.y = element_text(size = s))

m

# Save
ggsave(filename = paste("kl-tpd-", type, "-1D.pdf", sep = ""), plot = m,
       width = 10.5, height = 10.5)

}


## Computations 2D
{

# Load packages
library(sdetorus)
library(colorRamps)
library(doSNOW)
library(abind)
library(ggplot2)
library(reshape2)
library(Matrix)
packages <- c("Rcpp", "mvtnorm", "foreach", "colorRamps", "sdetorus",
              "abind", "Matrix")

# Type of process
# type <- "WN"
type <- "vM"

# Parameters for accuray of PDE solution
Mx <- 240 # If Mx^2 > 32767 apparently this slows down cluster computations
Mt <- 1500
sdInitial <- 0.1
xFullGrid <- seq(-pi, pi, l = Mx + 1)[-(Mx + 1)]
gxFullGrid <- as.matrix(expand.grid(xFullGrid, xFullGrid))

# State subgrids for approximating the integrals
MxSub <- Mx / 3
indAux <- rep(FALSE, Mx)
indAux[seq.int(1, Mx, by = round(Mx / MxSub))] <- TRUE
indSubGrid <- which(apply(as.matrix(expand.grid(indAux, indAux)), 1, all))
gxSubGrid <- gxFullGrid[indSubGrid, ]
x0SubGrid <- gxSubGrid
nx0SubGrid <- nrow(x0SubGrid)

# Times for measurement of divergences
tx <- c(seq(0.05, 0.30, by = 0.05), seq(0.4, 1.5, by = 0.1), 2, 2.5, 3)
indtx <- round(tx * Mt / max(tx))

# Tolerance and number of winding numbers used (2 * K + 1)
tol <- 1e-10
xMax <- sqrt(.Machine$double.xmax)
maxK <- 2

# Control over PDE
loadPde <- FALSE
onlyPde <- FALSE & !loadPde
savePde <- FALSE & !loadPde

# Drifts
if (type == "WN") {

  b <- function(x, pars) {
    driftWn2D(x = x, A = alphaToA(alpha = pars[1:3], sigma = pars[6:7]),
              mu = pars[4:5], sigma = pars[6:7], maxK = maxK)
  }
  jac.b <- function(x, pars, h = 1e-5) {
    l <- nrow(x)
    res <- b(x = rbind(cbind(x[, 1] + h, x[, 2]),
                       cbind(x[, 1] - h, x[, 2]),
                       cbind(x[, 1], x[, 2] + h),
                       cbind(x[, 1], x[, 2] - h)),
             pars = pars)
    cbind(res[1:l, ] - res[(l + 1):(2 * l), ],
          res[2 * l + 1:l, ] - res[2 * l + (l + 1):(2 * l), ]) / (2 * h)
  }
  sigma2 <- function(x, pars) matrix(pars[6:7]^2, nrow = length(x) / 2L, ncol = 2)

} else if (type == "vM") {

  b <- function(x, pars) {

    s1 <- sin(pars[4] - x[, 1])
    s2 <- sin(pars[5] - x[, 2])
    c1 <- cos(pars[4] - x[, 1])
    c2 <- cos(pars[5] - x[, 2])

    cbind(pars[1] * s1 - pars[3] * c1 * s2,
          pars[2] * s2 - pars[3] * s1 * c2)

  }
  jac.b <- function(x, pars) {

    s1 <- sin(pars[4] - x[1])
    s2 <- sin(pars[5] - x[2])
    c1 <- cos(pars[4] - x[1])
    c2 <- cos(pars[5] - x[2])

    matrix(c(-pars[1] * c1 - pars[3] * s1 * s2,  pars[3] * c1 * c2,
              pars[3] * c1 * c2, -pars[2] * c2 - pars[3] * s1 * s2),
           nrow = 2, ncol = 2, byrow = TRUE)

  }
  sigma2 <- function(x, pars) matrix(mean(pars[6:7])^2, nrow = length(x) / 2L,
                                     ncol = 2)

}

# Array of parameters
arrayPars <- rbind(c(0.5, 0.5, 0.25, pi/2, pi/2, 0.5, 0.5),
                   c(1, 1, 0.5, pi/2, pi/2, 0.5, 0.5),
                   c(2, 2, 1, pi/2, pi/2, 0.5, 0.5),
                   c(0.5, 0.5, 0.25, pi/2, pi/2, 0.75, 1.25),
                   c(1, 1, 0.5, pi/2, pi/2, 0.75, 1.25),
                   c(2, 2, 1, pi/2, pi/2, 0.75, 1.25),
                   c(0.5, 0.5, 0.25, pi/2, pi/2, 2, 2),
                   c(1, 1, 0.5, pi/2, pi/2, 2, 2),
                   c(2, 2, 1, pi/2, pi/2, 2, 2))

# Parallel execution. We can not use the strategy of the 1D case because
# computing the full PDE solution results in an object too large to be shared
# across clusters. We then use the following strategy:
# - Two sequential parallel loops, the first in the parameters and the second
#   in times.
# - The ouputs from the nested loops are merged in an array by columns and
#   slices, giving the same structure as in the 1D case.
# - Drifts, stationary distribution are computed inside the innermost loop,
#   redoing unnecesarily some computations for the scope of loop independence.
#   This gives a minor loss of performance.
# - The computation of the PDE and approximations is parallelized for different
#   starting points as an approach to achieve the best possible parallelization
#   with a number of cores available larger than length(tx).
maxCores <- 10 # min(20, detectCores() - 1) # 6 #detectCores() - 1 # 12
chunkSize <- 40 # (nx0SubGrid %/% rev(which((nx0SubGrid %% 1:maxCores) == 0))[1])
nTasks <- nx0SubGrid %/% chunkSize
nCores <- maxCores
cl <- makeSOCKcluster(maxCores)
registerDoSNOW(cl)
res2D <- array(dim = c(nrow(arrayPars), length(tx), 8))

# Construct smoothing matrix
S <- matrix(0, nx0SubGrid, nx0SubGrid)
pb <- txtProgressBar(max = nx0SubGrid, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
S <- foreach(i = 1:nx0SubGrid, .combine = rbind, .packages = packages,
             .options.snow = opts) %dopar% {
  s <- dWn1D(x = gxSubGrid[i, 1] - gxSubGrid[, 1], mu = 0, sigma = sdInitial) *
    dWn1D(x = gxSubGrid[i, 2] - gxSubGrid[, 2], mu = 0, sigma = sdInitial)
  s[s < tol] <- 0
  s
}
close(pb)
S <- S / apply(S, 1, periodicTrapRule2D)
S <- S / nx0SubGrid * 4 * pi^2

# Declare it as sparse
S <- Matrix(S, sparse = TRUE)

# Sequential loop over parameters
for (l in 1:nrow(arrayPars)) {

  # Selection of parameters
  pars <- arrayPars[l, ]

  # Drifts for PDE
  if (type == "WN") {

    bGrid <- driftWn2D(x = gxFullGrid, A = alphaToA(alpha = pars[1:3],
                                                    sigma = pars[6:7]),
                       mu = pars[4:5], sigma = pars[6:7])
    bxGrid <- matrix(bGrid[, 1], nrow = Mx, ncol = Mx)
    byGrid <- matrix(bGrid[, 2], nrow = Mx, ncol = Mx)
    sigma2xGrid <- matrix(pars[6]^2, nrow = Mx, ncol = Mx)
    sigma2yGrid <- matrix(pars[7]^2, nrow = Mx, ncol = Mx)
    sigmaxyGrid <- matrix(0, nrow = Mx, ncol = Mx)

  } else if (type == "vM") {

    bxGrid <- matrix(pars[1] * sin(pars[4] - gxFullGrid[, 1]) -
                       pars[3] * cos(pars[4] - gxFullGrid[, 1]) *
                       sin(pars[5] - gxFullGrid[, 2]), nrow = Mx, ncol = Mx)
    byGrid <- matrix(pars[2] * sin(pars[5] - gxFullGrid[, 2]) -
                       pars[3] * sin(pars[4] - gxFullGrid[, 1]) *
                       cos(pars[5] - gxFullGrid[, 2]), nrow = Mx, ncol = Mx)
    sigma2xGrid <- matrix(mean(pars[6:7])^2, nrow = Mx, ncol = Mx)
    sigma2yGrid <- sigma2xGrid
    sigmaxyGrid <- matrix(0, nrow = Mx, ncol = Mx)

  }

  # Stationary density
  if (type == "WN") {

    dStat <- drop(dStatWn2D(x = x0SubGrid, alpha = pars[1:3], mu = pars[4:5],
                            sigma = pars[6:7]))

  } else if (type == "vM") {

    # dBvm employs kappa[3] = alpha[3] / 2
    dStat <- drop(sdetorus:::dBvm(x = x0SubGrid, kappa = c(2, 2, 4) *
                                     pars[1:3] / (mean(pars[6:7])^2),
                                   mu = pars[4:5]))

  }
  dStat[dStat > xMax] <- xMax
  dStat[dStat < tol | !is.finite(dStat)] <- tol
  dStat <- dStat / periodicTrapRule2D(fxy = dStat)

  # Progress
  cat("\n PDE for", l, "/", nrow(arrayPars), "\n")

  # Load PDE?
  if (loadPde) {

    if (exists("p")) {

      rm(p)
      gc()

    }
    load(file = paste("pde-", type, "-", l, "-2D.RData", sep = ""))

  } else {

    # Parallel loop for PDE for all times with task chunking
    # Result stored as sparse matrix (initial + evaluation, ti)
    pb <- txtProgressBar(max = nTasks, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    p <- foreach(k = 1:nTasks, .combine = rBind, .packages = packages,
                 .noexport = "S", .options.snow = opts) %dopar% {

      # Chunk of tasks
      foreach(iChunk = 1:chunkSize, .combine = rBind) %do% {

        # Index in the complete task
        i <- (k - 1) * chunkSize + iChunk

        # Solve PDE from initial condition and get subsolution
        u0 <- dWn1D(x = xFullGrid, mu = x0SubGrid[i, 1], sigma = sdInitial) %*%
          t(dWn1D(x = xFullGrid, mu = x0SubGrid[i, 2], sigma = sdInitial))
        u0 <- matrix(u0 / periodicTrapRule2D(fxy = u0), ncol = 1)
        q <- crankNicolson2D(u0 = u0, bx = bxGrid, by = byGrid,
                             sigma2x = sigma2xGrid, sigma2y = sigma2yGrid,
                             sigmaxy = sigmaxyGrid, Mx = Mx,
                             deltax = xFullGrid[2] - xFullGrid[1], My = Mx,
                             deltay = xFullGrid[2] - xFullGrid[1], N = indtx,
                             deltat = max(tx) / Mt)[indSubGrid, ]

        # Result as a sparse matrix to save space - normalization comes in
        # later!
        q[q < tol] <- 0
        Matrix(q, sparse = TRUE)

      }

    }
    close(pb)

    # Save PDE?
    if (savePde) {

      save(list = "p", file = paste("pde-", type, "-", l, "-2D.RData",
                                    sep = ""))

    }

    # Free memory
    gc()

  }

  # Compute approximations?
  if (!onlyPde) {

    # Progress
    cat("\n Approxs for", l, "/", nrow(arrayPars), "\n")

    # Parallel/sequential loop over times
    res2D[l, , ] <- foreach(ti = 1:length(tx), .combine = rbind,
                            .packages = packages) %do% {

      # Progress
      cat("\n Time", ti, "/", length(tx))

      # Time
      t <- tx[ti]

      # Only the i-th time and as a dense matrix (eval, intial)
      pti <- matrix(p[, ti], nrow = nx0SubGrid, ncol = nx0SubGrid)

      # Normalize PDE: rounding to tolerance, set zeros to tolerance and
      # ensure it is a proper density
      if (any(pti < -tol)) {

        cat(paste("PDE negative values in i =", i, "with val =",
                  signif(min(pti), 2), "\n"))

      }
      pti <- round(pti, digits = -log10(tol))
      pti[pti > xMax] <- xMax
      pti[pti < tol | !is.finite(pti)] <- tol
      pti <- sweep(pti, 2, apply(pti, 2,
                                 function(x) periodicTrapRule2D(fxy = x)), "/")

      # Stationary
      klS <- repCol(dStat, nx0SubGrid)
      klS <- apply(log(pti / klS) * pti, 2, periodicTrapRule2D) * dStat
      gc()

      # Parallel KL computation for pseudo-likelihoods and approximation (with
      # smooth initial condition)
      kl <- function(method, circular, vmApprox) {

        # Progress
        cat("\n", paste(ifelse(circular, "", "U"), method,
                        ifelse(vmApprox, "vM", ""), sep = ""), "\n")

        # Pseudo-likelihood with task chunking (export objects not defined
        # inside kl)
        pb <- txtProgressBar(max = length(tx), style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        ap <- foreach(k = 1:nTasks, .combine = cBind,
                      .export = c("gxSubGrid", "x0SubGrid", "nx0SubGrid", "t",
                                  "pars", "maxK", "tol", "xMax", "b", "jac.b",
                                  "sigma2", "chunkSize", "nTasks"),
                      .noexport = c("p", "pti", "S"), .packages = packages,
                      .options.snow = opts) %dopar% {

          # Chunk of tasks
          foreach(iChunk = 1:chunkSize, .combine = cbind) %do% {

            # Index in the complete task
            i <- (k - 1) * chunkSize + iChunk

            # Result
            if (method == "WOU") {

              q <- tryCatch(dTpdWou2D(x = gxSubGrid,
                                     x0 = repRow(x0SubGrid[i, , drop = FALSE],
                                                 nx0SubGrid), t = t,
                                     alpha = pars[1:3], mu = pars[4:5],
                                     sigma = pars[6:7], maxK = maxK),
                            error = function(e) {
                message(paste("Error in method = ", method, ", circular = ",
                              circular, ", vmApprox = ", vmApprox, sep = ""))
                rep(tol, nx0SubGrid)
                })

            } else {

              q <- tryCatch(dPsTpd(x = gxSubGrid,
                                   x0 = x0SubGrid[i, , drop = FALSE], t = t,
                                   b = b, jac.b = jac.b, sigma2 = sigma2,
                                   method = method, pars = pars, maxK = maxK,
                                   circular = circular, vmApprox = vmApprox),
                            error = function(e) {
                message(paste("Error in method = ", method, ", circular = ",
                              circular, ", vmApprox = ", vmApprox, sep = ""))
                rep(tol, nx0SubGrid)
              })

            }

            # Return as a sparse matrix
            q[q > xMax] <- xMax
            q[q < tol | !is.finite(q)] <- 0
            Matrix(q, sparse = TRUE)

          }

        }
        close(pb)

        # Smoothing
        ap <- ap %*% S

        # Normalize solution: rounding to tolerance, set zeros to tolerance and
        # ensure it is a proper density
        ap <- round(ap, digits = -log10(tol))
        ap[ap < tol] <- tol
        ap <- sweep(ap, 2, apply(ap, 2,
                                 function(x) periodicTrapRule2D(fxy = x)), "/")

        # KL
        apply(log(pti / ap) * pti, 2, periodicTrapRule2D) * dStat

      }

      # Pseudo-likelihoods
      klEu <- kl(method = "E", circular = TRUE, vmApprox = FALSE)
      gc()
      klSo <- kl(method = "SO", circular = TRUE, vmApprox = FALSE)
      gc()
      klUEu <- kl(method = "E", circular = FALSE, vmApprox = FALSE)
      gc()
      klUSo <- kl(method = "SO", circular = FALSE, vmApprox = FALSE)
      gc()
      klEuVm <- kl(method = "E", circular = TRUE, vmApprox = TRUE)
      gc()
      klSoVm <- kl(method = "SO", circular = TRUE, vmApprox = TRUE)
      gc()
      klAp <- kl(method = "WOU", circular = TRUE, vmApprox = FALSE)
      gc()

      # Expectation
      eKlS <- periodicTrapRule2D(fxy = klS)
      eKlEu <- periodicTrapRule2D(fxy = klEu)
      eKlSo <- periodicTrapRule2D(fxy = klSo)
      eKlUEu <- periodicTrapRule2D(fxy = klUEu)
      eKlUSo <- periodicTrapRule2D(fxy = klUSo)
      eKlEuVm <- periodicTrapRule2D(fxy = klEuVm)
      eKlSoVm <- periodicTrapRule2D(fxy = klSoVm)
      eKlAp <- periodicTrapRule2D(fxy = klAp)

      # Result
      c(eKlS, eKlEu, eKlSo, eKlUEu, eKlUSo, eKlEuVm, eKlSoVm, eKlAp)

    }

    # Save data
    save(list = "res2D", file = paste("kl-tpd-", type, "-2D.RData", sep = ""))

  }

}

# Save data
save(list = "res2D", file = paste("kl-tpd-", type, "-2D.RData", sep = ""))

# Stop cluster
stopCluster(cl)

}


## Analysis 2D
{

# Load packages
library(ggplot2)
library(reshape2)

# Type of process
# type <- "WN"
type <- "vM"

# Load plot
load(file = paste("kl-tpd-", type, "-2D.RData", sep = ""))

# Format array into a data frame
dfRes2D <- melt(res2D)
colnames(dfRes2D) <- c("model", "t", "Approx", "KL")
dfRes2D$model <- factor(dfRes2D$model, levels = 1:9)
dfRes2D$t <- c(seq(0.05, 0.30, by = 0.05), seq(0.4, 1.5, by = 0.1),
               2, 2.5, 3)[dfRes2D$t]
dfRes2D$Approx <- factor(dfRes2D$Approx, levels = 1:8,
                         labels = c("S", "E", "SO", "UE", "USO", "EvM",
                                    "SOvM", "WOU"))

# Captions
if (type == "WN") {

  labels <- c(
    "1" = "list({alpha[1] == alpha[2]} == 0.5, alpha[3] == 0.25,
    {sigma[1] == sigma[2]} == 0.5)",
    "2" = "list({alpha[1] == alpha[2]} == 1, alpha[3] == 0.5,
    {sigma[1] == sigma[2]} == 0.5)",
    "3" = "list({alpha[1] == alpha[2]} == 2, alpha[3] == 1,
    {sigma[1] == sigma[2]} == 0.5)",
    "4" = "list({alpha[1] == alpha[2]} == 0.5, alpha[3] == 0.25,
    sigma[1] == 0.75, sigma[2] == 1.25)",
    "5" = "list({alpha[1] == alpha[2]} == 1, alpha[3] == 0.5,
    sigma[1] == 0.75, sigma[2] == 1.25)",
    "6" = "list({alpha[1] == alpha[2]} == 2, alpha[3] == 1,
    sigma[1] == 0.75, sigma[2] == 1.25)",
    "7" = "list({alpha[1] == alpha[2]} == 0.5, alpha[3] == 0.25,
    {sigma[1] == sigma[2]} == 2)",
    "8" = "list({alpha[1] == alpha[2]} == 1, alpha[3] == 0.5,
    {sigma[1] == sigma[2]} == 2)",
    "9" = "list({alpha[1] == alpha[2]} == 2, alpha[3] == 1,
    {sigma[1] == sigma[2]} == 2)"
    )

} else if (type == "vM") {

  labels <- c(
    "1" = "list({alpha[1] == alpha[2]} == 0.5, alpha[3] == 0.25, sigma == 0.5)",
    "2" = "list({alpha[1] == alpha[2]} == 1, alpha[3] == 0.5, sigma == 0.5)",
    "3" = "list({alpha[1] == alpha[2]} == 2, alpha[3] == 1, sigma == 0.5)",
    "4" = "list({alpha[1] == alpha[2]} == 0.5, alpha[3] == 0.25, sigma == 1)",
    "5" = "list({alpha[1] == alpha[2]} == 1, alpha[3] == 0.5, sigma == 1)",
    "6" = "list({alpha[1] == alpha[2]} == 2, alpha[3] == 1, sigma == 1)",
    "7" = "list({alpha[1] == alpha[2]} == 0.5, alpha[3] == 0.25, sigma == 2)",
    "8" = "list({alpha[1] == alpha[2]} == 1, alpha[3] == 0.5, sigma == 2)",
    "9" = "list({alpha[1] == alpha[2]} == 2, alpha[3] == 1, sigma == 2)"
  )

}

# Recreate ggplot2's colors manually
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Manual legends
breaks <- c("S", "E", "UE", "EvM", "SO", "USO", "SOvM", "WOU")
colours <- gg_color_hue(4)[c(1, rep(2, 3), rep(3, 3), 4)]
shapes <- c(16, rep(c(16, 1, 6), 2), 16)
names(colours) <- names(shapes) <- breaks

# Beautiful plot
s <- 10.5
m <- ggplot(data = dfRes2D, aes(x = t, y = KL, colour = Approx, shape = Approx))
m <- m + geom_line(alpha = 0.75) + geom_point(alpha = 0.75, size = 2) +
  scale_y_log10() +  facet_wrap(~ model, ncol = 3, nrow = 3,
                                labeller = as_labeller(labels, label_parsed)) +
  scale_color_manual(name = "Approx.", labels = breaks, breaks = breaks,
                     values = colours) +
  scale_shape_manual(name = "Approx.", labels = breaks, breaks = breaks,
                     values = shapes) +
  theme(aspect.ratio = 1, strip.text = element_text(size = s),
        legend.text = element_text(size = s),
        legend.title = element_text(size = s + 2),
        axis.text.x = element_text(size = s),
        axis.text.y = element_text(size = s),
        axis.title.x = element_text(size = s),
        axis.title.y = element_text(size = s))

m

# Save
ggsave(filename = paste("kl-tpd-", type, "-2D.pdf", sep = ""), plot = m,
       width = 10.5, height = 10.5)

}
