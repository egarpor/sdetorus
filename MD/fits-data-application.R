
library(circular)
library(ggplot2)
library(DirStats)
library(sdetorus)
library(gridExtra)
library(doParallel)
library(foreach)

## Load data
{

# Bottaro, S., Lindorff-Larsen, K. and Best, R. B. (2013).
# 100ns MD simulations on GB3 - one is with CHARMM36 + EEF1SB
# (the secondary structure is partially lost) - the second is performed using
# CHARMM22+GBSW (very stable)
load(file = "dih-bottaro.RData")

# Set NA's for first and last dihedrals
dih.c22$phi[dih.c22$n.aa == 1] <- NA
dih.c22$psi[dih.c22$n.aa == 56] <- NA
dih.c36$phi[dih.c36$n.aa == 1] <- NA
dih.c36$psi[dih.c36$n.aa == 56] <- NA

}


## Nonparametric fits
{

npFitDrift1D <- function(x = seq(-pi, pi, l = 201)[-201], data, delta, h = NULL,
                         p = 0) {

  # Map to an iid problem
  dataX <- DirStats::to.circle(data[-length(data)])
  dataY <- diffCirc(data) / delta

  # Automatic bandwidth
  if (is.null(h)) {

    h <- DirStats::bw.cv.loc(data.dir = dataX, data.lin = dataY, p = p)$h.opt

  }

  # Local linear estimate of the regression function
  npFit <- drop(DirStats::loc.directional.linear(
    x = DirStats::to.circle(x), data.dir = dataX,
    data.lin = dataY, h = h, p = p)$Yhat)

  return(list("dataX" = data[-length(data)], "dataY" = dataY, "npFit" = npFit,
              "h" = h))

}

npFitDiffusion1D <- function(x = seq(-pi, pi, l = 201)[-201], data, delta,
                             h = NULL, hDrift = h, hDens = h, p = 0,
                             type = c("simple", "drift", "integral")[1]) {

  # Map to an iid problem
  dataX <- DirStats::to.circle(data[-length(data)])
  dataY <- switch(type,
                  "simple" = diffCirc(data)^2 / delta,
                  "drift" = toPiInt(data[-1] - data[-length(data)] -
                              delta * npFitDrift1D(x = data[-length(data)],
                                                   data = data[-length(data)],
                                                   delta = delta, h = hDrift,
                                                   p = p)$npFit)^2 / delta)

  # Automatic bandwidth
  if (is.null(h)) {

    h <- DirStats::bw.cv.loc(data.dir = dataX, data.lin = dataY, p = p)$h.opt

  }

  # Local linear estimate of the regression function
  npFit <- drop(DirStats::loc.directional.linear(
    x = DirStats::to.circle(x), data.dir = dataX,
    data.lin = dataY, h = h, p = p)$Yhat)

  return(list("dataX" = data[-length(data)], "dataY" = dataY, "npFit" = npFit,
              "h" = h))

}

npFitDiffusion2D <- function(x = seq(-pi, pi, l = 201)[-201], data, delta,
                             h1 = NULL, h2 = NULL, target = c("b", "sigma2",
                                                              "sigmaxy")[1]) {

  # Map to an iid problem
  dataX1 <- DirStats::to.circle(data[-nrow(data), 1])
  dataX2 <- DirStats::to.circle(data[-nrow(data), 2])
  dataY <- switch(target,
                  b = diffCirc(data) / delta,
                  sigma2 = diffCirc(data)^2 / delta,
                  sigmaxy = repCol(apply(diffCirc(data), 1, prod) / delta, 2))

  # Automatic bandwidth
  if (is.null(h1) | is.null(h2)) {

    ha <- DirStats::bw.cv.loc(data.dir = dataX1, data.lin = dataY[, 1],
                              p = 0)$h.opt
    hb <- DirStats::bw.cv.loc(data.dir = dataX2, data.lin = dataY[, 1],
                              p = 0)$h.opt
    h1 <- c(ha, hb)

    ha <- DirStats::bw.cv.loc(data.dir = dataX1, data.lin = dataY[, 2],
                              p = 0)$h.opt
    hb <- DirStats::bw.cv.loc(data.dir = dataX2, data.lin = dataY[, 2],
                              p = 0)$h.opt
    h2 <- c(ha, hb)

  }

  # Local linear estimate of the regression function
  x1 <- DirStats::to.circle(x[, 1])
  x2 <- DirStats::to.circle(x[, 2])
  npFit1 <- drop(DirStats::nw.directional.directional.linear(
    x1 = x1, x2 = x2, data.dir.1 = dataX1, data.dir.2 = dataX2,
    data.lin = dataY[, 1], h1 = h1[1], h2 = h1[2]))
  npFit2 <- drop(DirStats::nw.directional.directional.linear(
    x1 = x1, x2 = x2, data.dir.1 = dataX1, data.dir.2 = dataX2,
    data.lin = dataY[, 2], h1 = h2[1], h2 = h2[2]))

  return(list("dataX" = data[-nrow(data), ], "dataY" = dataY,
              "npFit" = cbind(npFit1, npFit2), "h1" = h1, "h2" = h2))

}

smooth1D <- function(x, dataX, dataY, h, p = 0) {

  drop(DirStats::loc.directional.linear(
    x = DirStats::to.circle(x), data.dir = DirStats::to.circle(dataX),
    data.lin = dataY, h = h, p = p)$Yhat)

}

smooth2D <- function(x, dataX, dataY, h1, h2) {

  x1 <- DirStats::to.circle(x[, 1])
  x2 <- DirStats::to.circle(x[, 2])
  dataX1 <- DirStats::to.circle(dataX[ , 1])
  dataX2 <- DirStats::to.circle(dataX[ , 2])
  npFit1 <- drop(DirStats::nw.directional.directional.linear(
    x1 = x1, x2 = x2, data.dir.1 = dataX1, data.dir.2 = dataX2,
    data.lin = dataY[, 1], h1 = h1[1], h2 = h1[2]))
  npFit2 <- drop(DirStats::nw.directional.directional.linear(
    x1 = x1, x2 = x2, data.dir.1 = dataX1, data.dir.2 = dataX2,
    data.lin = dataY[, 2], h1 = h2[1], h2 = h2[2]))
  cbind(npFit1, npFit2)

}

npDensityDiffusion1D <- function(x, data, h = NULL) {

  dataX <- DirStats::to.circle(data)
  if (is.null(h)) {

    h <- DirStats::bw.mixt.pi.mise(x = dataX)

  }

  npFit <- kernel.directional(x = DirStats::to.circle(x), data = dataX, h = h)
  return(list("dataX" = data[-length(data)], "npFit" = npFit, "h" = h))

}

npDensityDiffusion2D <- function(x, data, h = NULL) {

  dataX1 <- DirStats::to.circle(data[, 1])
  dataX2 <- DirStats::to.circle(data[, 2])
  if (is.null(h)) {

    h1 <- DirStats::bw.mixt.pi.mise(x = dataX1)
    h2 <- DirStats::bw.mixt.pi.mise(x = dataX2)
    h <- c(h1, h2)

  }

  npFit <- DirStats::kernel.directional.directional(
    x1 = DirStats::to.circle(x[, 1]), x2 = DirStats::to.circle(x[, 2]),
    data.dir.1 = dataX1, data.dir.2 = dataX2, h1 = h[1], h2 = h[2],
    vec.grid = TRUE)
  return(list("dataX" = data[-nrow(data), ], "npFit" = npFit, "h" = h))

}

npCheck1D <- function(samp, delta, b, sigma2, rTraj, hdrift = NULL,
                      hdiff = NULL, hdens = NULL, nth = 200, type = "simple",
                      p = 0) {

  # Grid
  th <- seq(-pi, pi, l = nth)

  # Drifts NP and P smoothed
  driftNP <- npFitDrift1D(x = th, data = samp, delta = delta, h = hdrift,
                          p = p)
  driftP <- smooth1D(x = th, dataX = samp, dataY = b(samp), h = driftNP$h,
                     p = p)
  driftNP <- driftNP$npFit

  # Diffusions NP and P smoothed
  diffusionNP <- npFitDiffusion1D(x = th, data = samp, delta = delta, h = hdiff,
                                  hDrift = hdrift, p = p, type = type)
  diffusionP <- smooth1D(x = th, dataX = samp, dataY = sigma2(samp),
                         h = diffusionNP$h, p = p)
  diffusionNP <- diffusionNP$npFit

  # Density data
  kde <- npDensityDiffusion1D(x = th, data = samp, h = hdens)
  dens <- smooth1D(x = th, dataX = th, dataY = d(th), h = kde$h, p = 0)
  kde <- kde$npFit

  # Create data frames
  df <- data.frame(x = th, driftP = driftP, driftNP = driftNP,
                   diffusionP = diffusionP, diffusionNP = diffusionNP,
                   normDrift = abs(driftP - driftNP),
                   normDiff = abs(diffusionP - diffusionNP), kde = kde,
                   dens = dens)
  colnames(df) <- c("x", "driftP", "driftNP", "diffusionP", "diffusionNP",
                    "normDrift", "normDiff", "kde", "dens")
  fitSamp <- rTraj(x0 = samp[1], N = length(samp) - 1, delta = delta)
  dsamp <- data.frame(xData = samp, xSamp = fitSamp,
                      t = seq(0, T, l = length(samp)))
  q <- vector("list", 4)

  # ggplot2's colors
  gg_color_hue <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65,
                                  c = 100)[1:n]

  # Common settings
  the <- theme(aspect.ratio = 1, legend.position = "bottom",
               legend.box = "horizontal", legend.direction = "horizontal",
               legend.justification = "center", legend.key.size = unit(2, "cm"),
               legend.key.height = unit(0.5, "cm"),
               legend.key.width = unit(2, "cm"),
               legend.text = element_text(size = 16),
               legend.title = element_text(size = 16),
               axis.text.x = element_text(size = 16),
               axis.text.y = element_text(size = 16),
               axis.title.x = element_text(size = 16),
               axis.title.y = element_text(size = 16))
  gui1 <- guides(colour = guide_legend(order = 1),
                 fill = guide_colourbar(order = 2),
                 alpha = guide_legend(order = 3))
  gui2 <- guides(colour = guide_legend(override.aes = list(size = 1)))

  # Trajectory
  a <- 0.3
  siz <- 0.3
  q[[1]] <- ggplot(data = dsamp) +
    scale_x_continuous(name = expression(t)) +
    scale_y_continuous(name = expression(psi[t]), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    coord_cartesian(ylim = c(-pi, pi))

  q[[1]] <- q[[1]] +
    geom_path(mapping = aes(x = t, y = xData - 2 * pi, colour = "Data"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = t, y = xData, colour = "Data"), data = dsamp,
              size = siz, alpha = a) +
    geom_path(mapping = aes(x = t, y = xData + 2 * pi, colour = "Data"),
              data = dsamp, size = siz, alpha = a)

  q[[1]] <- q[[1]] +
    geom_path(mapping = aes(x = t, y = xSamp - 2 * pi, colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = t, y = xSamp, colour = "Fit"), data = dsamp,
              size = siz, alpha = a) +
    geom_path(mapping = aes(x = t, y = xSamp + 2 * pi, colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_point(mapping = aes(x = t, y = toPiInt(xData), colour = "Data"),
               size = siz * 1.5, alpha = a)  +
    geom_point(mapping = aes(x = t, y = toPiInt(xSamp), colour = "Fit"),
               size = siz * 1.5, alpha = a)

  q[[1]] <- q[[1]] +
    scale_color_manual(name = "Trajectory",
                       values = c("Fit" = gg_color_hue(5)[3],
                                  "Data" = gg_color_hue(5)[5])) +
    the + gui1 + gui2

  # Density
  q[[2]] <- ggplot(data = df) +
    geom_line(mapping = aes(x = x, y = dens, colour = "P",
                            alpha = kde), size = 1) +
    geom_line(mapping = aes(x = x, y = kde, colour = "NP",
                            alpha = kde), size = 1) +
    scale_x_continuous(name = expression(psi), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    scale_y_continuous(name = expression(f),
                       limits = c(0, max(c(df$kde, df$dens)))) +
    scale_color_manual(name = expression(hat(f)),
                       values = c("P" = gg_color_hue(5)[3],
                                  "NP" = gg_color_hue(5)[5])) +
    scale_alpha_continuous("KDE", range = c(0.1, 1),
                           breaks = pretty(c(0, max(df$kde)), n = 2),
                           limits = c(0, max(df$kde))) +
    the + gui1 + gui2

  # Drift
  q[[3]] <- ggplot(data = df) +
    geom_line(mapping = aes(x = x, y = driftP, colour = "P",
                            alpha = kde), size = 1) +
    geom_line(mapping = aes(x = x, y = driftNP, colour = "NP",
                            alpha = kde), size = 1) +
    scale_x_continuous(name = expression(psi), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    scale_y_continuous(name = expression(b),
                       limits = max(abs(c(df$driftP, df$driftNP))) * c(-1, 1)) +
    scale_color_manual(name = expression(hat(b)),
                       values = c("P" = gg_color_hue(5)[3],
                                  "NP" = gg_color_hue(5)[5])) +
    scale_alpha_continuous("KDE", range = c(0.1, 1),
                           breaks = pretty(c(0, max(df$kde)), n = 2),
                           limits = c(0, max(df$kde))) +
    the + gui1 + gui2

  # Diffusion
  q[[4]] <- ggplot(data = df) +
    geom_line(mapping = aes(x = x, y = sqrt(diffusionP), colour = "P",
                            alpha = kde), size = 1) +
    geom_line(mapping = aes(x = x, y = sqrt(diffusionNP), colour = "NP",
                            alpha = kde), size = 1) +
    scale_x_continuous(name = expression(theta), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    scale_y_continuous(name = expression(sigma),
                       limits = c(0, sqrt(max(c(df$diffusionP,
                                                df$diffusionNP))))) +
    scale_color_manual(name = expression(hat(sigma)),
                       values = c("P" = gg_color_hue(5)[3],
                                  "NP" = gg_color_hue(5)[5])) +
    scale_alpha_continuous("KDE", range = c(0.1, 1),
                           breaks = pretty(c(0, max(df$kde)), n = 2),
                           limits = c(0, max(df$kde))) +
    the + gui1 + gui2

  return(q)

}

npCheck2D <- function(samp, delta, b, sigma2, sigmaxy, rTraj, hdrift1 = NULL,
                      hdrift2 = NULL, hdiff1 = NULL, hdiff2 = NULL,
                      hdiff3 = NULL, hdens = NULL, nth = 100, p = 0) {

  # Grid
  sth <- 1.25 * seq(-pi, pi, l = nth)
  th <- as.matrix(expand.grid(sth, sth))
  sth2 <- 1.25 * seq(-pi, pi, l = 20)
  th2 <- as.matrix(expand.grid(sth2, sth2))

  # Drifts NP and P smoothed
  h1 <- hdrift1
  h2 <- hdrift2
  driftNP <- npFitDiffusion2D(x = th, data = samp, delta = delta, h1 = h1,
                              h2 = h2, target = "b")
  driftNP2 <- npFitDiffusion2D(x = th2, data = samp, delta = delta,
                               h1 = driftNP$h1, h2 = driftNP$h2,
                               target = "b")$npFit
  driftP <- smooth2D(x = th, dataX = samp, dataY = b(samp), h1 = driftNP$h1,
                     h2 = driftNP$h2)
  driftP2 <- smooth2D(x = th2, dataX = samp, dataY = b(samp), h1 = driftNP$h1,
                      h2 = driftNP$h2)
  driftNP <- driftNP$npFit

  # Diffusions NP and P smoothed
  h1 <- hdiff1
  h2 <- hdiff2
  diffusionNP <- npFitDiffusion2D(x = th, data = samp, delta = delta, h1 = h1,
                                  h2 = h2, target = "sigma2")
  diffusionNP2 <- npFitDiffusion2D(x = th2, data = samp, delta = delta,
                                   h1 = diffusionNP$h1, h2 = diffusionNP$h2,
                                   target = "sigma2")$npFit
  diffusionP <- smooth2D(x = th, dataX = samp, dataY = sigma2(samp),
                         h1 = diffusionNP$h1, h2 = diffusionNP$h2)
  diffusionP2 <- smooth2D(x = th2, dataX = samp, dataY = sigma2(samp),
                          h1 = diffusionNP$h1, h2 = diffusionNP$h2)
  diffusionNP <- diffusionNP$npFit

  # Cov diffusions NP and P smoothed
  h3 <- hdiff3
  diffusionCovNP <- npFitDiffusion2D(x = th, data = samp, delta = delta,
                                     h1 = h3, h2 = h3, target = "sigmaxy")
  diffusionCovNP2 <- npFitDiffusion2D(x = th2, data = samp, delta = delta,
                                      h1 = diffusionCovNP$h1,
                                      h2 = diffusionCovNP$h2,
                                      target = "sigmaxy")$npFit[, 1]
  diffusionCovP <- smooth2D(x = th, dataX = samp,
                            dataY = repCol(sigmaxy(samp), 2),
                            h1 = diffusionCovNP$h1, h2 = diffusionCovNP$h2)[, 1]
  diffusionCovP2 <- smooth2D(x = th2, dataX = samp,
                             dataY = repCol(sigmaxy(samp), 2),
                             h1 = diffusionCovNP$h1,
                             h2 = diffusionCovNP$h2)[, 1]
  diffusionCovNP <- diffusionCovNP$npFit[, 1]

  # Density data
  kde <- npDensityDiffusion2D(x = th, data = samp, h = hdens)
  kde2 <- npDensityDiffusion2D(x = th2, data = samp, h = kde$h)$npFit
  dens <- smooth2D(x = th, dataX = th, dataY = repCol(d(th), 2), h1 = kde$h,
                   h2 = kde$h)[, 1]
  dens2 <- smooth2D(x = th2, dataX = th2, dataY = repCol(d(th2), 2), h1 = kde$h,
                    h2 = kde$h)[, 1]
  kde <- kde$npFit

  # Create data frames
  normDrift <- sqrt(rowSums((driftP - driftNP)^2)) / sqrt(rowSums(driftNP^2))
  normDiff <- sqrt(rowSums((diffusionP - diffusionNP)^2 +
                             2 * (diffusionCovNP - diffusionCovP)^2)) /
    sqrt(rowSums(diffusionNP^2 + 2 * diffusionCovNP^2))
  df <- data.frame(x = th[, 1], y = th[, 2], driftP = driftP, driftNP = driftNP,
                   diffusionP = diffusionP, diffusionNP = diffusionNP,
                   diffusionCovP = diffusionCovP,
                   diffusionCovNP = diffusionCovNP, normDrift = normDrift,
                   normDiff = normDiff, kde = kde, dens = dens)
  df2 <- data.frame(x = th2[, 1], y = th2[, 2], driftP = driftP2,
                    driftNP = driftNP2, diffusionP = diffusionP2,
                    diffusionNP = diffusionNP2, diffusionCovP = diffusionCovP2,
                    diffusionCovNP = diffusionCovNP2,
                    normDrift = sqrt(rowSums((driftP2 - driftNP2)^2)),
                    normDiff = sqrt(rowSums((diffusionP2 - diffusionNP2)^2)),
                    kde = kde2, dens = dens2)
  colnames(df) <- colnames(df2) <- c("x", "y", "driftP.1", "driftP.2",
                                     "driftNP.1", "driftNP.2", "diffusionP.1",
                                     "diffusionP.2", "diffusionNP.1",
                                     "diffusionNP.2", "diffusionCovP",
                                     "diffusionCovNP", "normDrift",  "normDiff",
                                     "kde", "dens")
  fitSamp <- rTraj(x0 = samp[1, , drop = TRUE], N = nrow(samp) - 1,
                   delta = delta)
  dsamp <- data.frame(xData = samp[, 1], yData = samp[, 2],
                      xSamp = fitSamp[, 1], ySamp = fitSamp[, 2],
                      t = seq(0, T, l = nrow(samp)))
  q <- vector("list", 4)

  # ggplot2's colors
  gg_color_hue <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65,
                                  c = 100)[1:n]

  # Common settings
  the <- theme(aspect.ratio = 1, legend.position = "bottom",
               legend.box = "horizontal", legend.direction = "horizontal",
               legend.justification = "center", legend.key.size = unit(2, "cm"),
               legend.key.height = unit(0.5, "cm"),
               legend.key.width = unit(2, "cm"),
               legend.text = element_text(size = 16),
               legend.title = element_text(size = 16),
               axis.text.x = element_text(size = 16),
               axis.text.y = element_text(size = 16),
               axis.title.x = element_text(size = 16),
               axis.title.y = element_text(size = 16))
  gui1 <- guides(colour = guide_legend(order = 1),
                 fill = guide_colourbar(order = 2),
                 alpha = guide_legend(order = 3))
  gui2 <- guides(colour = guide_legend(override.aes = list(size = 1)))

  # Trajectory
  a <- 0.2
  siz <- 0.2
  q[[1]] <- ggplot(data = dsamp) +
    scale_x_continuous(name = expression(phi[t]), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    scale_y_continuous(name = expression(psi[t]), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi, pi))

  q[[1]] <- q[[1]] +
    geom_path(mapping = aes(x = xData - 2 * pi, y = yData - 2 * pi,
                            colour = "Data"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xData - 2 * pi, y = yData, colour = "Data"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xData - 2 * pi, y = yData + 2 * pi,
                            colour = "Data"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xData, y = yData - 2 * pi, colour = "Data"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xData, y = yData, colour = "Data"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xData, y = yData + 2 * pi, colour = "Data"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xData + 2 * pi, y = yData - 2 * pi,
                            colour = "Data"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xData + 2 * pi, y = yData, colour = "Data"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xData + 2 * pi, y = yData + 2 * pi,
                            colour = "Data"),
              data = dsamp, size = siz, alpha = a)

  q[[1]] <- q[[1]] +
    geom_path(mapping = aes(x = xSamp - 2 * pi, y = ySamp - 2 * pi,
                            colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xSamp - 2 * pi, y = ySamp, colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xSamp - 2 * pi, y = ySamp + 2 * pi,
                            colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xSamp, y = ySamp - 2 * pi, colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xSamp, y = ySamp, colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xSamp, y = ySamp + 2 * pi, colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xSamp + 2 * pi, y = ySamp - 2 * pi,
                            colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xSamp + 2 * pi, y = ySamp, colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_path(mapping = aes(x = xSamp + 2 * pi, y = ySamp + 2 * pi,
                            colour = "Fit"),
              data = dsamp, size = siz, alpha = a) +
    geom_point(mapping = aes(x = toPiInt(xData), y = toPiInt(yData),
                             colour = "Data"), size = siz * 1.5, alpha = a)  +
    geom_point(mapping = aes(x = toPiInt(xSamp), y = toPiInt(ySamp),
                             colour = "Fit"), size = siz * 1.5, alpha = a)

  q[[1]] <- q[[1]] +
    geom_hline(yintercept = c(-pi, pi), colour = "gray") +
    geom_vline(xintercept = c(-pi, pi), colour = "gray") +
    scale_color_manual(name = "Trajectory",
                       values = c("Fit" = gg_color_hue(5)[3],
                                  "Data" = gg_color_hue(5)[5])) +
    the + gui1 + gui2

  # Density
  siz <- 0.1
  maxDensDiff <- max(dens)
  oob <- function(x) {
    x[x > maxDensDiff] <- maxDensDiff
    x[x < -maxDensDiff] <- -maxDensDiff
    return(x)
  }
  q[[2]] <- ggplot(data = df, aes(x = x, y = y)) +
    geom_raster(aes(fill = oob(kde - dens), alpha = kde), interpolate = TRUE) +
    scale_x_continuous(name = expression(theta[1]),
                       breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    scale_y_continuous(name = expression(theta[2]),
                       breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi, pi)) +
    scale_fill_gradient2(guide = "colourbar", name = "Difference",
                         limits = c(-maxDensDiff, maxDensDiff)) +
    scale_alpha_continuous("KDE", range = c(0.1, 1),
                           breaks = pretty(c(0, max(df$kde)), n = 2),
                           limits = c(0, max(df$kde))) +
    geom_hline(yintercept = c(-pi, pi), colour = "gray") +
    geom_vline(xintercept = c(-pi, pi), colour = "gray") +
    the + gui1 + gui2 + theme(legend.key.width = unit(1, "cm"))

  # Drift
  siz <- 1
  rDrift <- 1.5 / max(c(sqrt(df2$driftNP.1^2 + df2$driftNP.2^2),
                 sqrt(df2$driftP.1^2 + df2$driftP.2^2)), na.rm = TRUE)
  maxNormDrift <- max(normDrift)
  oob <- function(x) {
    x[x > maxNormDrift] <- maxNormDrift
    return(x)
  }
  q[[3]] <- ggplot(data = df, aes(x = x, y = y)) +
    # geom_raster(aes(fill = normDrift, alpha = kde), interpolate = TRUE) +
    scale_x_continuous(name = expression(phi), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    scale_y_continuous(name = expression(psi), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi, pi)) +
    # scale_fill_gradientn(guide = "colourbar", name = "Distance",
    #                      colours = matlab.like(202),
    #                      limits = c(0, maxNormDrift),
    #                      oob = oob) +
    scale_alpha_continuous("KDE", range = c(0.1, 1),
                           breaks = pretty(c(0, max(df$kde)), n = 2),
                           limits = c(0, max(df$kde))) +
    geom_segment(data = df2, aes(xend = x + rDrift * driftP.1,
                                 yend = y + rDrift * driftP.2, col = "P",
                                 alpha = kde), size = siz,
                 arrow = arrow(length = unit(0.1, "cm"))) +
    geom_segment(data = df2, aes(xend = x + rDrift * driftNP.1,
                                 yend = y + rDrift * driftNP.2, col = "NP",
                                 alpha = kde), size = siz,
                 arrow = arrow(length = unit(0.1, "cm"))) +
    scale_color_manual(name = expression(group("(", list(hat(b)[1],
                                                         hat(b)[2]), ")")),
                       values = c("P" = gg_color_hue(5)[3],
                                  "NP" = gg_color_hue(5)[5])) +
    geom_hline(yintercept = c(-pi, pi), colour = "gray") +
    geom_vline(xintercept = c(-pi, pi), colour = "gray") +
    the + gui1 + gui2

  # Diffusion
  siz <- 1
  rDiff <- 1.5 / max(c(sqrt(df2$diffusionNP.1^2 + df2$diffusionNP.2^2),
                 sqrt(df2$diffusionP.1^2 + df2$diffusionP.2^2)), na.rm = TRUE)
  maxNormDiff <- max(normDiff)#3
  oob <- function(x) {
    x[x > maxNormDiff] <- maxNormDiff
    return(x)
  }
  q[[4]] <- ggplot(data = df, aes(x = x, y = y)) +
    # geom_raster(aes(fill = normDiff, alpha = kde), interpolate = TRUE) +
    scale_x_continuous(name = expression(phi), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    scale_y_continuous(name = expression(psi), breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi, pi)) +
    # scale_fill_gradientn(guide = "colourbar", name = "Distance",
    #                      colours = matlab.like(202),
    #                      limits = c(0, maxNormDiff), oob = oob) +
    scale_alpha_continuous("KDE", range = c(0.1, 1),
                           breaks = pretty(c(0, max(df$kde)), n = 2),
                           limits = c(0, max(df$kde))) +
    geom_segment(data = df2, aes(xend = x + rDiff * sqrt(diffusionP.1),
                                 yend = y + rDiff * sqrt(diffusionP.2),
                                 col = "P", alpha = kde), size = siz,
                 arrow = arrow(length = unit(0.1, "cm"))) +
    geom_segment(data = df2, aes(xend = x + rDiff * sqrt(diffusionNP.1),
                                 yend = y + rDiff * sqrt(diffusionNP.2),
                                 col = "NP", alpha = kde), size = siz,
                 arrow = arrow(length = unit(0.1, "cm"))) +
    scale_color_manual(name = expression(group("(", list(hat(sigma)[1],
                                                         hat(sigma)[2]), ")")),
                       values = c("P" = gg_color_hue(5)[3],
                                  "NP" = gg_color_hue(5)[5])) +
    geom_hline(yintercept = c(-pi, pi), colour = "gray") +
    geom_vline(xintercept = c(-pi, pi), colour = "gray") +
    the + gui1 + gui2

  # Covariance diffusion
  maxNormDiff <- 1
  oob <- function(x) {
    x[x > maxNormDiff] <- maxNormDiff
    return(x)
  }
  q[[5]] <- ggplot(data = df, aes(x = x, y = y)) +
    geom_raster(aes(fill = oob(diffusionCovNP - df$diffusionCovP), alpha = kde),
                interpolate = TRUE) +
    scale_x_continuous(name = expression(theta[1]),
                       breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    scale_y_continuous(name = expression(theta[2]),
                       breaks = seq(-pi, pi, l = 5),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi, pi)) +
    scale_fill_gradient2(guide = "colourbar", name = "Difference",
                         limits = c(-maxNormDiff, maxNormDiff)) +
    scale_alpha_continuous("KDE", range = c(0.1, 1),
                           breaks = pretty(seq(0, max(df$kde), l = 10), n = 2),
                           limits = c(0, max(df$kde))) +
    geom_hline(yintercept = c(-pi, pi), colour = "gray") +
    geom_vline(xintercept = c(-pi, pi), colour = "gray") +
    the + gui1 + gui2 + theme(legend.key.width = unit(1, "cm"))

  return(q)

}

functionNpPlot <- function(data, k, hdrift, hdiff, hdens, display = "",
                           info = "") {

  # Common objects
  obs <- data[data$n.aa == k, ]
  N <- nrow(obs)
  T <- 100

  # Process
  plot(0, 0, type = "n", xlab = "t", ylab = expression(list(phi[t],psi[t])),
       xlim = c(0, T), ylim = c(-pi, pi), axes = FALSE,
       main = paste("Processeses for", display, k, "/ 56", info))
  box()
  axis(1)
  axis(2, at = seq(-pi, pi, pi/2), labels = expression(-pi, -pi/2, 0, pi/2, pi))
  abline(h = c(-pi, pi), lty = 2, col = "gray")
  linesCirc(seq(0, T, l = N), obs$phi, col = rgb(0, 0.25, 0.75, alpha = 0.75))
  linesCirc(seq(0, T, l = N), obs$psi, col = rgb(0.75, 0.25, 0, alpha = 0.75))

  # Densities
  res1 <- kernel.directional(x = DirStats::to.circle(seq(-pi, pi, l = 200)),
                             data = DirStats::to.circle(obs$phi), h = hdens)
  plot(seq(-pi, pi, l = 200), res1, type = "l",
       col = rgb(0, 0.25, 0.75, alpha = 0.75), lwd = 2, ylim = c(0, 2),
       xlab = expression(list(phi[t],psi[t])), ylab = "Density",
       main = paste("Densities for", k, "/ 56", info), axes = FALSE)
  hist(obs$phi, breaks = seq(-pi - 1e-3, pi + 1e-3, l = 50), freq = F,
       add = T, border = rgb(0, 0.25, 0.75, alpha = 0.75))
  rug(obs$phi, col = rgb(0, 0.25, 0.75, alpha = 0.25))
  res2 <- kernel.directional(x = DirStats::to.circle(seq(-pi, pi, l = 200)),
                             data = DirStats::to.circle(obs$psi), h = hdens)
  lines(seq(-pi, pi, l = 200), res2, type = "l",
        col = rgb(0.75, 0.25, 0, alpha = 0.75), lwd = 2)
  hist(obs$psi, breaks = seq(-pi - 1e-3, pi + 1e-3, l = 50), freq = FALSE,
       add = TRUE, border = rgb(0.75, 0.25, 0, alpha = 0.75))
  rug(obs$psi, col = rgb(0.75, 0.25, 0, alpha = 0.25))
  box(); axis(2); axis(1, at = seq(-pi, pi, pi/2),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi))
  abline(v = c(-pi, pi), lty = 2, col = "gray")


}

}


## Analysis 1D
{

# Load results of this block
# load("analysis-1D.RData")

# Chosen one
k1 <- 9
ang <- "psi"

# nCores <- min(detectCores() - 1, 10)
# cl <- makeCluster(nCores, outfile = "")
# registerDoParallel(cl = cl)
# # Select dataset: k-th AA and subset of observations
# foreach(ang = c("psi")) %:%
#   foreach(k1 = c(9, 21),
#           .packages = c("circular", "DirStats", "sdetorus",
#                         "ggplot2", "gridExtra")) %dopar% {

  # Read data
  data <- as.matrix(subset(dih.c36, select = ang, subset = (n.it %% 10 == 0) &
                             (dih.c36$n.aa == k1)))[-1, ]
  N <- length(data)
  T <- 100
  delta <- T / N

  # Starting values
  sigmaStat <- sqrt(sigmaDiff(data = cbind(data), delta = delta,
                              isotropic = TRUE)[1, 1])
  nStarts <- 10
  nK <- 3
  ml <- lapply(1:nK, function(k) lapply(1:nStarts, function(i)
    emVmf(data = cbind(data), k = k, kappaMax = 200, maxIter = 5000,
                      tol = rep(1e-7, 3))))
  ind <- matrix(apply(expand.grid(1:nK, 1:nStarts), 1, function(ij)
    ml[[ij[1]]][[ij[2]]]$BIC), nrow = nK, ncol = nStarts)
  ind <- which(ind == min(ind), arr.ind = TRUE)[1, ]
  ml <- ml[[ind[1]]][[ind[2]]]
  muStat <- drop(ml$M)
  kappaStat <- drop(ml$K)
  pStat <- ml$alpha
  alphaStat <- kappaStat * sigmaStat^2 / 2
  start <- rbind(c(log(alphaStat), muStat, log(sigmaStat)))
  lower <- c(rep(-Inf, length(pStat)), rep(-pi, length(pStat)),
             rep(-Inf, length(pStat)))[-6]
  upper <- -lower

  # Diffusion-specific
  b <- function(x, pars) {
    pt <- pStat#(1 + tanh(pars[6])) / 2
    cbind(driftMixVm(x = drop(x), alpha = exp(pars[1:length(pStat)]),
                     mu = pars[(length(pStat) + 1):(2 * length(pStat))],
                     sigma = exp(pars[2 * length(pStat) + 1]), p = pt))
  }
  b1 <- function(x, pars, h = 1e-4) {
    l <- length(x)
    res <- b(x = c(x + h, x - h), pars = pars)
    drop(res[1:l] - res[(l + 1):(2 * l)])/(2 * h)
  }
  sigma2 <- function(x, pars) rep(exp(2 * pars[5]), length(x))

  # Fits

  # Pseudo-likelihoods do not seem to get properly the proportion!
  # system.time(fit1b <- psMle(data = cbind(data), delta = delta, b = b,
  #                            sigma2 = sigma2, method = "E", start = start,
  #                            lower = lower, upper = upper, maxit = 5000,
  #                            optMethod = "nlm", verbose = 2, tol = 1e-6,
  #                            selectSolution = "lowest", K = 2))
  # fit1 <- fit1b
  #
  # system.time(fit1c <- psMle(data = cbind(data), delta = delta, b = b,
  #                            b1 = b1, sigma2 = sigma2, method = "SO",
  #                            start = start, lower = lower, upper = upper,
  #                            maxit = 5000, optMethod = "Nelder-Mead",
  #                            verbose = 2, tol = 1e-6,
  #                            selectSolution = "lowest", K = 2))
  # fit1 <- fit1c

  system.time(fit1a <- mlePde1D(data = drop(data), delta = delta, b = b,
                                sigma2 = sigma2, Mx = 500, Mt = 20,
                                sdInitial = 0.01, linearBinning = TRUE,
                                start = start, lower = lower, upper = upper,
                                maxit = 5000, optMethod = "Nelder-Mead",
                                verbose = 2, tol = 1e-8,
                                selectSolution = "lowest"))
  fit1 <- fit1a

  # head(fit1a)
  # head(fit1b)
  # head(fit1c)

  # Fitted parameters
  alpha <- exp(fit1$par[1:length(pStat)])
  mu <- fit1$par[(length(pStat) + 1):(2 * length(pStat))]
  sigma <- exp(fit1$par[2 * length(pStat) + 1])
  Sigma <- diag(sigma^2, nrow = 1, ncol = 1)
  p <- pStat #(1 + tanh(fit1$par[6])) / 2

  # Smoothing parameters
  pLoc <- 0
  hdens <- DirStats::bw.cv.dir(data.dir = DirStats::to.circle(data), h = 1,
                               efic = TRUE)$h.opt
  hdrift <- DirStats::bw.cv.loc(
    data.dir = DirStats::to.circle(data[-length(data)]),
    data.lin = diffCirc(data) / delta, p = pLoc, efic = TRUE,
    plot.it = FALSE)$h.opt
  hdiff <- DirStats::bw.cv.loc(
    data.dir = DirStats::to.circle(data[-length(data)]),
    data.lin = diffCirc(data)^2 / delta, p = pLoc, efic = TRUE,
    plot.it = FALSE)$h.opt

  # Check
  bCheck <- function(x) b(x, pars = c(log(alpha), mu, log(sigma),
                                      atanh(2 * p - 1)))
  sigma2Check <- function(x) sigma2(x, pars = c(log(alpha), mu, log(sigma),
                                                atanh(2 * p - 1)))
  d <- function(x) dVmf(x = x, M = cbind(mu),
                                K = cbind(2 * alpha / sigma^2), alpha = p)
  rTraj <- function(x0, N, delta) rTrajLangevin(x0 = x0, drift = driftMixVm,
                                                SigDif = Sigma, alpha = alpha,
                                                mu = mu, sigma = sigma, p = p,
                                                N = N, delta = delta,
                                                NFine = N * 10,
                                                deltaFine = delta / 10,
                                                circular = FALSE)
  set.seed(623679851)
  q1 <- npCheck1D(samp = unwrapCircSeries(data), delta = delta, b = bCheck,
                  sigma2 = sigma2Check, rTraj = rTraj, nth = 500,
                  hdrift = hdrift, hdiff = hdiff, hdens = hdens,
                  type = "simple", p = pLoc)
  # grid.arrange(q1[[1]], q1[[2]], q1[[3]], q1[[4]], ncol = 2, nrow = 2)
  # The discretization has an important effect on the NP estimate of sigma^2,
  # producing the curvy patterns. The bandwidth needs to be larger for it.

  g1 <- NULL
  while (is.null(g1)) {
    cat("Saving...\n")
    g1 <- tryCatch(grid.arrange(q1[[1]], q1[[2]], q1[[3]], q1[[4]],
                                ncol = 2, nrow = 2), error = function(e) NULL)
  }
  # ggsave(filename = paste("1d-", k1, "-", ang, ".png", sep = ""), plot = g1,
  #        width = 14, height = 14)

  for (j in 1:4) ggsave(filename = paste("1d-", k1, "-", j, ".pdf", sep = ""),
                        plot = q1[[j]], width = 7, height = 7)

# }
#
# stopCluster(cl)

}


## Analysis 2D
{

# Load results of this block
# load("analysis-2D.RData")

# Chosen one
k2 <- 14

# nCores <- 2
# cl <- makeCluster(nCores, outfile = "")
# registerDoParallel(cl = cl)
#
# # Select dataset: k-th AA and subset of observations
# foreach(k2 = c(54:50, 45:42, 16:14, 8:4),
#         .packages = c("circular", "DirStats", "sdetorus",
#                       "ggplot2", "gridExtra")) %dopar% {

  # Read data
  data <- as.matrix(subset(dih.c36, select = c(phi, psi),
                           subset = (n.it %% 10 == 0) &
                             (dih.c36$n.aa == k2)))[-1, ]
  N <- nrow(data)
  T <- 100
  delta <- T / N

  # Starting values
  ml1 <- mle.wrappednormal(x = circular(data[, 1]))
  ml2 <- mle.wrappednormal(x = circular(data[, 2]))
  ml <- mleOptimWrapper(minusLogLik = function(x) {
    -sum(log(DirStats::dwntd(th1 = data[, 1], th2 = data[, 2],
                             mu1 = x[1], mu2 = x[2], sigma1 = x[3],
                             sigma2 = x[4], rho = x[5], K = 5)))
  }, start = c(toPiInt(as.numeric(ml1$mu)), toPiInt(as.numeric(ml2$mu)),
               ml1$sd, ml2$sd, 0), lower = c(-pi, -pi, 0.01, 0.01, -0.99),
  upper = c(pi, pi, 20, 20, 0.99), optMethod = "Nelder-Mead", tol = 1e-8,
  maxit = 2000, verbose = 4, selectSolution = "lowest")
  muStat <- toPiInt(ml$par[1:2])
  SWN <- matrix(c(ml$par[3]^2, ml$par[5] * ml$par[3] * ml$par[4],
                  ml$par[5] * ml$par[3] * ml$par[4], ml$par[4]^2),
                nrow = 2, ncol = 2, byrow = TRUE)
  SigmaStat <- sigmaDiff(data = data, delta = delta)
  sigmaStat <- sqrt(diag(SigmaStat))
  rhoStat <- SigmaStat[1, 2] / prod(sigmaStat)
  AStat <- 0.5 * SigmaStat %*% solve(SWN)
  alphaStat <- aToAlpha(A = AStat, Sigma = SigmaStat)
  start <- c(alphaStat, muStat, sigmaStat, rhoStat)
  lower <- c(0.1, 0.1, -250, -pi, -pi, 0.1, 0.1, -0.99)
  upper <- c(300, 300, 250, pi, pi, 25, 25, 0.99)

  # Fit
  system.time(fit2 <- approxMleWn2D(data = data, delta = delta, start = start,
                                    K = 3, lower = lower, upper = upper,
                                    maxit = 10000, optMethod = "Nelder-Mead",
                                    verbose = 2, tol = 1e-8,
                                    selectSolution = "lowest"))

  # Fitted parameters
  head(fit2)
  alpha <- fit2$par[1:3]
  mu <- fit2$par[4:5]
  sigma <- fit2$par[6:7]
  rho <- fit2$par[8]
  A <- alphaToA(alpha = alpha, sigma = sigma, rho = rho)
  Sigma <- diag(sigma^2)
  Sigma[2, 1] <- Sigma[1, 2] <- rho * prod(sigma)

  # Smoothing parameters
  iso <- FALSE
  hdens <- unlist(head(DirStats::bw.cv.dirdir(
    data.dir.1 = DirStats::to.circle(data[, 1]),
    data.dir.2 = DirStats::to.circle(data[, 2]), diag = iso), 2))
  hdrift1 <- DirStats::bw.cv.nw.dir.dir(
    data.dir.1 = DirStats::to.circle(data[-nrow(data), 1]),
    data.dir.2 = DirStats::to.circle(data[-nrow(data), 2]),
    data.lin = diffCirc(data[, 1]) / delta, diag = iso, plot.it = FALSE,
    efic = TRUE)$h.opt
  hdrift2 <- DirStats::bw.cv.nw.dir.dir(
    data.dir.1 = DirStats::to.circle(data[-nrow(data), 1]),
    data.dir.2 = DirStats::to.circle(data[-nrow(data), 2]),
    data.lin = diffCirc(data[, 2]) / delta, diag = iso, plot.it = FALSE,
    efic = TRUE)$h.opt
  hdiff1 <- DirStats::bw.cv.nw.dir.dir(
    data.dir.1 = DirStats::to.circle(data[-nrow(data), 1]),
    data.dir.2 = DirStats::to.circle(data[-nrow(data), 2]),
    data.lin = diffCirc(data[, 1])^2 / delta, diag = iso, plot.it = FALSE,
    efic = TRUE)$h.opt
  hdiff2 <- DirStats::bw.cv.nw.dir.dir(
    data.dir.1 = DirStats::to.circle(data[-nrow(data), 1]),
    data.dir.2 = DirStats::to.circle(data[-nrow(data), 2]),
    data.lin = diffCirc(data[, 2])^2 / delta, diag = iso, plot.it = FALSE,
    efic = TRUE)$h.opt
  hdiff3 <- DirStats::bw.cv.nw.dir.dir(
    data.dir.1 = DirStats::to.circle(data[-nrow(data), 1]),
    data.dir.2 = DirStats::to.circle(data[-nrow(data), 2]),
    data.lin = diffCirc(data[, 2]) * diffCirc(data[, 1]) / delta, diag = iso,
    plot.it = FALSE, efic = TRUE)$h.opt

  # Diffusion-specific
  b <- function(x) driftWn2D(x = toPiInt(x), A = A, mu = mu, sigma = sigma,
                             rho = rho)
  sigma2 <- function(x) repRow(sigma^2, nrow(x))
  sigmaxy <- function(x) rep(rho * prod(sigma), nrow(x))
  d <- function(x) dStatWn2D(x = toPiInt(x), alpha = alpha, mu = mu,
                             sigma = sigma, rho = rho)
  rTraj <- function(x0, N, delta) rTrajLangevin(x0 = x0, drift = driftWn2D,
                                                SigDif = Sigma, N = N,
                                                delta = delta, NFine = N * 10,
                                                deltaFine = delta / 10, A = A,
                                                mu = mu, sigma = sigma,
                                                rho = rho, circular = FALSE)

  # Check
  set.seed(623679851)
  q2 <- npCheck2D(samp = unwrapCircSeries(data), delta = delta, b = b,
                  sigma2 = sigma2, sigmaxy = sigmaxy, rTraj = rTraj,
                  hdrift1 = hdrift1, hdrift2 = hdrift2, hdiff1 = hdiff1,
                  hdiff2 = hdiff2, hdiff3 = hdiff3, hdens = hdens, nth = 100)

  g2 <- NULL
  while (is.null(g2)) {
    cat("Saving...\n")
    g2 <- tryCatch(grid.arrange(q2[[1]], q2[[2]], q2[[3]], q2[[4]], ncol = 2,
                                nrow = 2), error = function(e) NULL)
  }
  # ggsave(filename = paste("2d-", k2, ".png", sep = ""),
  #        plot = g2, width = 14, height = 14)

  for (j in 1:4) ggsave(filename = paste("2d-", k2, "-", j, ".pdf", sep = ""),
                        plot = q2[[j]], width = 7, height = 7)

# }
#
# stopCluster(cl)

}
