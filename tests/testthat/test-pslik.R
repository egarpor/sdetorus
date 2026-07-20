# Pseudo-likelihood block: dPsTpd (wrapped Euler / Shoji-Ozaki pseudo-tpds)
# and psMle (its maximizer).

test_that("dPsTpd densities are non-negative and integrate to one", {

  alpha <- 1
  sigma <- 1
  mu <- 0

  # Drift and its finite-difference derivatives (closures over the parameters)
  b <- function(x) driftWn1D(x = x, alpha = alpha, mu = mu, sigma = sigma)
  b1 <- function(x, h = 1e-4) {
    l <- length(x)
    res <- driftWn1D(c(x + h, x - h), alpha = alpha, mu = mu, sigma = sigma)
    drop(res[1:l] - res[(l + 1):(2 * l)]) / (2 * h)
  }
  b2 <- function(x, h = 1e-4) {
    l <- length(x)
    res <- driftWn1D(c(x + h, x, x - h), alpha = alpha, mu = mu, sigma = sigma)
    drop(res[1:l] - 2 * res[(l + 1):(2 * l)] + res[(2 * l + 1):(3 * l)]) / (h^2)
  }
  sigma2 <- function(x) rep(sigma^2, length(x))

  grid <- seq(-pi, pi, l = 501)[-501]
  x0 <- pi / 2

  for (m in c("E", "SO", "SO2")) {
    dens <- dPsTpd(x = grid, x0 = x0, t = 0.5, method = m, b = b, b1 = b1,
                   b2 = b2, sigma2 = sigma2)
    expect_true(all(dens >= 0))
    expect_equal(periodicTrapRule1D(fx = dens), 1, tolerance = 5e-2)
  }

})

test_that("dPsTpd (Euler) approaches dTpdWou1D as t decreases", {

  alpha <- 1
  sigma <- 1
  mu <- 0
  b <- function(x) driftWn1D(x = x, alpha = alpha, mu = mu, sigma = sigma)
  b1 <- function(x, h = 1e-4) {
    l <- length(x)
    res <- driftWn1D(c(x + h, x - h), alpha = alpha, mu = mu, sigma = sigma)
    drop(res[1:l] - res[(l + 1):(2 * l)]) / (2 * h)
  }
  b2 <- function(x, h = 1e-4) {
    l <- length(x)
    res <- driftWn1D(c(x + h, x, x - h), alpha = alpha, mu = mu, sigma = sigma)
    drop(res[1:l] - 2 * res[(l + 1):(2 * l)] + res[(2 * l + 1):(3 * l)]) / (h^2)
  }
  sigma2 <- function(x) rep(sigma^2, length(x))

  grid <- seq(-pi, pi, l = 501)[-501]
  x0 <- pi / 2

  dif <- sapply(c(0.5, 0.1), function(tt) {
    dE <- dPsTpd(x = grid, x0 = x0, t = tt, method = "E", b = b, b1 = b1,
                 b2 = b2, sigma2 = sigma2)
    dW <- dTpdWou1D(x = grid, x0 = rep(x0, length(grid)), t = tt, alpha = alpha,
                    mu = mu, sigma = sigma)
    max(abs(dE - dW))
  })
  # The Euler pseudo-density converges to the exact tpd as t -> 0
  expect_lt(dif[2], dif[1])
  expect_lt(dif[2], 0.15)

})

test_that("dPsTpd accepts vector and matrix inputs equivalently (p = 1)", {

  alpha <- 1
  sigma <- 1
  mu <- 0
  b <- function(x) driftWn1D(x = x, alpha = alpha, mu = mu, sigma = sigma)
  b1 <- function(x, h = 1e-4) {
    l <- length(x)
    res <- driftWn1D(c(x + h, x - h), alpha = alpha, mu = mu, sigma = sigma)
    drop(res[1:l] - res[(l + 1):(2 * l)]) / (2 * h)
  }
  b2 <- function(x, h = 1e-4) {
    l <- length(x)
    res <- driftWn1D(c(x + h, x, x - h), alpha = alpha, mu = mu, sigma = sigma)
    drop(res[1:l] - 2 * res[(l + 1):(2 * l)] + res[(2 * l + 1):(3 * l)]) / (h^2)
  }
  sigma2 <- function(x) rep(sigma^2, length(x))

  grid <- seq(-pi, pi, l = 200)[-200]
  x0 <- pi / 2
  dv <- dPsTpd(x = grid, x0 = x0, t = 0.5, method = "E", b = b, b1 = b1,
               b2 = b2, sigma2 = sigma2)
  dm <- dPsTpd(x = cbind(grid), x0 = cbind(x0), t = 0.5, method = "E", b = b,
               b1 = b1, b2 = b2, sigma2 = sigma2)
  expect_equal(as.numeric(dv), as.numeric(dm))

})

test_that("psMle returns a finite estimate (smoke)", {

  # Note: psMle uses `method == "SO"`, so a single method string must be passed
  # (the default length-3 vector errors under R >= 4.3).
  set.seed(1)
  samp <- rTrajWn1D(x0 = 0, alpha = 0.5, mu = 0, sigma = 2, N = 40, delta = 0.5)
  b <- function(x, pars) driftWn1D(x = x, alpha = pars[1], mu = pars[2],
                                   sigma = pars[3])
  sigma2 <- function(x, pars) rep(pars[3]^2, length(x))

  fit <- psMle(data = samp, delta = 0.5, method = "E", b = b, sigma2 = sigma2,
               start = c(0.5, 0, 2), lower = c(0.1, -pi, 0.1),
               upper = c(10, pi, 10), selectSolution = "lowest", maxit = 20)
  expect_length(fit$par, 3)
  expect_true(all(is.finite(fit$par)))

})

test_that("psMle recovers the diffusion coefficient of a WN diffusion", {

  skip_on_cran()
  set.seed(12345678)
  samp <- rTrajWn1D(x0 = 0, alpha = 0.25, mu = 0, sigma = 2, N = 100,
                    delta = 0.5)
  b <- function(x, pars) driftWn1D(x = x, alpha = pars[1], mu = pars[2],
                                   sigma = pars[3])
  sigma2 <- function(x, pars) rep(pars[3]^2, length(x))

  fit <- psMle(data = samp, delta = 0.5, method = "E", b = b, sigma2 = sigma2,
               start = c(0.25, 0, 2), lower = c(0.1, -pi, 0.1),
               upper = c(10, pi, 10), selectSolution = "lowest")
  # sigma is the best-identified parameter at this sampling frequency
  expect_equal(fit$par[3], 2, tolerance = 0.3)

})

test_that("dPsTpd von Mises approximation gives a proper 1D density", {

  alpha <- 1
  sigma <- 1
  mu <- 0
  b <- function(x) driftWn1D(x = x, alpha = alpha, mu = mu, sigma = sigma)
  b1 <- function(x, h = 1e-4) {
    l <- length(x)
    res <- driftWn1D(c(x + h, x - h), alpha = alpha, mu = mu, sigma = sigma)
    drop(res[1:l] - res[(l + 1):(2 * l)]) / (2 * h)
  }
  b2 <- function(x, h = 1e-4) {
    l <- length(x)
    res <- driftWn1D(c(x + h, x, x - h), alpha = alpha, mu = mu, sigma = sigma)
    drop(res[1:l] - 2 * res[(l + 1):(2 * l)] + res[(2 * l + 1):(3 * l)]) / (h^2)
  }
  sigma2 <- function(x) rep(sigma^2, length(x))
  grid <- seq(-pi, pi, l = 301)[-301]

  for (m in c("E", "SO", "SO2")) {
    dens <- dPsTpd(x = grid, x0 = pi / 2, t = 0.3, method = m, b = b, b1 = b1,
                   b2 = b2, sigma2 = sigma2, vmApprox = TRUE)
    expect_true(all(dens >= 0))
    expect_equal(periodicTrapRule1D(fx = dens), 1, tolerance = 5e-2)
  }

})

test_that("dPsTpd handles the bivariate case (p = 2)", {

  alpha <- c(2, 1, 0.5)
  sigma <- c(1, 2)
  x0 <- c(pi / 2, -pi / 2)
  t <- 0.3
  A <- alphaToA(alpha = alpha, sigma = sigma)
  b <- function(x) driftWn2D(x = x, A = A, mu = rep(0, 2), sigma = sigma)
  jac.b <- function(x, h = 1e-4) {
    l <- nrow(x)
    res <- driftWn2D(x = rbind(cbind(x[, 1] + h, x[, 2]),
                               cbind(x[, 1] - h, x[, 2]),
                               cbind(x[, 1], x[, 2] + h),
                               cbind(x[, 1], x[, 2] - h)),
                     A = A, mu = rep(0, 2), sigma = sigma)
    cbind(res[1:l, ] - res[(l + 1):(2 * l), ],
          res[2 * l + 1:l, ] - res[2 * l + (l + 1):(2 * l), ]) / (2 * h)
  }
  # In the SO p >= 2 branch, sigma2 is called with a single-row vector, hence
  # the length(x) / 2 (not nrow(x)) form
  sigma2 <- function(x) matrix(sigma^2, nrow = length(x) / 2L, ncol = 2)

  g <- seq(-pi, pi, l = 61)[-61]
  xy <- as.matrix(expand.grid(g, g))
  for (m in c("E", "SO")) {
    dens <- dPsTpd(x = xy, x0 = rbind(x0), t = t, method = m, b = b,
                   jac.b = jac.b, sigma2 = sigma2)
    expect_true(all(dens >= 0))
    expect_equal(periodicTrapRule2D(matrix(dens, length(g), length(g))), 1,
                 tolerance = 5e-2)
  }

  # p = 2 SO with the bivariate von Mises approximation
  densVm <- dPsTpd(x = xy, x0 = rbind(x0), t = t, method = "SO", b = b,
                   jac.b = jac.b, sigma2 = sigma2, vmApprox = TRUE)
  expect_equal(periodicTrapRule2D(matrix(densVm, length(g), length(g))), 1,
               tolerance = 5e-2)

  # x0 supplied with one row per evaluation point exercises the per-row branch
  x0Rep <- repRow(x0, nrow(xy))
  densRep <- dPsTpd(x = xy, x0 = x0Rep, t = t, method = "SO", b = b,
                    jac.b = jac.b, sigma2 = sigma2)
  expect_equal(periodicTrapRule2D(matrix(densRep, length(g), length(g))), 1,
               tolerance = 5e-2)
  densRepVm <- dPsTpd(x = xy, x0 = x0Rep, t = t, method = "SO", b = b,
                      jac.b = jac.b, sigma2 = sigma2, vmApprox = TRUE)
  expect_equal(periodicTrapRule2D(matrix(densRepVm, length(g), length(g))), 1,
               tolerance = 5e-2)

})

test_that("psMle accepts a vector data argument (p = 1)", {

  set.seed(1)
  samp <- rTrajWn1D(x0 = 0, alpha = 0.5, mu = 0, sigma = 2, N = 40, delta = 0.5)
  b <- function(x, pars) driftWn1D(x = x, alpha = pars[1], mu = pars[2],
                                   sigma = pars[3])
  sigma2 <- function(x, pars) rep(pars[3]^2, length(x))
  # Pass the trajectory as a plain vector; psMle reshapes it internally
  fit <- psMle(data = as.vector(samp), delta = 0.5, method = "E", b = b,
               sigma2 = sigma2, start = c(0.5, 0, 2), lower = c(0.1, -pi, 0.1),
               upper = c(10, pi, 10), selectSolution = "lowest", maxit = 20)
  expect_length(fit$par, 3)
  expect_true(all(is.finite(fit$par)))

})
