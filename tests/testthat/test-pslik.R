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
