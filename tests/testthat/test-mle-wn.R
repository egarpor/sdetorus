test_that("dTpdWou agrees with dTpdWou1D and dTpdWou2D", {

  # 1D
  x <- seq(-pi, pi, l = 10)
  d1 <- dTpdWou(x = cbind(x), x0 = pi, t = 0.5, A = 1, mu = 0, Sigma = 1)
  d1ref <- dTpdWou1D(x = cbind(x), x0 = rep(pi, 10), t = 0.5, alpha = 1,
                     mu = 0, sigma = 1)
  expect_equal(as.numeric(d1), as.numeric(d1ref), tolerance = 1e-8)

  # 2D
  alpha <- c(2, 1, -1)
  sigma <- c(1.5, 2)
  rho <- 0.5
  Sigma <- diag(sigma^2)
  Sigma[1, 2] <- Sigma[2, 1] <- rho * prod(sigma)
  A <- alphaToA(alpha = alpha, sigma = sigma, rho = rho)
  mu <- c(pi, 0)
  x0 <- c(0, 0)
  xx <- as.matrix(expand.grid(seq(-pi, pi, l = 5), seq(-pi, pi, l = 5)))
  d2 <- dTpdWou(x = xx, x0 = x0, t = 0.5, A = A, mu = mu, Sigma = Sigma)
  d2ref <- dTpdWou2D(x = xx, x0 = rbind(x0), t = 0.5, alpha = alpha, mu = mu,
                     sigma = sigma, rho = rho)
  expect_equal(as.numeric(d2), as.numeric(d2ref), tolerance = 1e-8)

})

test_that("approxMleWn1D returns a finite estimate", {

  set.seed(1)
  samp <- rTrajWn1D(x0 = 0, alpha = 0.5, mu = 0, sigma = 2, N = 500,
                    delta = 0.1)
  fit <- approxMleWn1D(data = samp, delta = 0.1, start = c(0.5, 0, 2),
                       selectSolution = "lowest")
  expect_length(fit$par, 3)
  expect_true(all(is.finite(fit$par)))

})

test_that("approxMleWn2D runs with the positive-definiteness region active", {

  # The region function is now passed to mleOptimWrapper (regression: it used to
  # be dead code containing an NA-producing term). It must run without error and
  # return a finite, positive-definite estimate.
  set.seed(20240720)
  samp <- rTrajWn2D(x0 = c(0, 0), alpha = c(2, 2, -0.5), mu = c(0, 0),
                    sigma = c(1, 1), rho = 0.2, N = 500, delta = 0.1)
  fit <- suppressWarnings(
    approxMleWn2D(data = samp, delta = 0.1,
                  start = c(2, 2, -0.5, 0, 0, 1, 1, 0.2),
                  selectSolution = "lowest"))
  expect_length(fit$par, 8)
  expect_true(all(is.finite(fit$par)))

  # The estimate lies in the feasibility region enforced by 'region'
  pars <- fit$par
  testPosDef <- 0.25 * (pars[8] * (pars[2] - pars[1]))^2 + pars[1] * pars[2] -
    pars[3]^2
  expect_gte(testPosDef, 0)

})

test_that("sigmaDiff estimates the diffusion matrix and honors constraints", {

  set.seed(1)
  x <- drop(euler1D(x0 = 0, alpha = 1, mu = 0, sigma = 1.5, N = 2000,
                    delta = 0.01))
  expect_equal(as.numeric(sigmaDiff(x, delta = 0.01)), 1.5^2, tolerance = 0.2)

  set.seed(2)
  x2 <- t(euler2D(x0 = rbind(c(pi, pi)), A = rbind(c(2, 1), c(1, 2)),
                  mu = c(pi, pi), sigma = c(1, 1), N = 2000,
                  delta = 0.01)[1, , ])
  siso <- sigmaDiff(x2, delta = 0.01, isotropic = TRUE)
  expect_equal(siso[1, 1], siso[2, 2])
  sdiag <- sigmaDiff(x2, delta = 0.01, diagonal = TRUE)
  expect_equal(sdiag[1, 2], 0)

})
