test_that("dTpdPde1D approximates the WN transition density", {

  Mx <- 100
  grid <- seq(-pi, pi, l = Mx + 1)[-(Mx + 1)]
  x0 <- pi
  t <- 0.5
  alpha <- 1
  sigma <- 1

  pde <- dTpdPde1D(Mx = Mx, x0 = x0, t = t, alpha = alpha, mu = 0,
                   sigma = sigma)
  wou <- dTpdWou1D(x = grid, x0 = rep(x0, Mx), t = t, alpha = alpha, mu = 0,
                   sigma = sigma)

  # Both integrate to one and are close to each other
  expect_equal(periodicTrapRule1D(fx = pde), 1, tolerance = 1e-2)
  expect_equal(max(abs(pde - wou)), 0, tolerance = 5e-2)

})

test_that("dTpdPde2D returns a non-negative density with the right shape", {

  # On a fine grid the Crank-Nicolson solution approximately conserves mass;
  # on a coarse grid it only does so loosely (see the sdInitial caveat in the
  # documentation), so the mass check is intentionally lenient.
  M <- 40
  pde <- dTpdPde2D(Mx = M, My = M, x0 = c(0, pi), t = 1, alpha = c(1, 1, 0.5),
                   mu = c(pi / 2, 0), sigma = 1:2, sdInitial = 0.25)
  expect_equal(dim(pde), c(M, M))
  expect_true(all(pde >= -1e-8))
  expect_equal(periodicTrapRule2D(fxy = pde), 1, tolerance = 0.15)

})

test_that("mlePde2D runs with Mx != My (regression for asymmetric grids)", {

  # The closest-bin computation used to wrap the y coordinate with Mx instead
  # of My; with Mx != My this exercises the corrected code path.
  skip_on_cran()
  set.seed(2334567)
  data <- rTrajWn2D(x0 = c(0, 0), alpha = c(1, 0.5, 0.25), mu = c(0, 0),
                    sigma = c(2, 1), N = 30, delta = 0.5)
  sigma <- c(2, 1)
  b <- function(x, pars) driftWn2D(x = x, A = alphaToA(alpha = pars[1:3],
                                                       sigma = sigma),
                                   mu = pars[4:5], sigma = sigma)
  sigma2 <- function(x, pars) repRow(sigma^2, nrow(x))

  fit <- mlePde2D(data = data, delta = 0.5, b = b, sigma2 = sigma2,
                  Mx = 12, My = 8, Mt = 3, start = rbind(c(1, 1, 0, 1, 1)),
                  lower = c(0.1, 0.1, -25, -25, -25),
                  upper = c(25, 25, 25, 25, 25), maxit = 2,
                  selectSolution = "lowest")
  expect_length(fit$par, 5)
  expect_true(all(is.finite(fit$par)))

})
