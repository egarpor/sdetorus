# WN diffusion in 2D: the stationary density and the tpd sampler (WN-2D.cpp).

test_that("rTpdWn2D has the right shape, support and conditional mean", {

  set.seed(1)
  x0 <- rbind(c(0, 0), c(1, -1))
  s <- rTpdWn2D(n = 500, x0 = x0, t = 0.5, mu = c(0, 0), alpha = c(1, 1, 0.3),
                sigma = c(1, 1))
  expect_equal(dim(s), c(500, 2, 2))
  expect_true(all(s >= -pi & s < pi))

  # For x0 = mu = (0, 0), the conditional mean stays near 0
  cmean <- c(atan2(mean(sin(s[, 1, 1])), mean(cos(s[, 1, 1]))),
             atan2(mean(sin(s[, 2, 1])), mean(cos(s[, 2, 1]))))
  expect_lt(max(abs(cmean)), 0.2)

})

test_that("dStatWn2D integrates to one and factorizes when independent", {

  grid <- seq(-pi, pi, l = 101)[-101]
  xy <- as.matrix(expand.grid(grid, grid))
  alpha <- c(1, 2, 0.5)
  mu <- c(0.5, -0.5)
  sigma <- c(1, 1.5)

  dens <- dStatWn2D(x = xy, alpha = alpha, mu = mu, sigma = sigma, rho = 0.3)
  expect_true(all(dens >= 0))
  expect_equal(periodicTrapRule2D(matrix(dens, length(grid), length(grid))), 1,
               tolerance = 1e-3)

  # alpha3 = 0 and rho = 0 => product of wrapped-normal marginals with
  # stationary variance sigma^2 / (2 * alpha)
  df <- dStatWn2D(x = xy, alpha = c(1, 2, 0), mu = mu, sigma = sigma, rho = 0)
  prodMarg <- dWn1D(xy[, 1], mu = mu[1], sigma = sigma[1] / sqrt(2 * 1)) *
              dWn1D(xy[, 2], mu = mu[2], sigma = sigma[2] / sqrt(2 * 2))
  expect_equal(as.numeric(df), as.numeric(prodMarg), tolerance = 1e-6)

})

test_that("rStatWn2D samples lie on the torus with the right circular mean", {

  set.seed(2)
  x <- rStatWn2D(n = 2000, mu = c(1, -1), alpha = c(2, 2, 0.5), sigma = c(1, 1))
  expect_equal(dim(x), c(2000, 2))
  expect_true(all(x >= -pi & x < pi))
  cm <- c(atan2(mean(sin(x[, 1])), mean(cos(x[, 1]))),
          atan2(mean(sin(x[, 2])), mean(cos(x[, 2]))))
  expect_equal(cm, c(1, -1), tolerance = 0.15)

})

test_that("driftWn2D and dTpdWou2D handle rho != 0 and larger maxK", {

  x <- rbind(c(0, 1), c(pi, -pi), c(1, 0.5))
  A <- alphaToA(alpha = c(2, 1, 0.5), sigma = c(1, 2), rho = 0.5)
  d3 <- driftWn2D(x = x, A = A, mu = c(0, 0), sigma = c(1, 2), rho = 0.5,
                  maxK = 3)
  expect_equal(dim(d3), c(3, 2))
  expect_true(all(is.finite(d3)))

  p <- dTpdWou2D(x = x, x0 = rbind(c(0, 0)), t = 0.5, alpha = c(2, 1, 0.5),
                 mu = c(0, 0), sigma = c(1, 2), rho = 0.5, maxK = 3)
  expect_true(all(p >= 0 & is.finite(p)))

})
