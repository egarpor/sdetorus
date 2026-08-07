# Internal toroidal von Mises mixtures (dTvm / emTvm) and the C++ they drive.

test_that("dTvm reduces to a von Mises for a single 1D component", {

  grid <- seq(-pi, pi, l = 501)[-501]
  d <- sdetorus:::dTvm(x = grid, M = 0.4, K = 3, alpha = 1)
  expect_equal(d, dVm(x = grid, mu = 0.4, kappa = 3))
  expect_equal(periodicTrapRule1D(fx = d), 1, tolerance = 1e-4)

})

test_that("dTvm densities are non-negative and integrate to one", {

  grid <- seq(-pi, pi, l = 501)[-501]

  # Two-component 1D mixture
  dm <- sdetorus:::dTvm(x = grid, M = rbind(-1, 1), K = rbind(2, 4),
                        alpha = c(0.3, 0.7))
  expect_true(all(dm >= 0))
  expect_equal(periodicTrapRule1D(fx = dm), 1, tolerance = 1e-3)

  # 2D single component (product of two von Mises)
  g <- seq(-pi, pi, l = 201)[-201]
  xy <- as.matrix(expand.grid(g, g))
  d2 <- sdetorus:::dTvm(x = xy, M = c(0.5, -0.5), K = c(2, 3), alpha = 1)
  expect_equal(periodicTrapRule2D(matrix(d2, length(g), length(g))), 1,
               tolerance = 1e-3)

})

test_that("dTvm handles vector/matrix input, besselInterp and size checks", {

  grid <- seq(-pi, pi, l = 200)[-200]
  dvec <- sdetorus:::dTvm(x = grid, M = 0.4, K = 3, alpha = 1)
  dmat <- sdetorus:::dTvm(x = cbind(grid), M = 0.4, K = 3, alpha = 1)
  expect_equal(as.numeric(dvec), as.numeric(dmat))

  # Spline Bessel interpolation gives a close result
  dbi <- sdetorus:::dTvm(x = grid, M = 0.4, K = 3, alpha = 1,
                         besselInterp = TRUE)
  expect_equal(dbi, dvec, tolerance = 1e-6)

  # Incompatible sizes of x, M, K, alpha
  expect_error(sdetorus:::dTvm(x = grid, M = rbind(0, 1), K = 3, alpha = 1),
               "Incompatible")

})

test_that("emTvm fits a two-component toroidal mixture", {

  set.seed(42)
  n <- 200
  z <- sample(1:2, n, replace = TRUE)
  centre <- ifelse(z == 1, -1.5, 1.5)
  data <- cbind(toPiInt(centre + rnorm(n, sd = 0.3)),
                toPiInt(centre + rnorm(n, sd = 0.3)))

  em <- sdetorus:::emTvm(data = data, k = 2, maxIter = 30)
  expect_named(em, c("M", "K", "alpha", "pih", "BIC", "convergence"))
  expect_equal(dim(em$M), c(2, 2))
  expect_equal(dim(em$K), c(2, 2))
  expect_length(em$alpha, 2)
  expect_equal(sum(em$alpha), 1, tolerance = 1e-6)
  expect_true(is.finite(em$BIC))

  # Isotropic variant also runs
  emi <- sdetorus:::emTvm(data = data, k = 2, maxIter = 30, isotropic = TRUE)
  expect_equal(dim(emi$K), c(2, 2))

})
