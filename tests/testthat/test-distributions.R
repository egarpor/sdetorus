test_that("dVm integrates to one and is symmetric around mu", {

  grid <- seq(-pi, pi, l = 1001)[-1001]
  for (kappa in c(0, 0.5, 2, 10)) {
    dens <- dVm(x = grid, mu = 0.4, kappa = kappa)
    expect_equal(periodicTrapRule1D(fx = dens), 1, tolerance = 1e-4)
    expect_true(all(dens >= 0))
  }
  # kappa = 0 is the uniform density
  expect_equal(dVm(x = grid, mu = 0, kappa = 0), rep(1 / (2 * pi), length(grid)),
               tolerance = 1e-10)

})

test_that("dJp integrates to one and matches dVm for psi = 0", {

  grid <- seq(-pi, pi, l = 1001)[-1001]
  expect_equal(dJp(x = grid, mu = 0, kappa = 1, psi = 0),
               dVm(x = grid, mu = 0, kappa = 1))
  for (psi in c(-1, 0.5, 1)) {
    dens <- dJp(x = grid, mu = 0, kappa = 1, psi = psi)
    expect_equal(periodicTrapRule1D(fx = dens), 1, tolerance = 1e-3)
  }

})

test_that("dWn1D integrates to one", {

  grid <- seq(-pi, pi, l = 1001)[-1001]
  for (sigma in c(0.3, 1, 2)) {
    dens <- dWn1D(x = grid, mu = 0.2, sigma = sigma)
    expect_equal(periodicTrapRule1D(fx = dens), 1, tolerance = 1e-4)
  }

})

test_that("dBvm integrates to one over the torus", {

  x <- seq(-pi, pi, l = 121)[-121]
  xy <- as.matrix(expand.grid(x, x))
  dens <- dBvm(x = xy, mu = c(0, pi / 2), kappa = c(2, 3, 1))
  fxy <- matrix(dens, nrow = length(x), ncol = length(x))
  expect_equal(periodicTrapRule2D(fxy = fxy), 1, tolerance = 1e-3)

})

test_that("constBvm reduces to the product of Bessel functions for lambda = 0", {

  kappa <- c(2, 3, 0)
  expect_equal(constBvm(kappa = kappa),
               4 * pi^2 * besselI(kappa[1], nu = 0) * besselI(kappa[2], nu = 0))

})
