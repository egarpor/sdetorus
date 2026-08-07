test_that("periodicTrapRule and integrateSimp match known 1D integral", {

  # True value of int_{-pi}^{pi} sin(x)^2 * exp(cos(x)) dx
  true <- 3.55099937
  N <- 21
  grid <- seq(-pi, pi, l = N)
  fx <- sin(grid)^2 * exp(cos(grid))

  expect_equal(periodicTrapRule1D(fx = fx, endsMatch = TRUE), true,
               tolerance = 1e-4)
  expect_equal(periodicTrapRule1D(fx = fx[-N], endsMatch = FALSE), true,
               tolerance = 1e-4)
  expect_equal(integrateSimp1D(fx = fx, lengthInterval = 2 * pi), true,
               tolerance = 1e-3)

})

test_that("periodicTrapRule and integrateSimp match known 2D integral", {

  true <- 22.31159
  N <- 21
  grid <- seq(-pi, pi, l = N)
  fxy <- outer(grid, grid, function(x, y) (sin(x)^2 * exp(cos(x)) +
                                             sin(y)^2 * exp(cos(y))) / 2)

  expect_equal(periodicTrapRule2D(fxy = fxy, endsMatch = TRUE), true,
               tolerance = 1e-3)
  expect_equal(integrateSimp2D(fxy = fxy), true, tolerance = 1e-2)

})

test_that("integrateSimp3D matches known 3D integral", {

  true <- 140.1878
  N <- 21
  grid <- seq(-pi, pi, l = N)
  fxy <- outer(grid, grid, function(x, y) (sin(x)^2 * exp(cos(x)) +
                                             sin(y)^2 * exp(cos(y))) / 2)
  fxyz <- array(0, dim = c(N, N, N))
  for (i in 1:N) fxyz[i, , ] <- fxy

  expect_equal(periodicTrapRule3D(fxyz = fxyz, endsMatch = TRUE), true,
               tolerance = 1e-2)
  expect_equal(integrateSimp3D(fxyz = fxyz), true, tolerance = 1e-1)

})

test_that("mcTorusIntegrate approximates a known torus integral", {

  set.seed(123)
  # int over the torus of sin(x1) * cos(x2) is 0
  expect_equal(mcTorusIntegrate(f = function(x) sin(x[, 1]) * cos(x[, 2]),
                                p = 2), 0, tolerance = 0.1)
  # int over the torus of a constant c is c * (2 * pi)^p
  expect_equal(mcTorusIntegrate(f = function(x) rep(2, nrow(x)), p = 2),
               2 * (2 * pi)^2, tolerance = 1e-8)

})
