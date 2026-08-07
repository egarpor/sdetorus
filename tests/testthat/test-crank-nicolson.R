# Crank-Nicolson Fokker-Planck solvers (crankNicolson.cpp).

test_that("crankNicolson1D conserves mass and preserves the initial condition", {

  Mx <- 100
  grid <- seq(-pi, pi, l = Mx + 1)[-(Mx + 1)]
  u0 <- matrix(dWn1D(x = grid, mu = 0, sigma = 0.5), ncol = 1)
  b <- driftWn1D(x = grid, alpha = 1, mu = 0, sigma = 1)
  s2 <- rep(1, Mx)
  deltax <- grid[2] - grid[1]

  u <- crankNicolson1D(u0 = u0, b = b, sigma2 = s2, N = 0:10, deltat = 0.05,
                       Mx = Mx, deltax = deltax)
  expect_equal(dim(u), c(Mx, 11))
  expect_lt(max(abs(colMeans(u) * 2 * pi - 1)), 1e-3)     # mass conserved
  expect_equal(u[, 1], u0[, 1])                            # IC preserved

  # A single requested time equals the last column of the trajectory
  v <- crankNicolson1D(u0 = u0, b = b, sigma2 = s2, N = 10, deltat = 0.05,
                       Mx = Mx, deltax = deltax)
  expect_equal(sum(abs(u[, 11] - v[, 1])), 0, tolerance = 1e-8)

  # imposePositive keeps the solution non-negative
  uP <- crankNicolson1D(u0 = u0, b = b, sigma2 = s2, N = 10, deltat = 0.05,
                        Mx = Mx, deltax = deltax, imposePositive = 1)
  expect_gte(min(uP), -1e-12)

})

test_that("crankNicolson2D conserves mass and preserves the initial condition", {

  M <- 24
  grid <- seq(-pi, pi, l = M + 1)[-(M + 1)]
  gg <- as.matrix(expand.grid(grid, grid))
  u0 <- matrix(c(outer(dWn1D(grid, 0, 0.5), dWn1D(grid, 0, 0.5))), ncol = 1)
  bGrid <- driftWn2D(x = gg, A = alphaToA(c(1, 1, 0.3), c(1, 1)), mu = c(0, 0),
                     sigma = c(1, 1))
  bx <- matrix(bGrid[, 1], M, M)
  by <- matrix(bGrid[, 2], M, M)
  s2 <- matrix(1, M, M)
  sxy <- matrix(0, M, M)
  deltax <- grid[2] - grid[1]

  u <- crankNicolson2D(u0 = u0, bx = bx, by = by, sigma2x = s2, sigma2y = s2,
                       sigmaxy = sxy, N = 0:5, deltat = 0.05, Mx = M,
                       deltax = deltax, My = M, deltay = deltax)
  expect_equal(dim(u), c(M * M, 6))
  expect_equal(u[, 1], u0[, 1])                              # IC preserved
  expect_lt(max(abs(colMeans(u) * 4 * pi^2 - 1)), 0.05)      # mass ~ conserved

  uP <- crankNicolson2D(u0 = u0, bx = bx, by = by, sigma2x = s2, sigma2y = s2,
                        sigmaxy = sxy, N = 5, deltat = 0.05, Mx = M,
                        deltax = deltax, My = M, deltay = deltax,
                        imposePositive = 1)
  expect_gte(min(uP), -1e-12)

})
