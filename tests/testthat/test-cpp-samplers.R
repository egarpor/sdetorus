# Euler samplers and their Monte Carlo step-ahead wrappers (euler.cpp).

test_that("stepAheadWn1D matches the final euler1D step (same RNG state)", {

  set.seed(7)
  sa <- stepAheadWn1D(x0 = c(0.5, -1), alpha = 1, mu = 0, sigma = 1, M = 1,
                      N = 20, delta = 0.05)
  set.seed(7)
  eu <- euler1D(x0 = c(0.5, -1), alpha = 1, mu = 0, sigma = 1, N = 20,
                delta = 0.05)
  expect_equal(dim(sa), c(2, 1))
  expect_equal(sa[, 1], eu[, 21], tolerance = 1e-10)
  expect_true(all(sa >= -pi & sa < pi))

})

test_that("stepAheadWn2D matches the final euler2D slice (same RNG state)", {

  A <- alphaToA(alpha = c(1, 1, 0.3), sigma = c(1, 1))
  set.seed(9)
  sa <- stepAheadWn2D(x0 = rbind(c(0.5, -0.5)), mu = c(0, 0), A = A,
                      sigma = c(1, 1), M = 1, N = 15, delta = 0.05)
  set.seed(9)
  eu <- euler2D(x0 = rbind(c(0.5, -0.5)), A = A, mu = c(0, 0), sigma = c(1, 1),
                N = 15, delta = 0.05)
  expect_equal(dim(sa), c(1, 2, 1))
  expect_equal(sa[, , 1], eu[, , 16], tolerance = 1e-10)

})

test_that("stepAheadWn returns M trajectory ends per starting value", {

  sa1 <- stepAheadWn1D(x0 = c(0, 1), alpha = 1, mu = 0, sigma = 1, M = 50,
                       N = 10, delta = 0.05)
  expect_equal(dim(sa1), c(2, 50))

  A <- alphaToA(alpha = c(1, 1, 0), sigma = c(1, 1))
  sa2 <- stepAheadWn2D(x0 = rbind(c(0, 0), c(1, 1)), mu = c(0, 0), A = A,
                       sigma = c(1, 1), M = 30, N = 10, delta = 0.05)
  expect_equal(dim(sa2), c(2, 2, 30))
  expect_true(all(sa2 >= -pi & sa2 < pi))

})

test_that("euler samplers run for the von Mises drift (type = 2)", {

  set.seed(3)
  e1 <- euler1D(x0 = 0, alpha = 1, mu = 0, sigma = 1, N = 20, delta = 0.05,
                type = 2)
  expect_equal(dim(e1), c(1, 21))
  expect_equal(e1[, 1], 0)

  A <- alphaToA(alpha = c(1, 1, 0), sigma = c(1, 1))
  e2 <- euler2D(x0 = rbind(c(0, 0)), A = A, mu = c(0, 0), sigma = c(1, 1),
                N = 20, delta = 0.05, type = 2)
  expect_equal(dim(e2), c(1, 2, 21))
  expect_equal(e2[, , 1], c(0, 0))

})
