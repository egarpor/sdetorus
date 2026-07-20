test_that("logBesselI0Scaled agrees with the direct computation", {

  x <- seq(0, 1e3, l = 500)
  expect_equal(logBesselI0Scaled(x = x, splineApprox = TRUE),
               logBesselI0Scaled(x = x, splineApprox = FALSE),
               tolerance = 1e-4)

})

test_that("a1Inv inverts A1 and flags values outside its image", {

  # A1(k) = I1(k) / I0(k)
  k <- c(0.5, 1, 2, 5)
  a1 <- besselI(k, nu = 1, expon.scaled = TRUE) /
    besselI(k, nu = 0, expon.scaled = TRUE)
  expect_equal(a1Inv(x = a1, splineApprox = FALSE), k, tolerance = 1e-4)

  # x >= 1 is outside the image of A1
  expect_message(res <- a1Inv(x = c(0.5, 1.5)),
                 "not on the image")
  expect_true(is.infinite(res[2]))

})

test_that("scoreMatchWnVm and momentMatchWnVm are positive and finite", {

  sigma <- c(0.25, 0.5, 1, 2)
  km <- momentMatchWnVm(sigma = sigma)
  ks <- scoreMatchWnVm(sigma = sigma)
  expect_true(all(km >= 0 & is.finite(km)))
  expect_true(all(ks >= 0 & is.finite(ks)))

})

test_that("scoreMatchWnBvm returns c(0, 0, 0) for a non-invertible invSigma", {

  invSigma <- matrix(1, nrow = 2, ncol = 2)
  expect_message(kappa <- scoreMatchWnBvm(invSigma = invSigma),
                 "non-invertible")
  expect_equal(kappa, c(0, 0, 0))

})

test_that("scoreMatchWnBvm returns c(0, 0, 0) when W is singular (regression)", {

  # Sigma with a zero variance yields a singular matrix W in the score-matching
  # linear system. Previously this returned a scalar NA instead of c(0, 0, 0).
  Sigma <- matrix(c(0, 0, 0, 1), nrow = 2, ncol = 2)
  expect_message(kappa <- scoreMatchWnBvm(Sigma = Sigma), "singular")
  expect_equal(kappa, c(0, 0, 0))

})

test_that("scoreMatchWnBvm returns a length-3 vector for a regular input", {

  Sigma <- rbind(c(0.5, 0.1), c(0.1, 0.4))
  kappa <- scoreMatchWnBvm(Sigma = Sigma)
  expect_length(kappa, 3)
  expect_true(all(is.finite(kappa)))

})
