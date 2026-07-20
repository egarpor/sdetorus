test_that("driftMixVm forwards its expTrc argument (regression)", {

  # Prior to the fix, driftMixVm ignored 'expTrc' (it was hard-coded to 30),
  # so the two calls below returned identical drifts.
  x <- seq(-pi, pi, l = 60)
  alpha <- c(2, 2)
  mu <- c(0, 3)
  sigma <- 1
  p <- c(0.5, 0.5)

  dFull <- driftMixVm(x = x, alpha = alpha, mu = mu, sigma = sigma, p = p,
                      expTrc = 30)
  dTrunc <- driftMixVm(x = x, alpha = alpha, mu = mu, sigma = sigma, p = p,
                       expTrc = 1)
  expect_false(isTRUE(all.equal(dFull, dTrunc)))

})

test_that("driftMixVm agrees with the general driftMixIndVm for any expTrc", {

  x <- seq(-pi, pi, l = 60)
  alpha <- c(2, 2)
  mu <- c(0, 3)
  sigma <- 1
  p <- c(0.5, 0.5)

  for (e in c(30, 5, 1)) {
    a <- driftMixVm(x = x, alpha = alpha, mu = mu, sigma = sigma, p = p,
                    expTrc = e)
    b <- drop(driftMixIndVm(x = cbind(x), A = cbind(alpha), M = cbind(mu),
                            sigma = sigma, p = p, expTrc = e))
    expect_equal(a, b, tolerance = 1e-10)
  }

})

test_that("driftWn1D agrees with the general driftWn in 1D", {

  x <- seq(-pi, pi, l = 50)
  alpha <- 1.5
  mu <- 0.5
  sigma <- 1.2
  d1 <- driftWn1D(x = x, alpha = alpha, mu = mu, sigma = sigma, maxK = 2)
  dG <- driftWn(x = cbind(x), A = alpha, mu = mu, Sigma = sigma^2, maxK = 2)
  expect_equal(as.numeric(d1), as.numeric(dG), tolerance = 1e-8)

})

test_that("driftWn2D agrees with the general driftWn in 2D", {

  alpha <- 3:1
  mu <- c(0, 0)
  sigma <- 1:2
  rho <- 0.5
  Sigma <- diag(sigma^2)
  Sigma[1, 2] <- Sigma[2, 1] <- rho * prod(sigma)
  A <- alphaToA(alpha = alpha, sigma = sigma, rho = rho)
  x <- rbind(c(0, 1), c(1, 0.1), c(pi, pi), c(-pi, -pi), c(pi / 2, 0))
  d2 <- driftWn2D(x = x, A = A, mu = mu, sigma = sigma, rho = rho)
  dG <- driftWn(x = x, A = A, mu = mu, Sigma = Sigma)
  expect_equal(unname(d2), unname(dG), tolerance = 1e-8)

})

test_that("driftJp reduces to the von Mises drift for psi = 0", {

  x <- seq(-pi, pi, l = 50)
  expect_equal(driftJp(x = x, alpha = 1, mu = 0.3, psi = 0), sin(0.3 - x))

})
