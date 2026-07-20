test_that("alphaToA and aToAlpha are inverse of each other", {

  alpha <- 3:1
  Sigma <- rbind(c(1, 0.5), c(0.5, 4))
  A <- alphaToA(alpha = alpha, Sigma = Sigma)

  expect_equal(aToAlpha(A = A, Sigma = Sigma), as.numeric(alpha))
  expect_equal(alphaToA(alpha = aToAlpha(A = A, Sigma = Sigma), Sigma = Sigma),
               A)

  # solve(A) %*% Sigma is symmetric by construction
  S <- solve(A) %*% Sigma
  expect_equal(S, t(S))

})

test_that("OU transition moments match their closed forms", {

  alpha <- 1.3
  mu <- -0.4
  sigma <- 0.8
  t <- 0.7
  x0 <- 2

  expect_equal(meantOu(x0 = x0, t = t, alpha = alpha, mu = mu),
               mu + (x0 - mu) * exp(-alpha * t))
  expect_equal(vartOu(t = t, alpha = alpha, sigma = sigma),
               sigma^2 / (2 * alpha) * (1 - exp(-2 * alpha * t)))
  # covstOu diagonal equals vartOu
  cov <- covstOu(s = t, t = t, alpha = alpha, sigma = sigma)
  expect_equal(as.numeric(cov), vartOu(t = t, alpha = alpha, sigma = sigma))

})

test_that("univariate and multivariate OU tpds agree", {

  # A diagonal bivariate OU factorizes into two independent univariate OUs.
  alpha <- c(1, 2)
  mu <- c(0.5, -0.5)
  sigma <- c(1, 1.5)
  t <- 0.4
  x0 <- c(1, -1)
  x <- rbind(c(0, 0), c(0.5, -0.5), c(1, 1))

  A <- diag(alpha)
  Sigma <- diag(sigma^2)
  dMulti <- dTpdMou(x = x, x0 = x0, t = t, A = A, mu = mu, Sigma = Sigma)
  dUni <- dTpdOu(x = x[, 1], x0 = x0[1], t = t, alpha = alpha[1], mu = mu[1],
                 sigma = sigma[1]) *
    dTpdOu(x = x[, 2], x0 = x0[2], t = t, alpha = alpha[2], mu = mu[2],
           sigma = sigma[2])
  expect_equal(as.numeric(dMulti), as.numeric(dUni), tolerance = 1e-8)

})

test_that("mleOu recovers OU parameters on simulated data", {

  set.seed(345678)
  data <- rTrajOu(x0 = 0, alpha = 1, mu = 0, sigma = 1, N = 500, delta = 0.1)
  fit <- mleOu(data = data, delta = 0.1, start = c(2, 1, 2),
               lower = c(0.1, -10, 0.1), upper = c(25, 10, 25))
  expect_length(fit$par, 3)
  expect_true(all(is.finite(fit$par)))
  # sigma is the best-identified parameter at this sampling frequency
  expect_equal(fit$par[3], 1, tolerance = 0.2)

})
