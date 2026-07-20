test_that("rTrajOu returns x0 and the right length", {

  set.seed(345678)
  traj <- rTrajOu(x0 = 0, alpha = 1, mu = 0, sigma = 1, N = 100, delta = 0.1)
  expect_length(traj, 101)
  expect_equal(traj[1], 0)

})

test_that("rTrajMou has the right shape and is reproducible", {

  A <- alphaToA(alpha = c(1, 2, 0.5), sigma = 1:2)
  args <- list(x0 = c(0, 0), A = A, mu = c(1, 1), Sigma = diag((1:2)^2),
               N = 50, delta = 0.1)

  set.seed(42)
  m1 <- do.call(rTrajMou, args)
  set.seed(42)
  m2 <- do.call(rTrajMou, args)

  expect_equal(dim(m1), c(51, 2))
  expect_equal(m1[1, ], c(0, 0))
  expect_identical(m1, m2)

})

test_that("rTrajMou matches an rmvnorm-based reference (efficiency rewrite)", {

  # Reference implementation that draws each step with mvtnorm::rmvnorm. The
  # optimized rTrajMou precomputes the covariance root once instead, but must
  # produce the identical trajectory for a given RNG state.
  refRTrajMou <- function(x0, A, mu, Sigma, N, delta) {
    p <- ncol(A)
    eigA <- eigen(A)
    covt <- covtMou(t = delta, Sigma = Sigma, eigA = eigA)
    samp <- matrix(x0, nrow = N + 1, ncol = p, byrow = TRUE)
    for (i in 2:(N + 1)) {
      samp[i, ] <- mvtnorm::rmvnorm(n = 1,
        mean = meantMou(t = delta, x0 = samp[i - 1, ], mu = mu, eigA = eigA),
        sigma = covt)
    }
    samp
  }

  A <- alphaToA(alpha = c(1, 2, 0.5), sigma = 1:2)
  args <- list(x0 = c(0, 0), A = A, mu = c(1, 1), Sigma = diag((1:2)^2),
               N = 100, delta = 0.1)

  set.seed(987658)
  fast <- do.call(rTrajMou, args)
  set.seed(987658)
  ref <- do.call(refRTrajMou, args)

  expect_equal(fast, ref, tolerance = 1e-12)

})

test_that("rTrajWn1D and rTrajWn2D return the right shapes", {

  set.seed(1)
  s1 <- rTrajWn1D(x0 = 0, alpha = 1, mu = 0, sigma = 1, N = 100, delta = 0.01)
  expect_length(s1, 101)
  expect_equal(s1[1], 0)

  s2 <- rTrajWn2D(x0 = c(0, 0), alpha = c(1, 1, -0.5), mu = c(pi, pi),
                  sigma = c(1, 1), N = 100, delta = 0.01)
  expect_equal(dim(s2), c(101, 2))
  expect_equal(s2[1, ], c(0, 0))

})

test_that("rTrajLangevin returns the documented shapes for p = 1 and p = 2", {

  set.seed(1)
  # p = 1: a vector of length N + 1
  s1 <- rTrajLangevin(x0 = 0, drift = driftJp, SigDif = 1, alpha = 1, mu = 0,
                      psi = 1, N = 50, delta = 0.1)
  expect_length(s1, 51)
  expect_equal(s1[1], 0)
  expect_true(all(s1 >= -pi & s1 < pi))

  # p = 2: a matrix of size c(N + 1, 2)
  s2 <- rTrajLangevin(x0 = c(0, 0), drift = driftMvm, alpha = c(1, 1),
                      mu = c(1, -1), A = diag(0, 2), SigDif = diag(2), N = 50,
                      delta = 0.1)
  expect_equal(dim(s2), c(51, 2))
  expect_equal(s2[1, ], c(0, 0))
  expect_true(all(s2 >= -pi & s2 < pi))

  # circular = FALSE leaves the trajectory unwrapped
  s3 <- rTrajLangevin(x0 = 0, drift = driftJp, SigDif = 4, alpha = 1, mu = 0,
                      psi = 1, N = 50, delta = 0.1, circular = FALSE)
  expect_length(s3, 51)

})
