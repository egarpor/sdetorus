test_that("solveTridiag solves tridiagonal systems", {

  set.seed(1)
  n <- 10
  a <- rnorm(n, 3, 1)
  b <- rnorm(n, 10, 1)
  cc <- rnorm(n, 0, 1)
  d <- rnorm(n, 0, 1)
  A <- matrix(0, nrow = n, ncol = n)
  diag(A) <- b
  for (i in 1:(n - 1)) {
    A[i + 1, i] <- a[i + 1]
    A[i, i + 1] <- cc[i]
  }

  expect_equal(drop(solveTridiag(a = a, b = b, c = cc, d = d)),
               drop(solve(a = A, b = d)), tolerance = 1e-8)

  # LU presaving path
  LU <- forwardSweepTridiag(a = a, b = b, c = cc)
  expect_equal(drop(solveTridiag(a = a, b = LU[, 1], c = LU[, 2], d = d,
                                 LU = 1)),
               drop(solve(a = A, b = d)), tolerance = 1e-8)

  # Several constant vectors at once
  expect_equal(unname(solveTridiagMatConsts(a = a, b = b, c = cc,
                                            d = cbind(d, d + 1))),
               unname(cbind(solve(A, d), solve(A, d + 1))), tolerance = 1e-8)

})

test_that("solvePeriodicTridiag solves circulant tridiagonal systems", {

  set.seed(2)
  n <- 10
  a <- rnorm(n, 3, 1)
  b <- rnorm(n, 10, 1)
  cc <- rnorm(n, 0, 1)
  d <- rnorm(n, 0, 1)
  A <- matrix(0, nrow = n, ncol = n)
  diag(A) <- b
  for (i in 1:(n - 1)) {
    A[i + 1, i] <- a[i + 1]
    A[i, i + 1] <- cc[i]
  }
  A[1, n] <- a[1]
  A[n, 1] <- cc[n]

  expect_equal(drop(solvePeriodicTridiag(a = a, b = b, c = cc, d = d)),
               drop(solve(a = A, b = d)), tolerance = 1e-8)

  LU <- forwardSweepPeriodicTridiag(a = a, b = b, c = cc)
  expect_equal(drop(solvePeriodicTridiag(a = a, b = LU[, 1], c = LU[, 2],
                                         d = d, LU = 1)),
               drop(solve(a = A, b = d)), tolerance = 1e-8)

})

test_that("safeSoftMax matches a direct softmax computation", {

  m <- rbind(1:10, 20:11)
  expect_equal(safeSoftMax(m),
               rbind(exp(1:10) / sum(exp(1:10)),
                     exp(20:11) / sum(exp(20:11))), tolerance = 1e-12)
  # Rows sum to one (no truncation for moderate ranges)
  expect_equal(rowSums(safeSoftMax(rbind(c(0, 1, 2)))), 1)

})
