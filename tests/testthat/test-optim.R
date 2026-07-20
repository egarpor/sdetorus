test_that("mleOptimWrapper finds the minimum of a convex quadratic", {

  res <- mleOptimWrapper(minusLogLik = function(x) sum((x - 1:4)^2),
                         start = rbind(10:13), optMethod = "BFGS")
  expect_equal(res$par, as.numeric(1:4), tolerance = 1e-3)
  expect_true(res$convergence)

})

test_that("mleOptimWrapper respects lower and upper box constraints", {

  # Unconstrained minimum is at (100, 100); it must be clamped to the upper
  # bound (1, 1). This exercises the upper-bound penalty branch (regression:
  # the penalty used to be computed against 'lower' instead of 'upper').
  res <- mleOptimWrapper(minusLogLik = function(x) sum((x - c(100, 100))^2),
                         start = rbind(c(5, 5)), lower = c(-10, -10),
                         upper = c(1, 1), optMethod = "Nelder-Mead",
                         selectSolution = "lowest")
  expect_equal(res$par, c(1, 1), tolerance = 0.05)

  # Unconstrained minimum at (-100, -100) must be clamped to the lower bound
  res2 <- mleOptimWrapper(minusLogLik = function(x) sum((x + 100)^2),
                          start = rbind(c(-5, -5)), lower = c(-1, -1),
                          upper = c(10, 10), optMethod = "Nelder-Mead",
                          selectSolution = "lowest")
  expect_equal(res2$par, c(-1, -1), tolerance = 0.05)

})

test_that("mleOptimWrapper honors selectSolution", {

  # A spurious start (1:2) has no local minimum for a concave function.
  starts <- rbind(10:13, 1:2)
  lowest <- mleOptimWrapper(minusLogLik = function(x) -sum((x - 1:4)^2),
                            start = starts, selectSolution = "lowest")
  locMin <- mleOptimWrapper(minusLogLik = function(x) -sum((x - 1:4)^2),
                            start = starts, selectSolution = "lowestLocMin")
  # 'lowest' always returns something; 'lowestLocMin' returns NA for a concave
  # objective (no local minimum exists).
  expect_true(all(!is.na(lowest$par)))
  expect_true(all(is.na(locMin$par)))

})

test_that("mleOptimWrapper applies a user region as a penalty", {

  # Constrain to y <= x^2 with the region argument.
  region <- function(pars) {
    x <- pars[1]
    y <- pars[2]
    if (y <= x^2) {
      list("pars" = pars, "penalty" = 0)
    } else {
      list("pars" = c(sqrt(y), y), "penalty" = y - x^2)
    }
  }
  res <- mleOptimWrapper(minusLogLik = function(x) sum((x - 1:2)^2),
                         start = rbind(10:11), region = region,
                         lower = c(0.5, 1), upper = c(Inf, Inf),
                         optMethod = "Nelder-Mead", selectSolution = "lowest")
  expect_true(all(is.finite(res$par)))
  expect_true(res$value >= 0)

})
