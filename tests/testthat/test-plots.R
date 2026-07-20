# Plotting helpers are called for their side effects; open a null device so the
# tests run headless. We assert they run without error and that the value-
# returning ones give the documented output.

test_that("linesCirc and linesTorus draw wrapped segments and arrows", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  x <- 1:50
  y <- toPiInt(pi * cos(2 * pi * x / 50))
  plot(x, y, ylim = c(-pi, pi))
  expect_no_error(linesCirc(x = x, y = y, col = rainbow(length(x)),
                            ltyCross = 2))
  # arrows() warns on zero-length segments; that is expected, not a failure
  expect_no_error(suppressWarnings(linesCirc(x = x, y = y, arrows = TRUE)))
  expect_error(linesCirc(x = x, y = y[-1]), "lengths differ")

  xx <- toPiInt(rnorm(50))
  yy <- toPiInt(rnorm(50))
  plot(xx, yy, xlim = c(-pi, pi), ylim = c(-pi, pi))
  expect_no_error(linesTorus(x = xx, y = yy, col = rainbow(50), ltyCross = 2))
  expect_no_error(suppressWarnings(linesTorus(x = xx, y = yy, arrows = TRUE)))

})

test_that("plotSurface2D returns the evaluated grid and torusAxis draws", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  grid <- seq(-pi, pi, l = 40)
  z <- plotSurface2D(grid, grid, f = function(x) sin(x[1]) * cos(x[2]),
                     nLev = 10)
  expect_equal(dim(z), c(length(grid), length(grid)))

  # Vectorized f and explicit levels branch
  z2 <- plotSurface2D(grid, grid, f = function(x) sin(x[, 1]) * cos(x[, 2]),
                      levels = seq(-1, 1, l = 8), fVect = TRUE)
  expect_equal(dim(z2), c(length(grid), length(grid)))

  plot(grid, grid, type = "n", axes = FALSE)
  expect_no_error(torusAxis())
  expect_no_error(torusAxis(twoPi = TRUE))

})

test_that("matlab.like.colorRamps returns the requested number of colors", {

  expect_length(matlab.like.colorRamps(10), 10)
  expect_length(matlab.like.colorRamps(10, two = TRUE), 10)

})

test_that("rgl-based helpers run under the null device", {

  skip_if_not_installed("rgl")
  options(rgl.useNULL = TRUE)

  n <- 15
  x <- toPiInt(rnorm(n))
  y <- toPiInt(rnorm(n))
  z <- toPiInt(rnorm(n))
  rgl::plot3d(x, y, z)
  on.exit(rgl::close3d())

  expect_no_error(linesTorus3d(x = x, y = y, z = z, col = rainbow(n)))
  expect_no_error(suppressWarnings(
    linesTorus3d(x = x, y = y, z = z, arrows = TRUE, col = 1)))
  expect_no_error(torusAxis3d())

  f <- function(x) 10 * (sin(x[, 1]) * cos(x[, 2]) - sin(x[, 3]))^2
  g <- seq(-pi, pi, l = 12)
  t <- plotSurface3D(g, g, g, size = 5, alpha = 0.1, fVect = TRUE, f = f)
  expect_length(t, length(g)^3)

})
