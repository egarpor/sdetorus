test_that("toPiInt, to2PiInt and toInt wrap to the correct intervals", {

  x <- seq(-3 * pi, 5 * pi, l = 200)
  yPi <- toPiInt(x)
  y2Pi <- to2PiInt(x)

  expect_true(all(yPi >= -pi & yPi < pi))
  expect_true(all(y2Pi >= 0 & y2Pi < 2 * pi))

  # Wrapping leaves the angle unchanged modulo 2 * pi
  expect_equal(cos(x), cos(yPi))
  expect_equal(sin(x), sin(y2Pi))

  # toInt to a generic [a, b)
  z <- toInt(1:10, a = 1, b = 5)
  expect_true(all(z >= 1 & z < 5))

})

test_that("diffCirc equals diff up to wrapping and to2PiInt of principal values", {

  x <- c(-pi, -pi / 2, pi - 0.1, -pi + 0.2)
  # Non-circular version is exactly diff()
  expect_equal(diffCirc(x, circular = FALSE), diff(x))
  # Circular differences lie in [-pi, pi)
  d <- diffCirc(x)
  expect_true(all(d >= -pi & d < pi))
  expect_equal(sin(d), sin(diff(x)))

})

test_that("unwrapCircSeries is a continuous representative of the wrapped series", {

  x <- c(-pi, -pi / 2, pi - 0.1, -pi + 0.2)
  u <- unwrapCircSeries(x)
  expect_equal(toPiInt(u), x)
  # Matrix input is unwrapped column-wise
  m <- cbind(x, rev(x))
  um <- unwrapCircSeries(m)
  expect_equal(toPiInt(um), m)

})

test_that("weightsLinearInterp1D reproduces linear functions exactly", {

  x <- seq(-4, 4, by = 0.5)
  g1 <- x - 0.25
  g2 <- x + 0.5
  w <- weightsLinearInterp1D(x = x, g1 = g1, g2 = g2)
  f <- function(x) 2 * x + 1
  expect_equal(rowSums(w * cbind(f(g1), f(g2))), f(x))
  # Weights add up to one
  expect_equal(rowSums(w), rep(1, length(x)))
  # Degenerate bins get 0.5-0.5
  wDeg <- weightsLinearInterp1D(x = 0, g1 = 1, g2 = 1)
  expect_equal(as.numeric(wDeg), c(0.5, 0.5))

})

test_that("weightsLinearInterp2D reproduces bilinear functions exactly", {

  x <- seq(-4, 4, by = 0.5)
  y <- 2 * x
  gx1 <- x - 0.25
  gx2 <- x + 0.5
  gy1 <- y - 0.75
  gy2 <- y + 0.25
  w <- weightsLinearInterp2D(x = x, y = y, gx1 = gx1, gx2 = gx2,
                             gy1 = gy1, gy2 = gy2)
  f <- function(x, y) 2 * x + 3 * y + 1
  # Column order is lower-lower, upper-lower, lower-upper, upper-upper
  expect_equal(rowSums(w * cbind(f(gx1, gy1), f(gx2, gy1),
                                 f(gx1, gy2), f(gx2, gy2))), f(x, y))
  expect_equal(rowSums(w), rep(1, length(x)))

})
