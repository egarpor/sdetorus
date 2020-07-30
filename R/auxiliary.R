

#' @title Lines and arrows with vertical wrapping
#'
#' @description Joins the corresponding points with line segments or arrows that exhibit wrapping in \eqn{[-\pi,\pi)} in the vertical axis.
#'
#' @param x vector with horizontal coordinates.
#' @param y vector with vertical coordinates, wrapped in \eqn{[-\pi,\pi)}.
#' @param col color vector of length \code{1} or the same length of \code{x} and \code{y}.
#' @param lty line type as in \code{\link[graphics]{par}}.
#' @param ltyCross specific line type for crossing segments.
#' @param arrows flag for drawing arrows instead of line segments.
#' @param ... further graphical parameters passed to \code{\link[graphics]{segments}} or \code{\link[graphics]{arrows}}.
#' @return Nothing. The functions are called for drawing wrapped lines.
#' @details \code{y} is wrapped to \eqn{[-\pi,\pi)} before plotting.
#' @examples
#' x <- 1:100
#' y <- toPiInt(pi * cos(2 * pi * x / 100) + 0.5 * runif(50, -pi, pi))
#' plot(x, y, ylim = c(-pi, pi))
#' linesCirc(x = x, y = y, col = rainbow(length(x)), ltyCross = 2)
#' plot(x, y, ylim = c(-pi, pi))
#' linesCirc(x = x, y = y, col = rainbow(length(x)), arrows = TRUE)
#' @export
linesCirc <- function(x = seq_along(y), y, col = 1, lty = 1, ltyCross = lty,
                      arrows = FALSE, ...) {

  # For determining crossings
  twoPi <- 2 * pi
  y <- toPiInt(y)
  dy <- diff(y)
  k <- -(dy < -pi) + (dy > pi)

  # Draw segments and avoid plotting two times the non-crossing ones
  l <- length(x)
  if (length(y) != l) stop("'x' and 'y' lengths differ")
  cross <- abs(k) + 1
  func <- ifelse(arrows, graphics::arrows, segments)
  func(x0 = x[-l], y0 = y[-l], x1 = x[-1], y1 = y[-1] - k * twoPi,
       lty = c(lty, ltyCross)[cross], col = col, ...)
  func(x0 = x[-l], y0 = y[-l] + k * twoPi, x1 = x[-1], y1 = y[-1],
       lty = c(0, ltyCross)[cross], col = col, ...)

}


#' @title Lines and arrows with wrapping in the torus
#'
#' @description Joins the corresponding points with line segments or arrows that exhibit wrapping in \eqn{[-\pi,\pi)} in the horizontal and vertical axes.
#'
#' @param x vector with horizontal coordinates, wrapped in \eqn{[-\pi,\pi)}.
#' @inheritParams linesCirc
#' @return Nothing. The functions are called for drawing wrapped lines.
#' @details \code{x} and \code{y} are wrapped to \eqn{[-\pi,\pi)} before plotting.
#' @examples
#' x <- toPiInt(rnorm(50, mean = seq(-pi, pi, l = 50), sd = 0.5))
#' y <- toPiInt(x + rnorm(50, mean = seq(-pi, pi, l = 50), sd = 0.5))
#' plot(x, y, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(length(x)),
#'      pch = 19)
#' linesTorus(x = x, y = y, col = rainbow(length(x)), ltyCross = 2)
#' plot(x, y, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(length(x)),
#'      pch = 19)
#' linesTorus(x = x, y = y, col = rainbow(length(x)), arrows = TRUE)
#' @export
linesTorus <- function(x, y, col = 1, lty = 1, ltyCross = lty, arrows = FALSE,
                       ...) {

  # For determining crossings
  twoPi <- 2 * pi
  x <- toPiInt(x)
  y <- toPiInt(y)
  dx <- diff(x)
  dy <- diff(y)
  kx <- -(dx < -pi) + (dx > pi)
  ky <- -(dy < -pi) + (dy > pi)

  # Draw segments and avoid plotting two times the non-crossing ones
  l <- length(x)
  if (length(y) != l) stop("'x' and 'y' lengths differ")
  cross <- (abs(kx) | abs(ky)) + 1
  func <- ifelse(arrows, graphics::arrows, segments)
  func(x0 = x[-l], y0 = y[-l], x1 = x[-1] - kx * twoPi, y1 = y[-1] - ky * twoPi,
       lty = c(lty, ltyCross)[cross], col = col, ...)
  func(x0 = x[-l] + kx * twoPi, y0 = y[-l] + ky * twoPi, x1 = x[-1], y1 = y[-1],
       lty = c(0, ltyCross)[cross], col = col, ...)

}


#' @title Lines and arrows with wrapping in the torus
#'
#' @description Joins the corresponding points with line segments or arrows that exhibit wrapping in \eqn{[-\pi,\pi)} in the horizontal and vertical axes.
#'
#' @param x,y vectors with horizontal coordinates, wrapped in \eqn{[-\pi,\pi)}.
#' @param z vector with vertical coordinates, wrapped in \eqn{[-\pi,\pi)}.
#' @param col color vector of length \code{1} or the same length of \code{x}, \code{y}, and \code{z}.
#' @inheritParams linesTorus
#' @return Nothing. The functions are called for drawing wrapped lines.
#' @details \code{x}, \code{y}, and \code{z} are wrapped to \eqn{[-\pi,\pi)} before plotting. \code{arrows = TRUE} makes sequential calls to \code{\link[rgl]{arrow3d}}, and is substantially slower than \code{arrows = FALSE}.
#' @examples
#' \dontrun{
#' library(rgl)
#' x <- toPiInt(rnorm(20, mean = seq(-pi, pi, l = 20), sd = 0.5))
#' y <- toPiInt(rnorm(20, mean = seq(-pi, pi, l = 20), sd = 0.5))
#' z <- toPiInt(x + y + rnorm(20, mean = seq(-pi, pi, l = 20), sd = 0.5))
#' plot3d(x, y, z, xlim = c(-pi, pi), ylim = c(-pi, pi), zlim = c(-pi, pi), 
#'        col = rainbow(length(x)), size = 2, box = FALSE, axes = FALSE)
#' linesTorus3d(x = x, y = y, z = z, col = rainbow(length(x)), lwd = 2)
#' torusAxis3d()
#' plot3d(x, y, z, xlim = c(-pi, pi), ylim = c(-pi, pi), zlim = c(-pi, pi), 
#'        col = rainbow(length(x)), size = 2, box = FALSE, axes = FALSE)
#' linesTorus3d(x = x, y = y, z = z, col = rainbow(length(x)), ltyCross = 2,
#'              arrows = TRUE, theta = 0.1 * pi / 180, barblen = 0.1)
#' torusAxis3d()
#' }
#' @export
linesTorus3d <- function(x, y, z, col = 1, arrows = FALSE, ...) {
  
  # For determining crossings
  twoPi <- 2 * pi
  x <- toPiInt(x)
  y <- toPiInt(y)
  z <- toPiInt(z)
  dx <- diff(x)
  dy <- diff(y)
  dz <- diff(z)
  kx <- -(dx < -pi) + (dx > pi)
  ky <- -(dy < -pi) + (dy > pi)
  kz <- -(dz < -pi) + (dz > pi)
  
  # Draw segments
  l <- length(x)
  if (length(y) != l | length(z) != l) stop("'x', 'y' or 'z' lengths differ")
  if (length(col) == 1) {
    
    col <- rep(col, l - 1)
    
  }
  if (arrows) {
    
    xyz1 <- cbind(x, y, z)
    k <- cbind(kx, ky, kz) * twoPi
    xyz2 <- xyz1[-1, , drop = FALSE]
    xyz1 <- xyz1[-l, , drop = FALSE]
    xyzk2 <- xyz2 - k
    xyzk1 <- xyz1 + k
    sapply(1:(l - 1), function(i) {
      rgl::arrow3d(p0 = xyz1[i, ], p1 = xyzk2[i, ], type = "lines", 
                   col = col[i], ...)
      rgl::arrow3d(p0 = xyzk1[i, ], p1 = xyz2[i, ], type = "lines", 
                   col = col[i], ...)
    })
    invisible()
    
  } else {
    
    rgl::segments3d(x = c(rbind(x[-l], x[-1] - kx * twoPi)),
                    y = c(rbind(y[-l], y[-1] - ky * twoPi)),
                    z = c(rbind(z[-l], z[-1] - kz * twoPi)),
                    col = col, ...)
    rgl::segments3d(x = c(rbind(x[-l] + kx * twoPi, x[-1])),
                    y = c(rbind(y[-l] + ky * twoPi, y[-1])),
                    z = c(rbind(z[-l] + kz * twoPi, z[-1])),
                    col = col, ...)
    
  }
  
}


#' @title Quadrature rules in 1D, 2D and 3D
#'
#' @description Quadrature rules for definite integrals over intervals in 1D, \eqn{\int_{x_1}^{x_2} f(x)dx}, rectangles in 2D, \eqn{\int_{x_1}^{x_2}\int_{y_1}^{y_2} f(x,y)dydx} and cubes in 3D, \eqn{\int_{x_1}^{x_2}\int_{y_1}^{y_2}\int_{z_1}^{z_2} f(x,y,z)dzdydx}. The trapezoidal rules assume that the function is periodic, whereas the Simpson rules work for arbitrary functions.
#'
#' @param fx vector containing the evaluation of the function to integrate over a uniform grid in \eqn{[x_1,x_2]}.
#' @param fxy matrix containing the evaluation of the function to integrate over a uniform grid in \eqn{[x_1,x_2]\times[y_1,y_2]}.
#' @param fxyz three dimensional array containing the evaluation of the function to integrate over a uniform grid in \eqn{[x_1,x_2]\times[y_1,y_2]\times[z_1,z_2]}.
#' @param endsMatch flag to indicate whether the values of the last entries of \code{fx}, \code{fxy} or \code{fxyz} are the ones in the first entries (elements, rows, columns, slices). See examples for usage.
#' @inheritParams base::sum
#' @param lengthInterval vector containing the lengths of the intervals of integration.
#' @return The value of the integral.
#' @details The simple trapezoidal rule has a very good performance for periodic functions in 1D and 2D(order of error ). The higher dimensional extensions are obtained by iterative usage of the 1D rules.
#' @references
#' Press, W. H., Teukolsky, S. A., Vetterling, W. T., Flannery, B. P. (1996). \emph{Numerical Recipes in Fortran 77: The Art of Scientific Computing (Vol. 1 of Fortran Numerical Recipes)}. Cambridge University Press, Cambridge.
#' @examples
#' # In 1D. True value: 3.55099937
#' N <- 21
#' grid <- seq(-pi, pi, l = N)
#' fx <- sin(grid)^2 * exp(cos(grid))
#' periodicTrapRule1D(fx = fx, endsMatch = TRUE)
#' periodicTrapRule1D(fx = fx[-N], endsMatch = FALSE)
#' integrateSimp1D(fx = fx, lengthInterval = 2 * pi)
#' integrateSimp1D(fx = fx[-N]) # Worse, of course
#'
#' # In 2D. True value: 22.31159
#' fxy <- outer(grid, grid, function(x, y) (sin(x)^2 * exp(cos(x)) +
#'                                          sin(y)^2 * exp(cos(y))) / 2)
#' periodicTrapRule2D(fxy = fxy, endsMatch = TRUE)
#' periodicTrapRule2D(fxy = fxy[-N, -N], endsMatch = FALSE)
#' periodicTrapRule1D(apply(fxy[-N, -N], 1, periodicTrapRule1D))
#' integrateSimp2D(fxy = fxy)
#' integrateSimp1D(apply(fxy, 1, integrateSimp1D))
#'
#' # In 3D. True value: 140.1878
#' fxyz <- array(fxy, dim = c(N, N, N))
#' for (i in 1:N) fxyz[i, , ] <- fxy
#' periodicTrapRule3D(fxyz = fxyz, endsMatch = TRUE)
#' integrateSimp3D(fxyz = fxyz)
#' @export
periodicTrapRule1D <- function(fx, endsMatch = FALSE, na.rm = TRUE,
                               lengthInterval = 2 * pi) {

  # If the intial and final points are the included (they are redundant)
  if (endsMatch) {

    fx <- fx[-1]

  }

  # Extension of equation (4.1.11) in Numerical Recipes in F77
  int <- mean(fx, na.rm = na.rm) * lengthInterval

  return(int)

}


#' @rdname periodicTrapRule1D
#' @export
periodicTrapRule2D <- function(fxy, endsMatch = FALSE, na.rm = TRUE,
                               lengthInterval = rep(2 * pi, 2)) {

  # If the intial and final points are the included (they are redundant)
  if (endsMatch) {

    fxy <- fxy[-1, -1]

  }

  # Extension of equation (4.1.11) in Numerical Recipes in F77
  int <- mean(fxy, na.rm = na.rm) * prod(lengthInterval)

  return(int)

}


#' @rdname periodicTrapRule1D
#' @export
periodicTrapRule3D <- function(fxyz, endsMatch = FALSE, na.rm = TRUE,
                               lengthInterval = rep(2 * pi, 3)) {

  # If the intial and final points are the included (they are redundant)
  if (endsMatch) {

    fxyz <- fxyz[-1, -1, -1]

  }

  # Extension of equation (4.1.11) in Numerical Recipes in F77
  int <- mean(fxyz, na.rm = na.rm) * prod(lengthInterval)

  return(int)

}


#' @rdname periodicTrapRule1D
#' @export
integrateSimp1D <- function(fx, lengthInterval = 2 * pi, na.rm = TRUE) {

  # Length and spacing of the grid
  n <- sum(!is.na(fx))
  h <- lengthInterval / (n - 1)

  # Equation (4.1.14) in Numerical Recipes in F77.
  # Much better than equation (4.1.13).
  int <- 3 / 8 * (fx[1] + fx[n]) +
         7 / 6 * (fx[2] + fx[n - 1]) +
         23 / 24 * (fx[3] + fx[n - 2]) +
         sum(fx[4:(n - 3)], na.rm = na.rm)
  int <- h * int

  return(int)

}


#' @rdname periodicTrapRule1D
#' @export
integrateSimp2D <- function(fxy, lengthInterval = rep(2 * pi, 2),
                            na.rm = TRUE) {

  # Length and spacing of the grid
  n <- dim(fxy)
  h <- lengthInterval / (n - 1)

  # Iteration of equation (4.1.14) in Numerical Recipes in F77
  # Integral on X
  intX <- 3 / 8 * (fxy[1, ] + fxy[n[1], ]) +
          7 / 6 * (fxy[2, ] + fxy[n[1] - 1, ]) +
          23 / 24 * (fxy[3, ] + fxy[n[1] - 2, ]) +
          colSums(fxy[4:(n[1] - 3), ], na.rm = na.rm)

  # Integral on Y
  int <- 3 / 8 * (intX[1] + intX[n[2]]) +
         7 / 6 * (intX[2] + intX[n[2] - 1]) +
         23 / 24 * (intX[3] + intX[n[2] - 2]) +
         sum(intX[4:(n[2] - 3)], na.rm = na.rm)
  int <- int * prod(h)

  return(int)

}


#' @rdname periodicTrapRule1D
#' @export
integrateSimp3D <- function(fxyz, lengthInterval = rep(2 * pi, 3),
                            na.rm = TRUE) {

  # Length and spacing of the grid
  n <- dim(fxyz)
  h <- lengthInterval / (n - 1)

  # Iteration of equation (4.1.14) in Numerical Recipes in F77
  # Integral on X
  intX <- 3 / 8 * (fxyz[1, , ] + fxyz[n[1], , ]) +
          7 / 6 * (fxyz[2, , ] + fxyz[n[1] - 1, , ]) +
          23 / 24 * (fxyz[3, , ] + fxyz[n[1] - 2, , ]) +
          colSums(fxyz[4:(n[1] - 3), , ], na.rm = na.rm, dims = 1)

  # Integral on Y
  intY <- 3 / 8 * (intX[1, ] + intX[n[2], ]) +
          7 / 6 * (intX[2, ] + intX[n[2] - 1, ]) +
          23 / 24 * (intX[3, ] + intX[n[2] - 2, ]) +
          colSums(intX[4:(n[2] - 3), ], na.rm = na.rm)

  # Integral on Z
  int <- 3 / 8 * (intY[1] + intY[n[3]]) +
         7 / 6 * (intY[2] + intY[n[3] - 1]) +
         23 / 24 * (intY[3] + intY[n[3] - 2]) +
         sum(intY[4:(n[3] - 3)], na.rm = na.rm)
  int <-  int * prod(h)

  return(int)

}


#' @title Wrapping of radians to its principal values
#'
#' @description Utilities for transforming a reals into \eqn{[-\pi, \pi)}, \eqn{[0, 2\pi)} or \eqn{[a, b)}.
#'
#' @param x a vector, matrix or object for whom \code{\link[base]{Arithmetic}} is defined.
#' @param a,b the lower and upper limits of \eqn{[a, b)}.
#' @return The wrapped vector in the chosen interval.
#' @details Note that \eqn{b} is \bold{excluded} from the result, see examples.
#' @examples
#' # Wrapping of angles
#' x <- seq(-3 * pi, 5 * pi, l = 100)
#' toPiInt(x)
#' to2PiInt(x)
#'
#' # Transformation to [1, 5)
#' x <- 1:10
#' toInt(x, 1, 5)
#' toInt(x + 1, 1, 5)
#'
#' # Transformation to [1, 5]
#' toInt(x, 1, 6)
#' toInt(x + 1, 1, 6)
#' @export
toPiInt <- function(x) {

  (x + pi) %% (2 * pi) - pi

}


#' @rdname toPiInt
#' @export
to2PiInt <- function(x) {

  x %% (2 * pi)

}


#' @rdname toPiInt
#' @export
toInt <- function(x, a, b) {

  (x - a) %% (b - a) + a

}


#' @title Lagged differences for circular time series
#'
#' @description Returns suitably lagged and iterated circular differences.
#'
#' @param x wrapped or unwrapped angles to be differenced. Must be a vector or a matrix, see details.
#' @param circular convenience flag to indicate whether wrapping should be done. If \code{FALSE}, the function is exactly \code{\link{diff}}.
#' @param ... parameters to be passed to \code{\link{diff}}.
#' @return The value of \code{diff(x, ...)}, circularly wrapped. Default parameters give an object of the kind of \code{x} with one less entry or row.
#' @details If \code{x} is a matrix then the difference operations are carried out row-wise, on each column separately.
#' @examples
#' # Vectors
#' x <- c(-pi, -pi/2, pi - 0.1, -pi + 0.2)
#' diffCirc(x) - diff(x)
#'
#' # Matrices
#' set.seed(234567)
#' N <- 100
#' x <- t(euler2D(x0 = rbind(c(0, 0)), A = diag(c(1, 1)), sigma = rep(2, 2),
#'                mu = c(pi, pi), N = N, delta = 1, type = 2)[1, , ])
#' diffCirc(x) - diff(x)
#' @export
diffCirc <- function(x, circular = TRUE, ...) {

  # Ordinary differences
  dif <- diff(x, ...)

  # Signed minimum circular differences
  if (circular) {

    dif <- toPiInt(dif)

  }

  return(dif)

}


#' @title Unwrapping of circular time series
#'
#' @description Completes a circular time series to a linear one by computing the closest wind numbers. Useful for plotting circular trajectories with crossing of boundaries.
#'
#' @param x wrapped angles. Must be a vector or a matrix, see details.
#' @return A vector or matrix containing a choice of unwrapped angles of \code{x} that maximizes linear continuity.
#' @details If \code{x} is a matrix then the unwrapping is carried out row-wise, on each column separately.
#' @examples
#' # Vectors
#' x <- c(-pi, -pi/2, pi - 0.1, -pi + 0.2)
#' u <- unwrapCircSeries(x)
#' max(abs(toPiInt(u) - x))
#' plot(toPiInt(x), ylim = c(-pi, pi))
#' for(k in -1:1) lines(u + 2 * k * pi)
#'
#' # Matrices
#' N <- 100
#' set.seed(234567)
#' x <- t(euler2D(x0 = rbind(c(0, 0)), A = diag(c(1, 1)), sigma = rep(1, 2),
#'                  mu = c(pi, pi), N = N, delta = 1, type = 2)[1, , ])
#' u <- unwrapCircSeries(x)
#' max(abs(toPiInt(u) - x))
#' plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi))
#' for(k1 in -3:3) for(k2 in -3:3) lines(u[, 1] + 2 * k1 * pi,
#'                                       u[, 2] + 2 * k2 * pi, col = gray(0.5))
#' @export
unwrapCircSeries <- function(x) {

  if (is.matrix(x)) {

    # Recursive call
    res <- apply(x, 2, unwrapCircSeries)

  } else {

    # Add + (-) one winding number if the difference with the previous point is
    # smaller (larger) than -pi (pi)
    dif <- diff(x)
    res <- x + c(0, 2 * pi * (cumsum(dif < -pi) - cumsum(dif > pi)))

  }

  return(res)

}


#' @title Weights for linear interpolation in 1D and 2D
#'
#' @description Computation of weights for linear interpolation in 1D and 2D.
#'
#' @param x,y vectors of length \code{n} containing the points where the weights should be computed.
#' @param g1,g2,gx1,gx2,gy1,gy2 vectors of length \code{n} such that \code{g1 <= x <= g2} for 1D and \code{gx1 <= x <= gx2} and \code{gy1 <= y <= gy2} for 2D.
#' @param circular flag to indicate whether the differences should be circularly wrapped.
#' @return \itemize{
#' \item For 1D, a matrix of size \code{c(n, 2)} containing the weights for lower (\code{g1}) and upper (\code{g1}) bins.
#' \item For 2D, a matrix of size \code{c(n, 4)} containing the weights for lower-lower (\code{g1x}, \code{g1y}), \emph{upper-lower} (\code{g2x}, \code{g1y}), \emph{lower-upper} (\code{g1x}, \code{g2y}) and upper-upper (\code{g2x}, \code{g2y}) bins.
#' \code{cbind(g1x, g1y)}, \code{cbind(g1x, g1y)}, \code{cbind(g1x, g1y)} and \code{cbind(g2x, g2y)}.
#' }
#' @details See the examples for how to use the weights for linear binning.
#' @examples
#' # 1D, linear
#' x <- seq(-4, 4, by = 0.5)
#' g1 <- x - 0.25
#' g2 <- x + 0.5
#' w <- weightsLinearInterp1D(x = x, g1 = g1, g2 = g2)
#' f <- function(x) 2 * x + 1
#' rowSums(w * cbind(f(g1), f(g2)))
#' f(x)
#'
#' # 1D, circular
#' x <- seq(-pi, pi, by = 0.5)
#' g1 <- toPiInt(x - 0.25)
#' g2 <- toPiInt(x + 0.5)
#' w <- weightsLinearInterp1D(x = x, g1 = g1, g2 = g2, circular = TRUE)
#' f <- function(x) 2 * sin(x) + 1
#' rowSums(w * cbind(f(g1), f(g2)))
#' f(x)
#'
#' # 2D, linear
#' x <- seq(-4, 4, by = 0.5)
#' y <- 2 * x
#' gx1 <- x - 0.25
#' gx2 <- x + 0.5
#' gy1 <- y - 0.75
#' gy2 <- y + 0.25
#' w <- weightsLinearInterp2D(x = x, y = y, gx1 = gx1, gx2 = gx2,
#'                            gy1 = gy1, gy2 = gy2)
#' f <- function(x, y) 2 * x + 3 * y + 1
#' rowSums(w * cbind(f(gx1, gy1), f(gx2, gy1), f(gx1, gy2), f(gx2, gy2)))
#' f(x, y)
#'
#' # 2D, circular
#' x <- seq(-pi, pi, by = 0.5)
#' y <- toPiInt(2 * x)
#' gx1 <- toPiInt(x - 0.25)
#' gx2 <- toPiInt(x + 0.5)
#' gy1 <- toPiInt(y - 0.75)
#' gy2 <- toPiInt(y + 0.25)
#' w <- weightsLinearInterp2D(x = x, y = y, gx1 = gx1, gx2 = gx2,
#'                            gy1 = gy1, gy2 = gy2, circular = TRUE)
#' f <- function(x, y) 2 * sin(x) + 3 * cos(y) + 1
#' rowSums(w * cbind(f(gx1, gy1), f(gx2, gy1), f(gx1, gy2), f(gx2, gy2)))
#' f(x, y)
#' @export
weightsLinearInterp1D <- function(x, g1, g2, circular = FALSE) {

  if (circular) {

    # Quotient
    dif <- toPiInt(g2 - g1)

    # Weights
    w <- toPiInt(cbind(g2 - x, x - g1))
    w <- w / dif

  } else {

    # Quotient
    dif <- g2 - g1

    # Weights
    w <- cbind(g2 - x, x - g1)
    w <- w / dif

  }

  # Avoid NaNs
  indDifZero <- which(dif == 0)
  w[indDifZero, ] <- c(0.5, 0.5)

  return(w)

}


#' @rdname weightsLinearInterp1D
#' @export
weightsLinearInterp2D <- function(x, y, gx1, gx2, gy1, gy2, circular = FALSE) {

  if (circular) {

    # Quotient
    dif <- toPiInt(gx2 - gx1) * toPiInt(gy2 - gy1)

    # Weights
    gx1 <- toPiInt(x - gx1)
    gx2 <- toPiInt(gx2 - x)
    gy1 <- toPiInt(y - gy1)
    gy2 <- toPiInt(gy2 - y)
    w <- cbind(gx2 * gy2, gx1 * gy2, gx2 * gy1, gx1 * gy1) / dif

  } else {

    # Quotient
    dif <- (gx2 - gx1) * (gy2 - gy1)

    # Weights
    gx1 <- x - gx1
    gx2 <- gx2 - x
    gy1 <- y - gy1
    gy2 <- gy2 - y
    w <- cbind(gx2 * gy2, gx1 * gy2, gx2 * gy1, gx1 * gy1) / dif

  }

  # Avoid NaNs
  indDifZero <- which(dif == 0)
  w[indDifZero, ] <- c(0.5, 0.5, 0.5, 0.5)

  return(w)

}


#' @title Contour plot of a 2D surface
#'
#' @description Convenient wrapper for plotting a contour plot of a real function of two variables.
#'
#' @param x,y numerical grids fore each dimension. They must be in ascending order.
#' @param f function to be plot. Must take a single argument (see examples).
#' @param z a vector of length \code{length(x) * length(y)} containing the evaluation of \code{f} in the bivariate grid. If not provided, it is computed internally.
#' @param nLev the number of levels the range of \code{z} will be divided into.
#' @param levels vector of contour levels. If not provided, it is set to \code{quantile(z, probs = seq(0, 1, l = nLev))}.
#' @param fVect flag to indicate whether \code{f} is a vectorized function (see examples).
#' @param ... further arguments passed to \code{\link[graphics]{image}}
#' @return The matrix \code{z}, invisible.
#' @examples
#' \dontrun{
#' grid <- seq(-pi, pi, l = 100)
#' plotSurface2D(grid, grid, f = function(x) sin(x[1]) * cos(x[2]), nLev = 20)
#' plotSurface2D(grid, grid, f = function(x) sin(x[, 1]) * cos(x[, 2]),
#'               levels = seq(-1, 1, l = 10), fVect = TRUE)
#' }
#' @export
plotSurface2D <- function(x = 1:nrow(z), y = 1:ncol(z), f, z, nLev = 20, levels,
                          fVect = FALSE, ...) {

  # Grid xy
  xy <- as.matrix(expand.grid(x = x, y = y))

  # z = f(xy)
  if (missing(z)) {

    if (fVect) {

      z <- f(xy)

    } else {

      z <- apply(xy, 1, f)

    }

    z <- matrix(z, nrow = length(x), ncol = length(y))

  }

  # Levels as quantiles
  if (missing(levels)) {

    levels <- quantile(c(z), prob = seq(0, 1, length.out = nLev), na.rm = TRUE)

  } else {

    nLev <- length(levels)

  }

  # 'Nice' plot
  image(x, y, z, breaks = levels, col = colorRamps::matlab.like(nLev - 1), ...)
  contour(x, y, z, levels = levels, add = TRUE,
          labels = format(levels, digits = 1))

  return(invisible(z))

}


#' @title Visualization of a 3D surface
#'
#' @description Convenient wrapper for visualizing a real function of three variables by means of a colour scale and alpha shading.
#'
#' @param x,y,z numerical grids for each dimension.
#' @param t a vector of length \code{length(x) * length(y) * length(z)} containing the evaluation of \code{f} in the trivariate grid. If not provided, it is computed internally.
#' @inheritParams plotSurface2D
#' @param nLev number of levels in the colour scale.
#' @param levels vector of breaks in the colour scale. If not provided, it is set to \code{quantile(z, probs = seq(0, 1, l = nLev))}.
#' @param size size of points in pixels.
#' @param alpha alpha value between \code{0} (fully transparent) and \code{1} (opaque).
#' @param ... further arguments passed to \code{\link[rgl]{plot3d}}
#' @return The vector \code{t}, invisible.
#' @examples
#' \dontrun{
#' grid <- seq(-pi, pi, l = 50)
#' t <- plotSurface3D(grid, grid, grid, size = 10, alpha = 0.01, fVect = TRUE,
#'                    f = function(x) 10 * (sin(x[, 1]) * cos(x[, 2]) - sin(x[, 3]))^2)
#' plotSurface3D(grid, grid, grid, t = t, size = 15, alpha = 0.1,
#'               levels = quantile(t, probs = seq(0, 0.1, l = 20)))
#' }
#' @export
plotSurface3D <- function(x = 1:nrow(t), y = 1:ncol(t), z = 1:dim(t)[3], f, t,
                          nLev = 20, levels, fVect = FALSE, size = 15,
                          alpha = 0.05, ...) {

  # Grids xy and xyz
  xyz <- as.matrix(expand.grid(x, y, z))

  # t = f(xyz)
  if (missing(t)) {

    if (fVect) {

      t <- f(xyz)

    } else {

      t <- apply(xyz, 1, f)

    }

  }

  # Levels as quantiles
  if (missing(levels)) {

    levels <- quantile(t, prob = seq(0, 1, length.out = nLev), na.rm = TRUE)

  } else {

    nLev <- length(levels)

  }

  # Subset of points with color
  tCut <- cut(x = t, breaks = levels)
  tCutNoNA <- !is.na(tCut)
  xyz <- subset(xyz, subset = tCutNoNA)
  col <- colorRamps::matlab.like(nLev)[tCut[tCutNoNA]]

  # Nice plot
  rgl::plot3d(xyz, col = col, size = size, alpha = alpha, lit = FALSE,
              point_antialias = FALSE, smooth = FALSE, aspect = FALSE, ...)

  return(invisible(t))

}


#' @title Visualization of a 3D surface using VMD
#'
#' @description Visualization of a real function of three variables by means of the 3D contour feature of VMD. The function creates an XPLOR map file that can be read by VMD and, optionally, establishes a connection for inspection of the surface.
#'
#' @param t a vector of length \code{nx * ny * nz} containing the evaluation of the function to be represented, i.e. \code{t = f(expand.grid(x, y, z))}.
#' @param nx,ny,nz length of the grid along each dimension.
#' @param outFile file in XPLOR format to be read by VMD.
#' @param call whether to invoke VMD for immediate visualization or not.
#' @param path path to VMD in the system. See details.
#' @return The function creates a XPLOR map file in the current working directory. This file contains \code{t} and other data that allows its readily interpretation by VMD. An optional call to VMD can also be done for immediate visualization.
#' @details If \code{call = TRUE}, then a call to VMD with \code{outFile} is produced. VMD needs to be installed (\url{http://www.ks.uiuc.edu/Research/vmd/}) and its \code{path} must be properly specified. For example, for Mac OS it should be something similar to \code{path = "/Applications/VMD\\ 1.9.2.app/Contents/MacOs/startup.command"}, with the version number depending on the installed one. The R session will wait until finalization of VMD.
#'
#' After launching VMD, a 3D scatterplot is shown. The options for exploration of the 3D contours can be accessed in "VMD Main -> Graphics -> Representations ...". This will open a new window with the available choices. The drawing methods 'Isosurface', 'VolumeSlice' and 'Fieldlines' are the most useful.
#'
#' Recall that the spacing is assumed to be the same in each of the grids.
#'
#' See \url{http://www.ks.uiuc.edu/Research/vmd/} for help on using VMD.
#' @examples
#' \dontrun{
#' gridx <- seq(-pi, pi, l = 50)
#' gridy <- seq(0, pi, l = 50)
#' gridz <- seq(0, pi, l = 20)
#' t <- plotSurface3D(gridx, gridy, gridz, size = 10, alpha = 0.01, fVect = TRUE,
#'                    f = function(x) 10 * (sin(x[, 1]) * cos(x[, 2]) - sin(x[, 3]))^2)
#'
#' # VMD has to be installed in the system (http://www.ks.uiuc.edu/Research/vmd/)
#' # Also, the path must be properly specified
#' plotVmdSurface3D(t = t, nx = length(gridx), ny = length(gridy), nz = length(gridz),
#'                  outFile = "test.map", call = TRUE, path =
#'                  "/Applications/VMD\\ 1.9.2.app/Contents/MacOs/startup.command")
#' }
#' @references
#' Humphrey, W., Dalke, A. and Schulten, K. (1996). VMD - Visual Molecular Dynamics, \emph{J. Molec. Graphics}, 14(1):33-38.
#' @export
plotVmdSurface3D <- function(t, nx, ny, nz, outFile = "test.map", call = FALSE,
                             path) {

  # Variables of little use since there is no scale on the axis
  spacing <- 1
  center <- c(0, 0, 0)

  # Preface
  utils::write.table(x = rbind("GRID_PARAMETER_FILE NONE",
                               "GRID_DATA_FILE NONE",
                               "MACROMOLECULE NONE",
                               paste("SPACING", spacing),
                               paste("NELEMENTS", nx - 1, ny - 1, nz - 1),
                               paste("CENTER", paste(center, collapse = " "))),
                     file = outFile, row.names = FALSE, col.names = FALSE,
                     quote = FALSE)

  # Data
  utils::write.table(x = t, file = outFile, append = TRUE, row.names = FALSE,
                     col.names = FALSE)

  # Call
  if (call) {

    system(paste(path, paste(getwd(), outFile, sep = "/")))

  }

}


#' @title Replication of rows and columns
#'
#' @description Wrapper for replicating a matrix/vector by rows or columns.
#'
#' @param x a numerical vector or matrix of dimension \code{c(nr, nc)}.
#' @param n the number of replicates of \code{x} by rows or columns.
#' @return A matrix of dimension \code{c(nr * n, nc)} for \code{repRow} or \code{c(nr, nc * n)} for \code{repCol}.
#' @examples
#' repRow(1:5, 2)
#' repCol(1:5, 2)
#' A <- rbind(1:5, 5:1)
#' A
#' repRow(A, 2)
#' repCol(A, 2)
#' @export
repRow <- function(x, n) {

  # If matrix, replicate blocks
  if (is.matrix(x)) {

    nr <- nrow(x) * n
    nc <- ncol(x)

  } else {

    nr <- n
    nc <- length(x)

  }

  matrix(t(x), nrow = nr, ncol = nc, byrow = TRUE)

}


#' @rdname repRow
#' @export
repCol <- function(x, n) {

  # If matrix, replicate blocks
  if (is.matrix(x)) {

    nr <- nrow(x)
    nc <- ncol(x) * n

  } else {

    nr <- length(x)
    nc <- n

  }

  matrix(x, nrow = nr, ncol = nc)

}


#' @title Utilities for conversion between row-column indexing and linear indexing of matrices
#'
#' @description Conversions between \code{cbind(i, j)} and \code{k} such that \code{A[i, j] == A[k]} for a matrix \code{A}. Either column or row ordering can be specified for the linear indexing, and also direct conversions between both types.
#'
#' @param i row index.
#' @param j column index.
#' @param k linear indexes for column-stacking or row-stacking ordering (if \code{byRows = TRUE}).
#' @param nr number of rows.
#' @param nc number of columns.
#' @param byRows whether to use row-ordering instead of the default column-ordering.
#' @return Depending on the function:
#' \itemize{
#' \item \code{kIndex}: a vector of length \code{nr * nc} with the linear indexes for \code{A}.
#' \item \code{ijIndex}: a matrix of dimension \code{c(length(k), 2)} giving \code{cbind(i, j)}.
#' \item \code{kColToRow} and \code{kRowToCol}: a vector of length \code{nr * nc} giving the permuting indexes to change the ordering of the linear indexes.
#' }
#' @examples
#' # Indexes of a 3 x 5 matrix
#' ij <- expand.grid(i = 1:3, j = 1:5)
#' kCols <- kIndex(i = ij[, 1], j = ij[, 2], nr = 3, nc = 5)
#' kRows <- kIndex(i = ij[, 1], j = ij[, 2], nr = 3, nc = 5, byRows = TRUE)
#'
#' # Checks
#' ijIndex(kCols, nr = 3, nc = 5)
#' ij
#' ijIndex(kRows, nr = 3, nc = 5, byRows = TRUE)
#' ij
#'
#' # Change column to row (and viceversa) ordering in the linear indexes
#' matrix(1:10, nr = 2, nc = 5)
#' kColToRow(1:10, nr = 2, nc = 5)
#' kRowToCol(kColToRow(1:10, nr = 2, nc = 5), nr = 2, nc = 5)
#' @export
kIndex <- function(i, j, nr, nc, byRows = FALSE) {

  if (byRows) {

    k <- (j - 1) * nr + i

  } else {

    k <- (i - 1) * nc + j

  }

  return(k)

}


#' @rdname kIndex
#' @export
ijIndex <- function(k, nr, nc, byRows = FALSE) {

  if (byRows) {

    i <- (k - 1) %% nr + 1
    j <- (k - i) / nr + 1

  } else {

    j <- (k - 1) %% nc + 1
    i <- (k - j) / nc + 1

  }

  cbind(i, j)

}


#' @rdname kIndex
#' @export
kColToRow <- function(k, nr, nc) {

  f <- floor((k - 1) / nc)
  (k - 1 - f * nc) * nr + f + 1

}


#' @rdname kIndex
#' @export
kRowToCol <- function(k, nr, nc) {

  f <- floor((k - 1) / nr)
  (k - 1 - f * nr) * nc + f + 1

}


#' @title Matching of matrices
#'
#' @description Wrapper for matching a matrix against another, by rows or columns.
#'
#' @param x matrix with the values to be matched.
#' @param mat matrix with the values to be matched against.
#' @param rows whether the match should be done by rows (\code{TRUE}) or columns (\code{FALSE}).
#' @param useMatch whether to rely on \code{\link[base]{match}} or not. Might give unexpected mismatches due to working with lists.
#' @param ... further parameters passed to \code{\link[base]{match}}.
#' @return An integer vector of length \code{nrow(x)} (or \code{ncol(x)}) giving the row (or col) position in table of the first match, if there is a match.
#' @examples
#' # By rows
#' A <- rbind(5:6, repRow(1:2, 3), 3:4)
#' B <- unique(A)
#' ind <- matMatch(x = A, mat = B)
#' A
#' B[ind, ]
#'
#' # By columns
#' A <- cbind(5:6, repCol(1:2, 3), 3:4)
#' B <- t(unique(t(A)))
#' ind <- matMatch(x = A, mat = B, rows = FALSE)
#' A
#' B[, ind]
#' @export
matMatch <- function(x, mat, rows = TRUE, useMatch = FALSE, ...) {

  if (useMatch) {

    # Match rows
    if (rows) {

      x <- t(x)
      mat <- t(mat)

    }

    # Call match using list matching
    match(x = data.frame(x), table = data.frame(mat), ...)

  } else {

    # Match columns
    if (!rows) {

      x <- t(x)
      mat <- t(mat)

    }

    # Naive approach
    sapply(1:nrow(x), function(i) which(apply(mat, 1,
                                              function(y) all(y == x[i, ])))[1])

  }

}


#' @title Draws pretty axis labels for circular variables
#'
#' @description Wrapper for drawing pretty axis labels for circular variables. To be invoked after \code{plot} with \code{axes = FALSE} has been called.
#'
#' @param sides an integer vector specifying which side of the plot the axes are to be drawn on. The axes are placed as follows: \code{1} = below, \code{2} = left, \code{3} = above, and \code{4} = right.
#' @param twoPi flag indicating that \eqn{[0,2\pi)} is the support, instead of \eqn{[-\pi,\pi)}.
#' @param ... further parameters passed to \code{\link[graphics]{axis}}.
#' @return This function is usually invoked for its side effect, which is to add axes to an already existing plot.
#' @details The function calls \code{\link[graphics]{box}}.
#' @examples
#' grid <- seq(-pi, pi, l = 100)
#' plotSurface2D(grid, grid, f = function(x) sin(x[1]) * cos(x[2]),
#'               nLev = 20, axes = FALSE)
#' torusAxis()
#' @export
torusAxis <- function(sides = 1:2, twoPi = FALSE, ...) {

  # Choose (-pi, pi) or (0, 2*pi)
  if (twoPi) {

    at <- seq(0, 2 * pi, l = 5)
    labels <- expression(0, pi/2, pi, 3*pi/2, 2*pi)

  } else {

    at <- seq(-pi, pi, l = 5)
    labels <- expression(-pi, -pi/2, 0, pi/2, pi)

  }

  # Draw box + axis
  box()
  for (i in sides) {

    axis(i, at = at, labels = labels, ...)

  }

}


#' @title Draws pretty axis labels for circular variables
#'
#' @description Wrapper for drawing pretty axis labels for circular variables. To be invoked after \code{plot3d} with \code{axes = FALSE} and \code{box = FALSE} has been called.
#'
#' @param sides an integer vector specifying which side of the plot the axes are to be drawn on. The axes are placed as follows: \code{1} = x, \code{2} = y, \code{3} = z.
#' @inheritParams torusAxis
#' @param ... further parameters passed to \code{\link[rgl:axes3d]{axis3d}}.
#' @return This function is usually invoked for its side effect, which is to add axes to an already existing plot.
#' @details The function calls \code{\link[rgl:axes3d]{box3d}}.
#' @examples
#' library(rgl)
#' x <- toPiInt(rnorm(50, mean = seq(-pi, pi, l = 50), sd = 0.5))
#' y <- toPiInt(rnorm(50, mean = seq(-pi, pi, l = 50), sd = 0.5))
#' z <- toPiInt(x + y + rnorm(50, mean = seq(-pi, pi, l = 50), sd = 0.5))
#' plot3d(x, y, z, xlim = c(-pi, pi), ylim = c(-pi, pi), zlim = c(-pi, pi), 
#'        col = rainbow(length(x)), size = 2, box = FALSE, axes = FALSE)
#' torusAxis3d()
#' @export
torusAxis3d <- function(sides = 1:3, twoPi = FALSE, ...) {
  
  # Choose (-pi, pi) or (0, 2*pi)
  if (twoPi) {
    
    at <- seq(0, 2 * pi, l = 5)
    labels <- expression(0, pi/2, pi, 3*pi/2, 2*pi)
    
  } else {
    
    at <- seq(-pi, pi, l = 5)
    labels <- expression(-pi, -pi/2, 0, pi/2, pi)
    
  }
  
  # Draw box + axis
  rgl::box3d()
  for (i in sides) {
    
    suppressWarnings(rgl::axis3d(c("x", "y", "z")[i], at = at,
                                 labels = labels, ...))
    
  }
  
}
