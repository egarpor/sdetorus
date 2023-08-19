
library(animation)
library(png)
library(sdetorus)
library(spatstat)
library(fields)
library(manipulate)
library(viridis)

# Read image
img <- readPNG("sdetorus-logo.png")

# Smooth image and put in appropriate layout
A <- as.matrix(Smooth(1 - as.im(img[, , 1]), sigma = 2))

# Avoid zeros by adding a sparse Wrapped Normal
Mx <- nrow(A)
My <- ncol(A)
x <- seq(-pi, pi, l = My + 1)[-c(My + 1)]
A <- A + 0.5 * matrix(dWn1D(x, mu = -0.25, sigma = 2), nrow = Mx, ncol = My)

# Normalize matrix as a density
A <- A / periodicTrapRule2D(fxy = A)

# Rotate
A <- t(A)[, Mx:1]

# Drifts
p1x <- c(2:Mx, 1)
m1x <- c(Mx, 1:(Mx - 1))
p1y <- c(2:My, 1)
m1y <- c(My, 1:(My - 1))
bx <- (log(A)[p1x, ] - log(A)[m1x, ]) / (2 * 2 * pi / Mx)
by <- (log(A)[, p1y] - log(A)[, m1y]) / (2 * 2 * pi / My)

# Diffusions
sigma2 <- matrix(1, Mx, My)
sigmaxy <- matrix(0, Mx, My)

# Crank-Nicolson
N <- 1e3
deltat <- 10 / N
u0 <- 0.5 * matrix(repRow(dWn1D(x, mu = pi, sigma = 0.5), Mx) +
                     repCol(dWn1D(x, mu = pi, sigma = 0.5), My), nrow = Mx * My)
u <- crankNicolson2D(u0 = u0, bx = bx, by = by, sigma2x = sigma2,
                     sigma2y = sigma2, sigmaxy = sigmaxy, N = 0:N,
                     deltat = deltat, Mx = Mx, deltax = 2 * pi / Mx, My = My,
                     deltay = 2 * pi / My)
m <- max(u)

# Visualization density
manipulate({
  graphics::image(x, x, matrix(u[, j], nrow = Mx, ncol = My),
                  breaks = seq(0, m, l = 31), col = matlab.like2(30))
}, j = slider(1, N))

# Gif
saveGIF({
  par(mar = rep(0, 4))
  for (i in seq(1, N, by = 5)) {
    graphics::image(x, x, matrix(u[, i], nrow = Mx, ncol = My),
                    breaks = seq(0, m, l = 31), col = matlab.like2(30),
                    axes = FALSE, xlab = "", ylab = "")
  }
}, movie.name = "sdetorus.gif", interval = 0.1)

# Video
saveVideo({
  par(mar = rep(0, 4))
  for (i in seq(1, N, by = 5)) {
    graphics::image(x, x, matrix(u[, i], nrow = Mx, ncol = My),
                    breaks = seq(0, m, l = 31), col = matlab.like2(30),
                    axes = FALSE, xlab = "", ylab = "")
  }
}, video.name = "sdetorus.mp4", other.opts = "-pix_fmt yuv420p -b 800k",
interval = 0.05)
