% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample-trajectories.R
\name{rTrajWn2D}
\alias{rTrajWn2D}
\title{Simulation of trajectories for the WN diffusion in 2D}
\usage{
rTrajWn2D(x0, alpha, mu, sigma, rho = 0, N = 100, delta = 0.01,
  NFine = ceiling(N * delta/deltaFine), deltaFine = min(delta/100, 0.001))
}
\arguments{
\item{x0}{vector of length \code{2} giving the initial point.}

\item{alpha}{vector of length \code{3} parametrizing the \code{A} matrix as in \code{\link{alphaToA}}.}

\item{mu}{a vector of length \code{2} giving the mean.}

\item{sigma}{vector of length \code{2} containing the \strong{square root} of the diagonal of \eqn{\Sigma}, the diffusion matrix.}

\item{rho}{correlation coefficient of \eqn{\Sigma}.}

\item{N}{number of discretization steps in the resulting trajectory.}

\item{delta}{discretization step.}

\item{NFine}{number of discretization steps for the fine trajectory. Must
be larger than \code{N}.}

\item{deltaFine}{discretization step for the fine trajectory. Must be
smaller than \code{delta}.}
}
\value{
A matrix of size \code{c(N + 1, 2)} containing \code{x0} in the
first entry and the discretized trajectory.
}
\description{
Simulation of the Wrapped Normal (WN) diffusion in 2D by
subsampling a fine trajectory obtained by the Euler discretization.
}
\details{
The fine trajectory is subsampled using the indexes
\code{seq(1, NFine + 1, by = NFine / N)}.
}
\examples{
samp <- rTrajWn2D(x0 = c(0, 0), alpha = c(1, 1, -0.5), mu = c(pi, pi),
                    sigma = c(1, 1), N = 1000, delta = 0.01)
plot(samp, xlim = c(-pi, pi), ylim = c(-pi, pi), pch = 19, cex = 0.25,
     xlab = expression(X[t]), ylab = expression(Y[t]), col = rainbow(1000))
linesTorus(samp[, 1], samp[, 2], col = rainbow(1000))
}
