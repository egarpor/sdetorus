% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pde.R
\name{dTpdPde2D}
\alias{dTpdPde2D}
\title{Transition probability density in 2D by PDE solving}
\usage{
dTpdPde2D(Mx = 50, My = 50, x0, t, alpha, mu, sigma, rho = 0,
  type = "WN", Mt = ceiling(100 * t), sdInitial = 0.1, ...)
}
\arguments{
\item{Mx, My}{sizes of the equispaced spatial grids in \eqn{[-\pi,\pi)} for each component.}

\item{x0}{point giving the mean of the initial circular density, an
isotropic WN with standard deviations equal to \code{sdInitial}.}

\item{t}{time separating \code{x0} and the evaluation of the tpd.}

\item{alpha}{for \code{"WN"}, a vector of length \code{3} parametrizing
the \code{A} matrix as in \code{\link{alphaToA}}. For \code{"vM"}, a vector
of length \code{3} containing \code{c(alpha[1:2], A[1, 2])}, from the
arguments \code{alpha} and \code{A} in \code{\link{driftMvm}}.}

\item{mu}{vector of length \code{2} giving the mean.}

\item{sigma}{for \code{"WN"}, a vector of length \code{2} containing the
\strong{square root} of the diagonal of the diffusion matrix. For
\code{"vM"}, the standard deviation giving the isotropic diffusion matrix.}

\item{rho}{for \code{"WN"}, the correlation of the diffusion matrix.}

\item{type}{either \code{"WN"} or \code{"vM"}.}

\item{Mt}{size of the time grid in \eqn{[0, t]}.}

\item{sdInitial}{standard deviations of the concentrated WN giving the
initial condition.}

\item{...}{Further parameters passed to \code{\link{crankNicolson2D}}.}
}
\value{
A matrix of size \code{c(Mx, My)} with the tpd evaluated at the
combinations of \code{seq(-pi, pi, l = Mx + 1)[-(Mx + 1)]} and
\code{seq(-pi, pi, l = My + 1)[-(My + 1)]}.
}
\description{
Computation of the transition probability density (tpd) of the
Wrapped Normal (WN) or Multivariate von Mises (MvM) diffusion, by solving
its associated Fokker--Planck Partial Differential Equation (PDE) in 2D.
}
\details{
A combination of small \code{sdInitial} and coarse space-time
discretization (small \code{Mx} and \code{Mt}) is prone to create numerical
instabilities. See Sections 3.4.2, 2.2.1 and 2.2.2 in García-Portugués et al.
(2019) for details.
}
\examples{
M <- 100
x <- seq(-pi, pi, l = M + 1)[-c(M + 1)]
image(x, x, dTpdPde2D(Mx = M, My = M, x0 = c(0, pi), t = 1,
                      alpha = c(1, 1, 0.5), mu = c(pi / 2, 0), sigma = 1:2),
      zlim = c(0, 0.25), col = matlab.like.colorRamps(20),
      xlab = "x", ylab = "y")
}
\references{
García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019)
Langevin diffusions on the torus: estimation and applications.
\emph{Statistics and Computing}, 29(2):1--22. \doi{10.1007/s11222-017-9790-2}
}
