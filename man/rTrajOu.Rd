% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OU.R
\name{rTrajOu}
\alias{rTrajOu}
\title{Simulation of trajectories for the univariate OU diffusion}
\usage{
rTrajOu(x0, alpha, mu, sigma, N = 100, delta = 0.001)
}
\arguments{
\item{x0}{initial point.}

\item{alpha}{strength of the drift.}

\item{mu}{unconditional mean of the diffusion.}

\item{sigma}{diffusion coefficient.}

\item{N}{number of discretization steps in the resulting trajectory.}

\item{delta}{time discretization step.}
}
\value{
A vector of length \code{N + 1} containing \code{x0} in the first
entry and the exact discretized trajectory on the remaining elements.
}
\description{
Simulation of trajectories of the \emph{univariate}
Ornstein--Uhlenbeck (OU) diffusion
\deqn{dX_t=\alpha(\mu - X_t)dt+\sigma dW_t, X_0=x_0} using the exact
transition probability density.
}
\details{
The law of the discretized trajectory is a multivariate normal
with mean \code{\link{meantOu}} and covariance matrix \code{\link{covstOu}}.
See \code{\link{rTrajMou}} for the multivariate case (less efficient for
dimension one).
}
\examples{
isRStudio <- identical(.Platform$GUI, "RStudio")
if (isRStudio) {
  manipulate::manipulate({
   set.seed(345678);
   plot(seq(0, N * delta, by = delta), rTrajOu(x0 = 0, alpha = alpha, mu = 0,
        sigma = sigma, N = N, delta = delta), ylim = c(-4, 4), type = "l",
        ylab = expression(X[t]), xlab = "t")
   }, delta = manipulate::slider(0.01, 5.01, step = 0.1),
   N = manipulate::slider(10, 500, step = 10, initial = 200),
   alpha = manipulate::slider(0.01, 5, step = 0.1, initial = 1),
   sigma = manipulate::slider(0.01, 5, step = 0.1, initial = 1))
}
}
