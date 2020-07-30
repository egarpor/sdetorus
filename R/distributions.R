

#' @title Density of the von Mises
#'
#' @description Computes the density of a von Mises in a numerically stable way.
#'
#' @param x evaluation angles, not necessary in \eqn{[\pi,\pi)}.
#' @param mu circular mean.
#' @param kappa non-negative concentration parameter.
#' @return A vector of the same length as \code{x} containing the density.
#' @references
#' Jammalamadaka, S. R. and SenGupta, A. (2001) \emph{Topics in Circular Statistics}. World Scientific Publishing, River Edge.
#' @examples
#' x <- seq(-pi, pi, l = 200)
#' plot(x, x, type = "n", ylab = "Density", ylim = c(0, 1))
#' for (i in 0:20) {
#'   lines(x, dVm(x = x, mu = 0, kappa = 5 * i / 20),
#'         col = rainbow(21)[i + 1])
#' }
#' @export
dVm <- function(x, mu, kappa) {

  exp(kappa * (cos(x - mu) - 1) - logBesselI0Scaled(x = kappa) - log(2 * pi))

}


#' @title Jones and Pewsey (2005)'s circular distribution
#'
#' @description Computes the circular density of Jones and Pewsey (2005).
#'
#' @inheritParams dVm
#' @param psi shape parameter, see details.
#' @param const normalizing constant, computed with \code{constJp} if not provided.
#' @param M grid size for computing the normalizing constant by numerical integration.
#' @return A vector of the same length as \code{x} containing the density.
#' @details Particular interesting choices for the shape parameter are:
#' \itemize{
#' \item \code{psi = -1}: gives the Wrapped Cauchy as stationary density.
#' \item \code{psi = 0}: is the sinusoidal drift of the vM diffusion.
#' \item \code{psi = 1}: gives the Cardioid as stationary density.
#' }
#' @references
#' Jammalamadaka, S. R. and SenGupta, A. (2001) \emph{Topics in Circular Statistics}. World Scientific Publishing, River Edge.
#'
#' Jones, M. C. and Pewsey, A. (2005). A family of symmetric distributions on the circle. \emph{J. Amer. Statist. Assoc.}, 100(472):1422-1428.
#' @examples
#' x <- seq(-pi, pi, l = 200)
#' plot(x, x, type = "n", ylab = "Density", ylim = c(0, 0.6))
#' for (i in 0:20) {
#'   lines(x, dJp(x = x, mu = 0, kappa = 1, psi = -2 + 4 * i / 20),
#'         col = rainbow(21)[i + 1])
#' }
#' @export
dJp <- function(x, mu, kappa, psi, const) {

  # Particular case von Mises?
  if (psi == 0) {

    dens <- dVm(x = x, mu = mu, kappa = kappa)

  } else {

    # Cosh and sinh
    e <- exp(kappa * psi)
    em <- 1 / e
    ch <- 0.5 * (e + em)
    sh <- 0.5 * (e - em)

    # Constant
    if (missing(const)) {

      const <- constJp(mu = mu, kappa = kappa, psi = psi)

    }

    dens <- (ch + sh * cos(x - mu))^(1 / psi) / const

  }

  return(dens)

}


#' @rdname dJp
#' @export
constJp <- function(mu, kappa, psi, M = 200) {

  # Cosh and sinh
  e <- exp(kappa * psi)
  em <- 1 / e
  ch <- 0.5 * (e + em)
  sh <- 0.5 * (e - em)

  # Density
  dens <- (ch + sh * cos(seq(-pi, pi, l = M + 1)[-(M + 1)]))^(1 / psi)
  periodicTrapRule1D(fx = dens, lengthInterval = 2 * pi)

}
