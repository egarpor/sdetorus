

#' @title sdetorus - Statistical Tools for Toroidal Diffusions
#'
#' @description Implementation of statistical methods for the estimation of
#' toroidal diffusions. Several diffusive models are provided, most of them
#' belonging to the Langevin family of diffusions on the torus. Specifically,
#' the wrapped normal and von Mises processes are included, which can be seen
#' as toroidal analogues of the Ornstein--Uhlenbeck diffusion. A collection of
#' methods for approximate maximum likelihood estimation, organized in four
#' blocks, is given: (i) based on the exact transition probability density,
#' obtained as the numerical solution to the Fokker-Plank equation;
#' (ii) based on wrapped pseudo-likelihoods; (iii) based on specific analytic
#' approximations by wrapped processes; (iv) based on maximum likelihood of
#' the stationary densities. The package allows the replicability of the
#' results in García-Portugués et al. (2019) <doi:10.1007/s11222-017-9790-2>.
#'
#' @author Eduardo García-Portugués.
#' @references
#' García-Portugués, E., Sørensen, M., Mardia, K. V. and Hamelryck, T. (2019)
#' Langevin diffusions on the torus: estimation and applications.
#' \emph{Statistics and Computing}, 29(2):1--22. \doi{10.1007/s11222-017-9790-2}
#' @docType package
#' @name sdetorus
#' @import graphics stats Rcpp mvtnorm colorRamps
#' @useDynLib sdetorus
#' @aliases sdetorus sdetorus-package
NULL
