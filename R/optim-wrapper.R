

#' @title Optimization wrapper for likelihood-based procedures
#'
#' @description A convenient wrapper to perform local optimization of the
#' likelihood function via \code{nlm} and \code{optim} including several
#' practical utilities.
#'
#' @param minusLogLik function computing the minus log-likelihood function.
#' Must have a single argument containing a vector of length \code{p}.
#' @param region function to impose a feasibility region via a penalty. See
#' details.
#' @param penalty imposed penalty if value is not finite.
#' @param optMethod one of the following strings: \code{"nlm"},
#' \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"},
#' \code{"SANN"}, or \code{"Brent"}.
#' @param start starting values, a matrix with \code{p} columns, with each
#' entry representing a different starting value.
#' @param lower,upper bound for box constraints as in method \code{"L-BFGS-B"}
#' of \code{\link[stats]{optim}}.
#' @param selectSolution which criterion is used for selecting a solution
#' among possible ones, either \code{"lowest"}, \code{"lowestConv"} or
#' \code{"lowestLocMin"}. \code{"lowest"} returns the solution with lowest
#' value in the \code{minusLogLik} function. \code{"lowestConv"} restricts
#' the search of the best solution among the ones for which the optimizer has
#' converged. \code{"lowestLocMin"} in addition imposes that the solution is
#' guaranteed to be a local minimum by examining the Hessian matrix.
#' @param checkCircular logical indicating whether to automatically treat the
#' variables with \code{lower} and \code{upper} entries equal to \code{-pi} and
#' \code{pi} as circular. See details.
#' @param maxit maximum number of iterations.
#' @param tol tolerance for convergence (passed to \code{reltol}, \code{pgtol}
#' or \code{gradtol}).
#' @param eigTol,condTol eigenvalue and condition number tolerance for the
#' Hessian in order to guarantee a local minimum. Used only if
#' \code{selectSolution = "lowestLocMin"}.
#' @param verbose an integer from \code{0} to \code{2} if
#' \code{optMethod = "Nelder-Mead"} or from \code{0} to \code{4} otherwise
#' giving the amount of information displayed.
#' @param ... further arguments passed to the \code{optMethod} selected. See
#' options in \code{\link[stats]{nlm}} or \code{\link[stats]{optim}}.
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{par}: estimated minimizing parameters
#'   \item \code{value}: value of \code{minusLogLik} at the minimum.
#'   \item \code{convergence}: if the optimizer has converged or not.
#'   \item \code{message}: a character string giving any additional information
#'   returned by the optimizer.
#'   \item \code{eigHessian}: eigenvalues of the Hessian at the minimum. Recall
#'   that for the same solution slightly different outputs may be obtained
#'   according to the different computations of the Hessian of \code{nlm} and
#'   \code{optim}.
#'   \item \code{localMinimumGuaranteed}: tests if the Hessian is positive
#'   definite (all eigenvalues larger than the tolerance \code{eigTol} and
#'   condition number smaller than \code{condTol}).
#'   \item \code{solutionsOutput}: a list containing the complete output of
#'   the selected method for the different starting values. It includes the
#'   extra objects \code{convergence} and \code{localMinimumGuaranteed}.
#' }
#' @details If \code{checkCircular = TRUE}, then the corresponding \code{lower}
#' and \code{upper} entries of the circular parameters are set to \code{-Inf}
#' and \code{Inf}, respectively, and \code{minusLogLik} is called with the
#' \emph{principal value} of the circular argument.
#'
#' If no solution is found satisfying the criterion in \code{selectSolution},
#' NAs are returned in the elements of the main solution.
#'
#' The Hessian is only computed if \code{selectSolution = "lowestLocMin"}.
#'
#' Region feasibility can be imposed by a function with the same arguments as
#' \code{minusLogLik} that resets \code{pars} in to the boundary of the
#' feasibility region and adds a penalty proportional to the violation of the
#' feasibility region. Note that this is \emph{not the best procedure at all}
#' to solve the constrained optimization problem, but just a relatively
#' flexible and quick approach (for a more advanced treatment of restrictions,
#' see \href{https://CRAN.R-project.org/view=Optimization}{
#' optimization-focused packages}). The value must be a list with objects
#' \code{pars} and \code{penalty}. By default no region is imposed, i.e.,
#' \code{region = function(pars) list("pars" = pars, "penalty" = 0)}. Note that
#'  the Hessian is computed from the unconstrained problem, hence
#'  \code{localMinimumGuaranteed} might be \code{FALSE} even if a local minimum
#'   to the constrained problem was found.
#' @examples
#' # No local minimum
#' head(mleOptimWrapper(minusLogLik = function(x) -sum((x - 1:4)^2),
#'                      start = rbind(10:13, 1:2), selectSolution = "lowest"))
#' head(mleOptimWrapper(minusLogLik = function(x) -sum((x - 1:4)^2),
#'                      start = rbind(10:13, 1:2),
#'                      selectSolution = "lowestConv"))
#' head(mleOptimWrapper(minusLogLik = function(x) -sum((x - 1:4)^2),
#'                      start = rbind(10:13, 1:2),
#'                      selectSolution = "lowestLocMin"))
#'
#' # Local minimum
#' head(mleOptimWrapper(minusLogLik = function(x) sum((x - 1:4)^2),
#'                      start = rbind(10:13), optMethod = "BFGS"))
#' head(mleOptimWrapper(minusLogLik = function(x) sum((x - 1:4)^2),
#'                      start = rbind(10:13), optMethod = "Nelder-Mead"))
#'
#' # Function with several local minimum and a 'spurious' one
#' f <- function(x)  0.75 * (x[1] - 1)^2 -
#'                   10 / (0.1 + 0.1 * ((x[1] - 15)^2 + (x[2] - 2)^2)) -
#'                   9.5 / (0.1 + 0.1 * ((x[1] - 15)^2 + (x[2] + 2)^2))
#' plotSurface2D(x = seq(0, 20, l = 100), y = seq(-3, 3, l = 100), f = f)
#' head(mleOptimWrapper(minusLogLik = f,
#'                      start = rbind(c(15, 2), c(15, -2), c(5, 0)),
#'                      selectSolution = "lowest"))
#' head(mleOptimWrapper(minusLogLik = f,
#'                      start = rbind(c(15, 2), c(15, -2), c(5, 0)),
#'                      selectSolution = "lowestConv"))
#' head(mleOptimWrapper(minusLogLik = f,
#'                      start = rbind(c(15, 2), c(15, -2), c(5, 0)),
#'                      selectSolution = "lowestLocMin", eigTol = 0.01))
#'
#' # With constraint region
#' head(mleOptimWrapper(minusLogLik = function(x) sum((x - 1:2)^2),
#'                      start = rbind(10:11),
#'                      region = function(pars) {
#'                        x <- pars[1]
#'                        y <- pars[2]
#'                        if (y <= x^2) {
#'                          return(list("pars" = pars, "penalty" = 0))
#'                        } else {
#'                         return(list("pars" = c(sqrt(y), y),
#'                                     "penalty" = y - x^2))
#'                        }
#'                      }, lower = c(0.5, 1), upper = c(Inf, Inf),
#'                 optMethod = "Nelder-Mead", selectSolution = "lowest"))
#' head(mleOptimWrapper(minusLogLik = function(x) sum((x - 1:2)^2),
#'                      start = rbind(10:11), lower = c(0.5, 1),
#'                      upper = c(Inf, Inf),optMethod = "Nelder-Mead"))
#' @export
mleOptimWrapper <- function(minusLogLik, region = function(pars)
  list("pars" = pars, "penalty" = 0), penalty = 1e10, optMethod = "Nelder-Mead",
  start, lower = rep(-Inf, ncol(start)), upper = rep(Inf, ncol(start)),
  selectSolution = "lowestLocMin", checkCircular = TRUE, maxit = 500,
  tol = 1e-5, verbose = 0, eigTol = 1e-4, condTol = 1e4, ...) {

  # Identify which parameters are circular to allow them vary in R but taking
  # its principal values. This avoids artificial local minimums created as a
  # combination of a bad starting value (for example -pi + 0.1) and a minimum
  # placed at the other side of the interval (for example at pi - 0.1)
  if (checkCircular) {

    circularPars <- which(lower == -pi & upper == pi)
    lower[circularPars] <- -Inf
    upper[circularPars] <- Inf
    start[circularPars] <- toPiInt(start[circularPars])

  }

  # Careful with overflows and boundaries
  minusLogLikFiniteBound <- function(pars) {

    # Wrap to [-pi, pi) the circular parameters
    if (checkCircular) {

      pars[circularPars] <- toPiInt(pars[circularPars])

    }

    # Check if all the parameters are in range
    # If not: reset them and create a continuous penalty

    # Lower
    belowLow <- lower > pars
    if (any(belowLow)) {

      penaltyLow <- sqrt(sum((pars[belowLow] - lower[belowLow])^2))
      pars[belowLow] <- lower[belowLow]

    } else {

      penaltyLow <- 0

    }

    # Upper
    aboveUp <- upper < pars
    if (any(aboveUp)) {

      penaltyUp <- sqrt(sum((pars[aboveUp] - lower[aboveUp])^2))
      pars[aboveUp] <- upper[aboveUp]

    } else {

      penaltyUp <- 0

    }

    # Check if parameters are in region
    # If not: reset them and create a penalty
    reg <- region(pars)

    # Call minusLogLik
    res <- minusLogLik(reg$pars) + reg$penalty + 10 * (penaltyLow + penaltyUp)

    # Check finiteness
    if (!is.finite(res)) {

      res <- penalty

    }

    # Return values
    return(res)

  }

  # Starting values
  nStart <- dim(start)[1]
  if (is.null(nStart)) {

    nStart <- 1
    start <- matrix(start, nrow = 1)

  }

  # Compute Hessian? Takes a significant ammount of time
  hessian <- selectSolution == "lowestLocMin"

  # List with solutions
  solutions <- vector(mode = "list", length = nStart)

  # Loop on different starting values
  for (i in 1:nStart) {

    if (optMethod == "nlm") {

      # Optimization
      solutions[[i]] <- tryCatch(nlm(f = minusLogLikFiniteBound, p = start[i, ],
                                     hessian = hessian, print.level = verbose,
                                     iterlim = maxit, gradtol = tol, ...),
                                 error = function(e) {
        show(e)
        list(minimum = NA, estimate = NA, gradient = NA, hessian = NA, code = 6,
             iterations = NA, message = paste(paste(e$call, collapse = " "),
                                              e$message, sep = ": "))
      })

      # Add Hessian slot if hessian = FALSE
      if (!hessian) {

        solutions[[i]]$hessian <- NA

      }

      # Convergence message
      solutions[[i]]$message <-
        switch(solutions[[i]]$code,
               paste("Relative gradient is close to zero,",
                     "current iterate is probably solution"),
               paste("Successive iterates within tolerance,",
                     "current iterate is probably solution"),
               paste("Last global step failed to locate a point",
                     "lower than estimate. Either estimate is an",
                     "approximate local minimum of the function or",
                     "steptol is too small"),
               "Iteration limit exceeded",
               paste("Maximum step size stepmax exceeded five consecutive",
                     "times. Either the function is unbounded below, becomes",
                     "asymptotic to a finite value from above in some",
                     "direction or stepmax is too small"),
               paste("Error in", solutions[[i]]$message))

      # Check if convergence happened
      solutions[[i]]$convergence <- solutions[[i]]$code %in% 1:3

      # Wrap circular parameters
      if (checkCircular) {

        solutions[[i]]$estimate[circularPars] <-
          toPiInt(solutions[[i]]$estimate[circularPars])

      }

    } else {

      # Optimization
      control <- list(maxit = maxit, trace = verbose, ...)
      if (optMethod == "L-BFGS-B") {

        control$pgtol <- tol
        solutions[[i]] <- tryCatch(optim(par = start[i, ],
                                         fn = minusLogLikFiniteBound,
                                         method = optMethod, lower = lower,
                                         upper = upper, hessian = hessian,
                                         control = control),
                                   error = function(e) {
                                     show(e)
                                     list(par = NA, value = NA, counts = NA,
                                          convergence = 6,
                                          message = paste(
                                            paste(e$call, collapse = " "),
                                            e$message, sep = ": "),
                                          hessian = NA)
                                     })

      } else {

        control$reltol <- tol
        solutions[[i]] <- tryCatch(optim(par = start[i, ],
                                         fn = minusLogLikFiniteBound,
                                         method = optMethod, hessian = hessian,
                                         control = control),
                                   error = function(e) {
                                     show(e)
                                     list(par = NA, value = NA, counts = NA,
                                          convergence = 6,
                                          message = paste(
                                            paste(e$call, collapse = " "),
                                            e$message, sep = ": "),
                                          hessian = NA)
                                     })

      }

      # Add Hessian slot if hessian = FALSE
      if (!hessian) {

        solutions[[i]]$hessian <- NA

      }

      # Convergence message
      solutions[[i]]$message <-
        paste(switch(as.character(solutions[[i]]$convergence),
                     "0" = "Successful completion",
                     "1" = "Iteration limit maxit had been reached",
                     "10" = "Degeneracy of the Nelder-Mead simplex",
                     "51" = "Warning from L-BFGS-B method",
                     "52" = "Error from L-BFGS-B method",
                     "6" = paste("Error in", solutions[[i]]$message)),
              ifelse(is.null(solutions[[i]]$message) |
                       solutions[[i]]$convergence == 6, "",
                     paste(".", solutions[[i]]$message)), sep = "")

      # Check if convergence happened
      solutions[[i]]$convergence <- solutions[[i]]$convergence == 0

      # Wrap circular parameters
      if (checkCircular) {

        solutions[[i]]$par[circularPars] <-
          toPiInt(solutions[[i]]$par[circularPars])

      }

    }

    # Add eigendecomposition of Hessian and flag indicating if the solution
    # is a local minimum
    if (anyNA(solutions[[i]]$hessian)) {

      solutions[[i]]$eigHessian <- NA
      solutions[[i]]$localMinimumGuaranteed <- FALSE

    } else {

      solutions[[i]]$eigHessian <- eigen(solutions[[i]]$hessian,
                                         symmetric = TRUE,
                                         only.values = TRUE)$values
      solutions[[i]]$localMinimumGuaranteed <-
        all(solutions[[i]]$eigHessian > eigTol) &
        (max(abs(solutions[[i]]$eigHessian)) <
           condTol * min(abs(solutions[[i]]$eigHessian)))

    }

  }

  # Which is the best found solution?
  if (selectSolution == "lowestLocMin") {

    indSelect <- sapply(1:nStart, function(i)
      solutions[[i]]$convergence & solutions[[i]]$localMinimumGuaranteed)

  } else if (selectSolution == "lowestConv") {

    indSelect <- sapply(1:nStart, function(i) solutions[[i]]$convergence)

  } else if (selectSolution == "lowest") {

    indSelect <- rep(TRUE, nStart)

  } else {

    stop("selectSolution must be 'lowestLocMin', 'lowestConv' or 'lowest'")

  }
  if (optMethod == "nlm") {

    ind <- which.min(sapply(1:nStart, function(i)
      solutions[[i]]$minimum)[indSelect])

  } else {

    ind <- which.min(sapply(1:nStart, function(i)
      solutions[[i]]$value)[indSelect])

  }
  ind <- (1:nStart)[which(indSelect)[ind]]

  # Add complete output of the solutions
  msgHessian <- ifelse(hessian, paste0(", eigTol = ", eigTol,
                                       " and condTol = ", condTol), "")
  msg <- paste0(paste0("No solution was found with selectSolution = '",
                       selectSolution, "'"), msgHessian)
  bestSolution <- list("par" = NA, "value" = NA, "convergence" = FALSE,
                       "message" = msg, "eigHessian" = NA,
                       "localMinimumGuaranteed" = FALSE,
                       "optMethod" = optMethod, "solutionsOutput" = solutions)

  # If found a solution
  if (length(ind) > 0) {

    if (optMethod == "nlm") {

      bestSolution$par <- solutions[[ind]]$estimate
      bestSolution$value <- solutions[[ind]]$minimum

    } else {

      bestSolution$par <- solutions[[ind]]$par
      bestSolution$value <- solutions[[ind]]$value

    }

    # Add method info
    bestSolution$convergence <- solutions[[ind]]$convergence
    bestSolution$message <- solutions[[ind]]$message
    bestSolution$eigHessian <- solutions[[ind]]$eigHessian
    bestSolution$localMinimumGuaranteed <-
      solutions[[ind]]$localMinimumGuaranteed

  }

  return(bestSolution)

}
