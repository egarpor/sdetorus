.onLoad <- function(libname = find.package("sdetorus"), pkgname = "sdetorus") {

  # Assign global variables
  assign(x = "logBesselI0ScaledSpline",
         value = splinefun(x = x1, y = logBesselI0ScaledEvalGrid),
         pos = environment(logBesselI0Scaled))
  assign(x = "a1InvSpline",
         value = splinefun(x = x2, y = a1InvEvalGrid),
         envir = environment(a1Inv))

  # CRAN NOTE avoidance
  if (getRversion() >= "2.15.1") {

    utils::globalVariables(c("a1InvSpline", "logBesselI0ScaledSpline"))

  }
  invisible()

}
