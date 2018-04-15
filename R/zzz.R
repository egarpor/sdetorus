.onLoad <- function(libname = find.package("sdetorus"), pkgname = "sdetorus") {

  # CRAN NOTE avoidance
  if (getRversion() >= "2.15.1") {
    
    utils::globalVariables(c("a1InvSpline", "logBesselI0ScaledSpline"))
    
  }
  invisible()
  
}
