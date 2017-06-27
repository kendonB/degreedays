#' @title Calculate degree days/degree-3 splines of temperature using the single sine method.
#' @description C++ function for calculating a data.frame
#' of degree days/degree-3 splines for vectors of tmin / tmax.
#' @inheritParams spl3_band_daily
#' @return \code{data.frame} with cubic spline terms.
#' @export
#' @useDynLib degreedays
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
spl3_band_ss <- function(t0, tmin, tmax, weights = NULL, parallel = FALSE){
  if(is.null(weights)){
    # rep fails with length(tmin) > .Machine$integer.max
    weights <- tmin*0 + 1
  }
  if(parallel){
    out <- data.frame(degreedays:::spl3_band_daily_par(t0, tmin, tmax, weights))
  } else {
    out <- data.frame(degreedays:::spl3_band_daily(t0, tmin, tmax, weights))
  }

  spl3_bottoms_names <- gsub("-", "n", as.character(t0))
  names(out) <- paste0("spl3_", spl3_bottoms_names, "_Inf")
  out
}
