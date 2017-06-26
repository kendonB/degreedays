#' @title Calculate degree days/degree-1 splines of temperature using the single sine method.
#' @description C++ function for calculating a data.frame
#' of degree days/degree-1 splines for vectors of tmin / tmax.
#' @inheritParams spl1_band_daily
#' @return \code{data.frame} with degree days.
#' @export
#' @useDynLib degreedays
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
spl1_band_ss <- function(t0, t1, tmin, tmax, weights = NULL, parallel = FALSE){
  if(is.null(weights)){
    # rep fails with length(tmin) > .Machine$integer.max
    weights <- tmin*0 + 1
  }
  if(parallel){
    out <- data.frame(degreedays:::spl1_band_daily_par(t0, t1, tmin, tmax, weights))
  } else {
    out <- data.frame(degreedays:::spl1_band_daily(t0, t1, tmin, tmax, weights))
  }

  spl1_bottoms_names <- gsub("-", "n", as.character(t0))
  spl1_tops_names <- gsub("-", "n", as.character(t1))
  names(out) <- paste0("spl1_", spl1_bottoms_names, "_", spl1_tops_names)
  out
}
