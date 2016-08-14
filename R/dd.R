#' @title Calculate degree days using the single sine method.
#' @description C++ function for calculating a data.frame
#' of degree days for vectors of tmin / tmax, and vectors of degree day bands.
#' @inheritParams degree_days_band_daily
#' @return \code{data.frame} with degree days.
#' @export
#' @useDynLib degreedays
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
dd_band_ss <- function(t0, t1, tmin, tmax, weights = NULL, parallel = FALSE){
  if(is.null(weights)){
    # rep fails with length(tmin) > .Machine$integer.max
    weights <- tmin*0 + 1
  }
  if(parallel){
    out <- data.frame(degreedays:::degree_days_band_daily_par(t0, t1, tmin, tmax, weights))
  } else {
    out <- data.frame(degreedays:::degree_days_band_daily(t0, t1, tmin, tmax, weights))
  }

  dd_bottoms_names <- gsub("-", "n", as.character(t0))
  dd_tops_names <- gsub("-", "n", as.character(t1))
  names(out) <- paste0("dd_", dd_bottoms_names, "_", dd_tops_names)
  out
}
