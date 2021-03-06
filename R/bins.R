#' @title Calculate degree bins using the single sine method.
#' @description C++ function for calculating a data.frame
#' of degree bins for vectors of tmin / tmax, and vectors of bins bounds.
#' @inheritParams spl1_band_daily
#' @return \code{data.frame} with degree days.
#' @export
#' @useDynLib degreedays
#' @importFrom Rcpp evalCpp
bins_ss <- function(t0, t1, tmin, tmax, weights = NULL, parallel = FALSE){
  if(is.null(weights)){
    # rep fails with length(tmin) > .Machine$integer.max
    weights <- tmin*0 + 1
  }
  if(parallel){
    out <- data.frame(degreedays:::days_in_bin_daily_par(t0, t1, tmin, tmax, weights))
  } else {
    out <- data.frame(degreedays:::days_in_bin_daily(t0, t1, tmin, tmax, weights))
  }

  bin_bottoms_names <- gsub("-", "n", as.character(t0))
  bin_tops_names <- gsub("-", "n", as.character(t1))
  names(out) <- paste0("bin_", bin_bottoms_names, "_", bin_tops_names)
  out
}