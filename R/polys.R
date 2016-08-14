#' @title Calculate degree bins using the single sine method.
#' @description C++ function for calculating a data.frame
#' of degree bins for vectors of tmin / tmax, and vectors of bins bounds.
#' @inheritParams degree_days_band_daily
#' @return \code{data.frame} with degree days.
#' @export
#' @useDynLib degreedays
#' @importFrom Rcpp evalCpp
poly_ss <- function(tmin, tmax, degree, weights = NULL, parallel = FALSE){
  if(is.null(weights)){
    # rep fails with length(tmin) > .Machine$integer.max
    weights <- tmin*0 + 1
  }
  if(parallel){
    out <- degreedays:::sin_poly_temp_par(tmin, tmax, 1:degree, weights)
  } else {
    out <- sapply(1:degree, function(x) degreedays:::sin_poly_temp(tmin, tmax, x, weights))
  }
  out <- data.frame(out)
  names(out) <- paste0("poly_", 1:degree)
  out
}
