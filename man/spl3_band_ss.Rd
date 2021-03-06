% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spl3.R
\name{spl3_band_ss}
\alias{spl3_band_ss}
\title{Calculate degree days/degree-3 splines of temperature using the single sine method.}
\usage{
spl3_band_ss(t0, tmin, tmax, weights = NULL, parallel = FALSE)
}
\arguments{
\item{t0}{vector of lower bounds}

\item{tmin}{vector of tmin values (1 per day)}

\item{tmax}{vector of tmax values (1 per day)}

\item{weights}{vector of optional weights to multiply output by. Default is 1.}
}
\value{
\code{data.frame} with cubic spline terms.
}
\description{
C++ function for calculating a data.frame
of degree days/degree-3 splines for vectors of tmin / tmax.
}
