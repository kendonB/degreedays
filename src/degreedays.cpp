#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;
using namespace Rcpp;

#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// From http://www.jstor.org/stable/pdf/1933912.pdf

// Potentially could implement this for general degree
// Right now, band derives from a degree 1 approximation (piecewise linear)
// and bin derives from a degree 0 approximation (bins or piecewise constant).

/// Unit is days - i.e. x == 0.5 -> 12pm
// Horizontal shift parameter is irrelevant - this assumes the day starts at the minimum.
//' @title Sine function with max and min given.
//' @description Sine function evaluated at x for single sine interpolation.
//' @param x Vector of x values in [0, 1] to evaluate at.
//' @param tmin \code{double} tmin value.
//' @param tmax \code{double} tmax value.
//' @return Vector of interpolated sine values.
// [[Rcpp::export]]
NumericVector sin_estimate_NumericVector(NumericVector x, double tmin, double tmax) {
  double shift = 6.0 / 24.0;
  int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = (tmax - tmin) / 2.0 * sin(2.0 * M_PI * (x[i] - shift)) + (tmax + tmin) / 2.0;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector sin_int_estimate_NumericVector(NumericVector x, double tmin, double tmax) {
  double shift = 6.0 / 24.0;
  int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = - (tmax - tmin) / 2.0 / (2.0*M_PI) * cos(2.0*M_PI * (x[i] - shift)) +
      (tmax + tmin) / 2.0 * x[i];
  }
  return out;
}

// [[Rcpp::export]]
double sin_estimate(double x, double tmin, double tmax) {
  double shift = 6.0 / 24.0;
  return (tmax - tmin) / 2.0 * sin(2.0 * M_PI * (x - shift)) + (tmax + tmin) / 2.0;
}

// [[Rcpp::export]]
double sin_int_estimate(double x, double tmin, double tmax) {
  double shift = 6.0 / 24.0;
  return - (tmax - tmin) / 2.0 / (2.0*M_PI) * cos(2.0*M_PI * (x - shift)) +
    (tmax + tmin) / 2.0 * x;
}

// [[Rcpp::export]]
double degree_days_one(double t0, double t1,
                               double tmin, double tmax) {
  if(tmin > tmax){
    throw std::invalid_argument("tmin > tmax");
  }
  double out;
  if (R_IsNA(tmin) || R_IsNA(tmax)){
    out = NA_REAL;
  } else if (tmax == tmin && t1 == tmax){
    // Strange case where tmax == tmin and they're on the top border of
    // a band.
    out = t1 - t0;
  } else if (tmax == tmin && t0 == tmax){
    // Strange case where tmax == tmin and they're on the bottom border of
    // a band.
    out = 0;
  } else if (t0 <= tmin && t1 >= tmax) {
    // Case A where the band bounds the entire day's temperature range.
    double m = (tmax + tmin) / 2.0;
    out = m - t0;
  } else if (t0 > tmin && t1 >= tmax && t0 < tmax) {
    // Case B where the band straddles the tmax value.
    double w = (tmax - tmin) / 2.0;
    double m = (tmax + tmin) / 2.0;
    double shift = 6.0 / 24.0;
    double theta_1 = asin((t0 - m) / w) / (2.0 * M_PI) + shift;
    double first_part = 2.0 * (sin_int_estimate(12.0 / 24.0, tmin, tmax) -
                               sin_int_estimate(theta_1, tmin, tmax));
    double second_part = t0 * 2.0 * (12.0 / 24.0 - theta_1);
    out = first_part - second_part;
  } else if (t1 < tmax && t1 >= tmin) {
    // Case C where the band is contained in the bounds of tmin and tmax
    // and the Case (not in Baskerville and Emin) where
    // the band straddles the tmin value.
    double tmpt0, tmpt1, outt0, outt1, infty;
    tmpt0 = t0;
    tmpt1 = t1;
    infty = std::numeric_limits<double>::infinity();
    outt0 = degree_days_one(tmpt0, infty, tmin, tmax);
    outt1 = degree_days_one(tmpt1, infty, tmin, tmax);
    out = outt0 - outt1;
  } else if (t1 < tmin) {
    // Band is below the minimum temperature.
    out = 1.0 * (t1 - t0);
  } else if (t0 >= tmax) {
    // Band is above the maximum temperature.
    out = 0.0;
  } else {
    throw std::invalid_argument("received incorrect t0/t1/tmax/tmin numbers");
  }
  return out;
}

//' @title Calculate degree days for daily data.
//' @description Calculate degree days for daily data using C++ code.
//' @param t0 vector of lower bounds
//' @param t1 vector of upper bounds
//' @param tmin vector of tmin values (1 per day)
//' @param tmax vector of tmax values (1 per day)
//' @return num_days x num_bins \code{matrix}
// [[Rcpp::export]]
NumericMatrix degree_days_band_daily(NumericVector t0, NumericVector t1,
                                     NumericVector tmin, NumericVector tmax){
  if (t0.size() != t1.size()) {
    throw std::invalid_argument("Lengths of t0 and t1 differ.");
  }
  if (tmin.size() != tmax.size()) {
    throw std::invalid_argument("Lengths of tmin and tmax differ.");
  }

  int nrow = tmin.size(), ncol = t0.size();
  NumericMatrix out(nrow, ncol);
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      out(i, j) = degree_days_one(t0[j], t1[j], tmin[i], tmax[i]);
    }
  }
  return out;
}

struct DegDay : public Worker {

  // input matrix to read from
  const RVector<double> t0;
  const RVector<double> t1;
  const RVector<double> tmin;
  const RVector<double> tmax;

  // output matrix to write to
  RMatrix<double> output;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  DegDay(const NumericVector t0, const NumericVector t1,
         const NumericVector tmin,const NumericVector tmax,
         NumericMatrix output)
    : t0(t0), t1(t1), tmin(tmin), tmax(tmax), output(output) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for(std::size_t j = 0; j < t0.size(); j++){
        // write to output matrix
        output(i,j) = degree_days_one(t0[j], t1[j], tmin[i], tmax[i]);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix degree_days_band_daily_par(NumericVector t0, NumericVector t1,
                             NumericVector tmin, NumericVector tmax) {

  // allocate the output matrix
  int nrow = tmin.size(), ncol = t0.size();
  NumericMatrix output(nrow, ncol);

  // SquareRoot functor (pass input and output matrixes)
  DegDay degDay(t0, t1, tmin, tmax, output);

  parallelFor(0, nrow, degDay);
  // call parallelFor to do the work

  // return the output matrix
  return output;
}

// [[Rcpp::export]]
double days_in_bin_one(double t0, double t1,
                          double tmin, double tmax) {
  if (tmin > tmax) {
    throw std::invalid_argument("tmin > tmax");
  }
  double out;
  if (R_IsNA(tmin) || R_IsNA(tmax)){
    out = NA_REAL;
  } else if (tmax == tmin && t1 == tmax) {
    // Strange case where tmax == tmin and they're on the top border of a bin
    out = 0.5;
  } else if (t0 >= t1) {
    throw std::invalid_argument("t0 < t1 is not TRUE");
  } else if (tmax == tmin && t0 == tmax) {
    // Strange case where tmax == tmin and they're on the bottom border of a bin
    out = 0.5;
  } else if (t0 <= tmin && t1 >= tmax) {
    // Case A where the band bounds the entire day's temperature range.
    out = 1.0;
  } else if (t0 > tmin && t1 >= tmax && t0 < tmax) {
    // Case B where the band straddles the tmax value.
    double w = (tmax - tmin) / 2.0;
    double m = (tmax + tmin) / 2.0;
    double shift = 6.0 / 24.0;
    double theta_1 = asin((t0 - m) / w) / (2.0 * M_PI) + shift;
    out = 2.0 * (0.5 - theta_1);
  } else if (t1 < tmax && t1 > tmin) {
    // Case C where the band is contained in the bounds of tmin and tmax
    // and the Case (not in Baskerville and Emin) where
    // the band straddles the tmin value.
    double infty = 99999.0;
    double outt0 = days_in_bin_one(t0, infty, tmin, tmax);
    double outt1 = days_in_bin_one(t1, infty, tmin, tmax);
    out = outt0 - outt1;
  } else if (t1 <= tmin) {
    // Band is below the minimum temperature.
    out = 0.0;
  } else if (t0 >= tmax) {
    // Band is above the maximum temperature.
    out = 0.0;
  } else {
    throw std::invalid_argument("received incorrect t0/t1/tmax/tmin numbers");
  }
  return out;
}

//' @title Calculate days in bin for daily data.
//' @description Calculate days in bin for daily data using C++ code.
//' @param t0 vector of lower bounds
//' @param t1 vector of upper bounds
//' @param tmin vector of tmin values (1 per day)
//' @param tmax vector of tmax values (1 per day)
//' @return num_days x num_bins \code{matrix}
// [[Rcpp::export]]
NumericMatrix days_in_bin_daily(NumericVector t0, NumericVector t1,
                                NumericVector tmin, NumericVector tmax){
  if (t0.size() != t1.size()) {
    throw std::invalid_argument("Lengths of t0 and t1 differ.");
  }
  if (tmin.size() != tmax.size()) {
    throw std::invalid_argument("Lengths of tmin and tmax differ.");
  }

  int nrow = tmin.size(), ncol = t0.size();
  NumericMatrix out(nrow, ncol);
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      out(i, j) = days_in_bin_one(t0[j], t1[j], tmin[i], tmax[i]);
    }
  }
  return out;
}

struct Bins : public Worker {

  // input matrix to read from
  const RVector<double> t0;
  const RVector<double> t1;
  const RVector<double> tmin;
  const RVector<double> tmax;

  // output matrix to write to
  RMatrix<double> output;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  Bins(const NumericVector t0, const NumericVector t1,
         const NumericVector tmin,const NumericVector tmax,
         NumericMatrix output)
    : t0(t0), t1(t1), tmin(tmin), tmax(tmax), output(output) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for(std::size_t j = 0; j < t0.size(); j++){
        // write to output matrix
        output(i,j) = days_in_bin_one(t0[j], t1[j], tmin[i], tmax[i]);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix days_in_bin_daily_par(NumericVector t0, NumericVector t1,
                                         NumericVector tmin, NumericVector tmax) {

  // allocate the output matrix
  int nrow = tmin.size(), ncol = t0.size();
  NumericMatrix output(nrow, ncol);

  // SquareRoot functor (pass input and output matrixes)
  Bins bins(t0, t1, tmin, tmax, output);

  parallelFor(0, nrow, bins);
  // call parallelFor to do the work

  // return the output matrix
  return output;
}

/*** R

*/
