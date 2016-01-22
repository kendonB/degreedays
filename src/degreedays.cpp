#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
using namespace Rcpp;

// From http://www.jstor.org/stable/pdf/1933912.pdf

// Potentially could implement this for general degree
// Right now, band derives from a degree 1 approximation (piecewise linear)
// and bin derives from a degree 0 approximation (bins or piecewise constant).

/// Unit is days - i.e. x == 0.5 -> 12pm
// Horizontal shift parameter is irrelevant - this assumes the day starts at the minimum.
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
NumericVector degree_days_band(NumericVector t0, NumericVector t1,
                               double tmin, double tmax) {
  int n0 = t0.size();
  int n1 = t1.size();
  if (n0 != n1) {
    throw std::invalid_argument("Lengths of t0 and t1 differ.");
  }
  NumericVector out(n0);
  for(int i = 0; i < n0; ++i) {
    if (R_IsNA(tmin) || R_IsNA(tmax)){
      out[i] = NA_REAL;
    } else if (t0[i] < tmin && t1[i] >= tmax) {
      // Case A where the band bounds the entire day's temperature range.
      double m = (tmax + tmin) / 2.0;
      out[i] = m - t0[i];
    } else if (t0[i] >= tmin && t1[i] >= tmax && t0[i] < tmax) {
      // Case B where the band straddles the tmax value.
      double w = (tmax - tmin) / 2.0;
      double m = (tmax + tmin) / 2.0;
      double shift = 6.0 / 24.0;
      double theta_1 = asin((t0[i] - m) / w) / (2.0 * M_PI) + shift;
      double first_part = 2.0 * (sin_int_estimate(12.0 / 24.0, tmin, tmax) -
                                 sin_int_estimate(theta_1, tmin, tmax));
      double second_part = t0[i] * 2.0 * (12.0 / 24.0 - theta_1);
      out[i] = first_part - second_part;
    } else if (t1[i] < tmax && t1[i] >= tmin) {
      // Case C where the band is contained in the bounds of tmin and tmax
      // and the Case (not in Baskerville and Emin) where
      // the band straddles the tmin value.
      NumericVector tmpt0(1), tmpt1(1), outt0(1), outt1(1), infty(1);
      tmpt0[0] = t0[i];
      tmpt1[0] = t1[i];
      infty[0] = std::numeric_limits<double>::infinity();
      outt0 = degree_days_band(tmpt0, infty, tmin, tmax);
      outt1 = degree_days_band(tmpt1, infty, tmin, tmax);
      out[i] = outt0[0] - outt1[0];
    } else if (t1[i] < tmin) {
      // Band is below the minimum temperature.
      out[i] = 1.0 * (t1[i] - t0[i]);
    } else if (t0[i] >= tmax) {
      // Band is above the maximum temperature.
      out[i] = 0.0;
    } else {
      throw std::invalid_argument("received incorrect t0[i]/t1[i]/tmax/tmin numbers");
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector days_in_bin(NumericVector t0, NumericVector t1,
                          double tmin, double tmax) {
  int n0 = t0.size();
  int n1 = t1.size();
  if (n0 != n1) {
    throw std::invalid_argument("Lengths of t0 and t1 differ.");
  }
  NumericVector out(n0);
  for(int i = 0; i < n0; ++i) {
    if (R_IsNA(tmin) || R_IsNA(tmax)){
      out[i] = NA_REAL;
    } else if (t0[i] < tmin && t1[i] >= tmax) {
      // Case A where the band bounds the entire day's temperature range.
      out[i] = 1.0;
    } else if (t0[i] >= tmin && t1[i] >= tmax && t0[i] < tmax) {
      // Case B where the band straddles the tmax value.
      double w = (tmax - tmin) / 2.0;
      double m = (tmax + tmin) / 2.0;
      double shift = 6.0 / 24.0;
      double theta_1 = asin((t0[i] - m) / w) / (2.0 * M_PI) + shift;
      out[i] = 2.0 * (0.5 - theta_1);
    } else if (t1[i] < tmax && t1[i] >= tmin) {
      // Case C where the band is contained in the bounds of tmin and tmax
      // and the Case (not in Baskerville and Emin) where
      // the band straddles the tmin value.
      NumericVector tmpt0(1), tmpt1(1), outt0(1), outt1(1), infty(1);
      tmpt0[0] = t0[i];
      tmpt1[0] = t1[i];
      infty[0] = 99999.0;
      outt0 = days_in_bin(tmpt0, infty, tmin, tmax);
      outt1 = days_in_bin(tmpt1, infty, tmin, tmax);
      out[i] = outt0[0] - outt1[0];
    } else if (t1[i] < tmin) {
      // Band is below the minimum temperature.
      out[i] = 0.0;
    } else if (t0[i] >= tmax) {
      // Band is above the maximum temperature.
      out[i] = 0.0;
    } else {
      throw std::invalid_argument("received incorrect t0[i]/t1[i]/tmax/tmin numbers");
    }
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
  NumericVector tmprow(t0.size() - 1);
  for (int i = 0; i < nrow; ++i) {
    tmprow = days_in_bin(t0, t1, tmin[i], tmax[i]);
    for (int j = 0; j < ncol; ++j) {
      out(i, j) = tmprow[j];
    }
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
  NumericVector tmprow(t0.size() - 1);
  for (int i = 0; i < nrow; ++i) {
    tmprow = degree_days_band(t0, t1, tmin[i], tmax[i]);
    for (int j = 0; j < ncol; ++j) {
      out(i, j) = tmprow[j];
    }
  }
  return out;
}
/*** R
#degree_days_band(14.9999999, 30.000001, 15, 30)
#  degree_days_band(15.0000001, 30.000001, 15, 30)
#  degree_days_band(t0 = 14.9999999, t1 = 29.999999, tmin = 15, tmax = 30)
#  3 degree_days_band(t0 = 15.0000001, t1 = 29.999999, tmin = 15, tmax = 30)
#   sum(degree_days_band(t0 = c(15, 16, 17, 18, 19, 20),
#                        t1 = c(16, 17, 18, 19, 20, 21), tmin = 15, tmax = 20))
#   degree_days_band(t0 = 30, t1 = 30.01, tmin = 15, tmax = 30)
#   degree_days_band(t0 = 31.0, t1 = 32.0, tmin = 15.0, tmax = 30.0)
#   degree_days_band(t0 = c(14.9999999, 15.0000001, 14.9999999, 15.0000001),
#                    t1 = c(30.000001, 30.000001, 29.999999, 29.999999),
#                    tmin = 15, tmax = 30)
# setwd("..")
*/
