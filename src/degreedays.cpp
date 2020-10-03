#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// From http://www.jstor.org/stable/pdf/1933912.pdf

// Potentially could implement this for general spline degree
// Right now, only linear splines, degree zero splines (labeled as bins),
// and polynomials have been implemented.

/// Unit is days of temperature - i.e. x == 0.5 -> 12pm
// Horizontal shift parameter is irrelevant for the single sine method
// - this assumes the day starts at the minimum.
//' @title Sine function with max and min given.
//' @description Sine function evaluated at x for single sine interpolation.
//' @param x Vector of x values in [0, 1] to evaluate at.
//' @param tmin \code{double} tmin value.
//' @param tmax \code{double} tmax value.
//' @return Vector of interpolated sine values.
// [[Rcpp::export]]
NumericVector sin_estimate_NumericVector(NumericVector x, double tmin, double tmax) {
  double shift = 6.0 / 24.0 + 0.0 + 0.0;
  std::size_t n = x.size();
  NumericVector out(n);
  for(std::size_t i = 0; i < n; ++i) {
    out[i] = (tmax - tmin) / 2.0 * sin(2.0 * M_PI * (x[i] - shift)) + (tmax + tmin) / 2.0;
  }
  return out;
}

/// Unit is days of temperature - i.e. x == 0.5 -> 12pm
// Horizontal shift parameter is irrelevant for the single sine method
// - this assumes the day starts at the minimum.
//' @title Sine function with max and min given.
//' @description Sine function evaluated at x for single sine interpolation.
//' @param x Vector of x values in [0, 1] to evaluate at.
//' @param tmin \code{double} tmin value.
//' @param tmax \code{double} tmax value.
//' @return Vector of interpolated sine values.
// [[Rcpp::export]]
NumericVector sin_estimate_NumericVectorxtmintmax(NumericVector x, NumericVector tmin, NumericVector tmax) {
  double shift = 6.0 / 24.0 + 0.0 + 0.0;
  std::size_t n = x.size();
  if (tmin.size() != x.size()) {
    throw std::invalid_argument("Lengths of tmin and x differ.");
  }
  if (tmax.size() != x.size()) {
    throw std::invalid_argument("Lengths of tmax and x differ.");
  }
  NumericVector out(n);
  for(std::size_t i = 0; i < n; ++i) {
    out[i] = (tmax[i] - tmin[i]) / 2.0 * sin(2.0 * M_PI * (x[i] - shift)) + (tmax[i] + tmin[i]) / 2.0;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector sin_int_estimate_NumericVector(NumericVector x, double tmin, double tmax) {
  double shift = 6.0 / 24.0;
  std::size_t n = x.size();
  NumericVector out(n);
  for(std::size_t i = 0; i < n; ++i) {
    out[i] = - (tmax - tmin) / 2.0 / (2.0*M_PI) * cos(2.0*M_PI * (x[i] - shift)) +
      (tmax + tmin) / 2.0 * x[i];
  }
  return out;
}

// [[Rcpp::export]]
double sin_poly_temp_one(double tmin, double tmax, double degree, double weight) {
  double out;
  if(degree == 1.0){
    out = (tmax + tmin)/2.0;
  } else if (degree == 2.0){
    // Copy pasted from integral-calculator.com
    out = (3.0*pow(tmax, 2.0)+2.0*tmin*tmax+3.0*pow(tmin,2.0))/8.0;
  } else if (degree == 3.0){
    // Copy pasted from integral-calculator.com
    out = (5.0*pow(tmax,3.0)+3.0*tmin*pow(tmax,2.0)+3.0*pow(tmin,2.0)*tmax+5.0*pow(tmin,3.0))/16.0;
  } else if (degree == 4.0){
    // Copy pasted from integral-calculator.com
    out = (35.0*pow(tmax, 4.0)+20.0*tmin*pow(tmax, 3.0)+18.0*pow(tmin, 2.0)*pow(tmax, 2.0)+20.0*pow(tmin, 3.0)*tmax+35.0*pow(tmin, 4.0))/128.0;
  } else if (degree == 5.0){
    // Copy pasted from integral-calculator.com
    out = (63.0*pow(tmax, 5.0)+35.0*tmin*pow(tmax, 4.0)+30.0*pow(tmin, 2.0)*pow(tmax, 3.0)+30.0*pow(tmin, 3.0)*pow(tmax, 2.0)+35.0*pow(tmin, 4.0)*tmax+63.0*pow(tmin, 5.0))/256.0;
  } else if (degree == 6.0){
    // Copy pasted from integral-calculator.com
    out = (231.0*pow(tmax, 6.0)+126.0*tmin*pow(tmax, 5.0)+105.0*pow(tmin, 2.0)*pow(tmax, 4.0)+100.0*pow(tmin, 3.0)*pow(tmax, 3.0)+105.0*pow(tmin, 4.0)*pow(tmax, 2.0)+126.0*pow(tmin, 5.0)*tmax+231.0*pow(tmin, 6.0))/1024.0;
  } else if (degree == 7.0){
    // Copy pasted from integral-calculator.com
    out = (429.0*pow(tmax, 7.0)+231.0*tmin*pow(tmax, 6.0)+189.0*pow(tmin, 2.0)*pow(tmax, 5.0)+175.0*pow(tmin, 3.0)*pow(tmax, 4.0)+175.0*pow(tmin, 4.0)*pow(tmax, 3.0)+189.0*pow(tmin, 5.0)*pow(tmax, 2.0)+231.0*pow(tmin, 6.0)*tmax+429.0*pow(tmin,7.0))/2048.0;
  } else if (degree == 8.0){
    // Copy pasted from integral-calculator.com
    double b = (tmax + tmin) / 2.0;
    double a = (tmax - tmin) / 2.0;
    out = (128.0*pow(b,8.0)+1792.0*pow(a,2.0)*pow(b,6.0)+3360.0*pow(a,4.0)*pow(b,4.0)+1120.0*pow(a, 6.0)*pow(b, 2.0)+35.0*pow(a,8.0))/128.0;
  } else {
    throw std::invalid_argument("Integer polynomial degrees up to 8 supported.");
  }
  out = out * weight;
  return out;
}

// [[Rcpp::export]]
NumericVector sin_poly_temp(NumericVector tmin, NumericVector tmax, double degree, NumericVector weights) {
  std::size_t n = tmin.size();
  if (tmin.size() != tmax.size()) {
    throw std::invalid_argument("Lengths of tmin and tmax differ.");
  }
  NumericVector out(n);
  if(degree == 1.0){
    for(std::size_t i = 0; i < n; ++i) {
      // Copy pasted from integral-calculator.com
      out[i] = (tmax[i] + tmin[i])/2.0*weights[i];
    }
  } else if (degree == 2.0){
    for(std::size_t i = 0; i < n; ++i) {
      // Copy pasted from integral-calculator.com
      out[i] = (3.0*pow(tmax[i], 2.0)+2.0*tmin[i]*tmax[i]+3.0*pow(tmin[i],2.0))/8.0*weights[i];
    }
  } else if (degree == 3.0){
    for(std::size_t i = 0; i < n; ++i) {
      // Copy pasted from integral-calculator.com
      out[i] = (5.0*pow(tmax[i],3.0)+3.0*tmin[i]*pow(tmax[i],2.0)+3.0*pow(tmin[i],2.0)*tmax[i]+5.0*pow(tmin[i],3.0))/16.0*weights[i];
    }
  } else if (degree == 4.0){
    for(std::size_t i = 0; i < n; ++i) {
      // Copy pasted from integral-calculator.com
      out[i] = (35.0*pow(tmax[i], 4.0)+20.0*tmin[i]*pow(tmax[i], 3.0)+18.0*pow(tmin[i], 2.0)*pow(tmax[i], 2.0)+20.0*pow(tmin[i], 3.0)*tmax[i]+35.0*pow(tmin[i], 4.0))/128.0*weights[i];
    }
  } else if (degree == 5.0){
    for(std::size_t i = 0; i < n; ++i) {
      // Copy pasted from integral-calculator.com
      out[i] = (63.0*pow(tmax[i], 5.0)+35.0*tmin[i]*pow(tmax[i], 4.0)+30.0*pow(tmin[i], 2.0)*pow(tmax[i], 3.0)+30.0*pow(tmin[i], 3.0)*pow(tmax[i], 2.0)+35.0*pow(tmin[i], 4.0)*tmax[i]+63.0*pow(tmin[i], 5.0))/256.0*weights[i];
    }
  } else if (degree == 6.0){
    for(std::size_t i = 0; i < n; ++i) {
      // Copy pasted from integral-calculator.com
      out[i] = (231.0*pow(tmax[i], 6.0)+126.0*tmin[i]*pow(tmax[i], 5.0)+105.0*pow(tmin[i], 2.0)*pow(tmax[i], 4.0)+100.0*pow(tmin[i], 3.0)*pow(tmax[i], 3.0)+105.0*pow(tmin[i], 4.0)*pow(tmax[i], 2.0)+126.0*pow(tmin[i], 5.0)*tmax[i]+231.0*pow(tmin[i], 6.0))/1024.0*weights[i];
    }
  } else if (degree == 7.0){
    for(std::size_t i = 0; i < n; ++i) {
      // Copy pasted from integral-calculator.com
      out[i] = (429.0*pow(tmax[i], 7.0)+231.0*tmin[i]*pow(tmax[i], 6.0)+189.0*pow(tmin[i], 2.0)*pow(tmax[i], 5.0)+175.0*pow(tmin[i], 3.0)*pow(tmax[i], 4.0)+175.0*pow(tmin[i], 4.0)*pow(tmax[i], 3.0)+189.0*pow(tmin[i], 5.0)*pow(tmax[i], 2.0)+231.0*pow(tmin[i], 6.0)*tmax[i]+429.0*pow(tmin[i],7.0))/2048.0*weights[i];
    }
  } else if (degree == 8.0){
    for(std::size_t i = 0; i < n; ++i) {
      // Copy pasted from integral-calculator.com
      double b = (tmax[i] + tmin[i]) / 2.0;
      double a = (tmax[i] - tmin[i]) / 2.0;
      out[i] = (128.0*pow(b,8.0)+1792.0*pow(a,2.0)*pow(b,6.0)+3360.0*pow(a,4.0)*pow(b,4.0)+1120.0*pow(a, 6.0)*pow(b, 2.0)+35.0*pow(a,8.0))/128.0*weights[i];
    }
  }
  return out;
}


struct Polys : public Worker {

  // input matrix to read from
  const RVector<double> tmin;
  const RVector<double> tmax;
  const RVector<double> degree;
  const RVector<double> weights;

  // output matrix to write to
  RMatrix<double> output;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  Polys(const NumericVector tmin, const NumericVector tmax,
        const NumericVector degree, const NumericVector weights,
        NumericMatrix output)
    : tmin(tmin), tmax(tmax), degree(degree), weights(weights), output(output) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for(std::size_t j = 0; j < degree.size(); j++){
        // write to output matrix
        output(i,j) = sin_poly_temp_one(tmin[i], tmax[i], degree[j], weights[i]);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix sin_poly_temp_par(NumericVector tmin, NumericVector tmax,
                                NumericVector degree, NumericVector weights) {

  // allocate the output matrix
  std::size_t nrow = tmin.size(), ncol = degree.size();
  NumericMatrix output(nrow, ncol);

  // SquareRoot functor (pass input and output matrixes)
  Polys polys(tmin, tmax, degree, weights, output);

  parallelFor(0, nrow, polys);
  // call parallelFor to do the work

  // return the output matrix
  return output;
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
double sin_cubed_int_estimate(double x, double tmin, double tmax) {
  double shift = 6.0 / 24.0;
  // Calculated at http://www.integral-calculator.com/
  double a = (tmax - tmax) / 2.0;
  double b = 2.0 * M_PI;
  double c = shift;
  double d = (tmax + tmin) / 2.0;
  return (pow(a, 3.0)*cos(3.0*b*(x-c))-
          9.0*pow(a, 2.0)*d*sin(2.0*b*(x-c))-
          9.0*a*(4.0*pow(d, 2.0)+pow(a,2.0))*cos(b*(x-c))+
          6.0*b*d*(2.0*pow(d, 2.0)+3*pow(a, 2.0))*x)/(12.0*b);
}

// [[Rcpp::export]]
double spl1_one(double t0, double t1,
                double tmin, double tmax, double weight) {
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
    out = 0.0;
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
    outt0 = spl1_one(tmpt0, infty, tmin, tmax, 1.0);
    outt1 = spl1_one(tmpt1, infty, tmin, tmax, 1.0);
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
  out = out * weight;
  return out;
}

//' @title Calculate degree days/degree-1 splines for daily data.
//' @description Calculate degree days/degree-1 splines for daily data using C++ code.
//' @param t0 vector of lower bounds
//' @param t1 vector of upper bounds
//' @param tmin vector of tmin values (1 per day)
//' @param tmax vector of tmax values (1 per day)
//' @param tmax vector of tmax values (1 per day)
//' @param weights vector of optional weights to multiply output by. Default is 1.
//' @return num_days x num_bins \code{matrix}
// [[Rcpp::export]]
NumericMatrix spl1_band_daily(NumericVector t0, NumericVector t1,
                              NumericVector tmin, NumericVector tmax,
                              NumericVector weights){
  if (t0.size() != t1.size()) {
    throw std::invalid_argument("Lengths of t0 and t1 differ.");
  }
  if (tmin.size() != tmax.size()) {
    throw std::invalid_argument("Lengths of tmin and tmax differ.");
  }

  std::size_t nrow = tmin.size(), ncol = t0.size();
  NumericMatrix out(nrow, ncol);
  for (std::size_t i = 0; i < nrow; ++i) {
    for (std::size_t j = 0; j < ncol; ++j) {
      out(i, j) = spl1_one(t0[j], t1[j], tmin[i], tmax[i], weights[i]);
    }
  }
  return out;
}

struct Spl1 : public Worker {

  // input matrix to read from
  const RVector<double> t0;
  const RVector<double> t1;
  const RVector<double> tmin;
  const RVector<double> tmax;
  const RVector<double> weights;

  // output matrix to write to
  RMatrix<double> output;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  Spl1(const NumericVector t0, const NumericVector t1,
       const NumericVector tmin, const NumericVector tmax,
       const NumericVector weights,
       NumericMatrix output)
    : t0(t0), t1(t1), tmin(tmin), tmax(tmax), weights(weights), output(output) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for(std::size_t j = 0; j < t0.size(); j++){
        // write to output matrix
        output(i,j) = spl1_one(t0[j], t1[j], tmin[i], tmax[i], weights[i]);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix spl1_band_daily_par(NumericVector t0, NumericVector t1,
                                  NumericVector tmin, NumericVector tmax,
                                  NumericVector weights) {

  // allocate the output matrix
  std::size_t nrow = tmin.size(), ncol = t0.size();
  NumericMatrix output(nrow, ncol);

  // SquareRoot functor (pass input and output matrixes)
  Spl1 spl1(t0, t1, tmin, tmax, weights, output);

  parallelFor(0, nrow, spl1);
  // call parallelFor to do the work

  // return the output matrix
  return output;
}

// [[Rcpp::export]]
double spl3_one(double t0, double tmin,
                double tmax, double weight) {
  if(tmin > tmax){
    throw std::invalid_argument("tmin > tmax");
  }
  double out;
  if (R_IsNA(tmin) || R_IsNA(tmax)){
    out = NA_REAL;
  } else if (tmax == tmin && t0 == tmax){
    // Strange case where tmax == tmin and they're on the bottom border of
    // the band.
    out = 0.0;
  } else if (t0 > tmin && t0 < tmax) {
    // Case B where the band straddles the tmax value.
    double w = (tmax - tmin) / 2.0;
    double m = (tmax + tmin) / 2.0;
    double shift = 6.0 / 24.0;
    double theta_1 = asin((t0 - m) / w) / (2.0 * M_PI) + shift;
    double first_part = 2.0 * (sin_cubed_int_estimate(12.0 / 24.0, tmin, tmax) -
                               sin_cubed_int_estimate(theta_1, tmin, tmax));
    double second_part = t0 * 2.0 * (12.0 / 24.0 - theta_1);
    out = first_part - second_part;
  } else if (t0 >= tmax) {
    // Band is above the maximum temperature.
    out = 0.0;
  } else if (t0 <= tmin) {
    // Case A where the band bounds the entire day's temperature range.
    double m = (tmax + tmin) / 2.0;
    double w = (tmax - tmin) / 2.0;
    double d = m - t0;
    return pow(d, 3.0)+3*pow(w, 2.0)*d/2.0;
  } else {
    throw std::invalid_argument("received incorrect t0/tmax/tmin numbers");
  }
  out = out * weight;
  return out;
}

//' @title Calculate degree days/degree-3 splines for daily data.
//' @description Calculate degree days/degree-3 splines for daily data using C++ code.
//' @param t0 vector of lower bounds
//' @param tmin vector of tmin values (1 per day)
//' @param tmax vector of tmax values (1 per day)
//' @param tmax vector of tmax values (1 per day)
//' @param weights vector of optional weights to multiply output by. Default is 1.
//' @return num_days x num_bands \code{matrix}
// [[Rcpp::export]]
NumericMatrix spl3_band_daily(NumericVector t0,
                              NumericVector tmin, NumericVector tmax,
                              NumericVector weights){
  if (tmin.size() != tmax.size()) {
    throw std::invalid_argument("Lengths of tmin and tmax differ.");
  }

  std::size_t nrow = tmin.size(), ncol = t0.size();
  NumericMatrix out(nrow, ncol);
  for (std::size_t i = 0; i < nrow; ++i) {
    for (std::size_t j = 0; j < ncol; ++j) {
      out(i, j) = spl3_one(t0[j], tmin[i], tmax[i], weights[i]);
    }
  }
  return out;
}

struct Spl3 : public Worker {

  // input matrix to read from
  const RVector<double> t0;
  const RVector<double> tmin;
  const RVector<double> tmax;
  const RVector<double> weights;

  // output matrix to write to
  RMatrix<double> output;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  Spl3(const NumericVector t0,
       const NumericVector tmin, const NumericVector tmax,
       const NumericVector weights,
       NumericMatrix output)
    : t0(t0), tmin(tmin), tmax(tmax), weights(weights), output(output) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for(std::size_t j = 0; j < t0.size(); j++){
        // write to output matrix
        output(i,j) = spl3_one(t0[j], tmin[i], tmax[i], weights[i]);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix spl3_band_daily_par(NumericVector t0,
                                  NumericVector tmin, NumericVector tmax,
                                  NumericVector weights) {

  // allocate the output matrix
  std::size_t nrow = tmin.size(), ncol = t0.size();
  NumericMatrix output(nrow, ncol);

  // SquareRoot functor (pass input and output matrixes)
  Spl3 spl3(t0, tmin, tmax, weights, output);

  parallelFor(0, nrow, spl3);
  // call parallelFor to do the work

  // return the output matrix
  return output;
}

// [[Rcpp::export]]
double days_in_bin_one(double t0, double t1,
                       double tmin, double tmax, double weight) {
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
    double outt0 = days_in_bin_one(t0, infty, tmin, tmax, 1.0);
    double outt1 = days_in_bin_one(t1, infty, tmin, tmax, 1.0);
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
  out = out * weight;
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
                                NumericVector tmin, NumericVector tmax,
                                NumericVector weights){
  if (t0.size() != t1.size()) {
    throw std::invalid_argument("Lengths of t0 and t1 differ.");
  }
  if (tmin.size() != tmax.size()) {
    throw std::invalid_argument("Lengths of tmin and tmax differ.");
  }

  std::size_t nrow = tmin.size(), ncol = t0.size();
  NumericMatrix out(nrow, ncol);
  for (std::size_t i = 0; i < nrow; ++i) {
    for (std::size_t j = 0; j < ncol; ++j) {
      out(i, j) = days_in_bin_one(t0[j], t1[j], tmin[i], tmax[i], weights[i]);
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
  const RVector<double> weights;

  // output matrix to write to
  RMatrix<double> output;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  Bins(const NumericVector t0, const NumericVector t1,
       const NumericVector tmin, const NumericVector tmax,
       const NumericVector weights,
       NumericMatrix output)
    : t0(t0), t1(t1), tmin(tmin), tmax(tmax), weights(weights), output(output) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for(std::size_t j = 0; j < t0.size(); j++){
        // write to output matrix
        output(i,j) = days_in_bin_one(t0[j], t1[j], tmin[i], tmax[i], weights[i]);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix days_in_bin_daily_par(NumericVector t0, NumericVector t1,
                                    NumericVector tmin, NumericVector tmax,
                                    NumericVector weights) {

  // allocate the output matrix
  std::size_t nrow = tmin.size(), ncol = t0.size();
  NumericMatrix output(nrow, ncol);

  // SquareRoot functor (pass input and output matrixes)
  Bins bins(t0, t1, tmin, tmax, weights, output);

  parallelFor(0, nrow, bins);
  // call parallelFor to do the work

  // return the output matrix
  return output;
}

struct InnerProduct : public Worker
{
  // source vectors
  const RVector<double> x;
  const RVector<double> y;

  // product that I have accumulated
  double product;

  // constructors
  InnerProduct(const NumericVector x, const NumericVector y)
    : x(x), y(y), product(0) {}
  InnerProduct(const InnerProduct& innerProduct, Split)
    : x(innerProduct.x), y(innerProduct.y), product(0) {}

  // process just the elements of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    product += std::inner_product(x.begin() + begin,
                                  x.begin() + end,
                                  y.begin() + begin,
                                  0.0);
  }

  // join my value with that of another InnerProduct
  void join(const InnerProduct& rhs) {
    product += rhs.product;
  }
};

//' @title Inner product calculator
//' @export
// [[Rcpp::export]]
double parInnProd(NumericVector x, NumericVector y) {

  // declare the InnerProduct instance that takes a pointer to the vector data
  InnerProduct innerProduct(x, y);

  // call paralleReduce to start the work
  parallelReduce(0, x.length(), innerProduct);

  // return the computed product
  return innerProduct.product;
}


/*** R

*/