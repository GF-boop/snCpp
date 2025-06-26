// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include "distributions.hpp"
#include "special_functions.hpp"

// --- Skew-Normal Family ---

/**
 * @title Skew-Normal Probability Density Function (C++ backend)
 * @description Internal C++ function for the skew-normal PDF. Not intended for direct user call.
 * @param x A numeric vector of quantiles.
 * @param xi Location parameter.
 * @param omega Scale parameter.
 * @param alpha Skewness parameter.
 * @param log_d Logical; if TRUE, probabilities are given as log(p).
 * @return A numeric vector of PDF values.
 * @keywords internal
 * @export
 */
// [[Rcpp::export]]
Rcpp::NumericVector dsn_cpp(Rcpp::NumericVector x, double xi, double omega, double alpha, bool log_d) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = sn::pdf(x[i], xi, omega, alpha, log_d);
  }
  return out;
}

/**
 * @title Skew-Normal Cumulative Distribution Function (C++ backend)
 * @description Internal C++ function for the skew-normal CDF. Not intended for direct user call.
 * @param x A numeric vector of quantiles.
 * @param xi Location parameter.
 * @param omega Scale parameter.
 * @param alpha Skewness parameter.
 * @param lower_tail Logical; if TRUE, probabilities are P(X <= x), otherwise P(X > x).
 * @param log_p Logical; if TRUE, probabilities p are given as log(p).
 * @return A numeric vector of CDF values.
 * @keywords internal
 * @export
 */
// [[Rcpp::export]]
Rcpp::NumericVector psn_cpp(Rcpp::NumericVector x, double xi, double omega, double alpha, bool lower_tail, bool log_p) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = sn::cdf(x[i], xi, omega, alpha, lower_tail, log_p);
  }
  return out;
}

/**
 * @title Skew-Normal Quantile Function (C++ backend)
 * @description Internal C++ function for the skew-normal quantile. Not intended for direct user call.
 * @param p A numeric vector of probabilities.
 * @param xi Location parameter.
 * @param omega Scale parameter.
 * @param alpha Skewness parameter.
 * @param tau Generalised skewness parameter (currently only tau=0 supported).
 * @param tol Tolerance for numerical approximation.
 * @return A numeric vector of quantile values.
 * @keywords internal
 * @export
 */
// [[Rcpp::export]]
Rcpp::NumericVector qsn_cpp(Rcpp::NumericVector p, double xi, double omega, double alpha, double tau, double tol) {
  int n = p.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = sn::quantile(p[i], xi, omega, alpha, tau, tol);
  }
  return out;
}

/**
 * @title Skew-Normal Random Variate Generation (C++ backend)
 * @description Internal C++ function for generating random variates from the skew-normal distribution. Not intended for direct user call.
 * @param n Number of observations.
 * @param xi Location parameter.
 * @param omega Scale parameter.
 * @param alpha Skewness parameter.
 * @param tau Generalised skewness parameter (currently only tau=0 supported).
 * @return A numeric vector of random variates.
 * @keywords internal
 * @export
 */
// [[Rcpp::export]]
Rcpp::NumericVector rsn_cpp(int n, double xi, double omega, double alpha, double tau) {
  return Rcpp::wrap(sn::random(n, xi, omega, alpha, tau));
}


// --- Skew-t Family ---

/**
 * @title Skew-t Probability Density Function (C++ backend)
 * @description Internal C++ function for the skew-t PDF. Not intended for direct user call.
 * @param x A numeric vector of quantiles.
 * @param xi Location parameter.
 * @param omega Scale parameter.
 * @param alpha Slant parameter.
 * @param nu Degrees of freedom.
 * @param log_d Logical; if TRUE, probabilities are given as log(p).
 * @return A numeric vector of PDF values.
 * @keywords internal
 * @export
 */
// [[Rcpp::export]]
Rcpp::NumericVector dst_cpp(Rcpp::NumericVector x, double xi, double omega, double alpha, double nu, bool log_d) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = st::pdf(x[i], xi, omega, alpha, nu, log_d);
  }
  return out;
}

/**
 * @title Skew-t Cumulative Distribution Function (C++ backend)
 * @description Internal C++ function for the skew-t CDF. Not intended for direct user call.
 * @param x A numeric vector of quantiles.
 * @param xi Location parameter.
 * @param omega Scale parameter.
 * @param alpha Slant parameter.
 * @param nu Degrees of freedom.
 * @param method Integer specifying the numerical method.
 * @param lower_tail Logical; if TRUE, probabilities are P(X <= x), otherwise P(X > x).
 * @param log_p Logical; if TRUE, probabilities p are given as log(p).
 * @return A numeric vector of CDF values.
 * @keywords internal
 * @export
 */
// [[Rcpp::export]]
Rcpp::NumericVector pst_cpp(Rcpp::NumericVector x, double xi, double omega, double alpha, double nu,
                            int method, bool lower_tail, bool log_p) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = st::cdf(x[i], xi, omega, alpha, nu, method, lower_tail, log_p);
  }
  return out;
}

/**
 * @title Skew-t Quantile Function (C++ backend)
 * @description Internal C++ function for the skew-t quantile. Not intended for direct user call.
 * @param p A numeric vector of probabilities.
 * @param xi Location parameter.
 * @param omega Scale parameter.
 * @param alpha Slant parameter.
 * @param nu Degrees of freedom.
 * @param tol Tolerance for numerical approximation.
 * @param method Integer specifying the numerical method.
 * @return A numeric vector of quantile values.
 * @keywords internal
 * @export
 */
// [[Rcpp::export]]
Rcpp::NumericVector qst_cpp(Rcpp::NumericVector p, double xi, double omega, double alpha, double nu,
                            double tol, int method) {
  int n = p.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = st::quantile(p[i], xi, omega, alpha, nu, tol, method);
  }
  return out;
}

/**
 * @title Skew-t Random Variate Generation (C++ backend)
 * @description Internal C++ function for generating random variates from the skew-t distribution. Not intended for direct user call.
 * @param n Number of observations.
 * @param xi Location parameter.
 * @param omega Scale parameter.
 * @param alpha Slant parameter.
 * @param nu Degrees of freedom.
 * @return A numeric vector of random variates.
 * @keywords internal
 * @export
 */
// [[Rcpp::export]]
Rcpp::NumericVector rst_cpp(int n, double xi, double omega, double alpha, double nu) {
  return Rcpp::wrap(st::random(n, xi, omega, alpha, nu));
}

/**
 * @title Owen's T function (C++ backend)
 * @description Internal C++ function for Owen's T function. Not intended for direct user call.
 * @param h A numeric vector.
 * @param a A numeric vector.
 * @return A numeric vector containing the values of T(h, a).
 * @keywords internal
 * @export
 */
// [[Rcpp::export]]
Rcpp::NumericVector owens_T(Rcpp::NumericVector h, Rcpp::NumericVector a) {
  int n = std::max(h.size(), a.size());
  Rcpp::NumericVector out(n);

  // Vectorize inputs
  Rcpp::NumericVector h_vec = Rcpp::rep_len(h, n);
  Rcpp::NumericVector a_vec = Rcpp::rep_len(a, n);

  for (int i = 0; i < n; ++i) {
    // Call the helper function from the utils namespace
    out[i] = utils::owens_t(h_vec[i], a_vec[i]);
  }
  return out;
}
