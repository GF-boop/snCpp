#' @useDynLib snCpp, .registration = TRUE
#' @importFrom Rcpp evalCpp

# --- Skew-Normal Family ---

#' Skew-Normal Probability Density Function
#'
#' Computes the probability density function of the skew-normal distribution.
#' This is an R wrapper for the C++ implementation.
#'
#' @param xi Location parameter (numeric, default: 0).
#' @param omega Scale parameter (numeric, default: 1). Must be positive.
#' @param alpha Skewness parameter (numeric, default: 0).
#' @param log_d Logical; if TRUE, probabilities p are given as log(p).
#' @return A numeric vector of PDF values.
#' @export
#' @examples
#' dsn(0)
#' dsn(1, alpha = 5)
dsn <- function(x, xi = 0, omega = 1, alpha = 0, log_d = FALSE) {
  if (!is.numeric(x)) stop("'x' must be numeric.")
  if (!is.numeric(xi) || length(xi) != 1) stop("'xi' must be a single numeric value.")
  if (!is.numeric(omega) || length(omega) != 1 || omega <= 0) stop("'omega' must be a single positive numeric value.")
  if (!is.numeric(alpha) || length(alpha) != 1) stop("'alpha' must be a single numeric value.")
  if (!is.logical(log_d) || length(log_d) != 1) stop("'log_d' must be a single logical value.")
  dsn_cpp(x, xi, omega, alpha, log_d)
}

#' Skew-Normal Cumulative Distribution Function
#'
#' Computes the cumulative distribution function of the skew-normal distribution.
#' This is an R wrapper for the C++ implementation.
#'
#' @param x A numeric vector of quantiles.
#' @param xi Location parameter (numeric, default: 0).
#' @param omega Scale parameter (numeric, default: 1). Must be positive.
#' @param alpha Skewness parameter (numeric, default: 0).
#' @param lower_tail Logical; if TRUE (default), probabilities are P(X <= x), otherwise P(X > x).
#' @param log_p Logical; if TRUE, probabilities p are given as log(p).
#' @return A numeric vector of CDF values.
#' @export
#' @examples
#' psn(0)
#' psn(1, alpha = 5)
psn <- function(x, xi = 0, omega = 1, alpha = 0, lower_tail = TRUE, log_p = FALSE) {
  if (!is.numeric(x)) stop("'x' must be numeric.")
  if (!is.numeric(xi) || length(xi) != 1) stop("'xi' must be a single numeric value.")
  if (!is.numeric(omega) || length(omega) != 1 || omega <= 0) stop("'omega' must be a single positive numeric value.")
  if (!is.numeric(alpha) || length(alpha) != 1) stop("'alpha' must be a single numeric value.")
  if (!is.logical(lower_tail) || length(lower_tail) != 1) stop("'lower_tail' must be a single logical value.")
  if (!is.logical(log_p) || length(log_p) != 1) stop("'log_p' must be a single logical value.")
  psn_cpp(x, xi, omega, alpha, lower_tail, log_p)
}

#' Skew-Normal Quantile Function
#'
#' Computes the quantile function of the skew-normal distribution.
#' This is an R wrapper for the C++ implementation.
#'
#' @param p A numeric vector of probabilities.
#' @param xi Location parameter (numeric, default: 0).
#' @param omega Scale parameter (numeric, default: 1). Must be positive.
#' @param alpha Skewness parameter (numeric, default: 0).
#' @param tau A parameter for generalized skew-normal (default: 0 for standard skew-normal).
#' @param tol Tolerance for numerical approximation (default: 1e-8).
#' @return A numeric vector of quantile values.
#' @export
#' @examples
#' qsn(0.5)
#' qsn(0.9, alpha = -3)
qsn <- function(p, xi = 0, omega = 1, alpha = 0, tau = 0, tol = 1e-8) {
  if (!is.numeric(p) || any(p < 0 | p > 1)) stop("'p' must be numeric and between 0 and 1.")
  if (!is.numeric(xi) || length(xi) != 1) stop("'xi' must be a single numeric value.")
  if (!is.numeric(omega) || length(omega) != 1 || omega <= 0) stop("'omega' must be a single positive numeric value.")
  if (!is.numeric(alpha) || length(alpha) != 1) stop("'alpha' must be a single numeric value.")
  if (!is.numeric(tau) || length(tau) != 1) stop("'tau' must be a single numeric value.")
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) stop("'tol' must be a single positive numeric value.")
  qsn_cpp(p, xi, omega, alpha, tau, tol)
}

#' Skew-Normal Random Variate Generation
#'
#' Generates random deviates from the skew-normal distribution.
#' This is an R wrapper for the C++ implementation.
#'
#' @param n Number of observations.
#' @param xi Location parameter (numeric, default: 0).
#' @param omega Scale parameter (numeric, default: 1). Must be positive.
#' @param alpha Skewness parameter (numeric, default: 0).
#' @param tau A parameter for generalized skew-normal (default: 0 for standard skew-normal).
#' @return A numeric vector of random variates.
#' @export
#' @examples
#' rsn(10)
#' rsn(5, alpha = 10)
rsn <- function(n, xi = 0, omega = 1, alpha = 0, tau = 0) {
  if (!is.numeric(n) || length(n) != 1 || n < 1 || floor(n) != n) stop("'n' must be a single positive integer.")
  if (!is.numeric(xi) || length(xi) != 1) stop("'xi' must be a single numeric value.")
  if (!is.numeric(omega) || length(omega) != 1 || omega <= 0) stop("'omega' must be a single positive numeric value.")
  if (!is.numeric(alpha) || length(alpha) != 1) stop("'alpha' must be a single numeric value.")
  if (!is.numeric(tau) || length(tau) != 1) stop("'tau' must be a single numeric value.")
  rsn_cpp(n, xi, omega, alpha, tau)
}


# --- Skew-t Family ---

#' Skew-t Probability Density Function
#'
#' Computes the probability density function of the skew-t distribution.
#' This is an R wrapper for the C++ implementation.
#'
#' @param x A numeric vector of quantiles.
#' @param xi Location parameter (numeric, default: 0).
#' @param omega Scale parameter (numeric, default: 1). Must be positive.
#' @param alpha Skewness parameter (numeric, default: 0).
#' @param nu Degrees of freedom (numeric, default: Inf for skew-normal limit). Must be positive.
#' @param log_d Logical; if TRUE, probabilities p are given as log(p).
#' @return A numeric vector of PDF values.
#' @export
#' @examples
#' dst(0, nu = 4)
#' dst(1, alpha = 5, nu = 10)
dst <- function(x, xi = 0, omega = 1, alpha = 0, nu = Inf, log_d = FALSE) {
  if (!is.numeric(x)) stop("'x' must be numeric.")
  if (!is.numeric(xi) || length(xi) != 1) stop("'xi' must be a single numeric value.")
  if (!is.numeric(omega) || length(omega) != 1 || omega <= 0) stop("'omega' must be a single positive numeric value.")
  if (!is.numeric(alpha) || length(alpha) != 1) stop("'alpha' must be a single numeric value.")
  if (!is.numeric(nu) || length(nu) != 1 || nu <= 0) stop("'nu' must be a single positive numeric value (or Inf).")
  if (!is.logical(log_d) || length(log_d) != 1) stop("'log_d' must be a single logical value.")
  dst_cpp(x, xi, omega, alpha, nu, log_d)
}

#' Skew-t Cumulative Distribution Function
#'
#' Computes the cumulative distribution function of the skew-t distribution.
#' This is an R wrapper for the C++ implementation.
#'
#' @param x A numeric vector of quantiles.
#' @param xi Location parameter (numeric, default: 0).
#' @param omega Scale parameter (numeric, default: 1). Must be positive.
#' @param alpha Skewness parameter (numeric, default: 0).
#' @param nu Degrees of freedom (numeric, default: Inf for skew-normal limit). Must be positive.
#' @param method Integer; method for CDF calculation (default: 0).
#' @param lower_tail Logical; if TRUE (default), probabilities are P(X <= x), otherwise P(X > x).
#' @param log_p Logical; if TRUE, probabilities p are given as log(p).
#' @return A numeric vector of CDF values.
#' @export
#' @examples
#' pst(0, nu = 4)
#' pst(1, alpha = 5, nu = 10)
pst <- function(x, xi = 0, omega = 1, alpha = 0, nu = Inf, method = 0L, lower_tail = TRUE, log_p = FALSE) {
  if (!is.numeric(x)) stop("'x' must be numeric.")
  if (!is.numeric(xi) || length(xi) != 1) stop("'xi' must be a single numeric value.")
  if (!is.numeric(omega) || length(omega) != 1 || omega <= 0) stop("'omega' must be a single positive numeric value.")
  if (!is.numeric(alpha) || length(alpha) != 1) stop("'alpha' must be a single numeric value.")
  if (!is.numeric(nu) || length(nu) != 1 || nu <= 0) stop("'nu' must be a single positive numeric value (or Inf).")
  if (!is.numeric(method) || length(method) != 1) stop("'method' must be a single integer.")
  if (!is.logical(lower_tail) || length(lower_tail) != 1) stop("'lower_tail' must be a single logical value.")
  if (!is.logical(log_p) || length(log_p) != 1) stop("'log_p' must be a single logical value.")
  pst_cpp(x, xi, omega, alpha, nu, method, lower_tail, log_p)
}

#' Skew-t Quantile Function
#'
#' Computes the quantile function of the skew-t distribution.
#' This is an R wrapper for the C++ implementation.
#'
#' @param p A numeric vector of probabilities.
#' @param xi Location parameter (numeric, default: 0).
#' @param omega Scale parameter (numeric, default: 1). Must be positive.
#' @param alpha Skewness parameter (numeric, default: 0).
#' @param nu Degrees of freedom (numeric, default: Inf for skew-normal limit). Must be positive.
#' @param tol Tolerance for numerical approximation (default: 1e-8).
#' @param method Integer; method for quantile calculation (default: 0).
#' @return A numeric vector of quantile values.
#' @export
#' @examples
#' qst(0.5, nu = 4)
#' qst(0.9, alpha = -3, nu = 10)
qst <- function(p, xi = 0, omega = 1, alpha = 0, nu = Inf, tol = 1e-8, method = 0L) {
  if (!is.numeric(p) || any(p < 0 | p > 1)) stop("'p' must be numeric and between 0 and 1.")
  if (!is.numeric(xi) || length(xi) != 1) stop("'xi' must be a single numeric value.")
  if (!is.numeric(omega) || length(omega) != 1 || omega <= 0) stop("'omega' must be a single positive numeric value.")
  if (!is.numeric(alpha) || length(alpha) != 1) stop("'alpha' must be a single numeric value.")
  if (!is.numeric(nu) || length(nu) != 1 || nu <= 0) stop("'nu' must be a single positive numeric value (or Inf).")
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) stop("'tol' must be a single positive numeric value.")
  if (!is.numeric(method) || length(method) != 1) stop("'method' must be a single integer.")
  qst_cpp(p, xi, omega, alpha, nu, tol, method)
}

#' Skew-t Random Variate Generation
#'
#' Generates random deviates from the skew-t distribution.
#' This is an R wrapper for the C++ implementation.
#'
#' @param n Number of observations.
#' @param xi Location parameter (numeric, default: 0).
#' @param omega Scale parameter (numeric, default: 1). Must be positive.
#' @param alpha Skewness parameter (numeric, default: 0).
#' @param nu Degrees of freedom (numeric, default: Inf for skew-normal limit). Must be positive.
#' @return A numeric vector of random variates.
#' @export
#' @examples
#' rst(10, nu = 4)
#' rst(5, alpha = 10, nu = 20)
rst <- function(n, xi = 0, omega = 1, alpha = 0, nu = Inf) {
  if (!is.numeric(n) || length(n) != 1 || n < 1 || floor(n) != n) stop("'n' must be a single positive integer.")
  if (!is.numeric(xi) || length(xi) != 1) stop("'xi' must be a single numeric value.")
  if (!is.numeric(omega) || length(omega) != 1 || omega <= 0) stop("'omega' must be a single positive numeric value.")
  if (!is.numeric(alpha) || length(alpha) != 1) stop("'alpha' must be a single numeric value.")
  if (!is.numeric(nu) || length(nu) != 1 || nu <= 0) stop("'nu' must be a single positive numeric value (or Inf).")
  rst_cpp(n, xi, omega, alpha, nu)
}

# --- Special Functions ---

#' Owen's T function
#'
#' Computes Owen's T function, T(h, a). This is a user-facing
#' wrapper for the C++ implementation.
#'
#' @param h A numeric vector.
#' @param a A numeric vector.
#' @return A numeric vector containing the values of T(h, a).
#' @export
#' @examples
#' owens.T(1, 1)
#' owens.T(c(1, 2, 3), a = 0.5)
owens.T <- function(h, a) {
  if (!is.numeric(h)) stop("'h' must be numeric.")
  if (!is.numeric(a)) stop("'a' must be numeric.")

  # CORRECT: Calling the C++ function with an underscore
  owens_T(h, a)
}
