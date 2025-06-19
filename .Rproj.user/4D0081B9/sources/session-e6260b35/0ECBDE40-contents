#include "skew_normal.hpp"
#include "special_functions.hpp" // For utils::owens_t

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#include <algorithm>

// =============================================================================
// ANONYMOUS NAMESPACE FOR ALL INTERNAL HELPERS
// These functions are private to this file and cannot be called elsewhere.
// =============================================================================
namespace {

// Helper for the PDF of the standard (xi=0, omega=1) SN distribution
double pdf_standard(double z, double alpha) {
  if (!R_finite(z)) return 0.0;
  return 2.0 * R::dnorm(z, 0.0, 1.0, 0) * R::pnorm(alpha * z, 0.0, 1.0, 1, 0);
}

  // Helper for the CDF of the standard (xi=0, omega=1) SN distribution
  double cdf_standard(double z, double alpha) {
    if (z == R_PosInf) return 1.0;
    if (z == R_NegInf) return 0.0;
    if (alpha == 0) return R::pnorm(z, 0.0, 1.0, 1, 0);

    // This is the definitive formula using a correct Owen's T function
    return R::pnorm(z, 0.0, 1.0, 1, 0) - 2.0 * utils::owens_t(z, alpha);
  }

  // Helper for the Cumulants of the standard (xi=0, omega=1) SN distribution
  std::vector<double> get_cumulants_standard(double alpha) {
    const double delta = alpha / std::sqrt(1.0 + alpha * alpha);
    const double b = std::sqrt(2.0 / M_PI);

    const double k1 = b * delta;                // Mean
    const double k2 = 1.0 - k1 * k1;            // Variance
    const double k3 = 0.5 * (4.0 - M_PI) * std::pow(k1, 3.0);
    const double k4 = 2.0 * (M_PI - 3.0) * std::pow(k1, 4.0);

    return {k1, k2, k3, k4};
  }

  // Newton-Raphson solver for the standard SN quantile
  double quantile_nr(double p, double alpha, double tol, bool& converged) {
    std::vector<double> cum = get_cumulants_standard(alpha);
    double k1 = cum[0], k2 = cum[1], k3 = cum[2], k4 = cum[3];

    if (k2 <= 0) { converged = false; return R_NaN; }

    double g1 = k3 / std::pow(k2, 1.5);
    double g2 = k4 / std::pow(k2, 2.0);

    double x = R::qnorm(p, 0.0, 1.0, 1, 0);
    x = x + (x*x - 1.0)*g1/6.0 + x*(x*x - 3.0)*g2/24.0 - x*(2.0*x*x - 5.0)*g1*g1/36.0;
    x = k1 + std::sqrt(k2) * x; // Initial guess is now correct

    for (int i = 0; i < 200; ++i) {
      double error = cdf_standard(x, alpha) - p; // Use correct helper
      if (std::abs(error) < tol) {
        converged = true;
        return x;
      }
      double density = pdf_standard(x, alpha); // Use correct helper
      if (density < 1e-15) break;
      x -= error / density;
    }
    converged = false;
    return x;
  }

  // Bisection solver for the standard SN quantile
  double quantile_bisection(double p, double alpha, double tol, bool& converged) {
    double x_low, x_high;
    // Start with a very wide but safe bracket
    x_low = R::qnorm(p, 0.0, 1.0, 1, 0) - 4.0;
    x_high = R::qnorm(p, 0.0, 1.0, 1, 0) + 4.0;

    auto f = [&](double x) { return cdf_standard(x, alpha) - p; };
    double f_low = f(x_low);
    double f_high = f(x_high);

    for (int i = 0; i < 10 && (f_low * f_high > 0); ++i) {
      x_low -= 2.0; x_high += 2.0;
      f_low = f(x_low); f_high = f(x_high);
    }

    if (f_low * f_high > 0) { converged = false; return R_NaN; }
    if (f_low > 0) std::swap(x_low, x_high);

    for (int i = 0; i < 100; ++i) {
      double x_mid = x_low + 0.5 * (x_high - x_low);
      if (x_high - x_low < tol) { converged = true; return x_mid; }
      if (f(x_mid) < 0) x_low = x_mid; else x_high = x_mid;
    }
    converged = true;
    return (x_low + x_high) / 2.0;
  }

} // End of unnamed namespace for helpers

// =============================================================================
// PUBLIC `sn` NAMESPACE IMPLEMENTATIONS
// =============================================================================

namespace sn {

double pdf(double x, double xi, double omega, double alpha, bool log_d) {
  if (omega <= 0) return R_NaN;
  double z = (x - xi) / omega;
  double log_pdf_val = std::log(2.0) + R::dnorm(z, 0.0, 1.0, true) + R::pnorm(alpha * z, 0.0, 1.0, true, true) - std::log(omega);
  return log_d ? log_pdf_val : std::exp(log_pdf_val);
}

double cdf(double x, double xi, double omega, double alpha, bool lower_tail, bool log_p) {
  if (omega <= 0) return R_NaN;
  double z = (x - xi) / omega;
  double p = cdf_standard(z, alpha); // Call the internal helper
  if (!lower_tail) p = 1.0 - p;
  if (log_p) { p = (p > 0) ? std::log(p) : R_NegInf; }
  return p;
}

std::vector<double> random(int n, double xi, double omega, double alpha, double tau) {
  if (n <= 0) return {};
  if (omega <= 0) return std::vector<double>(n, R_NaN);
  if (tau != 0) Rcpp::stop("sn::random only implemented for tau=0.");

  std::vector<double> y(n);
  Rcpp::RNGScope scope;
  double delta = alpha / std::sqrt(1.0 + alpha * alpha);
  double sqrt_1_minus_delta_sq = std::sqrt(1.0 - delta * delta);

  for (int i = 0; i < n; ++i) {
    double u0 = R::rnorm(0.0, 1.0);
    double v = R::rnorm(0.0, 1.0);
    y[i] = xi + omega * (delta * std::abs(u0) + sqrt_1_minus_delta_sq * v);
  }
  return y;
}

double quantile(double p, double xi, double omega, double alpha, double tau, double tol) {
  if (R_IsNA(p) || p < 0 || p > 1) return R_NaN;
  if (p == 0) return R_NegInf;
  if (p == 1) return R_PosInf;
  if (omega <= 0) return R_NaN;
  if (tau != 0) Rcpp::stop("sn::quantile only implemented for tau=0.");

  if (alpha == 0) return xi + omega * R::qnorm(p, 0.0, 1.0, 1, 0);
  if (alpha == R_PosInf) return xi + omega * std::sqrt(R::qchisq(p, 1.0, 1, 0));
  if (alpha == R_NegInf) return xi - omega * std::sqrt(R::qchisq(1.0 - p, 1.0, 1, 0));

  if (alpha < 0) {
    return xi - sn::quantile(1.0 - p, 0.0, omega, -alpha, 0.0, tol);
  }

  bool converged = false;
  double q_std = quantile_nr(p, alpha, tol, converged);

  // Final check on NR convergence before trying bisection
  if (!converged || std::abs(cdf_standard(q_std, alpha) - p) > tol) {
    q_std = quantile_bisection(p, alpha, tol, converged);
    if (!converged) {
      Rcpp::warning("Quantile solver failed to converge for p=%.4f, alpha=%.2f.", p, alpha);
    }
  }
  return xi + omega * q_std;
}

} // namespace sn
