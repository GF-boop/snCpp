// In file: core/skew_t.cpp

#include "skew_t.hpp"         // Self-header
#include "skew_normal.hpp"    // Dependency for nu = infinity
#include "integration.hpp"   // Dependency for cdf
#include "root_finding.hpp"  // Dependency for quantile
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h> 
#include <Rmath.h>   // For R's math functions (dt, pt, qf, etc.)
#include <cmath>     // For std::abs, std::sqrt, std::log, etc.
#include <algorithm> // For std::min, std::max

// Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// --- Internal Helper Functions ---
// These functions are only used within this file. We place them in an anonymous
// namespace to give them internal linkage, hiding them from other modules.
namespace {

// Corresponds to the original `st_tails_cpp`
double _cdf_tails_logp(double z, double alpha, double nu, bool lower_tail) {
    // This function remains complex, but its logic is unchanged.
    // We only call it for large |z|, and it computes the log-probability.
    if (alpha < 0) return _cdf_tails_logp(-z, -alpha, nu, !lower_tail);

    double lp; // log-probability
    if (z > 0) {
        double log_c2 = std::log(2.0) + R::lgammafn((nu + 1.0) / 2.0) + (nu / 2.0) * std::log(nu) +
                        R::pt(alpha * std::sqrt(nu + 1.0), nu + 1.0, 1, 1) -
                        R::lgammafn(nu / 2.0) - 0.5 * std::log(M_PI);
        if (!lower_tail) { // P(Z > z)
            lp = log_c2 - std::log(nu) - nu * std::log(z);
        } else { // P(Z < z) = 1 - P(Z > z)
            double upper_tail_log_prob = log_c2 - std::log(nu) - nu * std::log(z);
            lp = std::log(1.0 - std::exp(upper_tail_log_prob));
        }
    } else { // z < 0
        double log_c1 = std::log(2.0) + R::lgammafn((nu + 1.0) / 2.0) + (nu / 2.0) * std::log(nu) +
                        R::pt(-alpha * std::sqrt(nu + 1.0), nu + 1.0, 1, 1) -
                        R::lgammafn(nu / 2.0) - 0.5 * std::log(M_PI);
        if (lower_tail) { // P(Z < z)
            lp = log_c1 - std::log(nu) - nu * std::log(-z);
        } else { // P(Z > z) = 1 - P(Z < z)
            double lower_tail_log_prob = log_c1 - std::log(nu) - nu * std::log(-z);
            lp = std::log(1.0 - std::exp(lower_tail_log_prob));
        }
    }
    return lp;
}

// Corresponds to the original `pst_int_cpp`
double _cdf_integer_nu_recursive(double z, double alpha, double nu) {
    if (nu == 1) {
        return std::atan(z) / M_PI + std::acos(alpha / std::sqrt((1.0 + alpha * alpha) * (1.0 + z * z))) / M_PI;
    }
    if (nu == 2) {
        return 0.5 - std::atan(alpha) / M_PI + (0.5 + std::atan(z * alpha / std::sqrt(2.0 + z * z)) / M_PI) * z / std::sqrt(2.0 + z * z);
    }
    // The recursive part
    double term1 = _cdf_integer_nu_recursive(std::sqrt((nu - 2.0) / nu) * z, alpha, nu - 2.0);
    double term2_factor = z *
        std::exp(R::lgammafn((nu - 1.0) / 2.0) + (nu / 2.0 - 1.0) * std::log(nu) - 0.5 * std::log(M_PI) - R::lgammafn(nu / 2.0) - 0.5 * (nu - 1.0) * std::log(nu + z * z));
    double term2_pt = R::pt(std::sqrt(nu - 1.0) * alpha * z / std::sqrt(nu + z * z), nu - 1.0, 1, 0);

    return term1 + term2_pt * term2_factor;
}

} // end anonymous namespace


// --- Public Function Implementations ---
//' @export
// [[Rcpp::export]]
double st::pdf(double x, double xi, double omega, double alpha, double nu, bool log_d) {
    if (R_IsNA(x) || R_IsNA(xi) || R_IsNA(omega) || R_IsNA(alpha) || R_IsNA(nu)) return R_NaN;
    if (omega <= 0 || nu <= 0) return R_NaN;

    if (nu > 1e4) { // nu -> infinity case
        return sn::pdf(x, xi, omega, alpha, log_d);
    }

    double z = (x - xi) / omega;
    double pdf_dt = R::dt(z, nu, 0);
    double cdf_pt_arg = alpha * z * std::sqrt((nu + 1.0) / (z * z + nu));
    double cdf_pt = R::pt(cdf_pt_arg, nu + 1.0, 1, 0);

    double pdf_val = 2.0 * pdf_dt * cdf_pt / omega;

    return log_d ? std::log(pdf_val) : pdf_val;
}

//' @export
// [[Rcpp::export]]
std::vector<double> st::random(int n, double xi, double omega, double alpha, double nu) {
    if (n <= 0) return {};
    if (omega <= 0 || nu <= 0) return std::vector<double>(n, R_NaN);

    if (nu > 1e4) { // nu -> infinity case
        return sn::random(n, xi, omega, alpha, 0.0);
    }

    // Standard generation method
    std::vector<double> y(n);
    std::vector<double> z_std_val = sn::random(n, 0, 1, alpha, 0.0);
    
    GetRNGstate();
    for (int i = 0; i < n; ++i) {
        double v = R::rchisq(nu) / nu;
        y[i] = xi + omega * z_std_val[i] / std::sqrt(v);
    }
    PutRNGstate();

    return y;
}

//' @export
// [[Rcpp::export]]
double st::cdf(double x, double xi, double omega, double alpha, double nu, int method, bool lower_tail, bool log_p) {
    if (R_IsNA(x) || R_IsNA(xi) || R_IsNA(omega) || R_IsNA(alpha) || R_IsNA(nu)) return R_NaN;
    if (omega <= 0 || nu <= 0) return R_NaN;

    // --- Handle edge cases first ---
    if (x == R_PosInf) return log_p ? 0.0 : 1.0;
    if (x == R_NegInf) return log_p ? R_NegInf : 0.0;

    if (nu > 1e4) {
        return sn::cdf(x, xi, omega, alpha, lower_tail, log_p); // Pass all original parameters
    }
    if (alpha == 0) return R::pt((x - xi) / omega, nu, lower_tail, log_p);

    double z = (x - xi) / omega;
    double p_val;

    if (alpha == R_PosInf) {
        double z0 = (z < 0) ? 0 : z;
        p_val = R::pf(z0 * z0, 1.0, nu, 1, 0);
    } else if (alpha == R_NegInf) {
        double z0 = (z > 0) ? 0 : z;
        p_val = 1.0 - R::pf(z0 * z0, 1.0, nu, 1, 0);
    } else {
        // --- Method dispatcher ---
        bool use_int_nu_method = (nu == std::round(nu)) && (method == 4 || (method == 0 && nu <= 15));
        bool use_tails_method = method == 5 || (method == 0 && std::abs(z) > (30.0 + 1.0 / std::sqrt(nu)));
        
        if (use_int_nu_method) {
            p_val = _cdf_integer_nu_recursive(z, alpha, nu);
        } else if (use_tails_method) {
            double log_p_val = _cdf_tails_logp(z, alpha, nu, true); // This helper only does lower tail
            return log_p ? log_p_val : std::exp(log_p_val); // Return early
        } else { // Fallback to integration
            // NOTE: The original had methods 1, 2, 3. We simplify to one robust integration method.
            auto integrand = [&](double val) {
                return st::pdf(val, 0, 1, alpha, nu, false);
            };
            p_val = utils::integrate(integrand, R_NegInf, z);
        }
    }

    // --- Final processing for tail and log scale ---
    if (!lower_tail) p_val = 1.0 - p_val;
    if (log_p) p_val = std::log(p_val);

    return p_val;
}

//' @export
// [[Rcpp::export]]
double st::quantile(double p, double xi, double omega, double alpha, double nu, double tol, int method) {
    if (R_IsNA(p) || p < 0 || p > 1) return R_NaN;
    if (p == 0) return R_NegInf;
    if (p == 1) return R_PosInf;
    if (omega <= 0 || nu <= 0) return R_NaN;
    
    // --- Handle edge cases ---
    if (nu > 1e4) {
        return sn::quantile(p, xi, omega, alpha, 0.0, tol); // Pass tau = 0.0 for classic SN
    }
    if (nu == 1) { /* Rcpp::warning for skew-Cauchy placeholder */ return R_NaN; }

    if (alpha == 0) return xi + omega * R::qt(p, nu, 1, 0);

    // For infinite alpha, the distribution is a half-t
    if (alpha == R_PosInf) return xi + omega * std::sqrt(R::qf(p, 1.0, nu, 1, 0));
    if (alpha == R_NegInf) return xi - omega * std::sqrt(R::qf(1.0 - p, 1.0, nu, 1, 0));

    // For alpha < 0, use the symmetry property: Q(p; alpha) = -Q(1-p; -alpha)
    if (alpha < 0) {
        return xi - st::quantile(1.0 - p, 0, omega, -alpha, nu, tol, method);
    }

    // --- Main root-finding logic for alpha > 0 ---
    
    // 1. Define the objective function f(x) = cdf(x) - p = 0
    auto objective_function = [&](double x) {
        // Here we pass the `method` parameter down to the cdf
        return st::cdf(x, 0, 1, alpha, nu, method, true, false) - p;
    };

    // 2. Find good initial brackets [a, b] for the root
    double lower_bound = R::qt(p, nu, 1, 0) - 2.0; // A reasonable starting guess
    double upper_bound = std::sqrt(R::qf(p, 1.0, nu, 1, 0)) + 2.0; // Another guess

    // Refine the bracket to ensure the root is contained, similar to original code
    for(int i=0; i<10; ++i) { // Limit iterations to prevent infinite loop
        if (objective_function(lower_bound) > 0) { // a is too high
            lower_bound -= (2.0 + i*i);
        }
        if (objective_function(upper_bound) < 0) { // b is too low
            upper_bound += (2.0 + i*i);
        }
    }

    // 3. Call the generic root-finder!
    double q_std = utils::find_root_alternating(objective_function, lower_bound, upper_bound, tol);

    // 4. Rescale and return
    return xi + omega * q_std;
}
