// In file: core/skew_t.hpp

#ifndef CORE_SKEW_T_HPP
#define CORE_SKEW_T_HPP

#include <vector>

/**
 * @brief A namespace for core mathematical functions of the Skew-t distribution.
 */
namespace st {

/**
 * @brief The probability density function (PDF) for the Skew-t distribution.
 */
double pdf(double x, double xi = 0, double omega = 1, double alpha = 0, double nu = 1e100, bool log_d = false);

/**
 * @brief The cumulative distribution function (CDF) for the Skew-t distribution.
 * @param method An integer specifying the numerical method (0 for auto-selection).
 */
double cdf(double x, double xi = 0, double omega = 1, double alpha = 0, double nu = 1e100,
           int method = 0, bool lower_tail = true, bool log_p = false);

/**
 * @brief The quantile function (inverse CDF) for the Skew-t distribution.
 * @param method The method passed down to the CDF calculation during root finding.
 */
double quantile(double p, double xi = 0, double omega = 1, double alpha = 0, double nu = 1e100,
                double tol = 1e-8, int method = 0);

/**
 * @brief Generate random variates from the Skew-t distribution.
 */
std::vector<double> random(int n, double xi = 0, double omega = 1, double alpha = 0, double nu = 1e100);

} // namespace st

#endif // CORE_SKEW_T_HPP
