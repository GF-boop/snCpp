// In file: core/skew_normal.hpp

#ifndef CORE_SKEW_NORMAL_HPP
#define CORE_SKEW_NORMAL_HPP

#include <vector>

/**
 * @brief A namespace for core mathematical functions of the Skew-Normal distribution.
 */
namespace sn {

/**
 * @brief The probability density function (PDF) for the classic Skew-Normal distribution.
 */
double pdf(double x, double xi = 0, double omega = 1, double alpha = 0, bool log_d = false);

/**
 * @brief The cumulative distribution function (CDF) for the classic Skew-Normal distribution.
 * This is a correct implementation using Owen's T function.
 */
double cdf(double x, double xi = 0, double omega = 1, double alpha = 0, bool lower_tail = true, bool log_p = false);

/**
 * @brief The quantile function (inverse CDF) for the Skew-Normal distribution.
 * @note This implementation handles the classic (`tau=0`) and extended (`tau!=0`) cases for random number
 * generation, but the quantile function is only implemented for the classic case (`tau=0`).
 */
double quantile(double p, double xi = 0, double omega = 1, double alpha = 0, double tau = 0, double tol = 1e-8);

/**
 * @brief Generate random variates from the Skew-Normal distribution (classic or extended).
 * @param tau The extension parameter. `tau=0` gives the classic SN distribution.
 */
std::vector<double> random(int n, double xi = 0, double omega = 1, double alpha = 0, double tau = 0);

} // namespace sn

#endif // CORE_SKEW_NORMAL_HPP
