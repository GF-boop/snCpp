#ifndef UTILS_INTEGRATION_HPP
#define UTILS_INTEGRATION_HPP

#include <functional> // For std::function, the modern C++ way to handle callable objects

/**
 * @brief A namespace for common utility functions used across the library.
 */
namespace utils {

/**
 * @brief A robust wrapper for GSL numerical integration routines.
 *
 * This function takes a C++ standard function object and integrates it over
 * a specified interval. It automatically handles finite, semi-infinite, and
 * fully infinite integration ranges by selecting the appropriate GSL algorithm.
 *
 * @param f The integrand, a std::function that takes a double and returns a double.
 * @param lower The lower limit of integration. Can be R_NegInf.
 * @param upper The upper limit of integration. Can be R_PosInf.
 * @param abs_tol The absolute error tolerance for the integration result.
 * @param rel_tol The relative error tolerance for the integration result.
 *
 * @return The numerical result of the integration. Returns R_NaN if the
 * integration fails to converge or if an error occurs.
 */
double integrate(
    std::function<double(double)> f,
    double lower,
    double upper,
    double abs_tol = 1e-9,
    double rel_tol = 1e-9
);

} // namespace utils

#endif // UTILS_INTEGRATION_HPP