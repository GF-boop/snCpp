#ifndef UTILS_ROOT_FINDING_HPP
#define UTILS_ROOT_FINDING_HPP

#include <functional> // For std::function

namespace utils {

/**
 * @brief Finds the root of a function using an alternating Bisection and Regula Falsi method.
 *
 * This function searches for a value `x` such that `f(x) = 0` within a given
 * interval `[lower_bound, upper_bound]`. The root must be bracketed by the
 * initial bounds (i.e., f(lower_bound) and f(upper_bound) must have opposite signs).
 *
 * The algorithm alternates between the safe but slow Bisection method and the
 * faster but potentially unreliable Regula Falsi (False Position) method to
 * ensure both speed and convergence.
 *
 * @param f The function for which to find a root. It takes a double and returns a double.
 * @param lower_bound The lower boundary of the search interval.
 * @param upper_bound The upper boundary of the search interval.
 * @param tol The desired tolerance for convergence. The algorithm stops when |f(x)| < tol.
 * @param max_iter The maximum number of iterations to perform before giving up.
 *
 * @return The estimated root `x`. Returns `R_NaN` if the root is not bracketed,
 * if the algorithm fails to converge within `max_iter`, or if a numerical error occurs.
 */
double find_root_alternating(
    std::function<double(double)> f,
    double lower_bound,
    double upper_bound,
    double tol = 1e-9,
    int max_iter = 100
);

} // namespace utils

#endif // UTILS_ROOT_FINDING_HPP