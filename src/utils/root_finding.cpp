#include "root_finding.hpp" // Include the header for the function we are implementing.
#include <Rcpp.h> 
#include <cmath>     // For std::abs
#include <utility>   // For std::swap
#include <stdexcept> // For std::runtime_error (optional, for internal checks)
#include <Rmath.h>   // For R_NaN

namespace utils {

double find_root_alternating(
    std::function<double(double)> f,
    double lower_bound,
    double upper_bound,
    double tol,
    int max_iter
) {
    // --- Initial Sanity Checks ---
    if (R_IsNA(lower_bound) || R_IsNA(upper_bound)) {
        return R_NaN;
    }

    // Ensure lower_bound < upper_bound
    if (lower_bound > upper_bound) {
        std::swap(lower_bound, upper_bound);
    }

    double fa = f(lower_bound);
    double fb = f(upper_bound);

    // Check if the root is bracketed. This is a precondition for the algorithm.
    if (fa * fb > 0) {
        // Rcpp::warning("Root not bracketed in find_root_alternating. f(a)=%f, f(b)=%f", fa, fb);
        return R_NaN;
    }

    // Check if one of the bounds is already the root
    if (std::abs(fa) < tol) {
        return lower_bound;
    }
    if (std::abs(fb) < tol) {
        return upper_bound;
    }

    // --- Main Iteration Loop ---
    bool use_regula_falsi = false; // Start with the safer bisection method
    double xc = lower_bound;      // Initialize the current guess
    double fc = fa;

    for (int i = 0; i < max_iter; ++i) {
        // 1. Calculate the next guess for the root
        if (use_regula_falsi) {
            // Regula Falsi (False Position) method
            // Check for division by zero before calculating
            if (std::abs(fb - fa) < 1e-15) { // Denominator is too small, fallback to bisection
                 xc = (lower_bound + upper_bound) / 2.0;
            } else {
                 xc = upper_bound - fb * (upper_bound - lower_bound) / (fb - fa);
            }

            // A key problem with Regula Falsi is that it can get "stuck" if one
            // endpoint doesn't move. If the new guess `xc` is outside the
            // current bounds, it's a sign of numerical instability.
            // In such a case, we should fallback to a bisection step.
            if (xc < lower_bound || xc > upper_bound) {
                 xc = (lower_bound + upper_bound) / 2.0;
            }

        } else {
            // Bisection method
            xc = (lower_bound + upper_bound) / 2.0;
        }

        fc = f(xc);

        // 2. Check for convergence
        if (std::abs(fc) < tol) {
            return xc;
        }

        // 3. Update the bracket for the next iteration
        if (fc * fa < 0) {
            // The root is in the lower sub-interval [lower_bound, xc]
            upper_bound = xc;
            fb = fc;
        } else {
            // The root is in the upper sub-interval [xc, upper_bound]
            lower_bound = xc;
            fa = fc;
        }

        // 4. Toggle the method for the next iteration
        use_regula_falsi = !use_regula_falsi;
    }

    // If the loop finishes, the method failed to converge
    // Rcpp::warning("find_root_alternating failed to converge in %d iterations.", max_iter);
    return R_NaN;
}

} // namespace utils
