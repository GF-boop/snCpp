#include "integration.hpp" // Always include the header for the code you are implementing first.

#include <RcppGSL.h>             // For GSL integration with Rcpp
#include <gsl/gsl_integration.h> // GSL's numerical integration header
#include <Rmath.h>               // For R_IsNA, R_NegInf, R_PosInf, R_NaN
#include <Rcpp.h>                // For Rcpp::warning (optional, for debugging)

/**
 * @brief A C-style "trampoline" function to adapt a C++ std::function for GSL.
 *
 * The GSL integration routines are written in C and expect a plain C function
 * pointer of the form `double (* f) (double x, void * params)`.
 *
 * This adapter acts as a bridge. We pass a pointer to our C++ std::function object
 * via the `void * params` argument. This function then casts the void pointer
 * back to its true C++ type and calls it.
 *
 * This is an internal implementation detail and is not exposed in the header.
 */
double gsl_integrand_adapter(double x, void *params) {
    // Cast the void pointer back to a pointer to our std::function object
    auto* f_ptr = static_cast<std::function<double(double)>*>(params);
    // Call the C++ function
    return (*f_ptr)(x);
}


// Implementation of the integrate function declared in the header.
double utils::integrate(std::function<double(double)> f, double lower, double upper, double abs_tol, double rel_tol) {
    // --- Pre-computation checks ---
    if (R_IsNA(lower) || R_IsNA(upper)) {
        return R_NaN;
    }
    if (lower == upper) {
        return 0.0;
    }
    // GSL requires lower < upper, so we reverse the sign if the user provides inverted limits.
    if (lower > upper) {
        return -utils::integrate(f, upper, lower, abs_tol, rel_tol);
    }

    // --- GSL Setup ---
    // GSL workspace for storing intermediate results.
    // The limit of 1000 subintervals is generally sufficient.
    const size_t limit = 1000;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);

    double result, error;

    // The GSL function struct that GSL will operate on.
    gsl_function F;
    F.function = &gsl_integrand_adapter; // Point to our C-style adapter
    F.params = &f;                       // Pass the address of our C++ std::function as a void pointer

    // --- GSL Execution ---
    int status = 1; // Default to an error status

    // Choose the correct GSL routine based on the integration limits
    if (lower == R_NegInf && upper == R_PosInf) {
        // Integrate over (-inf, +inf)
        status = gsl_integration_qagi(&F, abs_tol, rel_tol, limit, w, &result, &error);
    } else if (lower == R_NegInf) {
        // Integrate over (-inf, b)
        status = gsl_integration_qagil(&F, upper, abs_tol, rel_tol, limit, w, &result, &error);
    } else if (upper == R_PosInf) {
        // Integrate over (a, +inf)
        status = gsl_integration_qagiu(&F, lower, abs_tol, rel_tol, limit, w, &result, &error);
    } else {
        // Integrate over a finite interval [a, b]
        status = gsl_integration_qags(&F, lower, upper, abs_tol, rel_tol, limit, w, &result, &error);
    }

    // --- Cleanup and Return ---
    gsl_integration_workspace_free(w);

    if (status) {
        // A non-zero status indicates a GSL error (e.g., failed to converge, max iterations reached).
        // Returning NaN is the standard way to signal a computation failure.
        // For debugging, you can uncomment the following line to see GSL error codes.
        // Rcpp::warning("GSL integration failed with status: %d", status);
        return R_NaN;
    }

    return result;
}
