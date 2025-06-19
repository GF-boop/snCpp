#include "special_functions.hpp" // For the declaration of utils::owens_t

#include <Rcpp.h> // For R::pnorm, M_2PI
#include <cmath>  // For std::abs, std::atan, std::exp, std::pow

namespace { // Start of unnamed namespace for private helpers

  // This is a file-local helper function. It is invisible to other .cpp files,
  // preventing any linker errors. It contains the core computational logic.
  // It is a direct C++ port of the R sn:::T.int function's algorithm.
  // It expects non-negative h and a, with a <= 1.
  double owens_t_engine(double h, double a) {
    
    // Define M_2PI if not already available
    #ifndef M_2PI
    #define M_2PI 6.28318530717958647692
    #endif

    if (a == 0.0) return 0.0;
    if (h == 0.0) return std::atan(a) / M_2PI;
    
    // Approximation for large h, from sn:::T.int
    if (h > 8.0) { // cut.point = 8 in R version
        return std::atan(a) * std::exp(-0.5 * h * h * a / std::atan(a)) * (1.0 + 0.00868 * std::pow(h * a, 4)) / M_2PI;
    }

    // Series expansion for small h, ported from R's sn:::T.int
    double h_sq = h * h;
    double sum_fui = 1.0;
    double fui = 1.0;
    double sum_series = 0.0;
    
    for (int i = 0; i <= 50; ++i) { // jmax = 50 in R version
        if (i > 0) {
            fui = fui * (h_sq / 2.0) / i; // calculates (h^2/2)^i / i!
            sum_fui += fui;
        }
        
        double b1 = std::exp(-0.5 * h_sq) * sum_fui;
        double term_i = (1.0 - b1) * std::pow(a, 2.0 * i + 1.0) / (2.0 * i + 1.0);
        
        if (i % 2 == 1) { // if i is odd, sign is negative
            sum_series -= term_i;
        } else { // if i is even, sign is positive
            sum_series += term_i;
        }
    }
    return (std::atan(a) - sum_series) / M_2PI;
  }

} // End of unnamed namespace

// --- Public Function Definition ---

namespace utils {
  // This is the implementation of the function declared in the header.
  // It's the only function from this file that other C++ files can see.
  double owens_t(double h, double a) {
    if (R_IsNA(h) || R_IsNA(a)) {
        return R_NaN;
    }
    
    if (!R_finite(a)) {
        if (a == R_PosInf) {
            return 0.5 * R::pnorm(-std::abs(h), 0.0, 1.0, 1, 0);
        }
        if (a == R_NegInf) {
            return -0.5 * R::pnorm(-std::abs(h), 0.0, 1.0, 1, 0);
        }
        // If a is NaN, it will fall through and be handled by the NA check on h
    }
    
    
    // Use symmetry properties to simplify the problem for the engine
    double sign_a = (a > 0) - (a < 0);
    double abs_a = std::abs(a);
    double abs_h = std::abs(h);

    if (abs_a <= 1.0) {
        return sign_a * owens_t_engine(abs_h, abs_a);
    } 
    
    // else, for |a| > 1, use the stable recursive identity
    double p_h = R::pnorm(abs_h, 0.0, 1.0, 1, 0);
    double p_ah = R::pnorm(abs_a * abs_h, 0.0, 1.0, 1, 0);
    double t_rec = owens_t_engine(abs_a * abs_h, 1.0 / abs_a);
    
    return sign_a * (0.5 * p_h + 0.5 * p_ah - p_h * p_ah - t_rec);
  }

} // namespace utils
