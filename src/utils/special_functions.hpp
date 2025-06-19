#ifndef UTILS_SPECIAL_FUNCTIONS_HPP
#define UTILS_SPECIAL_FUNCTIONS_HPP

// This header file declares the functions available in special_functions.cpp

namespace utils {

    /**
     * @brief Computes Owen's T function, T(h, a).
     *
     * This function provides a robust C++ implementation equivalent to the one
     * found in the R 'sn' package. It handles symmetry and uses a recursive
     * identity for numerical stability.
     *
     * @param h A numeric scalar for the 'h' parameter.
     * @param a A numeric scalar for the 'a' parameter.
     * @return The value of Owen's T function as a double.
     */
    double owens_t(double h, double a);

} // namespace utils

#endif // UTILS_SPECIAL_FUNCTIONS_HPP
