
// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include "distributions.hpp"
#include "special_functions.hpp"

// --- Skew-Normal Family ---

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dsn_cpp(Rcpp::NumericVector x, double xi = 0, double omega = 1, double alpha = 0, bool log_d = false) {
    int n = x.size();
    Rcpp::NumericVector out(n);
    for(int i = 0; i < n; ++i) {
        out[i] = sn::pdf(x[i], xi, omega, alpha, log_d);
    }
    return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector psn_cpp(Rcpp::NumericVector x, double xi = 0, double omega = 1, double alpha = 0, bool lower_tail = true, bool log_p = false) {
    int n = x.size();
    Rcpp::NumericVector out(n);
    for(int i = 0; i < n; ++i) {
        out[i] = sn::cdf(x[i], xi, omega, alpha, lower_tail, log_p);
    }
    return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qsn_cpp(Rcpp::NumericVector p, double xi = 0, double omega = 1, double alpha = 0, double tau = 0, double tol = 1e-8) {
    int n = p.size();
    Rcpp::NumericVector out(n);
    for(int i = 0; i < n; ++i) {
        out[i] = sn::quantile(p[i], xi, omega, alpha, tau, tol);
    }
    return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rsn_cpp(int n, double xi = 0, double omega = 1, double alpha = 0, double tau = 0) {
    return Rcpp::wrap(sn::random(n, xi, omega, alpha, tau));
}


// --- Skew-t Family ---

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dst_cpp(Rcpp::NumericVector x, double xi = 0, double omega = 1, double alpha = 0, double nu = R_PosInf, bool log_d = false) {
    int n = x.size();
    Rcpp::NumericVector out(n);
    for(int i = 0; i < n; ++i) {
        out[i] = st::pdf(x[i], xi, omega, alpha, nu, log_d);
    }
    return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pst_cpp(Rcpp::NumericVector x, double xi = 0, double omega = 1, double alpha = 0, double nu = R_PosInf,
                            int method = 0, bool lower_tail = true, bool log_p = false) {
    int n = x.size();
    Rcpp::NumericVector out(n);
    for(int i = 0; i < n; ++i) {
        out[i] = st::cdf(x[i], xi, omega, alpha, nu, method, lower_tail, log_p);
    }
    return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qst_cpp(Rcpp::NumericVector p, double xi = 0, double omega = 1, double alpha = 0, double nu = R_PosInf,
                            double tol = 1e-8, int method = 0) {
    int n = p.size();
    Rcpp::NumericVector out(n);
    for(int i = 0; i < n; ++i) {
        out[i] = st::quantile(p[i], xi, omega, alpha, nu, tol, method);
    }
    return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rst_cpp(int n, double xi = 0, double omega = 1, double alpha = 0, double nu = R_PosInf) {
    return Rcpp::wrap(st::random(n, xi, omega, alpha, nu));
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector owens_T(Rcpp::NumericVector h, Rcpp::NumericVector a) {
    int n = std::max(h.size(), a.size());
    Rcpp::NumericVector out(n);

    // Vectorize inputs
    Rcpp::NumericVector h_vec = Rcpp::rep_len(h, n);
    Rcpp::NumericVector a_vec = Rcpp::rep_len(a, n);

    for (int i = 0; i < n; ++i) {
        // Call the helper function from the utils namespace
        out[i] = utils::owens_t(h_vec[i], a_vec[i]);
    }
    return out;
}
