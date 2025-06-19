
# Load necessary libraries
library(sn) # For the reference implementation
library(snCpp) # Your re-implementation
library(ggplot2) # For plotting
library(dplyr) # For data manipulation
library(microbenchmark) # For performance comparison


h_vals <- c(0, 0.5, 1, 8, 10, -0.5, -8, 1e3)
a_vals <- c(0, 0.1, 0.9, 1, 1.1, 10, -0.5, -10, Inf, 4)

for (h in h_vals) {
  for (a in a_vals) {
    r_val <- sn:::T.Owen(h, a)
    cpp_val <- owens.T(h, a)
    cat(sprintf("h=%.2f, a=%.2f, R: %.10f, C++: %.10f, Diff: %.2e\n",
                h, a, r_val, cpp_val, r_val - cpp_val))
  }
}
cat("--- Comparing snCpp (C++ re-implementation) with sn (reference) ---\n")

# --- 2. Define Common Parameters ---
# Choose a set of parameters for comparison
common_xi <- 0.5
common_omega <- 1.5
common_alpha <- 3.0
common_nu <- 8.0 # Degrees of freedom for skew-t

# For skew-normal limit (nu = Inf)
sn_limit_nu <- 10e6

# --- 3. Probability Density Function (d...) ---
cat("\n--- Comparing PDF (d...) ---\n")
x_values <- seq(-5, 5, length.out = 500)

# Skew-t distribution
pdf_sn_st <- sn::dst(x_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu))
pdf_sncpp_st <- snCpp::dst_cpp(x_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu)

cat("Skew-t PDF Comparison (nu = ", common_nu, "):\n")
print(all.equal(pdf_sn_st, pdf_sncpp_st, tolerance = 1e-6)) # Adjust tolerance if needed

# Skew-student limit (nu = Inf)
pdf_sn_sn <- sn::dst(x_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = sn_limit_nu))
pdf_sncpp_sn <- snCpp::dst_cpp(x_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = sn_limit_nu)

cat("Skew-Normal PDF Comparison (nu = Inf):\n")
print(all.equal(pdf_sn_sn, pdf_sncpp_sn, tolerance = 1e-6))

# Plot PDF comparison
df_pdf <- data.frame(
  x = rep(x_values, 2),
  density = c(pdf_sncpp_st, pdf_sn_st),
  source = rep(c("snCpp (Skew-t)", "sn (Skew-t)"), each = length(x_values))
)
ggplot(df_pdf, aes(x = x, y = density, color = source)) +
  geom_line(lwd = 1) +
  ggtitle(paste0("Skew-t PDF Comparison (xi=", common_xi, ", omega=", common_omega, ", alpha=", common_alpha, ", nu=", common_nu, ")")) +
  theme_minimal() +
  ylab("Density") +
  xlab("x") +
  theme(legend.position = "bottom")
# ggsave("skew_t_pdf_comparison.png", width = 8, height = 6) # Uncomment to save plot

# --- 4. Cumulative Distribution Function (p...) ---
cat("\n--- Comparing CDF (p...) ---\n")

# Skew-t distribution
cdf_sn_st <- sn::pst(x_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu))
cdf_sncpp_st <- snCpp::pst_cpp(x_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu)

cat("Skew-t CDF Comparison (nu = ", common_nu, "):\n")
print(all.equal(cdf_sn_st, cdf_sncpp_st, tolerance = 1e-6))

# Skew-normal limit (nu = Inf)
cdf_sn_sn <- sn::pst(x_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = sn_limit_nu))
cdf_sncpp_sn <- snCpp::pst_cpp(x_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = sn_limit_nu)

cat("Skew-Normal CDF Comparison (nu = Inf):\n")
print(all.equal(cdf_sn_sn, cdf_sncpp_sn, tolerance = 1e-6))

# Plot CDF comparison
df_cdf <- data.frame(
  x = rep(x_values, 2),
  cdf = c(cdf_sncpp_sn, cdf_sn_sn),
  source = rep(c("snCpp (Skew-t)", "sn (Skew-t)"), each = length(x_values))
)
ggplot(df_cdf, aes(x = x, y = cdf, color = source)) +
  geom_line(lwd = 1) +
  ggtitle(paste0("Skew-t CDF Comparison (xi=", common_xi, ", omega=", common_omega, ", alpha=", common_alpha, ", nu=", common_nu, ")")) +
  theme_minimal() +
  ylab("CDF") +
  xlab("x") +
  theme(legend.position = "bottom")
# ggsave("skew_t_cdf_comparison.png", width = 8, height = 6) # Uncomment to save plot


# --- 5. Quantile Function (q...) ---
cat("\n--- Comparing Quantile (q...) ---\n")
p_values <- seq(0.01, 0.99, length.out = 100) # Probabilities

# Skew-t distribution
quantile_sn_st <- sn::qst(p_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu))
quantile_sncpp_st <- snCpp::qst_cpp(p_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu)

cat("Skew-t Quantile Comparison (nu = ", common_nu, "):\n")
print(all.equal(quantile_sn_st, quantile_sncpp_st, tolerance = 1e-5)) # Quantiles might need slightly higher tolerance

# Skew-student limit (nu = Inf)
quantile_sn_sn <- sn::qst(p_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = sn_limit_nu))
quantile_sncpp_sn <- snCpp::qst_cpp(p_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = sn_limit_nu)

cat("Skew-Normal Quantile Comparison (nu = Inf):\n")
print(all.equal(quantile_sn_sn, quantile_sncpp_sn, tolerance = 1e-5))

# Plot Quantile comparison
df_quantile <- data.frame(
  p = rep(p_values, 2),
  quantile = c(quantile_sncpp_sn, quantile_sn_sn),
  source = rep(c("snCpp (Skew-t)", "sn (Skew-t)"), each = length(p_values))
)
ggplot(df_quantile, aes(x = p, y = quantile, color = source)) +
  geom_line(lwd = 1) +
  ggtitle(paste0("Skew-t Quantile Comparison (xi=", common_xi, ", omega=", common_omega, ", alpha=", common_alpha, ", nu=", common_nu, ")")) +
  theme_minimal() +
  ylab("Quantile") +
  xlab("Probability") +
  theme(legend.position = "bottom")
# ggsave("skew_t_quantile_comparison.png", width = 8, height = 6) # Uncomment to save plot

# --- 6. Random Variate Generation (r...) ---
cat("\n--- Comparing Random Variate Generation (r...) ---\n")
n_samples <- 100000

# Skew-t distribution
set.seed(123) # For reproducibility
samples_sn_st <- sn::rst(n_samples, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu))
set.seed(123) # Use same seed for snCpp
samples_sncpp_st <- snCpp::rst_cpp(n_samples, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu)

cat("Skew-t Random Variate Summary (nu = ", common_nu, "):\n")
cat("sn::rpst:\n"); print(summary(samples_sn_st))
cat("snCpp::rst:\n"); print(summary(samples_sncpp_st))

# (Note: Random variates won't be identical, but their statistical properties should be similar)
# Plot histograms for visual comparison
df_random <- data.frame(
  value = c(samples_sncpp_st, samples_sn_st),
  source = rep(c("snCpp (Skew-t)", "sn (Skew-t)"), each = n_samples)
)
ggplot(df_random, aes(x = value, fill = source)) +
  geom_histogram(binwidth = 0.2, alpha = 0.6, position = "identity") +
  ggtitle(paste0("Skew-t Random Variate Comparison (xi=", common_xi, ", omega=", common_omega, ", alpha=", common_alpha, ", nu=", common_nu, ")")) +
  theme_minimal() +
  theme(legend.position = "bottom")
# ggsave("skew_t_random_comparison.png", width = 8, height = 6) # Uncomment to save plot


# --- 7. Performance Comparison (Optional) ---
cat("\n--- Performance Comparison ---\n")
cat("Running microbenchmarks for d/p/q/r functions...\n")

mb_d <- microbenchmark(
  sn_dst = sn::dst(x_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu)),
  snCpp_dst = snCpp::dst_cpp(x_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu),
  times = 100
)
cat("\nPDF (d) Benchmark:\n"); print(mb_d)

mb_p <- microbenchmark(
  sn_pst = sn::pst(x_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu)),
  snCpp_pst = snCpp::pst_cpp(x_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu),
  times = 100
)
cat("\nCDF (p) Benchmark:\n"); print(mb_p)

mb_q_inf <- microbenchmark(
  sn_qst = sn::pst(p_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = sn_limit_nu)),
  snCpp_qst = snCpp::pst_cpp(p_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = sn_limit_nu),
  times = 10
)
cat("\nCDF (p), nu = inf, Benchmark:\n"); print(mb_q_inf)

mb_q <- microbenchmark(
  sn_qst = sn::qst(p_values, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu)),
  snCpp_qst = snCpp::qst_cpp(p_values, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu),
  times = 100
)
cat("\nQuantile (q) Benchmark:\n"); print(mb_q)


mb_r <- microbenchmark(
  sn_rst = sn::rst(n_samples, dp = c(xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu)),
  snCpp_rst = snCpp::rst_cpp(n_samples, xi = common_xi, omega = common_omega, alpha = common_alpha, nu = common_nu),
  times = 10
)
cat("\nRandom Variate (r) Benchmark (fewer times for large n):\n"); print(mb_r)

cat("\n--- Comparison Complete ---\n")
