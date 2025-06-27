# This file goes in tests/testthat/test-distributions.R

# Load the package we are testing
library(snCpp)

# It's good practice to also have a known-good implementation to test against.
# We'll use the 'sn' package as a reference.
# Make sure to add 'sn' to the "Suggests" field in your DESCRIPTION file
# if you use this comparison test.
has_sn_package <- requireNamespace("sn", quietly = TRUE)

context("Skew-Normal Distribution Functions")

test_that("dsn works for standard normal case (alpha = 0)", {
  # For alpha=0, dsn should be identical to dnorm
  x <- c(-1, 0, 1)
  expect_equal(dsn(x, alpha = 0), dnorm(x))
})

test_that("rsn returns correct number of variates", {
  n <- 10
  expect_length(rsn(n, alpha = 5), n)
  expect_is(rsn(n, alpha = 5), "numeric")
})

test_that("qsn has correct symmetry for negative alpha", {
  p <- 0.8
  alpha <- 5
  # The property is Q(p; -alpha) = -Q(1-p; alpha)
  expect_equal(qsn(p, alpha = -alpha), -qsn(1 - p, alpha = alpha))
})


context("Skew-t Distribution Functions")

test_that("dst works for standard t case (alpha = 0)", {
  # For alpha=0, dst should be identical to dt
  x <- c(-1, 0, 1)
  nu <- 5
  expect_equal(dst(x, alpha = 0, nu = nu), dt(x, df = nu))
})

test_that("Parameter checks throw errors", {
  # Your functions should throw an error for invalid parameters
  expect_error(dsn(0, omega = -1), "'omega' must be a single positive numeric value.")
  expect_error(rst(5, nu = 0), "'nu' must be a single positive numeric value")
})


context("Special Functions (Owen's T)")

test_that("owens.T gives known values for specific cases", {
  # For a=1, T(h, 1) has a known closed-form solution
  h <- 1.5
  expect_equal(owens.T(h, 1), 0.5 * (pnorm(h) - pnorm(h)^2), tolerance = 1e-8)
})

test_that("owens.T handles boundary conditions correctly", {
  # T(h, 0) should be 0 for any h
  expect_equal(owens.T(h = 1.23, a = 0), 0)
  
  # T(0, a) should be atan(a) / (2*pi)
  a <- 0.75
  expect_equal(owens.T(h = 0, a = a), atan(a) / (2 * pi))
})

test_that("owens.T exhibits correct symmetry", {
  h <- 1.5
  a <- 0.5
  
  # Symmetry property 1: T(h, a) is even in h, so T(h, a) = T(-h, a)
  expect_equal(owens.T(h, a), owens.T(-h, a))
  
  # Symmetry property 2: T(h, a) is odd in a, so T(h, a) = -T(h, -a)
  expect_equal(owens.T(h, a), -owens.T(h, -a))
})


# Advanced: Compare against the reference 'sn' package
if (has_sn_package) {
  test_that("snCpp output matches sn package output", {
    # Define a set of parameters to test
    xi <- 1
    omega <- 0.3
    alpha <- 4
    nu <- 5
    x <- seq(-3, 3, length.out=10)
    p <- seq(0.01, 0.99, length.out=10)
    
    # Set a tolerance for numerical comparisons
    tol <- 1e-8
    
    # Compare dsn vs sn::dsn
    expect_equal(
      dsn(x, xi, omega, alpha),
      sn::dsn(x, xi, omega, alpha),
      tolerance = tol
    )
    
    # Compare pst vs sn::pst
    expect_equal(
      pst(x, xi, omega, alpha, nu),
      sn::pst(x, dp = c(xi, omega, alpha, nu)),
      tolerance = tol
    )
    
    # Compare qst vs sn::qst
    expect_equal(
      qst(p, xi, omega, alpha, nu),
      sn::qst(p, dp = c(xi, omega, alpha, nu)),
      tolerance = tol
    )
    
    # Compare owens.T vs sn::owensT
    expect_equal(
      owens.T(x, alpha),
      sn::T.Owen(x, alpha),
      tolerance = tol
    )
  })
}
