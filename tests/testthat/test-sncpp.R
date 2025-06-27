# This file goes in tests/testthat/test-distributions.R

# Load the package we are testing
library(snCpp)

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

test_that("owens.T gives known values", {
  # T(h, 1) = 0.5 * (pnorm(h) - pnorm(h)^2)
  h <- 1.5
  expect_equal(owens.T(h, 1), 0.5 * (pnorm(h) - pnorm(h)^2), tolerance = 1e-7)
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


#Compare against the reference 'sn' package
if (has_sn_package) {
  test_that("snCpp output matches sn package output", {
    # Define a set of parameters to test
    xi <- 1
    omega <- 2
    alpha <- 4
    nu <- 5
    x <- 1.5
    p <- 0.75

    # Set a tolerance for numerical comparisons
    tol <- 1e-7

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
  })
}

