test_that("Link functions work correctly", {
  # Test logit
  expect_equal(logit(0.5), 0)
  expect_equal(logit(0.75), log(3))
  expect_equal(logit(0.25), -log(3))

  # Test log_complement
  expect_equal(log_complement(0), 0)
  expect_equal(log_complement(0.5), log(0.5))

  # Test hawkins
  expect_equal(hawkins(0), 0)
  expect_equal(hawkins(1), log(1 + sqrt(2)))
})

test_that("Inverse link functions work correctly", {
  # Test inv_logit
  expect_equal(inv_logit(0), 0.5)
  expect_equal(inv_logit(log(3)), 0.75, tolerance = 1e-7)

  # Test inv_log_complement
  expect_equal(inv_log_complement(0), 0)
  expect_equal(inv_log_complement(log(0.5)), 0.5)

  # Test inv_hawkins
  expect_equal(inv_hawkins(0), 0)
  expect_equal(inv_hawkins(log(1 + sqrt(2))), 1, tolerance = 1e-7)
})

test_that("Derivatives of link functions work correctly", {
  # Test deriv_logit
  expect_equal(deriv_logit(0.5), 4)
  expect_equal(deriv_logit(0.25), 16/3, tolerance = 1e-7)

  # Test deriv_log_complement
  expect_equal(deriv_log_complement(0.5), -2)
  expect_equal(deriv_log_complement(0.75), -4)

  # Test deriv_hawkins
  expect_equal(deriv_hawkins(0), 1)
  expect_equal(deriv_hawkins(1), 1/sqrt(2), tolerance = 1e-7)
})

test_that("Link parsers work correctly", {
  # Test parse_link
  expect_identical(parse_link("logit"), logit)
  expect_identical(parse_link("log-complement"), log_complement)
  expect_identical(parse_link("hawkins"), hawkins)
  expect_identical(parse_link("identity"), identity)
  expect_identical(parse_link("exponential"), exp)

  # Test parse_inv_link
  expect_identical(parse_inv_link("logit"), inv_logit)
  expect_identical(parse_inv_link("log-complement"), inv_log_complement)
  expect_identical(parse_inv_link("hawkins"), inv_hawkins)
  expect_identical(parse_inv_link("identity"), identity)
  expect_identical(parse_inv_link("exponential"), log)

  # Test function passthrough
  custom_fn <- function(x) x^2
  expect_identical(parse_link(custom_fn), custom_fn)
  expect_identical(parse_inv_link(custom_fn), custom_fn)

  # Test error for unknown link
  expect_error(parse_link("unknown_link"))
  expect_error(parse_inv_link("unknown_link"))
})

test_that("Link and inverse functions are inverses of each other", {
  test_values <- c(0.1, 0.25, 0.5, 0.75, 0.9)

  # Test logit and inv_logit
  for (x in test_values) {
    expect_equal(inv_logit(logit(x)), x, tolerance = 1e-7)
    expect_equal(logit(inv_logit(x)), x, tolerance = 1e-7)
  }

  # Test log_complement and inv_log_complement
  for (x in test_values) {
    expect_equal(inv_log_complement(log_complement(x)), x, tolerance = 1e-7)
  }

  # Test hawkins and inv_hawkins
  for (x in test_values) {
    expect_equal(inv_hawkins(hawkins(x)), x, tolerance = 1e-7)
    expect_equal(hawkins(inv_hawkins(x)), x, tolerance = 1e-7)
  }
})

test_that("Edge cases are handled properly", {
  # Test logit at boundaries (should produce Inf/-Inf)
  expect_equal(logit(1), Inf)
  expect_equal(logit(0), -Inf)

  # Test log_complement at boundaries
  expect_equal(log_complement(1), -Inf)

  # Test inv_log_complement with large positive input
  expect_equal(inv_log_complement(100), 1 - exp(100))

  # Test derivatives at boundaries
  expect_equal(abs(deriv_logit(1)), Inf)
  expect_equal(abs(deriv_logit(0)), Inf)
  expect_equal(abs(deriv_log_complement(1)), Inf)

  # Test invalid inputs
  expect_error(logit(2)) # pif > 1
  expect_error(logit(-1)) # pif < 0
  expect_error(log_complement(2))
  expect_error(deriv_logit(2))
})
