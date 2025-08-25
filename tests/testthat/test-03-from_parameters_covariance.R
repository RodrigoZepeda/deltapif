test_that("from_parameters_covariance_p_component works correctly", {
  # Simple test case
  p1 <- c(0.3, 0.7)
  p2 <- c(0.4, 0.6)
  p1_cft <- c(0.2, 0.8)
  p2_cft <- c(0.3, 0.7)
  rr1 <- c(1.5, 2.0)
  rr2 <- c(1.8, 1.2)
  mu_obs1 <- sum(p1 * rr1)
  mu_obs2 <- sum(p2 * rr2)
  mu_cft1 <- sum(p1_cft * rr1)
  mu_cft2 <- sum(p2_cft * rr2)
  var_p <- matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2)

  # Calculate expected value manually
  deriv1 <- deriv_pif_p(p1, p1_cft, rr1, mu_obs1, mu_cft1)
  deriv2 <- deriv_pif_p(p2, p2_cft, rr2, mu_obs2, mu_cft2)
  expected <- t(deriv1) %*% var_p %*% deriv2

  # Test function
  result <- from_parameters_covariance_p_component(
    p1, p2, p1_cft, p2_cft, rr1, rr2,
    var_p, FALSE
  )

  expect_equal(as.numeric(result), as.numeric(expected), tolerance = 1e-6)

  # Test upper bound case
  result_upper <- from_parameters_covariance_p_component(
    p1, p2, p1_cft, p2_cft, rr1, rr2,
    var_p, TRUE
  )
  expect_true(as.numeric(result_upper) >= abs(as.numeric(result)))
})

test_that("from_parameters_covariance_beta_component works correctly", {
  # Simple test case
  p1 <- c(0.3, 0.7)
  p2 <- c(0.4, 0.6)
  p1_cft <- c(0.2, 0.8)
  p2_cft <- c(0.3, 0.7)
  rr1 <- c(1.5, 2.0)
  rr2 <- c(1.8, 1.2)
  rr_deriv1 <- c(0.5, 0.3) # arbitrary derivative values
  rr_deriv2 <- c(0.4, 0.2)
  mu_obs1 <- sum(p1 * rr1)
  mu_obs2 <- sum(p2 * rr2)
  mu_cft1 <- sum(p1_cft * rr1)
  mu_cft2 <- sum(p2_cft * rr2)
  var_beta <- matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2)

  # Calculate expected value manually
  deriv1 <- deriv_pif_beta(p1, p1_cft, rr1, rr_deriv1)
  deriv2 <- deriv_pif_beta(p2, p2_cft, rr2, rr_deriv2)
  expected <- t(deriv1) %*% var_beta %*% deriv2

  # Test function
  result <- from_parameters_covariance_beta_component(
    p1, p2, p1_cft, p2_cft, rr1, rr2, rr_deriv1, rr_deriv2,
    var_beta, FALSE
  )

  expect_equal(as.numeric(result), as.numeric(expected), tolerance = 1e-6)

  # Test upper bound case
  result_upper <- from_parameters_covariance_beta_component(
    p1, p2, p1_cft, p2_cft, rr1, rr2, rr_deriv1, rr_deriv2,
    var_beta, TRUE
  )
  expect_true(as.numeric(result_upper) >= abs(as.numeric(result)))
})

test_that("from_parameters_pif_covariance combines components correctly", {
  # Use same test case as above
  p1 <- c(0.3, 0.7)
  p2 <- c(0.4, 0.6)
  p1_cft <- c(0.2, 0.8)
  p2_cft <- c(0.3, 0.7)
  rr1 <- c(1.5, 2.0)
  rr2 <- c(1.8, 1.2)
  rr_deriv1 <- c(0.5, 0.3)
  rr_deriv2 <- c(0.4, 0.2)
  mu_obs1 <- sum(p1 * rr1)
  mu_obs2 <- sum(p2 * rr2)
  mu_cft1 <- sum(p1_cft * rr1)
  mu_cft2 <- sum(p2_cft * rr2)
  var_p <- matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2)
  var_beta <- matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2)

  # Calculate expected values
  p_comp <- from_parameters_covariance_p_component(
    p1, p2, p1_cft, p2_cft, rr1, rr2,
    var_p, FALSE
  )

  beta_comp <- from_parameters_covariance_beta_component(
    p1, p2, p1_cft, p2_cft, rr1, rr2, rr_deriv1, rr_deriv2,
    var_beta, FALSE
  )

  expected <- as.numeric(p_comp + beta_comp)

  # Test function
  result <- from_parameters_pif_covariance(
    p1, p2, p1_cft, p2_cft, rr1, rr2, rr_deriv1, rr_deriv2,
    var_p, var_beta,
    FALSE, FALSE
  )

  expect_equal(result, expected, tolerance = 1e-6)

  # Test with upper bounds
  result_upper <- from_parameters_pif_covariance(
    p1, p2, p1_cft, p2_cft, rr1, rr2, rr_deriv1, rr_deriv2,
    var_p, var_beta,
    TRUE, TRUE
  )
  expect_true(result_upper >= abs(result))
})

test_that("from_parameters_pif_variance works correctly", {
  # Test case where p1 = p2, etc (variance case)
  p <- c(0.3, 0.7)
  p_cft <- c(0.2, 0.8)
  rr <- c(1.5, 2.0)
  rr_deriv <- c(0.5, 0.3)
  mu_obs <- sum(p * rr)
  mu_cft <- sum(p_cft * rr)
  var_p <- matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2)
  var_beta <- matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2)

  # Should be equivalent to covariance with itself
  expected <- from_parameters_pif_covariance(
    p, p, p_cft, p_cft, rr, rr, rr_deriv, rr_deriv,
    var_p, var_beta,
    FALSE, FALSE
  )

  result <- from_parameters_pif_variance(
    p, p_cft, rr, rr_deriv, var_p, var_beta,
    FALSE, FALSE
  )

  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("Edge cases are handled properly", {
  # Single-element vectors
  p <- 0.5
  p_cft <- 0.4
  rr <- 1.5
  rr_deriv <- 0.2
  mu_obs <- p * rr
  mu_cft <- p_cft * rr
  var_p <- matrix(0.01)
  var_beta <- matrix(0.02)

  expect_silent(
    from_parameters_pif_variance(
      p, p_cft, rr, rr_deriv, var_p, var_beta
    )
  )

  # Zero variance case
  zero_var_p <- matrix(0, nrow = 2, ncol = 2)
  zero_var_beta <- matrix(0, nrow = 2, ncol = 2)
  p <- c(0.5, 0.5)
  p_cft <- c(0.4, 0.6)
  rr <- c(1.0, 1.0) # RR = 1 should give zero derivative

  result <- from_parameters_pif_variance(
    p, p_cft, rr, c(0.5, 0.5),
    zero_var_p, zero_var_beta
  )
  expect_equal(result, 0)
})

test_that("Input validation works", {
  # Length mismatches
  p <- c(0.3, 0.7)
  p_cft <- c(0.2, 0.8)
  rr <- c(1.5) # Too short

  expect_error(
    from_parameters_pif_variance(
      p, p_cft, rr, c(0.5, 0.5),
      sum(p * c(1.5, 2.0)), sum(p_cft * c(1.5, 2.0)),
      matrix(0.01, 2, 2), matrix(0.02, 2, 2)
    )
  )

  # Invalid probabilities
  expect_error(
    from_parameters_pif_variance(
      c(-0.1, 1.1), p_cft, c(1.5, 2.0), c(0.5, 0.5),
      sum(p * c(1.5, 2.0)), sum(p_cft * c(1.5, 2.0)),
      matrix(0.01, 2, 2), matrix(0.02, 2, 2)
    ),
    "Invalid probability"
  )

  # Invalid variance matrices
  expect_error(
    from_parameters_pif_variance(
      p, p_cft, c(1.5, 2.0), c(0.5, 0.5),
      sum(p * c(1.5, 2.0)), sum(p_cft * c(1.5, 2.0)),
      matrix(0.01, 1, 1), # Wrong dimension
      matrix(0.02, 2, 2)
    ),
  )
})
