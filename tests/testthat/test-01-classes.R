# Helper function to create valid pif_atomic_class instance
create_valid_pif_atomic <- function() {
  pif_atomic_class(
    p = c(0.3, 0.7),
    p_cft = c(0.2, 0.8),
    beta = c(0.5, 1.0),
    var_p = matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2),
    var_beta = matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = identity,
    link_deriv = function(x) 1,
    link_inv = identity,
    conf_level = 0.95,
    type = "PIF",
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE
  )
}

test_that("pif_class construction and validation works", {
  # Valid construction
  expect_silent(
    pif_class(
      pif = 0.3,
      variance = 0.01,
      conf_level = 0.95,
      type = "PIF",
      link = identity,
      link_inv = identity,
      link_deriv = function(x) 1
    )
  )

  # Invalid confidence level
  expect_error(
    pif_class(
      pif = 0.3,
      variance = 0.01,
      conf_level = 1.1,
      type = "PIF",
      link = identity,
      link_inv = identity,
      link_deriv = function(x) 1
    ),
    "Invalid confidence level"
  )

  # Invalid PIF value
  expect_error(
    pif_class(
      pif = 1.1,
      variance = 0.01,
      conf_level = 0.95,
      type = "PIF",
      link = identity,
      link_inv = identity,
      link_deriv = function(x) 1
    ),
    "PIF > 1"
  )

  # Invalid type
  expect_error(
    pif_class(
      pif = 0.3,
      variance = 0.01,
      conf_level = 0.95,
      type = "INVALID",
      link = identity,
      link_inv = identity,
      link_deriv = function(x) 1
    ),
    "should be either `PIF` or `PAF`"
  )

  # Test property access
  pif <- pif_class(
    pif = 0.3,
    variance = 0.01,
    conf_level = 0.95,
    type = "PIF",
    link = logit,
    link_inv = inv_logit,
    link_deriv = deriv_logit
  )

  expect_equal(pif@pif, 0.3)
  expect_equal(pif@type, "PIF")
  expect_equal(pif@link_vals, logit(0.3))
  expect_length(pif@ci, 2)
})

test_that("pif_atomic_class construction and validation works", {
  # Valid construction
  expect_silent(create_valid_pif_atomic())

  # Length mismatches
  expect_error(
    pif_atomic_class(
      p = c(0.3, 0.7),
      p_cft = c(0.2), # Different length
      beta = c(0.5, 1.0),
      var_p = matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2),
      var_beta = matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2),
      rr_link = function(x) exp(x),
      rr_link_deriv = function(x) exp(x),
      link = identity,
      link_deriv = function(x) 1,
      link_inv = identity,
      conf_level = 0.95,
      type = "PIF",
      upper_bound_p = FALSE,
      upper_bound_beta = FALSE
    ),
    "must be of the same length"
  )

  # Invalid probabilities
  expect_error(
    pif_atomic_class(
      p = c(-0.1, 1.1), # Invalid values
      p_cft = c(0.2, 0.8),
      beta = c(0.5, 1.0),
      var_p = matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2),
      var_beta = matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2),
      rr_link = function(x) exp(x),
      rr_link_deriv = function(x) exp(x),
      link = identity,
      link_deriv = function(x) 1,
      link_inv = identity,
      conf_level = 0.95,
      type = "PIF",
      upper_bound_p = FALSE,
      upper_bound_beta = FALSE
    ),
    "values < 0"
  )

  # Invalid variance matrices
  expect_error(
    pif_atomic_class(
      p = c(0.3, 0.7),
      p_cft = c(0.2, 0.8),
      beta = c(0.5, 1.0),
      var_p = matrix(c(0.01, 0.005, 0.005), nrow = 1), # Wrong dim
      var_beta = matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2),
      rr_link = function(x) exp(x),
      rr_link_deriv = function(x) exp(x),
      link = identity,
      link_deriv = function(x) 1,
      link_inv = identity,
      conf_level = 0.95,
      type = "PIF",
      upper_bound_p = FALSE,
      upper_bound_beta = FALSE
    ),
    "is not symmetric"
  )

  # Test property access
  pif_atomic <- create_valid_pif_atomic()
  expect_equal(pif_atomic@rr, exp(pif_atomic@beta))
  expect_equal(pif_atomic@mu_obs, sum(pif_atomic@p * pif_atomic@rr))
  expect_equal(pif_atomic@pif, 1 - (sum(pif_atomic@p_cft * pif_atomic@rr) / sum(pif_atomic@p * pif_atomic@rr)))
})

test_that("pif_total_class construction and validation works", {
  # Create some atomic PIFs first
  pif1 <- create_valid_pif_atomic()
  pif2 <- create_valid_pif_atomic()

  # Valid construction
  expect_silent(
    pif_total_class(
      pif_list = list(pif1, pif2),
      weights = c(0.5, 0.5),
      sigma_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
      conf_level = 0.95,
      link = identity,
      link_inv = identity,
      link_deriv = function(x) 1
    )
  )

  # Invalid pif_list
  expect_error(
    pif_total_class(
      pif_list = list(pif1, "not_a_pif"),
      weights = c(0.5, 0.5),
      sigma_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
      conf_level = 0.95,
      link = identity,
      link_inv = identity,
      link_deriv = function(x) 1
    ),
    "must be a 'pif_class'"
  )

  # Length mismatch
  expect_error(
    pif_total_class(
      pif_list = list(pif1, pif2),
      weights = c(0.5), # Too short
      sigma_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
      conf_level = 0.95,
      link = identity,
      link_inv = identity,
      link_deriv = function(x) 1
    ),
    "weights provided have length 1"
  )

  # Test property access
  pif_total <- pif_total_class(
    pif_list = list(pif1, pif2),
    weights = c(0.5, 0.5),
    sigma_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
    conf_level = 0.95,
    link = identity,
    link_inv = identity,
    link_deriv = function(x) 1
  )

  expect_equal(pif_total@coefs, c(pif1@pif, pif2@pif))
  expect_equal(pif_total@pif, as.numeric(t(c(0.5, 0.5)) %*% c(pif1@pif, pif2@pif)))
  expect_equal(suppressWarnings(dim(pif_total@covariance)), c(2, 2))
})

test_that("pif_ensemble_class construction and validation works", {
  # Create some atomic PIFs first
  pif1 <- create_valid_pif_atomic()
  pif2 <- create_valid_pif_atomic()

  # Valid construction
  expect_silent(
    pif_ensemble_class(
      pif_list = list(pif1, pif2),
      conf_level = 0.95,
      link = log_complement,
      link_inv = inv_log_complement,
      link_deriv = deriv_log_complement,
      weights = rep(1, 2),
      sigma_weights = matrix(0, 2, 2)
    )
  )

  # Invalid pif_list
  expect_error(
    pif_ensemble_class(
      pif_list = list(pif1, "not_a_pif"),
      conf_level = 0.95,
      link = log_complement,
      link_inv = inv_log_complement,
      link_deriv = deriv_log_complement,
      weights = rep(1, 2),
      sigma_weights = matrix(0, 2, 2)
    ),
    "must be a 'pif_class'"
  )

  # Test property access
  pif_ensemble <- pif_ensemble_class(
    pif_list = list(pif1, pif2),
    conf_level = 0.95,
    link = log_complement,
    link_inv = inv_log_complement,
    link_deriv = deriv_log_complement,
    weights = rep(1, 2),
    sigma_weights = matrix(0, 2, 2)
  )

  expect_equal(pif_ensemble@coefs, c(pif1@pif, pif2@pif))
  expect_equal(pif_ensemble@pif, 1 - (1 - pif1@pif) * (1 - pif2@pif))
  expect_equal(suppressWarnings(dim(pif_ensemble@covariance)), c(2, 2))

  # Verify link functions are set correctly
  expect_equal(pif_ensemble@link, log_complement)
  expect_equal(pif_ensemble@link_inv, inv_log_complement)
  expect_equal(pif_ensemble@link_deriv, deriv_log_complement)
})

test_that("Class inheritance works correctly", {
  pif_atomic <- create_valid_pif_atomic()
  pif1 <- create_valid_pif_atomic()
  pif2 <- create_valid_pif_atomic()
  pif_total <- pif_total_class(
    pif_list = list(pif1, pif2),
    weights = c(0.5, 0.5),
    sigma_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
    conf_level = 0.95,
    link = identity,
    link_inv = identity,
    link_deriv = function(x) 1
  )
  pif_ensemble <- pif_ensemble_class(
    pif_list = list(pif1, pif2),
    conf_level = 0.95,
    link = log_complement,
    link_inv = inv_log_complement,
    link_deriv = deriv_log_complement,
    weights = rep(1, 2),
    sigma_weights = matrix(0, 2, 2)
  )

  pif_global_ensemble <- pif_global_ensemble_class(
    pif_list = list(pif1, pif2),
    conf_level = 0.95,
    link = log_complement,
    link_inv = inv_log_complement,
    link_deriv = deriv_log_complement,
    weights = rep(1, 2),
    sigma_weights = matrix(0, 2, 2),
    pif_transform = identity,
    pif_deriv_transform = function(x) 1,
    pif_inverse_transform = identity
  )

  # Verify inheritance
  expect_true(S7::S7_inherits(pif_atomic, pif_class))
  expect_true(S7::S7_inherits(pif_global_ensemble, pif_class))
  expect_true(S7::S7_inherits(pif_total, pif_global_ensemble_class))
  expect_true(S7::S7_inherits(pif_ensemble, pif_global_ensemble_class))
  expect_true(S7::S7_inherits(pif_ensemble, pif_class))
  expect_true(S7::S7_inherits(pif_total, pif_class))

  # Verify methods work through inheritance
  expect_length(pif_atomic@ci, 2)
  expect_length(suppressWarnings(pif_total@ci), 2)
  expect_length(suppressWarnings(pif_ensemble@ci), 2)
})

test_that("Edge cases are handled properly", {
  # Single exposure category
  expect_silent(
    pif_atomic_class(
      p = 1,
      p_cft = 1,
      beta = 0.5,
      var_p = matrix(0.01),
      var_beta = matrix(0.02),
      rr_link = function(x) exp(x),
      rr_link_deriv = function(x) exp(x),
      link = identity,
      link_deriv = function(x) 1,
      link_inv = identity,
      conf_level = 0.95,
      type = "PIF",
      upper_bound_p = FALSE,
      upper_bound_beta = FALSE
    )
  )

  # Zero variance case
  zero_pif <- pif_atomic_class(
    p = c(0.3, 0.7),
    p_cft = c(0.2, 0.8),
    beta = c(0.5, 1.0),
    var_p = matrix(0, nrow = 2, ncol = 2),
    var_beta = matrix(0, nrow = 2, ncol = 2),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = identity,
    link_deriv = function(x) 1,
    link_inv = identity,
    conf_level = 0.95,
    type = "PIF",
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE
  )
  expect_equal(zero_pif@variance, 0)

  # Perfect correlation (upper bound)
  upper_pif <- pif_atomic_class(
    p = c(0.3, 0.7),
    p_cft = c(0.2, 0.8),
    beta = c(0.5, 1.0),
    var_p = matrix(c(0.01, 0.01, 0.01, 0.01), nrow = 2),
    var_beta = matrix(c(0.02, 0.02, 0.02, 0.02), nrow = 2),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = identity,
    link_deriv = function(x) 1,
    link_inv = identity,
    conf_level = 0.95,
    type = "PIF",
    upper_bound_p = TRUE,
    upper_bound_beta = TRUE
  )
  expect_true(upper_pif@variance > 0)
})
