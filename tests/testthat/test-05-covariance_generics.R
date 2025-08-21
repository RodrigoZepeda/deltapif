# Helper function to create consistent mock PIF objects
create_mock_pif <- function(p = 0.5, p_cft = 0.25, beta = 1.5,
                            var_p = 0.01, var_beta = 0.02) {
  pif_atomic_class(
    p = p,
    p_cft = p_cft,
    beta = beta,
    var_p = as.matrix(var_p),
    var_beta = as.matrix(var_beta),
    rr_link = identity,
    rr_link_deriv = function(x) 1,
    link = identity,
    link_inv = identity,
    link_deriv = function(x) 1,
    conf_level = 0.95,
    type = "PIF",
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE
  )
}

test_that("cov_atomic_pif validates inputs correctly", {
  pif1 <- create_mock_pif()
  pif2 <- create_mock_pif(p = 0.4, p_cft = 0.2, beta = 1.8)

  # Valid case
  expect_silent(cov_atomic_pif(pif1, pif2))

  # Invalid class inputs
  expect_error(
    cov_atomic_pif(list(), pif2),
    "not a `deltapif::pif_atomic_class` object"
  )


})

test_that("cov_atomic_pif handles independence parameters correctly", {
  pif1 <- create_mock_pif()
  pif2 <- create_mock_pif()

  # Test guessing when parameters are identical
  expect_message(
    cov_atomic_pif(pif1, pif2, uncorrelated_p = "guess", uncorrelated_beta = TRUE),
    "appear to be the same"
  )
  expect_message(
    cov_atomic_pif(pif1, pif2, uncorrelated_p = TRUE, uncorrelated_beta = "guess"),
    "appear to be the same"
  )
  expect_no_warning(
    cov_atomic_pif(pif1, pif2, uncorrelated_p = TRUE, uncorrelated_beta = TRUE)
  )

  # Test explicit independence settings
  expect_silent(cov_atomic_pif(pif1, pif2, uncorrelated_p = TRUE, uncorrelated_beta = TRUE))
  expect_silent(cov_atomic_pif(pif1, pif2, uncorrelated_p = FALSE, uncorrelated_beta = FALSE))
})

test_that("cov_atomic_pif calculates correctly", {
  pif1 <- create_mock_pif(p = c(0.4, 0.2), p_cft = c(0.1, 0.1), beta = c(1.1, 1.7),
                          var_p = matrix(0, 2, 2), var_beta = diag(c(0.01, 0.02), 2))
  pif2 <- create_mock_pif(p = c(0.5, 0.3), p_cft = c(0.1, 0.2), beta = c(1.1, 1.7),
                          var_p = matrix(0, 2, 2), var_beta = diag(c(0.01, 0.02), 2))

  # Test basic calculation
  result <- cov_atomic_pif(pif1, pif2, uncorrelated_p = FALSE, uncorrelated_beta = FALSE)
  expect_type(result, "double")

  # Test with custom variance matrices
  var_p <- matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2)
  var_beta <- matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2)
  expect_silent(
    cov_atomic_pif(pif1, pif2, var_p = var_p, var_beta = var_beta)
  )

  # Test with identical PIFs (should equal variance)
  expect_equal(
    cov_atomic_pif(pif1, pif1, uncorrelated_p = FALSE, uncorrelated_beta = FALSE),
    variance(pif1)
  )

})

test_that("cov_total_pif handles different PIF types", {
  atomic1 <- create_mock_pif()
  atomic2 <- create_mock_pif(p = 0.4, p_cft = 0.2, beta = 1.8)

  # Create a total PIF with two atomic PIFs
  total_pif <- pif_total_class(
    pif_list = list(atomic1, atomic2),
    weights = c(0.5, 0.5),
    sigma_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
    conf_level = 0.95,
    link = identity,
    link_inv = identity,
    link_deriv = function(x) 1
  )

  # Test atomic vs atomic
  expect_silent(suppressMessages(
    cov_total_pif(atomic1, atomic2,  uncorrelated_p = FALSE, uncorrelated_beta = FALSE))
  )

  # Test atomic vs total
  expect_silent(suppressMessages(
    cov_total_pif(atomic1, total_pif,  uncorrelated_p = FALSE, uncorrelated_beta = FALSE))
  )

  expect_silent(suppressMessages(
    cov_total_pif(total_pif, atomic1,  uncorrelated_p = FALSE, uncorrelated_beta = FALSE))
  )

  # Test total vs total
  expect_silent(suppressMessages(
    cov_total_pif(total_pif, total_pif,  uncorrelated_p = FALSE, uncorrelated_beta = FALSE))
  )

  # Test unsupported types
  expect_error(
    cov_total_pif(atomic1, list(),  uncorrelated_p = FALSE, uncorrelated_beta = FALSE),
    "Unsupported types"
  )
})

test_that("covariance generic works correctly", {
  pif1 <- create_mock_pif()
  pif2 <- create_mock_pif(p = 0.4, p_cft = 0.2, beta = 1.8)
  pif3 <- create_mock_pif(p = 0.3, p_cft = 0.1, beta = 2.0)

  # Single PIF case
  expect_equal(dim(covariance(pif1)), c(1, 1))

  # Multiple PIF case
  cov_mat <- covariance(pif1, pif2, pif3)
  expect_equal(dim(cov_mat), c(3, 3))
  expect_true(isSymmetric(cov_mat))

  # Diagonal should equal variances
  expect_equal(diag(cov_mat), c(variance(pif1), variance(pif2), variance(pif3)))

  # Test with custom independence matrices
  ind_p <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3)
  ind_beta <- matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1), nrow = 3)
  expect_silent(
    covariance(pif1, pif2, pif3, uncorrelated_p = ind_p, uncorrelated_beta = ind_beta)
  )

  # Test variance matrix input validation
  expect_error(
    covariance(pif1, pif2, var_p = matrix(1, nrow = 2)),
    "should be 2 x 2"
  )
})

test_that("variance generic works correctly", {
  pif1 <- create_mock_pif()

  # Basic variance calculation
  expect_type(variance(pif1), "double")
  expect_true(variance(pif1) > 0)

  # Should equal covariance with itself
  expect_equal(variance(pif1), covariance(pif1, pif1,  uncorrelated_p = FALSE, uncorrelated_beta = FALSE)[1,1])

  # Test warning for extra arguments
  expect_warning(
    variance(pif1, "extra_arg"),
    "does not support more than 1 argument"
  )
})

test_that("standard_deviation generic works correctly", {
  pif1 <- create_mock_pif()

  # Basic calculation
  expect_type(standard_deviation(pif1), "double")
  expect_equal(standard_deviation(pif1), sqrt(variance(pif1)))
})

test_that("correlation generic works correctly", {
  pif1 <- create_mock_pif()
  pif2 <- create_mock_pif(p = 0.4, p_cft = 0.2, beta = 1.8)

  # Single PIF case
  expect_equal(dim(correlation(pif1)), c(1, 1))
  expect_equal(correlation(pif1)[1,1], 1)

  # Multiple PIF case
  cor_mat <- correlation(pif1, pif2)
  expect_equal(dim(cor_mat), c(2, 2))
  expect_true(isSymmetric(cor_mat))
  expect_true(all(diag(cor_mat) == 1))
  expect_true(all(cor_mat >= -1 & cor_mat <= 1))

  # Test with shared parameters
  pif_shared_beta <- create_mock_pif(p = 0.3, p_cft = 0.1, beta = pif1@beta)
  cor_mat <- correlation(pif1, pif_shared_beta, uncorrelated_beta = FALSE)
  expect_true(cor_mat[1,2] != 0) # Should be correlated
})

test_that("edge cases are handled properly", {
  # Zero variance case
  pif_zero_var <- create_mock_pif(var_p = 0, var_beta = 0)
  expect_equal(variance(pif_zero_var), 0)

  # Correlation with zero variance
  expect_warning(correlation(pif_zero_var, pif_zero_var,  uncorrelated_p = FALSE, uncorrelated_beta = FALSE))

  # Single exposure category
  pif_single <- create_mock_pif(p = 1, p_cft = 1)
  expect_silent(variance(pif_single))
})
