# Helper function to create consistent mock PIF objects
create_mock_pif <- function(p = 0.5, p_cft = 0.25, beta = 1.5,
                            var_p = 0.01, var_beta = 0.02, label = paste0("mock", rnorm(1))) {
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
    upper_bound_beta = FALSE,
    label = label
  )
}

test_that("cov_atomic_pif validates inputs correctly", {
  pif1 <- create_mock_pif()
  pif2 <- create_mock_pif(p = 0.4, p_cft = 0.2, beta = 1.8)

  # Valid case
  expect_silent(cov_atomic_pif(pif1, pif2,
                               var_p = default_parameter_covariance_structure2(pif1, pif2),
                               var_beta = default_parameter_covariance_structure2(pif1, pif2)))

  # Invalid class inputs
  expect_error(
    cov_atomic_pif(list(), pif2),
    "pif_atomic_class"
  )

  expect_error(
    cov_atomic_pif(pif1, list()),
    "pif_atomic_class"
  )

})

test_that("cov_atomic_pif calculates correctly", {
  pif1 <- create_mock_pif(p = c(0.4, 0.2), p_cft = c(0.1, 0.1), beta = c(1.1, 1.7),
                          var_p = matrix(0, 2, 2), var_beta = diag(c(0.01, 0.02), 2),
                          label = "pif1")
  pif2 <- create_mock_pif(p = c(0.5, 0.3), p_cft = c(0.1, 0.2), beta = c(1.1, 1.7),
                          var_p = matrix(0, 2, 2), var_beta = diag(c(0.01, 0.02), 2),
                          label = "pif2")

  var_p    <- default_parameter_covariance_structure2(pif1, pif2, parameter = "p")
  var_p    <- subset(var_p, cols = "pif1", rows = "pif2")
  var_beta <- default_parameter_covariance_structure2(pif1, pif2, parameter = "beta")
  var_beta <- subset(var_beta, cols = "pif1", rows = "pif2")

  # Test basic calculation
  result <- cov_atomic_pif(pif1, pif2, as.matrix(var_p), as.matrix(var_beta))
  expect_type(result, "double")

  # Test with custom variance matrices
  var_p <- matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2)
  var_beta <- matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2)
  expect_silent(
    cov_atomic_pif(pif1, pif2, var_p = var_p, var_beta = var_beta)
  )

  # Test with identical PIFs (should equal variance)
  expect_equal(
     cov_atomic_pif(pif1, pif1, var_p = pif1@var_p, var_beta = pif1@var_beta),
     variance(pif1)
  )

})

test_that("cov_total_pif handles different PIF types", {
  atomic1 <- create_mock_pif(label = "1")
  atomic2 <- create_mock_pif(p = 0.4, p_cft = 0.2, beta = 1.8, label = "2")
  atomic3 <- create_mock_pif(p = 0.4, p_cft = 0.2, beta = 1.8, label = "3")

  # Create a total PIF with two atomic PIFs
  pif_list <- list(atomic2, atomic3)
  names(pif_list) <- c(atomic2@label, atomic3@label)
  total_pif <- pif_total_class(
    pif_list = pif_list,
    weights = c(0.5, 0.5),
    var_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
    conf_level = 0.95,
    link = identity,
    link_inv = identity,
    label = "total",
    var_pif_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
    link_deriv = function(x) 1
  )

  # Test atomic vs atomic
  expect_silent(
    cov_total_pif(atomic1, atomic2)
  )

  expect_silent(
    cov_total_pif(atomic1, total_pif)
  )

  expect_silent(
    cov_total_pif(total_pif, atomic1)
  )

  # Test total vs total
  expect_silent(
    cov_total_pif(total_pif, total_pif)
  )

  # # Test unsupported types
  expect_error(
    cov_total_pif(atomic1, list()),
    " must be of `pif_class`"
  )
})

test_that("covariance generic works correctly", {
  pif1 <- create_mock_pif(label = "pif1")
  pif2 <- create_mock_pif(p = 0.4, p_cft = 0.2, beta = 1.8, label = "pif2")
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
    covariance(pif1, pif2, pif3)
  )

  # Test variance matrix input validation
  expect_error(
    covariance(pif1, pif2,
               var_p = matrix(1, nrow = 2, ncol = 2, dimnames = list(list("pif", "pif2"), list("pif1", "pif2")))),
    "was not found"
  )
})

test_that("variance generic works correctly", {
  pif1 <- create_mock_pif()

  # Basic variance calculation
  expect_type(variance(pif1), "double")
  expect_true(variance(pif1) > 0)

  # Should equal covariance with itself
  expect_equal(variance(pif1), covariance(pif1, pif1)[1,1])

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
  cor_mat <- correlation(pif1, pif_shared_beta)
  expect_true(cor_mat[1,2] != 0) # Should be correlated
})

test_that("edge cases are handled properly", {
  # Zero variance case
  pif_zero_var <- create_mock_pif(var_p = 0, var_beta = 0)
  expect_equal(variance(pif_zero_var), 0)

  # Correlation with zero variance
  expect_warning(correlation(pif_zero_var, pif_zero_var))

  # Single exposure category
  pif_single <- create_mock_pif(p = 1, p_cft = 1)
  expect_silent(variance(pif_single))
})

# Extended tests for covariance generics

# Test cov_ensemble_weights function
test_that("cov_ensemble_weights works correctly", {
  # Create atomic PIFs
  pif1 <- create_mock_pif(label = "atomic1")
  pif2 <- create_mock_pif(label = "atomic2")
  pif3 <- create_mock_pif(label = "atomic3")
  pif4 <- create_mock_pif(label = "atomic4")

  # Create ensemble PIFs
  ensemble1 <- pif_ensemble(pif1, pif2, weights = c(0.6, 0.4), label = "ens1")
  ensemble2 <- pif_total(pif3, pif4, weights = c(0.3, 0.7), label = "total1")

  # Should return 0 for atomic PIFs (no weights)
  expect_equal(cov_ensemble_weights(pif1, pif2), 0)
  expect_equal(cov_ensemble_weights(ensemble1, pif1), 0)

  # Should work for ensemble with ensemble
  expect_silent(cov_ensemble_weights(ensemble1, ensemble2))

  # Test with custom variance matrices
  var_weights <- matrix(c(0.1, 0.05, 0.05, 0.2), nrow = 2)
  var_pif_weights <- matrix(c(0.01, 0.02, 0.03, 0.04), nrow = 2)

  expect_silent(
    cov_ensemble_weights(ensemble1, ensemble2,
                         var_weights = var_weights,
                         var_pif_weights = var_pif_weights)
  )

  # Test error handling
  expect_error(
    cov_ensemble_weights(list(), ensemble1),
    "must be a `pif_global_ensemble_class`"
  )
})

# Test cov_ensemble_atomic function
test_that("cov_ensemble_atomic works correctly", {
  atomic1 <- create_mock_pif(label = "atomic1")
  atomic2 <- create_mock_pif(label = "atomic2")
  atomic3 <- create_mock_pif(label = "atomic3")

  # Create ensemble
  ensemble1 <- pif_ensemble(atomic1, atomic2, weights = c(0.5, 0.5), label = "ens1")

  # Test ensemble with atomic
  expect_silent(cov_ensemble_atomic(ensemble1, atomic3))

  # Test atomic with atomic (should call cov_atomic_pif)
  result_atomic <- cov_ensemble_atomic(atomic1, atomic2)
  result_direct <- cov_atomic_pif(atomic1, atomic2)
  expect_equal(result_atomic, result_direct)

  # Test with custom variance structures
  var_pifs <- matrix(0.01, nrow = 2, ncol = 1)
  var_pif_weights <- matrix(0.005, nrow = 2, ncol = 1)

  expect_silent(
    cov_ensemble_atomic(ensemble1, atomic3,
                        var_pifs = var_pifs,
                        var_pif_weights = var_pif_weights)
  )

  # Test error handling
  expect_error(
    cov_ensemble_atomic(list(), atomic1),
    "should be a `pif_global_ensemble_class` or a `pif_atomic`"
  )

  expect_error(
    cov_ensemble_atomic(ensemble1, list()),
    "should be a `pif_atomic_class`"
  )
})

# Test covariance structures integration
test_that("covariance functions work with covariance_structure_class", {
  pif1 <- create_mock_pif(label = "pif1")
  pif2 <- create_mock_pif(label = "pif2")

  # Create covariance structures
  var_p_struct <- default_parameter_covariance_structure2(pif1, pif2, parameter = "p")
  var_beta_struct <- default_parameter_covariance_structure2(pif1, pif2, parameter = "beta")

  # Test with structures
  expect_silent(
    cov_atomic_pif(pif1, pif2, var_p = var_p_struct, var_beta = var_beta_struct)
  )

  # Test structure validation
  expect_error(
    cov_atomic_pif(pif1, pif2, var_p = "invalid"),
    "should be a `covariance_structure_class`"
  )
})

# Test dimension validation thoroughly
test_that("dimension validation works correctly", {
  pif1 <- pif(p = c(0.3, 0.4), beta = c(1.1, 1.5), label = "pif1", quiet = T)
  pif2 <- pif(p = c(0.2, 0.5, 0.1), beta = c(1.2, 1.8, 2.0), label = "pif2", quiet = T)

  # Correct dimensions
  var_p_correct <- matrix(0.01, nrow = 2, ncol = 3)
  var_beta_correct <- matrix(0.02, nrow = 2, ncol = 3)

  expect_silent(
    cov_atomic_pif(pif1, pif2, var_p = var_p_correct, var_beta = var_beta_correct)
  )

  # Wrong p dimensions
  expect_error(
    cov_atomic_pif(pif1, pif2, var_p = matrix(0.01, nrow = 3, ncol = 2)),
    "should be a matrix of dimensions 2 x 3"
  )

  # Wrong beta dimensions
  expect_error(
    cov_atomic_pif(pif1, pif2, var_beta = matrix(0.02, nrow = 3, ncol = 3)),
    "should be a matrix of dimensions 2 x 3"
  )
})

# Test cases_class variance methods
test_that("cases_class variance methods work", {
  base_pif <- create_mock_pif()

  cases_obj <- cases_class(
    overall_cases = 1000,
    pif_obj = base_pif,
    variance_cases = 100,
    link = identity,
    link_deriv = function(x) rep(1, length(x)),
    link_inv = identity,
    conf_level = 0.95
  )

  # Test variance method
  expect_equal(variance(cases_obj), cases_obj@variance)

  # Test standard deviation
  expect_equal(standard_deviation(cases_obj), sqrt(cases_obj@variance))

  # Test warning for extra args
  expect_warning(
    variance(cases_obj, "extra"),
    "does not support more than 1 argument"
  )
})

# Test recursive covariance calculations
test_that("recursive covariance calculations work", {
  # Create nested structure
  atomic1 <- create_mock_pif(label = "atomic1")
  atomic2 <- create_mock_pif(label = "atomic2")
  atomic3 <- create_mock_pif(label = "atomic3")
  atomic4 <- create_mock_pif(label = "atomic4")

  # First level ensembles
  ensemble1 <- pif_ensemble(atomic1, atomic2, weights = c(0.5, 0.5), label = "ens1")
  ensemble2 <- pif_total(atomic3, atomic4, weights = c(0.3, 0.7), label = "total1")

  # Second level ensemble
  meta_ensemble <- pif_total(ensemble1, ensemble2, weights = c(0.4, 0.6), label = "meta")

  # Test recursive calculation
  #expect_silent(cov_total_pif(meta_ensemble, meta_ensemble))
  #expect_silent(cov_total_pif(meta_ensemble, atomic1))
  expect_silent(cov_total_pif(ensemble1, ensemble2))
})

# Test covariance with identical PIFs
test_that("covariance with identical PIFs equals variance", {
  pif1 <- create_mock_pif()

  # Single PIF covariance should equal variance
  expect_equal(
    covariance(pif1)[1,1],
    variance(pif1)
  )

  # Atomic PIF covariance with itself
  expect_equal(
    cov_atomic_pif(pif1, pif1, var_p = pif1@var_p, var_beta = pif1@var_beta),
    variance(pif1)
  )

  # Total PIF covariance with itself
  pif2 <- create_mock_pif(label = "pif2")
  total_pif <- pif_total(pif1, pif2, weights = c(0.4, 0.6))
  expect_equal(
    cov_total_pif(total_pif, total_pif),
    variance(total_pif)
  )
})

# Test covariance matrix properties
test_that("covariance matrices have correct properties", {
  pif1 <- create_mock_pif(label = "pif1")
  pif2 <- create_mock_pif(label = "pif2")
  pif3 <- create_mock_pif(label = "pif3")

  # Test matrix properties
  cov_mat <- covariance(pif1, pif2, pif3)

  # Should be symmetric
  expect_true(isSymmetric(cov_mat))

  # Diagonal should be positive (variances)
  expect_true(all(diag(cov_mat) > 0))

  # Should be positive semi-definite (all eigenvalues >= 0)
  eigenvals <- eigen(cov_mat)$values
  expect_true(all(eigenvals >= -1e-10)) # Allow for numerical precision

  # Off-diagonal elements should be covariances
  expect_equal(cov_mat[1,2], cov_mat[2,1])
  expect_equal(cov_mat[1,3], cov_mat[3,1])
  expect_equal(cov_mat[2,3], cov_mat[3,2])
})

# Test correlation matrix properties
test_that("correlation matrices have correct properties", {
  pif1 <- create_mock_pif(label = "pif1")
  pif2 <- create_mock_pif(label = "pif2")
  pif3 <- create_mock_pif(label = "pif3")

  cor_mat <- correlation(pif1, pif2, pif3)

  # Should be symmetric
  expect_true(isSymmetric(cor_mat))

  # Diagonal should be 1
  expect_true(all(abs(diag(cor_mat) - 1) < 1e-10))

  # All values should be between -1 and 1
  expect_true(all(cor_mat >= -1 & cor_mat <= 1))

  # Should be positive semi-definite
  eigenvals <- eigen(cor_mat)$values
  expect_true(all(eigenvals >= -1e-10))
})

# Test with shared parameters
test_that("covariance correctly handles shared parameters", {
  # Create PIFs with same beta
  shared_beta <- 1.5
  pif1 <- create_mock_pif(p = 0.3, beta = shared_beta, label = "pif1")
  pif2 <- create_mock_pif(p = 0.5, beta = shared_beta, label = "pif2")

  # Should have positive covariance due to shared beta
  cov_val <- cov_atomic_pif(pif1, pif2,
                            var_p = matrix(0, nrow = 1, ncol = 1),
                            var_beta = matrix(0.01, nrow = 1, ncol = 1))
  expect_true(cov_val > 0)

  # Test with shared p (less realistic but mathematically valid)
  shared_p <- 0.4
  pif3 <- create_mock_pif(p = shared_p, beta = 1.2, label = "pif3")
  pif4 <- create_mock_pif(p = shared_p, beta = 1.8, label = "pif4")

  cov_val2 <- cov_atomic_pif(pif3, pif4,
                             var_p = matrix(0.01, nrow = 1, ncol = 1),
                             var_beta = matrix(0, nrow = 1, ncol = 1))
  expect_true(cov_val2 > 0)
})

# Test zero covariance cases
test_that("zero covariance is handled correctly", {
  pif1 <- create_mock_pif(var_p = 0, var_beta = 0, label = "pif1")
  pif2 <- create_mock_pif(var_p = 0, var_beta = 0, label = "pif2")

  # Should be zero when no shared parameters and no variance
  expect_equal(
    cov_atomic_pif(pif1, pif2, var_p = 0, var_beta = 0),
    0
  )

  # Variance should also be zero
  expect_equal(variance(pif1), 0)
})

# Test edge cases with different exposure categories
test_that("different exposure categories work correctly", {
  # Single exposure vs multiple exposures
  pif_single <- create_mock_pif(p = 0.5, p_cft = 0.2, beta = 1.5, label = "single")
  pif_multi <- create_mock_pif(p = c(0.3, 0.2), p_cft = c(0.1, 0.1),
                               var_p = diag(c(0.1, 0.2)),
                               var_beta = diag(c(0.1, 0.2)),
                               beta = c(1.2, 1.8), label = "multi")

  # Should handle different dimensions
  var_p <- matrix(0.01, nrow = 1, ncol = 2)
  var_beta <- matrix(0.02, nrow = 1, ncol = 2)

  expect_silent(
    cov_atomic_pif(pif_single, pif_multi, var_p = var_p, var_beta = var_beta)
  )
})

# Test input validation for covariance generic
test_that("covariance generic validates inputs properly", {
  pif1 <- create_mock_pif(label = "pif1")
  pif2 <- create_mock_pif(label = "pif2")
  pif3 <- create_mock_pif(label = "pif3")

  # Test with named variance matrices
  var_beta_named <- matrix(c(0.01, 0.005, 0.002,
                             0.005, 0.02, 0.001,
                             0.002, 0.001, 0.015), nrow = 3)
  rownames(var_beta_named) <- c("pif1", "pif2", "pif3")
  colnames(var_beta_named) <- c("pif1", "pif2", "pif3")


  expect_silent(
    covariance(pif1, pif2, pif3, var_beta = var_beta_named)
  )

  # Test error when PIF name not in matrix
  var_beta_missing <- var_beta_named
  colnames(var_beta_missing)[1] <- "missing_pif"

  expect_error(
    covariance(pif1, pif2, pif3, var_beta = var_beta_missing),
    "was not found in `var_beta`"
  )

  # Test automatic column naming
  var_beta_unnamed <- var_beta_named
  colnames(var_beta_unnamed) <- NULL
  rownames(var_beta_unnamed) <- NULL

  expect_silent(
    covariance(pif1, pif2, pif3, var_beta = var_beta_unnamed)
  )


})

# Test numerical stability
test_that("covariance calculations are numerically stable", {
  # Test with very small variances
  pif1 <- create_mock_pif(var_p = 1e-10, var_beta = 1e-10, label = "pif1")
  pif2 <- create_mock_pif(var_p = 1e-10, var_beta = 1e-10, label = "pif2")

  result <- covariance(pif1, pif2)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))

  # Test with large variances
  pif3 <- create_mock_pif(var_p = 100, var_beta = 100, label = "pif3")
  pif4 <- create_mock_pif(var_p = 100, var_beta = 100, label = "pif4")

  result2 <- covariance(pif3, pif4)
  expect_true(all(is.finite(result2)))
})

# Test ensemble weight covariance edge cases
test_that("ensemble weight covariance handles edge cases", {
  atomic1 <- create_mock_pif(label = "atomic1")
  atomic2 <- create_mock_pif(label = "atomic2")

  # Single component ensemble
  single_ensemble <- pif_total(atomic1, weights = 1, label = "single")

  # Should handle single component
  expect_equal(cov_ensemble_weights(single_ensemble, single_ensemble), 0)

  # Test with zero weights
  zero_weight_ensemble <- pif_total(atomic1, atomic2, weights = c(0, 1), label = "zero")
  expect_silent(cov_ensemble_weights(zero_weight_ensemble, zero_weight_ensemble))
})

# Test complex nested structures
test_that("complex nested covariance structures work", {
  # Create a complex hierarchy
  base1 <- create_mock_pif(label = "base1")
  base2 <- create_mock_pif(label = "base2")
  base3 <- create_mock_pif(label = "base3")
  base4 <- create_mock_pif(label = "base4")
  base5 <- create_mock_pif(label = "base5")
  base6 <- create_mock_pif(label = "base6")

  # Level 1: ensembles
  ens1 <- pif_ensemble(base1, base2, weights = c(0.6, 0.4), label = "ens1")
  ens2 <- pif_ensemble(base3, base4, weights = c(0.5, 0.5), label = "ens2")
  ens3 <- pif_ensemble(base5, base6, weights = c(0.5, 0.5), label = "ens3")

  # Level 2: total of ensembles
  total1 <- pif_total(ens1, ens2, weights = c(0.3, 0.7), label = "total1")

  # Should handle all levels
  expect_silent(variance(total1))
  expect_silent(covariance(total1, ens3))
  expect_silent(covariance(total1, base5))

})
