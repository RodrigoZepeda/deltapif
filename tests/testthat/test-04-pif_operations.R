# Helper function to create mock PIF objects
create_mock_pif_atomic <- function(type = "PIF", label = "a") {
  pif_atomic_class(
    p = c(0.3, 0.7),
    p_cft = c(0.2, 0.5),
    beta = c(0.5, 1.0),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = logit,
    link_deriv = deriv_logit,
    link_inv = inv_logit,
    var_p = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
    var_beta = matrix(c(0.02, 0, 0, 0.02), nrow = 2),
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE,
    conf_level = 0.95,
    label = label,
    type = type
  )
}

create_var_pif_weights_example <- function(n, include_null_case = TRUE) {
  # Create a non-NULL example
  if (n <= 0) {
    return(NULL)
  }

  sigma_list <- vector("list", n)

  for (i in 1:n) {
    sigma_list[[i]] <- vector("list", n)
    for (j in 1:n) {
      # Create a positive semi-definite covariance matrix
      # Using a simple approach: A %*% t(A) + small diagonal for stability
      A <- matrix(stats::rnorm(n * n), nrow = n)
      cov_matrix <- A %*% t(A) + diag(0.1, n)

      # Ensure symmetry (numerical precision might cause tiny asymmetries)
      cov_matrix <- (cov_matrix + t(cov_matrix)) / 2

      sigma_list[[i]][[j]] <- cov_matrix
    }
  }

  # If requested, randomly return NULL sometimes
  if (include_null_case && runif(1) < 0.3) {
    return(NULL)
  }

  return(sigma_list)
}

test_that("paf_total validates inputs correctly", {
  # Create mock PAFs
  paf1 <- create_mock_pif_atomic(type = "PAF", label = "1")
  paf2 <- create_mock_pif_atomic(type = "PAF", label = "2")

  # Valid case
  expect_silent(paf_total(paf1, paf2, weights = c(0.5, 0.5)))

})

test_that("pif_total validates inputs correctly", {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic(label = "b")

  # Valid cases
  expect_silent(pif_total(pif1, pif2, weights = c(0.5, 0.5)))
  expect_silent(pif_total(pif1, weights = 1)) # Single PIF case

  # Invalid weights
  expect_error(
    pif_total(pif1, pif2, weights = c(0.5, 0.6)), # Doesn't sum to 1
    "should sum to 1"
  )

  expect_error(
    pif_total(pif1, pif2, weights = c(1)), # Wrong length
    "provided have length 1"
  )

  #This one used to give an error leaving it here just in case
  expect_silent(
    pif_total(pif1, pif2, weights = c(0.5, 0.5), var_weights = matrix(1:4, nrow = 2))
  )


  # Invalid var_weights
  expect_error(
     pif_total(pif1, pif1, weights = c(0.5, 0.5), var_weights = matrix(1:4, nrow = 2)),
     "duplicated labels"
  )

  expect_error(
    pif_total(pif1, pif2, weights = c(0.5, 0.5), var_weights = "a"),
    "should be a number"
  )

  expect_error(
    pif_total(pif1, pif2, weights = c(0.5, 0.5), var_weights = 1:3)
    ,
    "dimensions"
  )


  # Test link function handling
  expect_silent(pif_total(pif1, weights = 1, link = "logit"))
  expect_silent(pif_total(pif1, weights = 1, link = logit, link_inv = inv_logit, link_deriv = deriv_logit))

})

test_that("pif_total calculates correctly", {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic(label = "2")

  # Test weights
  result <- pif_total(pif1, pif2, weights = c(0.7, 0.3))
  expect_equal(result@pif, 0.7*pif1@pif + 0.3*pif2@pif)

  # Test single PIF
  result <- pif_total(pif1, weights = 1)
  expect_equal(result@pif, pif1@pif)

  # Test variance calculation
  result <- pif_total(pif1, pif2, weights = c(0.5, 0.5), var_weights = diag(0.01, 2))
  expect_type(result@variance, "double")
  expect_true(result@variance >= 0)
})

test_that("pif_ensemble validates inputs correctly", {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic(label = "b")
  paf <- create_mock_pif_atomic(label = "paf") # Should work with PAFs too

  # Valid cases
  expect_silent(pif_ensemble(pif1, pif2))
  expect_silent(pif_ensemble(pif1, paf))
  expect_silent(pif_ensemble(pif1)) # Single PIF case


})

test_that("pif_ensemble calculates correctly", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")

  # Test two PIFs
  result <- pif_ensemble(pif1, pif2)
  expect_equal(result@pif, 1 - (1 - pif1@pif)*(1 - pif2@pif))

  # Test single PIF
  result <- pif_ensemble(pif1)
  expect_equal(result@pif, pif1@pif)

  # Test variance calculation
  result <- pif_ensemble(pif1, pif2)
  expect_type(result@variance, "double")
  expect_true(result@variance > 0)
})


test_that("paf_total validates inputs correctly", {
  pif1 <- create_mock_pif_atomic("PIF", label = "a")
  paf1 <- create_mock_pif_atomic("PAF", label = "paf1")
  paf2 <- create_mock_pif_atomic("PAF", label = "paf2")

  # Valid cases
  expect_silent(paf_total(paf1, paf2, weights = c(0.5, 0.5)))
  expect_error(
    paf_total(paf1, pif1, weights = c(0.4, 0.6)),
    "is not a Population Attributable Fraction"
  )
  expect_error(
    paf_total(pif1, paf1, weights = c(0.4, 0.6)),
    "is not a Population Attributable Fraction"
  )


})

test_that("logit link for negative pif", {

  pif1 <- pif(0.4, 0.2, 1.1, var_p = 0, var_beta = 0)
  pif2 <- pif(0.4, 0.8, 1.1, var_p = 0, var_beta = 0)

  # Valid cases
  expect_warning(
    pif_total(pif1, pif2, weights = c(0.1, 0.9), link = "logit"),
    "is only valid for strictly positive PIF"
  )

})

# Additional tests for pif_operations.R

# Test paf_ensemble function
test_that("paf_ensemble validates inputs correctly", {
  paf1 <- create_mock_pif_atomic(type = "PAF", label = "paf1")
  paf2 <- create_mock_pif_atomic(type = "PAF", label = "paf2")
  pif1 <- create_mock_pif_atomic(type = "PIF", label = "pif1")

  # Valid cases
  expect_silent(paf_ensemble(paf1, paf2))
  expect_silent(paf_ensemble(paf1)) # Single PAF case

  # Should error when mixing PIF with PAF
  expect_error(
    paf_ensemble(paf1, pif1),
    "is not a Population Attributable Fraction"
  )
})

test_that("paf_ensemble calculates correctly", {
  paf1 <- create_mock_pif_atomic(type = "PAF", label = "paf1")
  paf2 <- create_mock_pif_atomic(type = "PAF", label = "paf2")

  # Test ensemble calculation
  result <- paf_ensemble(paf1, paf2)
  expected <- 1 - (1 - paf1@pif) * (1 - paf2@pif)
  expect_equal(result@pif, expected)

  # Test single PAF
  result <- paf_ensemble(paf1)
  expect_equal(result@pif, paf1@pif)
})

# Test custom weights handling
test_that("custom weights work correctly", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")

  # Test with custom weights for ensemble
  result <- pif_ensemble(pif1, pif2, weights = c(0.7, 0.3))
  expected <- 1 - (1 - 0.7 * pif1@pif) * (1 - 0.3 * pif2@pif)
  expect_equal(result@pif, expected)

  # Test with NULL weights (should default to 1)
  result <- pif_ensemble(pif1, pif2, weights = NULL)
  expected <- 1 - (1 - pif1@pif) * (1 - pif2@pif)
  expect_equal(result@pif, expected)
})

# Test var_pif_weights parameter
test_that("var_pif_weights parameter works correctly", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")

  # Test with matrix var_pif_weights
  var_pif_weights <- matrix(c(0.01, 0.005, 0.005, 0.02), nrow = 2)
  expect_silent(
    pif_total(pif1, pif2, weights = c(0.5, 0.5), var_pif_weights = var_pif_weights)
  )

  # Test with wrong dimensions
  expect_error(
    pif_total(pif1, pif2, weights = c(0.5, 0.5),
              var_pif_weights = matrix(c(0.01, 0.005), nrow = 1)),
    "incorrect dimensions"
  )

  # Test with non-symmetric matrix
  expect_error(
    pif_total(pif1, pif2, weights = c(0.5, 0.5),
              var_pif_weights = matrix(c(0.01, 0.005, 0.003, 0.02), nrow = 2)),
    "is not symmetric"
  )

  # Test with vector (should error)
  expect_error(
    pif_total(pif1, pif2, weights = c(0.5, 0.5),
              var_pif_weights = c(0.01, 0.02)),
    "should be a matrix"
  )
})

# Test different link functions
test_that("different link functions work correctly", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")

  # Test log-complement link (default)
  result_log <- pif_total(pif1, pif2, weights = c(0.5, 0.5), link = "log-complement")
  expect_type(result_log@pif, "double")

  # Test identity link
  result_identity <- pif_total(pif1, pif2, weights = c(0.5, 0.5), link = "identity")
  expect_type(result_identity@pif, "double")

  # Test logit link
  result_logit <- pif_total(pif1, pif2, weights = c(0.5, 0.5), link = "logit")
  expect_type(result_logit@pif, "double")

  # Results should be different for different links (due to variance calculation)
  expect_false(isTRUE(all.equal(result_log@ci, result_identity@ci)))
})

# Test custom link functions
test_that("custom link functions work", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")

  # Custom square root link
  custom_link <- sqrt
  custom_link_inv <- function(x) x^2
  custom_link_deriv <- function(x) 1/(2*sqrt(x))

  expect_silent(
    pif_total(pif1, pif2, weights = c(0.5, 0.5),
              link = custom_link,
              link_inv = custom_link_inv,
              link_deriv = custom_link_deriv)
  )
})

# Test edge cases
test_that("edge cases are handled correctly", {
  pif1 <- create_mock_pif_atomic(label = "a")

  # Single PIF/PAF cases
  result_total <- pif_total(pif1, weights = 1)
  expect_equal(result_total@pif, pif1@pif)

  result_ensemble <- pif_ensemble(pif1)
  expect_equal(result_ensemble@pif, pif1@pif)

  # Zero weights case
  pif2 <- create_mock_pif_atomic(label = "b")
  result <- pif_total(pif1, pif2, weights = c(1, 0))
  expect_equal(result@pif, pif1@pif)
})

# Test label handling and duplication
test_that("label handling works correctly", {
  pif1 <- create_mock_pif_atomic(label = "same")
  pif2 <- create_mock_pif_atomic(label = "same") # Duplicate label

  # Should error on duplicate labels
  expect_error(
    pif_total(pif1, pif2, weights = c(0.5, 0.5)),
    "duplicated"
  )

  # Should work with custom label
  pif3 <- create_mock_pif_atomic(label = "different")
  expect_silent(
    pif_total(pif1, pif3, weights = c(0.5, 0.5), label = "custom_total")
  )
})

# Test variance correlation assumptions
test_that("variance correlation warnings work", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")

  # Should warn about correlation assumption with vector var_weights
  expect_warning(
    pif_total(pif1, pif2, weights = c(0.5, 0.5), var_weights = c(0.01, 0.02)),
    "are correlated but correlation is unknown"
  )
})

# Test quiet parameter
test_that("quiet parameter suppresses warnings", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")

  # Should not warn when quiet = TRUE
  expect_silent(
    pif_total(pif1, pif2, weights = c(0.5, 0.5),
              var_weights = c(0.01, 0.02), quiet = TRUE)
  )

  # Should warn when quiet = FALSE (default)
  expect_warning(
    pif_total(pif1, pif2, weights = c(0.5, 0.5),
              var_weights = c(0.01, 0.02), quiet = FALSE),
    "correlated but correlation is unknown"
  )
})

# Test confidence level parameter
test_that("confidence level parameter works", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")

  # Test different confidence levels
  result_95 <- pif_total(pif1, pif2, weights = c(0.5, 0.5), conf_level = 0.95)
  result_90 <- pif_total(pif1, pif2, weights = c(0.5, 0.5), conf_level = 0.90)

  # 90% CI should be narrower than 95% CI
  ci_width_95 <- result_95@ci[2] - result_95@ci[1]
  ci_width_90 <- result_90@ci[2] - result_90@ci[1]
  expect_true(ci_width_90 < ci_width_95)
})

# Test combination of totals and ensembles
test_that("combining totals and ensembles works", {
  # Create some base PIFs
  pif1 <- create_mock_pif_atomic(label = "base1")
  pif2 <- create_mock_pif_atomic(label = "base2")
  pif3 <- create_mock_pif_atomic(label = "base3")
  pif4 <- create_mock_pif_atomic(label = "base4")

  # Create ensembles
  ensemble1 <- pif_ensemble(pif1, pif2)
  ensemble2 <- pif_ensemble(pif3, pif4)

  # Should be able to combine ensembles into totals
  expect_silent(
    pif_total(ensemble1, ensemble2, weights = c(0.6, 0.4))
  )

  # Should be able to combine totals into ensembles
  total1 <- pif_total(pif1, pif2, weights = c(0.5, 0.5))
  total2 <- pif_total(pif3, pif4, weights = c(0.3, 0.7))

  expect_silent(
    pif_ensemble(total1, total2)
  )
})

# Test matrix var_weights handling
test_that("matrix var_weights are handled correctly", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")
  pif3 <- create_mock_pif_atomic(label = "c")

  # Create a proper covariance matrix
  var_weights <- matrix(c(0.01, 0.005, 0.002,
                          0.005, 0.02, 0.001,
                          0.002, 0.001, 0.015), nrow = 3)

  expect_silent(
    pif_total(pif1, pif2, pif3, weights = c(0.3, 0.4, 0.3),
              var_weights = var_weights)
  )
})

# Test that ensemble and total return correct class types
test_that("correct class types are returned", {
  pif1 <- create_mock_pif_atomic(label = "a")
  pif2 <- create_mock_pif_atomic(label = "b")

  result_total <- pif_total(pif1, pif2, weights = c(0.5, 0.5))
  result_ensemble <- pif_ensemble(pif1, pif2)

  expect_true(S7::S7_inherits(result_total, pif_total_class))
  expect_true(S7::S7_inherits(result_ensemble, pif_ensemble_class))

  # Both should also inherit from pif_global_ensemble_class and pif_class
  expect_true(S7::S7_inherits(result_total, pif_global_ensemble_class))
  expect_true(S7::S7_inherits(result_ensemble, pif_global_ensemble_class))
  expect_true(S7::S7_inherits(result_total, pif_class))
  expect_true(S7::S7_inherits(result_ensemble, pif_class))
})

# Test error handling for non-pif objects
test_that("non-pif objects are rejected", {
  pif1 <- create_mock_pif_atomic(label = "a")
  not_pif <- list(pif = 0.5) # Not a proper pif object

  expect_error(
    pif_total(pif1, not_pif, weights = c(0.5, 0.5)),
    "is not of.*pif_atomic_class.*or.*pif_global_ensemble_class"
  )

  expect_error(
    pif_ensemble(pif1, not_pif),
    "is not of.*pif_atomic_class.*or.*pif_global_ensemble_class"
  )
})
