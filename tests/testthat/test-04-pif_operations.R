# Helper function to create mock PIF objects
create_mock_pif_atomic <- function(type = "PIF") {
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
  paf1 <- create_mock_pif_atomic(type = "PAF")
  paf2 <- create_mock_pif_atomic(type = "PAF")

  # Valid case
  expect_silent(paf_total(paf1, paf2, weights = c(0.5, 0.5)))

})

test_that("pif_total validates inputs correctly", {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic()

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

  # Invalid var_weights
  expect_error(
    pif_total(pif1, pif2, weights = c(0.5, 0.5), var_weights = matrix(1:4, nrow = 2)),
    "is not symmetric"
  )

  expect_error(
    pif_total(pif1, pif2, weights = c(0.5, 0.5), var_weights = "a"),
    "should be a number"
  )

  expect_error(
    suppressWarnings(
      pif_total(pif1, pif2, weights = c(0.5, 0.5), var_weights = 1:3)
    ),
    "has incorrect dimensions"
  )

  # Test link function handling
  expect_silent(pif_total(pif1, weights = 1, link = "logit"))
  expect_silent(pif_total(pif1, weights = 1, link = logit, link_inv = inv_logit, link_deriv = deriv_logit))

})

test_that("pif_total calculates correctly", {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic()

  # Test weights
  result <- pif_total(pif1, pif2, weights = c(0.7, 0.3))
  expect_equal(result@pif, 0.7*pif1@pif + 0.3*pif2@pif)

  # Test single PIF
  result <- pif_total(pif1, weights = 1)
  expect_equal(result@pif, pif1@pif)

  # Test variance calculation
  result <- pif_total(pif1, pif2, weights = c(0.5, 0.5), var_weights = diag(0.01, 2))
  expect_type(result@variance, "double")
  expect_true(result@variance > 0)
})

test_that("pif_ensemble validates inputs correctly", {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic()
  paf <- create_mock_pif_atomic() # Should work with PAFs too

  # Valid cases
  expect_silent(pif_ensemble(pif1, pif2))
  expect_silent(pif_ensemble(pif1, paf))
  expect_silent(pif_ensemble(pif1)) # Single PIF case


})

test_that("pif_ensemble calculates correctly", {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic()

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
  pif1 <- create_mock_pif_atomic("PIF")
  paf1 <- create_mock_pif_atomic("PAF")
  paf2 <- create_mock_pif_atomic("PAF")

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
