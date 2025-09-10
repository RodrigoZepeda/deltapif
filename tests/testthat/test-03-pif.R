test_that("paf correctly delegates to pif with zero counterfactual", {
  # Test that paf calls pif with p_cft = 0
  expect_equal(
    paf(p = 0.5, beta = 1.5, quiet = TRUE, label = "1"),
    pif(p = 0.5, p_cft = 0, beta = 1.5, type = "PAF", quiet = TRUE, label = "1")
  )

  # Test variance parameters are passed through
  expect_equal(
    paf(p = 0.5, beta = 1.5, var_p = 0.01, var_beta = 0.02, quiet = TRUE, label = "paf"),
    pif(p = 0.5, p_cft = 0, beta = 1.5, var_p = 0.01, var_beta = 0.02, type = "PAF", quiet = TRUE, label = "paf")
  )

  # Test link functions are passed through
  expect_equal(
    paf(p = 0.5, beta = 1.5, link = "logit", quiet = TRUE, label = "test"),
    pif(p = 0.5, p_cft = 0, beta = 1.5, link = "logit", type = "PAF", quiet = TRUE, label = "test")
  )

})

test_that("pif validates inputs correctly", {
  # Valid cases
  expect_silent(pif(p = 0.5, p_cft = 0.2, beta = 1.5, quiet = TRUE))
  expect_silent(pif(p = c(0.3, 0.7), p_cft = c(0.2, 0.8), beta = c(1.5, 2.0),
                    var_p = matrix(0, 2, 2), var_beta = matrix(0, 2, 2)))

  # Invalid probabilities
  expect_error(pif(p = -0.1, p_cft = 0.2, beta = 1.5, quiet = TRUE), "has values < 0")
  expect_error(pif(p = 1.1, p_cft = 0.2, beta = 1.5, quiet = TRUE), "values > 1")
  expect_error(pif(p = 0.5, p_cft = -0.1, beta = 1.5, quiet = TRUE), "has values < 0")
  expect_error(pif(p = 0.5, p_cft = 1.2, beta = 1.5, quiet = TRUE), "values > 1")
  expect_error(pif(p = c(0.3, 0.8), p_cft = c(0.2, 0.8), beta = c(1.5, 2.0), quiet = TRUE), "sum.* > 1")

  # Length mismatches
  expect_error(
    pif(p = c(0.3, 0.7), p_cft = c(0.2), beta = c(1.5, 2.0), quiet = TRUE)
  )
  expect_error(
    pif(p = c(0.3, 0.7), p_cft = c(0.2, 0.8), beta = 1.5, quiet = TRUE),
    "different lengths"
  )

  # Invalid type specification
  expect_warning(
     pif(p = 0.5, p_cft = 0.2, beta = 1.5, type = "PAF", var_p = 0, var_beta = 0),
     "Are you sure you are not"
  )
})

test_that("pif handles variance inputs correctly", {
  # No variance case
  expect_silent(pif(p = 0.5, p_cft = 0.2, beta = 1.5, quiet = TRUE))
  expect_message(pif(p = 0.5, p_cft = 0.2, beta = 1.5, var_p = 0), "have no variance")
  expect_message(pif(p = 0.5, p_cft = 0.2, beta = 1.5, var_beta = 0), "have no variance")

  #Has two messages
  expect_message(
    expect_message(pif(p = 0.5, p_cft = 0.2, beta = 1.5), "have no variance"),
    "have no variance"
  )

  # Scalar variance
  expect_silent(pif(p = 0.5, p_cft = 0.2, beta = 1.5, var_p = 0.01, var_beta = 0.02))

  # Vector variance (should trigger upper bound warning)
  expect_error(
    pif(p = 0.5, p_cft = 0.2, beta = 1.5, var_p = c(0.01, 0.02), quiet = TRUE),
    "has different length than its covariance matrix"
  )

  # Matrix variance
  expect_silent(
    pif(p = c(0.3, 0.7), p_cft = c(0.2, 0.8), beta = c(1.5, 2.0),
        var_p = matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2),
        var_beta = matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2))
  )

  # Invalid variance matrices
  expect_message(
    pif(p = c(0.3, 0.7), p_cft = c(0.2, 0.8), beta = c(1.5, 2.0),
        var_p = matrix(1:4, nrow = 2), quiet = TRUE),
    "is not symmetric"
  )
  expect_message(
    pif(p = c(0.3, 0.7), p_cft = c(0.2, 0.8), beta = c(1.5, 2.0),
        var_beta = matrix(1:4, nrow = 2), quiet = TRUE),
    "is not symmetric"
  )
})

test_that("pif handles link functions correctly", {
  # Default log-complement link
  result <- pif(p = 0.5, p_cft = 0.2, beta = 1.5, quiet = TRUE)
  expect_equal(result@link, log_complement)
  expect_equal(result@link_inv, inv_log_complement)
  expect_equal(result@link_deriv, deriv_log_complement)

  # Character-specified links
  expect_equal(pif(p = 0.5, p_cft = 0.2, beta = 1.5, link = "logit", quiet = TRUE)@link, logit)
  expect_equal(pif(p = 0.5, p_cft = 0.2, beta = 1.5, link = "identity", quiet = TRUE)@link, identity)
  expect_equal(pif(p = 0.5, p_cft = 0.2, beta = 1.5, link = "hawkins", quiet = TRUE)@link, hawkins)

  # Custom function links
  custom_link <- function(x) x^2
  custom_inv <- function(x) sqrt(x)
  custom_deriv <- function(x) 2*x
  result <- pif(p = 0.5, p_cft = 0.2, beta = 1.5, quiet = TRUE,
                link = custom_link, link_inv = custom_inv, link_deriv = custom_deriv)
  expect_equal(result@link, custom_link)
  expect_equal(result@link_inv, custom_inv)
  expect_equal(result@link_deriv, custom_deriv)
  expect_equal(result@link_vals, result@pif^2)

})

test_that("pif handles rr_link functions correctly", {
  # Default identity link
  result <- pif(p = 0.5, p_cft = 0.2, beta = 1.5, quiet = TRUE, rr_link = "identity")
  expect_equal(result@rr_link, identity)

  result <- pif(p = 0.5, p_cft = 0.2, beta = 1.5, quiet = TRUE, rr_link = "exp")
  expect_equal(result@rr_link, exp)

  # Custom function links
  custom_rr_link <- function(x) exp(x)
  custom_rr_deriv <- function(x) exp(x)
  result <- pif(p = 0.5, p_cft = 0.2, beta = 1.5, quiet = TRUE,
                rr_link = custom_rr_link, rr_link_deriv = custom_rr_deriv)
  expect_equal(result@rr_link, custom_rr_link)
  expect_equal(result@rr_link_deriv, custom_rr_deriv)

})


test_that("pif warns about logit with non-positive PIF", {
  expect_warning(
    pif(p = 0.5, p_cft = 0.6, beta = 1.0, link = "logit", quiet = TRUE),
    "<= 0"
  )
})

test_that("pif handles quiet parameter correctly", {
  # Should show warnings when quiet = FALSE
  expect_message(
    pif(p = 0.5, p_cft = 0.2, beta = 1.5, var_p = NULL, var_beta = 1, quiet = FALSE),
    "Assuming parameters `p` have no variance"
  )

  expect_message(
    pif(p = 0.5, p_cft = 0.2, beta = 1.5, var_p = 1, var_beta = NULL, quiet = FALSE),
    "Assuming parameters `beta` have no variance"
  )

  # Should suppress warnings when quiet = TRUE
  expect_silent(
    pif(p = 0.5, p_cft = 0.2, beta = 1.5, var_p = NULL, quiet = TRUE)
  )
})

test_that("pif works for var_p options",{
  expect_equal(
    pif(p = 0.5, p_cft = 0.2, beta = 1.5,
        var_p = as_covariance_structure(matrix(0.21)), var_beta = NULL, quiet = TRUE,
        label = "pif"),
    pif(p = 0.5, p_cft = 0.2, beta = 1.5,
        var_p = matrix(0.21), var_beta = NULL, quiet = TRUE, label = "pif")
  )

  expect_equal(
    pif(p = 0.5, p_cft = 0.2, beta = 1.5,
        var_p = NULL, var_beta = as_covariance_structure(matrix(0.21)), quiet = TRUE,
        label = "pif"),
    pif(p = 0.5, p_cft = 0.2, beta = 1.5,
        var_p = NULL, var_beta = as_covariance_structure(matrix(0.21)), quiet = TRUE, label = "pif")
  )

  expect_equal(
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_p = NULL, var_beta = matrix(0, ncol = 2, nrow = 2), quiet = TRUE,
        label = "pif"),
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_p = NULL, var_beta = 0, quiet = TRUE, label = "pif")
  )

  expect_equal(
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_beta = NULL, var_p = matrix(0, ncol = 2, nrow = 2), quiet = TRUE,
        label = "pif"),
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_beta = NULL, var_p = 0, quiet = TRUE, label = "pif")
  )

  expect_equal(
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_p = NULL, var_beta = matrix(0, ncol = 2, nrow = 2), quiet = TRUE,
        label = "pif"),
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_p = NULL, var_beta = c(0, 0), quiet = TRUE, label = "pif")
  )

  expect_equal(
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_beta = NULL, var_p = matrix(0, ncol = 2, nrow = 2), quiet = TRUE,
        label = "pif"),
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_beta = NULL, var_p = c(0, 0), quiet = TRUE, label = "pif")
  )

  var_p <- c(0.1, 0.3)
  expect_equal(
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_beta = NULL, var_p = var_p, quiet = TRUE,
        label = "pif"),
    pif(p = c(0.5, 0.3), p_cft = c(0.2, 0.2), beta = c(1.5, 1.8),
        var_beta = NULL, var_p = sqrt(var_p %*% t(var_p)),
        quiet = TRUE, label = "pif")
  )

})

