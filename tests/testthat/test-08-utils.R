test_that("mu_obs_fun calculates correctly and validates inputs", {
  # Correct calculation
  expect_equal(
    mu_obs_fun(c(0.3, 0.7), c(1.1, 1.3)),
    0.3*1.1 + 0.7*1.3
  )

  expect_equal(
    mu_cft_fun(c(0.5, 0.5), c(1.5, 2.0)),
    1 + 0.5*(1.5-1) + 0.5*(2.0-1)
  )

 expect_error(mu_obs_fun(c(0.3, NA), c(1.5, 2.0)), "missing values")
 expect_error(mu_obs_fun(c(0.3, 0.7), c(1.5, NA)), "missing values")
 expect_error(mu_cft_fun(c(0.3, NA), c(1.5, 2.0)), "missing values")
 expect_error(mu_cft_fun(c(0.3, 0.7), c(1.5, NA)), "missing values")

 # Empty inputs
 expect_error(mu_obs_fun(numeric(0), c(1.5, 2.0)), "length < 0")
 expect_error(mu_obs_fun(c(0.3, 0.7), numeric(0)), "length < 0")
 expect_error(mu_cft_fun(numeric(0), c(1.5, 2.0)), "length < 0")
 expect_error(mu_cft_fun(c(0.3, 0.7), numeric(0)), "length < 0")

 # Invalid probabilities
 expect_error(mu_obs_fun(c(0.3, 0.8), c(1.5, 2.0)), "sum > 1")
 expect_error(mu_obs_fun(c(-0.1, 1.1), c(1.5, 2.0)), "Invalid probability")
 expect_error(mu_cft_fun(c(0.3, 0.8), c(1.5, 2.0)), "sum > 1")
 expect_error(mu_cft_fun(c(-0.1, 1.1), c(1.5, 2.0)), "Invalid probability")

 # Length mismatch
 expect_error(mu_obs_fun(c(0.3, 0.7), c(1.5)), "different lengths")
 expect_error(mu_cft_fun(c(0.3, 0.7), c(1.5)), "different lengths")

 # Invalid relative risks
 expect_error(mu_obs_fun(c(0.3, 0.7), c(1.5, 0)), "which is <= 0")
 expect_error(mu_obs_fun(c(0.3, 0.7), c(1.5, -1)), "which is <= 0")
 expect_error(mu_cft_fun(c(0.3, 0.7), c(1.5, 0)), "which is <= 0")
 expect_error(mu_cft_fun(c(0.3, 0.7), c(1.5, -1)), "which is <= 0")
})



 test_that("pif_fun calculates correctly", {
   # Test basic calculation
   expect_equal(
     pif_fun(c(0.3, 0.7), c(0.2, 0.8), c(1.5, 2.0)),
     1 - (mu_cft_fun(c(0.2, 0.8), c(1.5, 2.0)) / mu_obs_fun(c(0.3, 0.7), c(1.5, 2.0)))
   )

   # Edge case where observed = counterfactual (PIF = 0)
   expect_equal(
     pif_fun(c(0.3, 0.7), c(0.3, 0.7), c(1.5, 2.0)),
     0
   )


 })

 test_that("pif_fun2 calculates correctly and validates inputs", {
   # Correct calculation
   expect_equal(pif_fun2(1.5, 1.2), 1 - (1.2/1.5))
   expect_equal(pif_fun2(2.0, 1.0), 0.5)

   # Input validation
   # NA values
   expect_error(pif_fun2(NA, 1.5), "Missing input")
   expect_error(pif_fun2(1.5, NA), "Missing input")

   # Length > 1
   expect_error(pif_fun2(c(1.5, 2.0), 1.2), "Invalid length")
   expect_error(pif_fun2(1.5, c(1.2, 1.0)), "Invalid length")

   # Edge cases
   expect_equal(pif_fun2(1.0, 1.0), 0) # No difference
   expect_equal(pif_fun2(1.0, 0.0), 1) # Complete elimination
 })

 test_that("pif_atomic_ci calculates correctly and validates inputs", {
   # Correct calculation with identity link
   link_val <- 0.2
   link_var <- 0.01
   conf <- 0.95
   z <- qnorm((1 - conf)/2)
   expected <- c(
     link_val - z * sqrt(link_var),
     link_val + z * sqrt(link_var)
   )

   expect_equal(
     pif_atomic_ci(link_val, link_var, conf, identity),
     sort(expected)
   )

   # With logit link
   pif_val <- 0.3
   link_val <- logit(pif_val)
   ci <- pif_atomic_ci(link_val, link_var, conf, inv_logit)
   expect_true(ci[1] < ci[2])
   expect_true(ci[1] > 0 && ci[2] < 1) # Should be within (0,1)

   # Input validation
   # NA values
   expect_error(pif_atomic_ci(NA, 0.01, 0.95, identity), "Missing values")
   expect_error(pif_atomic_ci(0.2, NA, 0.95, identity), "Missing values")

   # Length > 1
   expect_error(pif_atomic_ci(c(0.2, 0.3), 0.01, 0.95, identity), "Invalid length")
   expect_error(pif_atomic_ci(0.2, c(0.01, 0.02), 0.95, identity), "Invalid length")

   # Invalid confidence level
   expect_error(pif_atomic_ci(0.2, 0.01, 1.1, identity), "Invalid confidence level")
   expect_error(pif_atomic_ci(0.2, 0.01, -0.1, identity), "Invalid confidence level")
 })

 test_that("pif_class_apply_1st works correctly", {
   # Create mock S7 objects
   pif_atomic <- pif_atomic_class(p = c(0.1, 0.2), p_cft = rep(0, 2), beta = c(1.1, 1.2),
                                  conf_level = 0.95, type = "PIF", link = logit,
                                  link_inv = inv_logit, link_deriv = deriv_logit,
                                  var_p = matrix(0, 2, 2), var_beta = matrix(0, 2, 2),
                                  rr_link = exp, rr_link_deriv = exp, upper_bound_p = FALSE,
                                  upper_bound_beta = FALSE)
   pif_total <- pif_total_class(pif_list = list(pif_atomic, pif_atomic), link = logit,
                                link_inv = inv_logit, link_deriv = deriv_logit,
                                weights = c(0.5, 0.5), sigma_weights = matrix(0, 2, 2))

   # Test with atomic pif
   expect_equal(
     pif_class_apply_1st(pif_atomic, identity, "pif"),
     pif_atomic@pif
   )

   # Test with total pif (should use first element)
   expect_equal(
     pif_class_apply_1st(pif_total, identity, "pif"),
     pif_total@pif_list[[1]]@pif
   )

   # Test invalid class
   expect_error(
     pif_class_apply_1st(list(), identity, "pif"),
     "Invalid class"
   )
 })

 test_that("link_deriv_vals works correctly", {
   pif_obj <- pif_atomic_class(p = c(0.1, 0.2), p_cft = rep(0, 2), beta = c(1.1, 1.2),
                                  conf_level = 0.95, type = "PIF", link = exp,
                                  link_inv = exp, link_deriv = exp,
                                  var_p = matrix(0, 2, 2), var_beta = matrix(0, 2, 2),
                                  rr_link = exp, rr_link_deriv = exp, upper_bound_p = FALSE,
                                  upper_bound_beta = FALSE)

   expect_equal(link_deriv_vals(pif_obj), exp(pif_obj@pif))

   # Test invalid input
   expect_error(link_deriv_vals(list()), "Invalid class")
 })

 test_that("fraction_type works correctly", {
   pif_obj <- pif_atomic_class(p = c(0.1, 0.2), p_cft = rep(0, 2), beta = c(1.1, 1.2),
                                  conf_level = 0.95, type = "PIF", link = logit,
                                  link_inv = inv_logit, link_deriv = deriv_logit,
                                  var_p = matrix(0, 2, 2), var_beta = matrix(0, 2, 2),
                                  rr_link = exp, rr_link_deriv = exp, upper_bound_p = FALSE,
                                  upper_bound_beta = FALSE)

   paf_obj <- pif_atomic_class(p = c(0.1, 0.2), p_cft = rep(0, 2), beta = c(1.1, 1.2),
                                  conf_level = 0.95, type = "PAF", link = logit,
                                  link_inv = inv_logit, link_deriv = deriv_logit,
                                  var_p = matrix(0, 2, 2), var_beta = matrix(0, 2, 2),
                                  rr_link = exp, rr_link_deriv = exp, upper_bound_p = FALSE,
                                  upper_bound_beta = FALSE)


   expect_equal(fraction_type(pif_obj), "PIF")
   expect_equal(fraction_type(paf_obj), "PAF")

   # Test invalid input
   expect_error(fraction_type(list()), "Invalid class")
 })

 test_that("change_link works correctly", {
   pif_obj <- pif_atomic_class(p = c(0.1, 0.2), p_cft = rep(0, 2), beta = c(1.1, 1.2),
                               conf_level = 0.95, type = "PIF", link = identity,
                               link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                               var_p = matrix(0, 2, 2), var_beta = matrix(0, 2, 2),
                               rr_link = exp, rr_link_deriv = exp, upper_bound_p = FALSE,
                               upper_bound_beta = FALSE)

   # Change to logit
   new_pif <- change_link(pif_obj, "logit")
   expect_equal(new_pif@link, logit)
   expect_equal(new_pif@link_inv, inv_logit)
   expect_equal(new_pif@link_deriv, deriv_logit)

   # Change with custom functions
   custom_link <- function(x) x^2
   custom_inv <- function(x) sqrt(x)
   custom_deriv <- function(x) 2*x
   new_pif <- change_link(pif_obj, custom_link, custom_inv, custom_deriv)
   expect_equal(new_pif@link, custom_link)
   expect_equal(new_pif@link_inv, custom_inv)
   expect_equal(new_pif@link_deriv, custom_deriv)

   # Test invalid inputs
   expect_error(change_link(list(), "logit"), "Invalid class")
   expect_error(change_link(pif_obj, "logit", inv_logit), "`link_inv` was given")

   # Test warning for logit with PIF <= 0
   neg_pif <- pif_atomic_class(p = c(0.1, 0.2), p_cft = rep(0, 2), beta = c(0.1, 0.2),
                    conf_level = 0.95, type = "PIF", link = identity,
                    link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                    var_p = matrix(0, 2, 2), var_beta = matrix(0, 2, 2),
                    rr_link = identity, rr_link_deriv = function(x) rep(1, length(x)), upper_bound_p = FALSE,
                    upper_bound_beta = FALSE)
   expect_warning(change_link(neg_pif, "logit"), "<= 0")
 })

