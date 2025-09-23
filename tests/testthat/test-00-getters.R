# Helper to create mock pif_atomic_class object
create_mock_pif_atomic <- function(type = "PIF", label = paste0("atomic_pif", stats::rnorm(1))) {
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
    type = type,
    label = label
  )
}

# Helper to create mock pif_class object
create_mock_pif <- function() {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic()
  piflist <- list(pif1, pif2)
  names(piflist) <- c(pif1@label, pif2@label)
  pif_total_class(
    pif_list = piflist,
    weights = c(0.5, 0.5),
    var_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
    link = logit,
    link_deriv = deriv_logit,
    link_inv = inv_logit,
    var_pif_weights = matrix(c(0.01, -0.2, -0.2, 0.01), nrow = 2),
    label = paste0("total", stats::rnorm(1))
  )
}

test_that("Basic getters work correctly with S7 objects", {

  mock_pif <- create_mock_pif_atomic()

  # Test get_rr
  expect_equal(get_rr(mock_pif), exp(mock_pif@beta))

  # Test get_mu_obs
  expect_equal(get_mu_obs(mock_pif), mock_pif@mu_obs)
  expect_equal(get_mu_obs(mock_pif),
               mu_obs_fun(mock_pif@p, mock_pif@rr))

  # Test get_mu_cft
  expect_equal(get_mu_cft(mock_pif), mock_pif@mu_cft)
  expect_equal(get_mu_cft(mock_pif),
               mu_cft_fun(mock_pif@p_cft, mock_pif@rr))

  # Test get_pif
  expect_equal(get_pif(mock_pif), mock_pif@pif)
  expect_equal(get_pif(mock_pif),
               pif_fun2(mock_pif@mu_obs, mock_pif@mu_cft))

  # Test get_link_vals
  expect_equal(get_link_vals(mock_pif), mock_pif@link_vals)
  expect_equal(get_link_vals(mock_pif), logit(mock_pif@pif))

  # Test get_link_deriv_vals
  expect_equal(get_link_deriv_vals(mock_pif), mock_pif@link_deriv_vals)
  expect_equal(get_link_deriv_vals(mock_pif), mock_pif@link_deriv(mock_pif@pif))
  expect_equal(get_link_deriv_vals(mock_pif),
               deriv_logit(mock_pif@pif))

  # Test get_rr_link_deriv_vals
  expect_equal(get_rr_link_deriv_vals(mock_pif), mock_pif@rr_link_deriv_vals)
  expect_equal(get_rr_link_deriv_vals(mock_pif), exp(mock_pif@beta))
  expect_equal(get_rr_link_deriv_vals(mock_pif), mock_pif@rr_link_deriv(mock_pif@beta))
  expect_equal(get_rr_link_deriv_vals(mock_pif), exp(mock_pif@beta))

})

test_that("Variance and CI getters work with S7 objects", {
  mock_pif <- create_mock_pif_atomic()

  # Test get_link_variance
  expected_var <- (mock_pif@link_deriv_vals)^2 * mock_pif@variance
  expect_equal(get_link_variance(mock_pif), expected_var)

  # Test get_variance_atomic
  # This will depend on your from_parameters_pif_variance implementation
  expect_type(get_variance_atomic(mock_pif), "double")

  # Test get_ci
  ci <- get_ci(mock_pif)
  expect_length(ci, 2)
  expect_true(ci[1] < ci[2])
})

test_that("Composite mock_pif getters work with S7 objects", {
  mock_pif <- create_mock_pif()

  # Test get_total_pif
  expect_equal(get_total_pif(mock_pif),
               as.numeric(t(mock_pif@weights) %*% mock_pif@coefs))

  # Test get_ensemble_pif
  expect_equal(get_ensemble_pif(mock_pif),
               1 - exp(sum(log(1 - 0.5*mock_pif@coefs))))
})

test_that("Covariance getters work with S7 objects", {
  mock_pif <- pif_total(
    paf(0.1, 1.2, var_p = 0.01, var_beta = 0, link = logit,
        link_deriv = deriv_logit, link_inv = inv_logit,
        label = "1",
        rr_link = identity, rr_link_deriv = function(x) {1}), #Github actions crashes if not specifying compeltely this
    paf(0.1, 1.2, var_p = 0.01, var_beta = 0, link = logit,
        link_deriv = deriv_logit, link_inv = inv_logit,
        label = "2",
        rr_link = identity, rr_link_deriv = function(x) {1}),
    weights = c(0.5, 0.5),
    link = logit,
    link_inv = inv_logit,
    link_deriv = deriv_logit,
    label = "total"
  )

  mock_ensemble <- pif_ensemble(
    paf(0.1, 1.2, var_p = 0.01, var_beta = 0, link = logit,
        link_deriv = deriv_logit, link_inv = inv_logit,
        label = "1",
        rr_link = identity, rr_link_deriv = function(x) {1}), #Github actions crashes if not specifying compeltely this
    paf(0.1, 1.2, var_p = 0.01, var_beta = 0, link = logit,
        link_deriv = deriv_logit, link_inv = inv_logit,
        label = "2",
        rr_link = identity, rr_link_deriv = function(x) {1})
  )

  # Test only work with total pif not with mock
  expect_equal(get_ensemble_coefs(mock_pif),
               sapply(mock_ensemble@pif_list, function(x) x@pif))


  # Test get_covariance_total
  cov_mat <- get_covariance(mock_pif)
  expect_equal(dim(cov_mat), c(2, 2))
  expect_true(isSymmetric(cov_mat))

  # Test get_ensemble_covariance
  ens_cov <- get_covariance(mock_ensemble)
  expect_equal(dim(ens_cov), c(2, 2))
})


test_that("Input validation works with S7 objects", {
  expect_error(get_rr("not a mock_pif object"))
  expect_error(get_pif(list()))
})


test_that("Types for totals", {
  pif1 <- create_mock_pif_atomic("PIF")
  pif2 <- create_mock_pif_atomic("PIF")
  paf1 <- create_mock_pif_atomic("PAF")
  paf2 <- create_mock_pif_atomic("PAF")

  piflist <- list(pif1, pif2)
  names(piflist) <- c(pif1@label, pif2@label)

  tp <- pif_total_class(piflist, link = identity, label = "total_class",
                        var_pif_weights = matrix(0, 2, 2),
                        link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                        weights = c(0.5, 0.5), var_weights = matrix(0, 2, 2))

  ep <- pif_ensemble_class(piflist, link = identity, label = "ensemble_class",
                           var_pif_weights = matrix(0, 2, 2),
                           link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                           weights = rep(1, 2), var_weights = matrix(0, 2, 2))

  expect_equal(get_ensemble_type(tp), "PIF")
  expect_equal(get_ensemble_type(ep), "PIF")

  paflist <- list(paf1, pif2)
  names(paflist) <- c(paf1@label, pif2@label)


  tp <- pif_total_class(paflist, link = identity,
                        label = "total_class",
                        var_pif_weights = matrix(0, 2, 2),
                        link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                        weights = c(0.5, 0.5), var_weights = matrix(0, 2, 2))

  ep <- pif_ensemble_class(paflist, link = identity,
                           label = "ensemble_class",
                           var_pif_weights = matrix(0, 2, 2),
                           link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                           weights = rep(1, 2), var_weights = matrix(0, 2, 2))


  expect_equal(get_ensemble_type(tp), "PIF")
  expect_equal(get_ensemble_type(ep), "PIF")

  papaflist <- list(paf1, paf2)
  names(papaflist) <- c(paf1@label, paf2@label)


  tp <- pif_total_class(papaflist, link = identity,
                        label = "total_class",
                        var_pif_weights = matrix(0, 2, 2),
                        link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                        weights = c(0.5, 0.5), var_weights = matrix(0, 2, 2))

  ep <- pif_ensemble_class(papaflist, link = identity,
                           label = "ensemble_class",
                           var_pif_weights = matrix(0, 2, 2),
                           link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                           weights = rep(1, 2), var_weights = matrix(0, 2, 2))


  expect_equal(get_ensemble_type(tp), "PAF")
  expect_equal(get_ensemble_type(ep), "PAF")
})



# Test edge cases and error conditions for basic getters
test_that("getter functions handle edge cases correctly", {
  # Test with minimal variance
  mock_pif_zero_var <- pif_atomic_class(
    p = c(0.5),
    p_cft = c(0.2),
    beta = c(0.5),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = identity,
    link_deriv = function(x) rep(1, length(x)),
    link_inv = identity,
    var_p = matrix(0, nrow = 1),
    var_beta = matrix(0, nrow = 1),
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE,
    conf_level = 0.95,
    type = "PIF",
    label = "zero_var_test"
  )

  expect_equal(get_variance_atomic(mock_pif_zero_var), 0)
  expect_length(get_ci(mock_pif_zero_var), 2)

  # Test with single exposure category
  single_exposure <- create_mock_pif_atomic()
  expect_type(get_rr(single_exposure), "double")
  expect_true(all(get_rr(single_exposure) > 0))
})

# Test getters with different link functions
test_that("getters work correctly with different link functions", {
  # Test with logit link
  mock_pif_logit <- pif_atomic_class(
    p = c(0.3, 0.4),
    p_cft = c(0.1, 0.2),
    beta = c(0.2, 0.5),
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
    type = "PIF",
    label = "logit_test"
  )

  # Test link transformation
  link_vals <- get_link_vals(mock_pif_logit)
  expect_equal(link_vals, logit(mock_pif_logit@pif))

  # Test link derivative
  link_deriv_vals <- get_link_deriv_vals(mock_pif_logit)
  expect_equal(link_deriv_vals, deriv_logit(mock_pif_logit@pif))

  # Test with log-complement link
  mock_pif_log_comp <- mock_pif_logit
  mock_pif_log_comp@link <- log_complement
  mock_pif_log_comp@link_deriv <- deriv_log_complement
  mock_pif_log_comp@link_inv <- inv_log_complement

  expect_equal(get_link_vals(mock_pif_log_comp), log_complement(mock_pif_log_comp@pif))
})

# Test getter functions with multiple exposure categories
test_that("getters handle multiple exposure categories correctly", {
  multi_exposure <- pif_atomic_class(
    p = c(0.2, 0.3, 0.1),
    p_cft = c(0.1, 0.15, 0.05),
    beta = c(0.3, 0.6, 1.2),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = identity,
    link_deriv = function(x) rep(1, length(x)),
    link_inv = identity,
    var_p = diag(0.01, 3),
    var_beta = diag(0.02, 3),
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE,
    conf_level = 0.95,
    type = "PIF",
    label = "multi_exposure"
  )

  # Test that RR vector has correct length
  rr_vals <- get_rr(multi_exposure)
  expect_length(rr_vals, 3)
  expect_true(all(rr_vals > 0))

  # Test mu calculations
  mu_obs <- get_mu_obs(multi_exposure)
  mu_cft <- get_mu_cft(multi_exposure)
  expect_type(mu_obs, "double")
  expect_type(mu_cft, "double")
  expect_true(mu_obs > 0)
  expect_true(mu_cft > 0)

  # Test PIF calculation consistency
  pif_val <- get_pif(multi_exposure)
  expected_pif <- (mu_obs - mu_cft) / mu_obs
  expect_equal(pif_val, expected_pif)
})

# Test cases_class getters
test_that("cases_class getters work correctly", {
  base_pif <- create_mock_pif_atomic()

  cases_obj <- cases_class(
    overall_cases = 1000,
    pif_obj = base_pif,
    variance_cases = 100,
    link = identity,
    link_deriv = function(x) rep(1, length(x)),
    link_inv = identity,
    conf_level = 0.95
  )

  # Test cases calculation
  expected_cases <- 1000 * base_pif@pif
  expect_equal(get_cases(cases_obj), expected_cases)

  # Test variance calculation
  variance_val <- get_variance_cases(cases_obj)
  expect_type(variance_val, "double")
  expect_true(variance_val >= 0)

  # Test link functions for cases
  link_vals_cases <- get_link_vals_cases(cases_obj)
  expect_equal(link_vals_cases, identity(cases_obj@cases))

  link_deriv_vals_cases <- get_link_deriv_vals_cases(cases_obj)
  expect_equal(link_deriv_vals_cases, 1)  # identity derivative
})

# Test ensemble and global getters more thoroughly
test_that("ensemble getters handle complex structures", {
  # Create a more complex ensemble structure
  pif1 <- create_mock_pif_atomic(label = "component1")
  pif2 <- create_mock_pif_atomic(label = "component2")
  pif3 <- create_mock_pif_atomic(label = "component3")

  # Test total with three components
  total_pif <- pif_total(pif1, pif2, pif3, weights = c(0.2, 0.5, 0.3))

  coefs <- get_ensemble_coefs(total_pif)
  expect_length(coefs, 3)
  expect_equal(as.vector(coefs), c(pif1@pif, pif2@pif, pif3@pif))

  # Test type determination with mixed types
  paf1 <- create_mock_pif_atomic(type = "PAF", label = "paf_comp")
  mixed_total <- pif_total(pif1, paf1, weights = c(0.6, 0.4))
  expect_equal(get_ensemble_type(mixed_total), "PIF")  # Should be PIF if any component is PIF

  # Test pure PAF ensemble
  paf2 <- create_mock_pif_atomic(type = "PAF", label = "paf_comp2")
  paf_total <- pif_total(paf1, paf2, weights = c(0.5, 0.5))
  expect_equal(get_ensemble_type(paf_total), "PAF")
})

# Test variance and covariance calculations
test_that("variance and covariance getters work correctly", {
  # Create PIFs with known variance structures
  pif1 <- paf(0.2, 1.5, var_p = 0.01, var_beta = 0.05, quiet = TRUE, label = "var_test1")
  pif2 <- paf(0.3, 1.8, var_p = 0.02, var_beta = 0.08, quiet = TRUE, label = "var_test2")

  total_pif <- pif_total(pif1, pif2, weights = c(0.4, 0.6),
                         var_weights = matrix(c(0.001, 0, 0, 0.001), nrow = 2))

  # Test covariance matrix
  cov_matrix <- get_covariance(total_pif)
  expect_true(is.matrix(cov_matrix))
  expect_equal(dim(cov_matrix), c(2, 2))
  expect_true(isSymmetric(cov_matrix))

  # Test variance calculation
  variance_val <- get_variance(total_pif)
  expect_type(variance_val, "double")
  expect_true(variance_val >= 0)

  # Test single component case
  single_total <- pif_total(pif1, weights = 1)
  single_cov <- get_covariance(single_total)
  expect_equal(single_cov, variance(pif1))
})

# Test transform functions for ensembles
test_that("transform getters work correctly", {
  pif1 <- create_mock_pif_atomic(label = "trans1")
  pif2 <- create_mock_pif_atomic(label = "trans2")

  # Test ensemble with log-complement transform
  ensemble_pif <- pif_ensemble(pif1, pif2, weights = c(0.7, 0.3))

  sum_transformed <- get_sum_transformed_weighted_coefs(ensemble_pif)
  expected_sum <- sum(sapply(c(0.7 * pif1@pif, 0.3 * pif2@pif), log_complement))
  expect_equal(sum_transformed, expected_sum)

  # Test global ensemble PIF
  global_pif <- get_global_ensemble_pif(ensemble_pif)
  expected_global <- inv_log_complement(sum_transformed)
  expect_equal(global_pif, expected_global)

  # Test direct ensemble calculation
  ensemble_calc <- get_ensemble_pif(ensemble_pif)
  expected_ensemble <- 1 - (1 - 0.7 * pif1@pif) * (1 - 0.3 * pif2@pif)
  expect_equal(ensemble_calc, expected_ensemble)
})

# Test getter robustness with extreme values
test_that("getters handle extreme values appropriately", {
  # Test with very small PIF values
  small_pif <- pif_atomic_class(
    p = c(0.01),
    p_cft = c(0.005),
    beta = c(0.01),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = identity,
    link_deriv = function(x) rep(1, length(x)),
    link_inv = identity,
    var_p = matrix(0.0001, nrow = 1),
    var_beta = matrix(0.0001, nrow = 1),
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE,
    conf_level = 0.95,
    type = "PIF",
    label = "small_pif"
  )

  expect_true(get_pif(small_pif) >= 0)
  expect_true(get_mu_obs(small_pif) > 0)
  expect_true(get_mu_cft(small_pif) > 0)

  # Test confidence intervals are ordered correctly
  ci <- get_ci(small_pif)
  expect_true(ci[1] <= ci[2])
})

# Test error handling in getters
test_that("getters provide appropriate error messages", {
  # Test with invalid objects
  expect_error(get_rr("not_a_pif"))
  expect_error(get_pif(list()))
  expect_error(get_mu_obs(42))

  # Test with NULL objects
  expect_error(get_rr(NULL))
  expect_error(get_pif(NULL))
})

# Test consistency between different getter approaches
test_that("getter calculations are internally consistent", {
  mock_pif <- create_mock_pif_atomic()

  # Test that manual calculation matches getter
  manual_mu_obs <- mu_obs_fun(mock_pif@p, get_rr(mock_pif))
  getter_mu_obs <- get_mu_obs(mock_pif)
  expect_equal(manual_mu_obs, getter_mu_obs)

  manual_mu_cft <- mu_cft_fun(mock_pif@p_cft, get_rr(mock_pif))
  getter_mu_cft <- get_mu_cft(mock_pif)
  expect_equal(manual_mu_cft, getter_mu_cft)

  manual_pif <- pif_fun2(getter_mu_obs, getter_mu_cft)
  getter_pif <- get_pif(mock_pif)
  expect_equal(manual_pif, getter_pif)

  # Test link consistency
  manual_link_vals <- mock_pif@link(getter_pif)
  getter_link_vals <- get_link_vals(mock_pif)
  expect_equal(manual_link_vals, getter_link_vals)
})

# Test getters with different confidence levels
test_that("confidence interval getters respect conf_level", {
  mock_pif_90 <- create_mock_pif_atomic()
  mock_pif_90@conf_level <- 0.90

  mock_pif_99 <- create_mock_pif_atomic()
  mock_pif_99@conf_level <- 0.99

  ci_90 <- get_ci(mock_pif_90)
  ci_99 <- get_ci(mock_pif_99)

  # 99% CI should be wider than 90% CI
  width_90 <- ci_90[2] - ci_90[1]
  width_99 <- ci_99[2] - ci_99[1]
  expect_true(width_99 > width_90)
})

# Test getters with custom RR link functions
test_that("getters work with custom RR link functions", {
  # Test with identity RR link (RR = beta directly)
  custom_pif <- pif_atomic_class(
    p = c(0.3),
    p_cft = c(0.1),
    beta = c(1.5),  # This IS the relative risk
    rr_link = identity,
    rr_link_deriv = function(x) rep(1, length(x)),
    link = identity,
    link_deriv = function(x) rep(1, length(x)),
    link_inv = identity,
    var_p = matrix(0.01, nrow = 1),
    var_beta = matrix(0.02, nrow = 1),
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE,
    conf_level = 0.95,
    type = "PIF",
    label = "custom_rr"
  )

  expect_equal(get_rr(custom_pif), 1.5)
  expect_equal(get_rr_link_deriv_vals(custom_pif), 1)

  # Test with square link
  square_link_pif <- custom_pif
  square_link_pif@rr_link <- function(x) x^2
  square_link_pif@rr_link_deriv <- function(x) 2*x
  square_link_pif@beta <- c(sqrt(1.5))

  expect_equal(get_rr(square_link_pif), 1.5)
  expect_equal(get_rr_link_deriv_vals(square_link_pif), 2*sqrt(1.5))
})

# Test boundary conditions
test_that("getters handle boundary conditions correctly", {
  # Test with zero counterfactual prevalence (PAF case)
  paf_case <- pif_atomic_class(
    p = c(0.4),
    p_cft = c(0),  # Zero counterfactual
    beta = c(0.5),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = identity,
    link_deriv = function(x) rep(1, length(x)),
    link_inv = identity,
    var_p = matrix(0.01, nrow = 1),
    var_beta = matrix(0.02, nrow = 1),
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE,
    conf_level = 0.95,
    type = "PAF",
    label = "paf_boundary"
  )

  expect_equal(get_mu_cft(paf_case), 1.0)  # Should be 1 when p_cft = 0
  expect_true(get_pif(paf_case) > 0)

  # Test with prevalence summing to 1
  full_prev <- pif_atomic_class(
    p = c(0.6, 0.4),
    p_cft = c(0.3, 0.2),
    beta = c(0.2, 0.8),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = identity,
    link_deriv = function(x) rep(1, length(x)),
    link_inv = identity,
    var_p = diag(0.01, 2),
    var_beta = diag(0.02, 2),
    upper_bound_p = FALSE,
    upper_bound_beta = FALSE,
    conf_level = 0.95,
    type = "PIF",
    label = "full_prev"
  )

  expect_true(sum(full_prev@p) <= 1)
  expect_true(get_mu_obs(full_prev) > 0)
  expect_true(get_pif(full_prev) >= 0)
})
