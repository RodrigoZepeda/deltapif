# Mock S7 class definitions for testing
mock_pif_atomic_class <- S7::new_class(
  "deltapif::pif_atomic_class",
  properties = list(
    p = S7::class_numeric,
    p_cft = S7::class_numeric,
    beta = S7::class_numeric,
    rr = S7::class_numeric,
    rr_link = S7::class_function,
    rr_link_deriv = S7::class_function,
    mu_obs = S7::class_numeric,
    mu_cft = S7::class_numeric,
    pif = S7::class_numeric,
    link = S7::class_function,
    link_deriv = S7::class_function,
    link_inv = S7::class_function,
    variance = S7::class_numeric,
    var_p = S7::class_numeric,
    var_beta = S7::class_numeric,
    upper_bound_p = S7::class_numeric,
    upper_bound_beta = S7::class_numeric,
    conf_level = S7::class_numeric,
    link_vals = S7::class_numeric,
    link_deriv_vals = S7::class_numeric,
    rr_link_deriv_vals = S7::class_numeric,
    link_variance = S7::class_numeric
  )
)

mock_pif_class <- S7::new_class(
  "deltapif::pif_class",
  properties = list(
    pif_list = S7::class_list,
    weights = S7::class_numeric,
    sigma_weights = S7::class_numeric,
    coefs = S7::class_numeric,
    covariance = S7::class_numeric
  )
)

# Helper to create mock pif_atomic_class object
create_mock_pif_atomic <- function() {
  mock_pif_atomic_class(
    p = c(0.3, 0.7),
    p_cft = c(0.2, 0.5),
    beta = c(0.5, 1.0),
    rr = c(1.5, 2.0),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    mu_obs = mu_obs_fun(c(0.3, 0.7), c(1.5, 2.0)),
    mu_cft = mu_cft_fun(c(0.2, 0.5), c(1.5, 2.0)),
    pif = pif_fun2(
      mu_obs_fun(c(0.3, 0.7), c(1.5, 2.0)),
      mu_cft_fun(c(0.2, 0.5), c(1.5, 2.0))
    ),
    link = logit,
    link_deriv = deriv_logit,
    link_inv = inv_logit,
    variance = 0.01,
    var_p = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
    var_beta = matrix(c(0.02, 0, 0, 0.02), nrow = 2),
    upper_bound_p = 1,
    upper_bound_beta = Inf,
    conf_level = 0.95,
    link_vals = logit(pif_fun2(
      mu_obs_fun(c(0.3, 0.7), c(1.5, 2.0)),
      mu_cft_fun(c(0.2, 0.5), c(1.5, 2.0))
    )),
    link_deriv_vals = deriv_logit(pif_fun2(
      mu_obs_fun(c(0.3, 0.7), c(1.5, 2.0)),
      mu_cft_fun(c(0.2, 0.5), c(1.5, 2.0))
    )),
    rr_link_deriv_vals = exp(c(0.5, 1.0)),
    link_variance = (deriv_logit(pif_fun2(
      mu_obs_fun(c(0.3, 0.7), c(1.5, 2.0)),
      mu_cft_fun(c(0.2, 0.5), c(1.5, 2.0))
    )))^2 * 0.01
  )
}

# Helper to create mock pif_class object
create_mock_pif <- function() {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic()
  mock_pif_class(
    pif_list = list(pif1, pif2),
    weights = c(0.5, 0.5),
    sigma_weights = matrix(c(0.01, 0, 0, 0.01), nrow = 2),
    coefs = c(pif1@pif, pif2@pif),
    covariance = matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2)
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
               exp(sum(log(mock_pif@coefs))))
})

test_that("Covariance getters work with S7 objects", {
  mock_pif <- pif_total(
    paf(0.1, 1.2, var_p = 0.01, var_beta = 0, link = logit,
        link_deriv = deriv_logit, link_inv = inv_logit,
        rr_link = identity, rr_link_deriv = function(x) {1}), #Github actions crashes if not specifying compeltely this
    paf(0.1, 1.2, var_p = 0.01, var_beta = 0, link = logit,
        link_deriv = deriv_logit, link_inv = inv_logit,
        rr_link = identity, rr_link_deriv = function(x) {1}),
    weights = c(0.5, 0.5)
  )

  # Test only work with total pif not with mock
  expect_equal(get_total_coefs(mock_pif),
               sapply(mock_pif@pif_list, function(x) x@pif))


  # Test get_covariance_total
  cov_mat <- get_covariance_total(mock_pif)
  expect_equal(dim(cov_mat), c(2, 2))
  expect_true(isSymmetric(cov_mat))

  # Test get_ensemble_covariance
  ens_cov <- get_ensemble_covariance(mock_pif)
  expect_equal(dim(ens_cov), c(2, 2))
})

test_that("Edge cases are handled properly", {
  # Test with empty mock_pif object
  empty_pif <- mock_pif_atomic_class(
    p = numeric(0),
    p_cft = numeric(0),
    beta = numeric(0),
    rr = numeric(0),
    rr_link = function(x) numeric(0),
    rr_link_deriv = function(x) numeric(0),
    mu_obs = NA_real_,
    mu_cft = NA_real_,
    pif = NA_real_,
    link = identity,
    link_deriv = function(x) NA_real_,
    link_inv = identity,
    variance = NA_real_,
    var_p = matrix(numeric(0), nrow = 0),
    var_beta = matrix(numeric(0), nrow = 0),
    upper_bound_p = NA_real_,
    upper_bound_beta = NA_real_,
    conf_level = 0.95,
    link_vals = NA_real_,
    link_deriv_vals = NA_real_,
    rr_link_deriv_vals = numeric(0),
    link_variance = NA_real_
  )

  expect_equal(get_rr(empty_pif), numeric(0))
  expect_error(get_mu_obs(empty_pif))
  expect_error(get_mu_cft(empty_pif))
  expect_error(get_pif(empty_pif))
  expect_error(get_variance_atomic(empty_pif))
  expect_error(get_ci(empty_pif))
})

test_that("Input validation works with S7 objects", {
  expect_error(get_rr("not a mock_pif object"))
  expect_error(get_pif(list()))
})
