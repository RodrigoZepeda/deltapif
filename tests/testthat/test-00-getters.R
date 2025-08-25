# Helper to create mock pif_atomic_class object
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
    type = type,
    label = paste0("atomic_pif", stats::rnorm(1))
  )
}

# Helper to create mock pif_class object
create_mock_pif <- function() {
  pif1 <- create_mock_pif_atomic()
  pif2 <- create_mock_pif_atomic()
  pif_total_class(
    pif_list = list(pif1, pif2),
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
               1 - exp(sum(log(1 - mock_pif@coefs))))
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

  tp <- pif_total_class(list(pif1, pif2), link = identity, label = "total_class",
                        var_pif_weights = matrix(0, 2, 2),
                        link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                        weights = c(0.5, 0.5), var_weights = matrix(0, 2, 2))

  ep <- pif_ensemble_class(list(pif1, pif2), link = identity, label = "ensemble_class",
                           var_pif_weights = matrix(0, 2, 2),
                           link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                           weights = rep(1, 2), var_weights = matrix(0, 2, 2))

  expect_equal(get_ensemble_type(tp), "PIF")
  expect_equal(get_ensemble_type(ep), "PIF")

  tp <- pif_total_class(list(paf1, pif2), link = identity,
                        label = "total_class",
                        var_pif_weights = matrix(0, 2, 2),
                        link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                        weights = c(0.5, 0.5), var_weights = matrix(0, 2, 2))

  ep <- pif_ensemble_class(list(pif1, paf2), link = identity,
                           label = "ensemble_class",
                           var_pif_weights = matrix(0, 2, 2),
                           link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                           weights = rep(1, 2), var_weights = matrix(0, 2, 2))


  expect_equal(get_ensemble_type(tp), "PIF")
  expect_equal(get_ensemble_type(ep), "PIF")

  tp <- pif_total_class(list(paf1, paf2), link = identity,
                        label = "total_class",
                        var_pif_weights = matrix(0, 2, 2),
                        link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                        weights = c(0.5, 0.5), var_weights = matrix(0, 2, 2))

  ep <- pif_ensemble_class(list(paf1, paf2), link = identity,
                           label = "ensemble_class",
                           var_pif_weights = matrix(0, 2, 2),
                           link_inv = identity, link_deriv = function(x) rep(1, length(x)),
                           weights = rep(1, 2), var_weights = matrix(0, 2, 2))


  expect_equal(get_ensemble_type(tp), "PAF")
  expect_equal(get_ensemble_type(ep), "PAF")
})

