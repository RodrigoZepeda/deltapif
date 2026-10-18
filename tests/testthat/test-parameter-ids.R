test_that("Parameters given as named lists keep their identifier", {
  frac <- pif(p = list("447" = c(0.12, 0.30, 0.10)), p_cft = c(0.05, 0.2, 0.1),
              beta = list("beta_1" = log(c(1.2, 1.5, 2))), var_p = 0, var_beta = 0,
              quiet = TRUE)
  expect_equal(frac@p_id, "447")
  expect_equal(frac@beta_id, "beta_1")
  expect_equal(frac@p, c(0.12, 0.30, 0.10))

  #Equivalent to the numeric version
  frac_num <- pif(p = c(0.12, 0.30, 0.10), p_cft = c(0.05, 0.2, 0.1),
                  beta = log(c(1.2, 1.5, 2)), var_p = 0, var_beta = 0, quiet = TRUE)
  expect_equal(coef(frac), coef(frac_num))
  expect_true(is.na(frac_num@p_id))
  expect_true(is.na(frac_num@beta_id))

  #paf uses the length of the vector inside the list
  frac_paf <- paf(p = list("447" = c(0.12, 0.30, 0.10)), beta = log(c(1.2, 1.5, 2)),
                  var_p = 0, var_beta = 0, quiet = TRUE)
  expect_equal(frac_paf@p_cft, rep(0, 3))
})

test_that("Invalid named lists are rejected", {
  expect_error(pif(p = list(0.1), beta = 0.2, quiet = TRUE), "named list")
  expect_error(pif(p = list(a = 0.1, b = 0.2), beta = 0.2, quiet = TRUE), "named list")
  expect_error(pif(p = 0.1, beta = list(a = "x"), quiet = TRUE), "numeric")
})

test_that("Fractions are correlated only when they share an identifier", {
  shared_a <- paf(p = 0.085, beta = list("rr" = log(1.59)), var_beta = 0.165^2,
                  var_p = 0, label = "a", quiet = TRUE)
  shared_b <- paf(p = 0.12, beta = list("rr" = log(1.59)), var_beta = 0.165^2,
                  var_p = 0, label = "b", quiet = TRUE)
  indep_a  <- paf(p = 0.085, beta = log(1.59), var_beta = 0.165^2,
                  var_p = 0, label = "a", quiet = TRUE)
  indep_b  <- paf(p = 0.12, beta = log(1.59), var_beta = 0.165^2,
                  var_p = 0, label = "b", quiet = TRUE)

  w <- c(0.4, 0.6)
  shared <- paf_total(shared_a, shared_b, weights = w, quiet = TRUE)
  indep  <- paf_total(indep_a, indep_b, weights = w, quiet = TRUE)

  #Same relative risk: perfectly correlated fractions
  expect_equal(variance(shared),
               (w[1] * standard_deviation(shared_a) + w[2] * standard_deviation(shared_b))^2)

  #No identifier: independent fractions
  expect_equal(variance(indep), w[1]^2 * variance(indep_a) + w[2]^2 * variance(indep_b))
})

test_that("Using the same identifier for different estimates is an error", {
  a <- paf(p = 0.085, beta = list("rr" = log(1.59)), var_beta = 0.027,
           var_p = 0, label = "a", quiet = TRUE)
  b_value <- paf(p = 0.12, beta = list("rr" = log(1.60)), var_beta = 0.027,
                 var_p = 0, label = "b", quiet = TRUE)
  b_var   <- paf(p = 0.12, beta = list("rr" = log(1.59)), var_beta = 0.030,
                 var_p = 0, label = "b", quiet = TRUE)

  expect_error(paf_total(a, b_value, weights = c(0.5, 0.5)), "same identifier")
  expect_error(paf_total(a, b_var, weights = c(0.5, 0.5)), "same identifier")
  expect_error(covariance(a, b_value), "same identifier")
})

test_that("Covariances between different identifiers enter totals and ensembles", {
  a <- paf(p = 0.3, beta = list("lead" = log(2.2)), var_beta = 0.04, var_p = 0,
           label = "Lead", quiet = TRUE)
  b <- paf(p = 0.1, beta = list("radiation" = log(1.2)), var_beta = 0.02, var_p = 0,
           label = "Radiation", quiet = TRUE)
  ids    <- c("lead", "radiation")
  cov_rr <- matrix(c(0.04, 0.01, 0.01, 0.02), ncol = 2, dimnames = list(ids, ids))

  #Gradient of each PAF with respect to its beta
  grad_a <- deriv_pif_beta(p = a@p, p_cft = a@p_cft, rr = a@rr, rr_link_deriv_vals = a@rr)
  grad_b <- deriv_pif_beta(p = b@p, p_cft = b@p_cft, rr = b@rr, rr_link_deriv_vals = b@rr)

  w <- c(0.7, 0.3)
  total <- paf_total(a, b, weights = w, var_beta = cov_rr, quiet = TRUE)
  expect_equal(
    variance(total),
    w[1]^2 * variance(a) + w[2]^2 * variance(b) + 2 * w[1] * w[2] * grad_a * grad_b * 0.01
  )

  ens       <- paf_ensemble(a, b, var_beta = cov_rr, quiet = TRUE)
  ens_indep <- paf_ensemble(a, b, quiet = TRUE)
  cov_ab    <- grad_a * grad_b * 0.01
  expect_equal(
    variance(ens) - variance(ens_indep),
    as.numeric(2 * (1 - coef(ens))^2 * cov_ab / ((1 - coef(a)) * (1 - coef(b))))
  )
})

test_that("Vector-valued identifiers repeat their name in the covariance matrix", {
  a <- pif(p = list("survey_a" = c(0.2, 0.1)), p_cft = c(0.1, 0.05), beta = log(c(1.5, 2)),
           var_p = diag(c(1e-4, 5e-5)), var_beta = 0, label = "a", quiet = TRUE)
  b <- pif(p = list("survey_b" = 0.3), p_cft = 0.1, beta = log(1.4),
           var_p = 2e-4, var_beta = 0, label = "b", quiet = TRUE)

  ids <- c("survey_a", "survey_a", "survey_b")
  cov_p <- matrix(c(1e-4, 0, 3e-5,
                    0, 5e-5, 1e-5,
                    3e-5, 1e-5, 2e-4), ncol = 3, dimnames = list(ids, ids))

  grad_a <- deriv_pif_p(p = a@p, p_cft = a@p_cft, rr = a@rr)
  grad_b <- deriv_pif_p(p = b@p, p_cft = b@p_cft, rr = b@rr)

  w <- c(0.5, 0.5)
  total <- pif_total(a, b, weights = w, var_p = cov_p, quiet = TRUE)
  expect_equal(
    variance(total),
    as.numeric(w[1]^2 * variance(a) + w[2]^2 * variance(b) +
                 2 * w[1] * w[2] * t(grad_a) %*% cov_p[1:2, 3, drop = FALSE] %*% grad_b)
  )

  #Wrong number of repetitions
  bad_ids <- c("survey_a", "survey_b", "survey_b")
  bad_cov <- cov_p
  dimnames(bad_cov) <- list(bad_ids, bad_ids)
  expect_error(pif_total(a, b, weights = w, var_p = bad_cov, quiet = TRUE))
})

test_that("Covariance matrices between identifiers are validated", {
  a <- paf(p = 0.3, beta = list("lead" = log(2.2)), var_beta = 0.04, var_p = 0,
           label = "Lead", quiet = TRUE)
  b <- paf(p = 0.1, beta = list("radiation" = log(1.2)), var_beta = 0.02, var_p = 0,
           label = "Radiation", quiet = TRUE)
  ids <- c("lead", "radiation")

  #No names
  expect_error(paf_ensemble(a, b, var_beta = matrix(c(0.04, 0.01, 0.01, 0.02), 2)),
               "row and column names")

  #Unknown identifier
  other <- c("lead", "asbestos")
  expect_error(
    paf_ensemble(a, b, var_beta = matrix(c(0.04, 0.01, 0.01, 0.02), 2,
                                         dimnames = list(other, other))),
    "asbestos"
  )

  #Diagonal differs from the variance of the fraction
  expect_error(
    paf_ensemble(a, b, var_beta = matrix(c(0.05, 0.01, 0.01, 0.02), 2,
                                         dimnames = list(ids, ids))),
    "differs"
  )

  #Not symmetric
  expect_error(
    paf_ensemble(a, b, var_beta = matrix(c(0.04, 0.01, 0.02, 0.02), 2,
                                         dimnames = list(ids, ids))),
    "symmetric"
  )
})

test_that("Weighted adjusted fractions use the identifiers", {
  a <- paf(p = 0.3, beta = list("rr" = log(1.8)), var_beta = 0.02, var_p = 0.001,
           label = "a", quiet = TRUE)
  b <- paf(p = 0.2, beta = list("rr" = log(1.8)), var_beta = 0.02, var_p = 0.001,
           label = "b", quiet = TRUE)
  a0 <- paf(p = 0.3, beta = log(1.8), var_beta = 0.02, var_p = 0.001,
            label = "a", quiet = TRUE)
  b0 <- paf(p = 0.2, beta = log(1.8), var_beta = 0.02, var_p = 0.001,
            label = "b", quiet = TRUE)

  adj_shared <- weighted_adjusted_paf(a, b, quiet = TRUE)
  adj_indep  <- weighted_adjusted_paf(a0, b0, quiet = TRUE)
  expect_equal(coef(adj_shared[["a"]]), coef(adj_indep[["a"]]))
  expect_false(isTRUE(all.equal(variance(adj_shared[["a"]]), variance(adj_indep[["a"]]))))
})
