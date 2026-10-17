test_that("variance_components splits the variance of an atomic fraction", {
  sd_beta <- (log(2.20) - log(1.15)) / (2 * qnorm(0.975))
  pif_obj <- pif(p = 0.085, p_cft = 0.07225, beta = log(1.59),
                 var_beta = sd_beta^2, var_p = 0.01^2, quiet = TRUE)

  mu_obs <- 1 + 0.085 * 0.59
  mu_cft <- 1 + 0.07225 * 0.59
  expected_p    <- (mu_cft * 0.59 / mu_obs^2)^2 * 0.01^2
  expected_beta <- (1.59 * (mu_obs * 0.07225 - mu_cft * 0.085) / mu_obs^2)^2 * sd_beta^2

  vc <- variance_components(pif_obj)
  expect_named(vc, c("p", "beta", "weights"))
  expect_equal(unname(vc["p"]), expected_p)
  expect_equal(unname(vc["beta"]), expected_beta)
  expect_equal(unname(vc["weights"]), 0)
  expect_equal(sum(vc), variance(pif_obj))
})

test_that("variance_components of a total separates parameters and weights", {
  a <- paf(p = 0.085, beta = log(1.59), var_beta = 0.165^2, var_p = 0.01^2,
           label = "a", quiet = TRUE)
  b <- paf(p = 0.12, beta = log(1.59), var_beta = 0.165^2, var_p = 0.01^2,
           label = "b", quiet = TRUE)
  total <- paf_total(a, b, weights = c(0.4, 0.6), var_weights = diag(c(0.001, 0.002)),
                     quiet = TRUE)

  va <- variance_components(a)
  vb <- variance_components(b)
  vc <- variance_components(total)

  #Independent prevalences, shared relative risk and independent weights
  expect_equal(unname(vc["p"]), unname(0.4^2 * va["p"] + 0.6^2 * vb["p"]))
  expect_equal(unname(vc["beta"]), unname((0.4 * sqrt(va["beta"]) + 0.6 * sqrt(vb["beta"]))^2))
  expect_equal(unname(vc["weights"]), coef(a)^2 * 0.001 + coef(b)^2 * 0.002)
  expect_equal(sum(vc), variance(total))
})

test_that("variance_components adds up for ensembles", {
  a <- paf(p = 0.3, beta = log(1.4), var_beta = 0.02, var_p = 0.001, label = "a", quiet = TRUE)
  b <- paf(p = 0.2, beta = log(1.9), var_beta = 0.03, var_p = 0.002, label = "b", quiet = TRUE)
  ens <- paf_ensemble(a, b, weights = c(0.8, 0.9), quiet = TRUE)
  vc  <- variance_components(ens)
  expect_equal(sum(vc), variance(ens))
  expect_true(all(vc >= 0))
})
