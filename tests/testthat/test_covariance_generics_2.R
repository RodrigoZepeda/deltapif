# Helper -----------------------------------------------------------------------

make_pif <- function(p = 0.3, p_cft = 0.1, beta = 1.5,
                     var_p = 0.01, var_beta = 0.05,
                     label = paste0("pif_", round(abs(stats::rnorm(1)), 4))) {
  pif(p = p, p_cft = p_cft, beta = beta,
      var_p = var_p, var_beta = var_beta,
      quiet = TRUE, label = label)
}

make_paf <- function(p = 0.3, beta = 1.5,
                     var_p = 0.01, var_beta = 0.05,
                     label = paste0("paf_", round(abs(stats::rnorm(1)), 4))) {
  paf(p = p, beta = beta,
      var_p = var_p, var_beta = var_beta,
      quiet = TRUE, label = label)
}

# cov_atomic_pif ---------------------------------------------------------------

test_that("cov_atomic_pif returns a numeric scalar for independent pifs", {
  p1 <- make_paf(label = "a1")
  p2 <- make_paf(p = 0.4, beta = 1.2, label = "a2")

  result <- cov_atomic_pif(p1, p2,
                           var_p    = matrix(0, 1, 1),
                           var_beta = matrix(0, 1, 1))
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_equal(result, 0)
})

test_that("cov_atomic_pif with same pif equals variance", {
  p1 <- make_paf(label = "same")
  expect_equal(
    cov_atomic_pif(p1, p1, var_p = p1@var_p, var_beta = p1@var_beta),
    variance(p1)
  )
})

test_that("cov_atomic_pif auto-detects shared beta and p", {
  p1 <- make_paf(label = "shared1")
  p2 <- make_paf(label = "shared2")   # same p/beta values as p1
  result <- cov_atomic_pif(p1, p2)
  expect_type(result, "double")
  expect_true(result > 0)             # correlated because parameters match
})

test_that("cov_atomic_pif errors on non-pif_atomic inputs", {
  p1 <- make_paf(label = "err1")
  expect_error(cov_atomic_pif(list(), p1), "pif_atomic_class")
  expect_error(cov_atomic_pif(p1, list()), "pif_atomic_class")
})

test_that("cov_atomic_pif errors on wrong var_p dimensions", {
  p1 <- make_paf(label = "dim1")
  p2 <- make_paf(label = "dim2")
  expect_error(
    cov_atomic_pif(p1, p2, var_p = matrix(c(1, 2, 3, 4), nrow = 2)),
    "dimensions"
  )
})

test_that("cov_atomic_pif errors on invalid var_p type", {
  p1 <- make_paf(label = "vp1")
  p2 <- make_paf(label = "vp2")
  expect_error(
    cov_atomic_pif(p1, p2, var_p = "bad"),
    "covariance_structure_class"
  )
})

test_that("cov_atomic_pif accepts covariance_structure_class inputs", {
  p1 <- make_paf(label = "cs1")
  p2 <- make_paf(label = "cs2")

  cs_p    <- default_parameter_covariance_structure2(p1, p2, parameter = "p")
  cs_beta <- default_parameter_covariance_structure2(p1, p2, parameter = "beta")

  expect_silent(cov_atomic_pif(p1, p2, var_p = cs_p, var_beta = cs_beta))
})

test_that("cov_atomic_pif works with multi-dimensional p and beta", {
  p1 <- paf(p = c(0.3, 0.2), beta = c(1.2, 0.8),
            var_p = diag(c(0.01, 0.01)),
            var_beta = diag(c(0.05, 0.05)),
            quiet = TRUE, label = "md1")
  p2 <- paf(p = c(0.4, 0.1), beta = c(0.9, 1.1),
            var_p = diag(c(0.02, 0.02)),
            var_beta = diag(c(0.04, 0.04)),
            quiet = TRUE, label = "md2")

  result <- cov_atomic_pif(p1, p2,
                           var_p    = matrix(0, nrow = 2, ncol = 2),
                           var_beta = matrix(0, nrow = 2, ncol = 2))
  expect_type(result, "double")
  expect_length(result, 1L)
})

# variance generic -------------------------------------------------------------

test_that("variance returns a non-negative double for atomic pif", {
  p1 <- make_paf(label = "var1")
  v  <- variance(p1)
  expect_type(v, "double")
  expect_true(v >= 0)
})

test_that("variance of atomic equals cov_atomic_pif with itself", {
  p1 <- make_paf(label = "var_eq")
  expect_equal(
    variance(p1),
    cov_atomic_pif(p1, p1, var_p = p1@var_p, var_beta = p1@var_beta)
  )
})

test_that("variance of zero-variance pif is zero", {
  p0 <- make_paf(var_p = 0, var_beta = 0, label = "var0")
  expect_equal(variance(p0), 0)
})

test_that("variance of ensemble pif is non-negative", {
  p1  <- make_paf(label = "ve1")
  p2  <- make_paf(p = 0.4, beta = 1.2, label = "ve2")
  ens <- pif_ensemble(p1, p2, weights = c(0.6, 0.4),
                      var_weights = matrix(c(0.1, 0.02, 0.02, 0.1), ncol = 2),
                      label = "ens_v")
  v <- variance(ens)
  expect_type(v, "double")
  expect_true(v >= 0)
})

test_that("variance of total pif is non-negative", {
  p1  <- make_paf(label = "vt1")
  p2  <- make_paf(p = 0.4, beta = 1.2, label = "vt2")
  tot <- pif_total(p1, p2, weights = c(0.5, 0.5),
                   var_weights = matrix(c(0.1, 0, 0, 0.1), ncol = 2),
                   label = "tot_v")
  expect_true(variance(tot) >= 0)
})

# standard_deviation generic ---------------------------------------------------

test_that("standard_deviation equals sqrt(variance)", {
  p1 <- make_paf(label = "sd1")
  expect_equal(standard_deviation(p1), sqrt(variance(p1)))
})

# covariance generic -----------------------------------------------------------

test_that("covariance of a single pif is a 1x1 matrix equal to variance", {
  p1 <- make_paf(label = "cov_single")
  cm <- covariance(p1)
  expect_true(is.matrix(cm))
  expect_equal(dim(cm), c(1L, 1L))
  expect_equal(cm[1, 1], variance(p1))
})

test_that("covariance matrix is symmetric with correct diagonal", {
  p1 <- make_paf(label = "cm1")
  p2 <- make_paf(p = 0.4, beta = 1.2, label = "cm2")
  p3 <- make_paf(p = 0.2, beta = 0.9, label = "cm3")

  cm <- covariance(p1, p2, p3)
  expect_equal(dim(cm), c(3L, 3L))
  expect_true(isSymmetric(cm))
  expect_equal(diag(cm), c(variance(p1), variance(p2), variance(p3)))
})

test_that("covariance of independent pifs has near-zero off-diagonal", {
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05, quiet = TRUE, label = "ind1")
  p2 <- paf(0.4, 0.8, var_p = 0.02, var_beta = 0.03, quiet = TRUE, label = "ind2")

  cm <- covariance(p1, p2,
                   var_p    = matrix(c(1, 0, 0, 1), nrow = 2,
                                     dimnames = list(c("ind1", "ind2"), c("ind1", "ind2"))),
                   var_beta = matrix(c(1, 0, 0, 1), nrow = 2,
                                     dimnames = list(c("ind1", "ind2"), c("ind1", "ind2"))))
  expect_equal(cm[1, 2], 0)
})

test_that("covariance errors on bad var_p row names", {
  p1 <- make_paf(label = "cn1")
  p2 <- make_paf(label = "cn2")
  bad_vp <- matrix(c(1, 0, 0, 1), nrow = 2,
                   dimnames = list(c("wrong", "names"), c("cn1", "cn2")))
  expect_error(covariance(p1, p2, var_p = bad_vp), "was not found")
})

# correlation generic ----------------------------------------------------------

test_that("correlation of a single pif is 1x1 matrix of 1", {
  p1  <- make_paf(label = "cor_s")
  cor <- correlation(p1)
  expect_equal(dim(cor), c(1L, 1L))
  expect_equal(cor[1, 1], 1)
})

test_that("correlation matrix has 1 on diagonal and values in [-1, 1]", {
  p1  <- make_paf(label = "cor1")
  p2  <- make_paf(p = 0.4, beta = 1.2, label = "cor2")
  cor <- correlation(p1, p2)
  expect_true(isSymmetric(cor))
  expect_true(all(diag(cor) == 1))
  expect_true(all(cor >= -1 & cor <= 1))
})

test_that("correlation warns on zero-variance pif", {
  p0 <- make_paf(var_p = 0, var_beta = 0, label = "cor0")
  p1 <- make_paf(label = "cor1b")
  expect_warning(correlation(p0, p1))
})

# cov_ensemble_weights ---------------------------------------------------------

test_that("cov_ensemble_weights returns 0 for two atomic pifs", {
  p1 <- make_paf(label = "ew1")
  p2 <- make_paf(label = "ew2")
  expect_equal(cov_ensemble_weights(p1, p2), 0)
})

test_that("cov_ensemble_weights returns 0 for ensemble vs atomic", {
  p1  <- make_paf(label = "ewa1")
  p2  <- make_paf(label = "ewa2")
  p3  <- make_paf(label = "ewa3")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "ens_ew")
  expect_equal(cov_ensemble_weights(ens, p3), 0)
})

test_that("cov_ensemble_weights runs silently for two ensembles", {
  p1   <- make_paf(label = "ee1")
  p2   <- make_paf(label = "ee2")
  p3   <- make_paf(label = "ee3")
  p4   <- make_paf(label = "ee4")
  ens1 <- pif_ensemble(p1, p2, weights = c(0.6, 0.4),
                       var_weights = matrix(c(0.1, 0.02, 0.02, 0.1), ncol = 2),
                       label = "ens_ee1")
  ens2 <- pif_total(p3, p4, weights = c(0.3, 0.7),
                    var_weights = matrix(c(0.2, 0, 0, 0.2), ncol = 2),
                    label = "tot_ee2")
  expect_silent(cov_ensemble_weights(ens1, ens2))
})

test_that("cov_ensemble_weights errors on non-pif inputs", {
  expect_error(cov_ensemble_weights(list(), make_paf(label = "err_ew")),
               "must be.*pif_class")
})

# cov_ensemble_atomic ----------------------------------------------------------

test_that("cov_ensemble_atomic with two atomics equals cov_atomic_pif", {
  p1 <- make_paf(label = "ea_at1")
  p2 <- make_paf(label = "ea_at2")
  expect_equal(
    cov_ensemble_atomic(p1, p2),
    cov_atomic_pif(p1, p2)
  )
})

test_that("cov_ensemble_atomic runs silently for ensemble vs atomic", {
  p1  <- make_paf(label = "ea1")
  p2  <- make_paf(label = "ea2")
  p3  <- make_paf(p = 0.2, beta = 0.9, label = "ea3")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "ens_ea")
  expect_silent(cov_ensemble_atomic(ens, p3))
})

test_that("cov_ensemble_atomic errors on bad first argument", {
  p1 <- make_paf(label = "ea_err1")
  expect_error(cov_ensemble_atomic(list(), p1),
               "pif_global_ensemble_class|pif_atomic")
})

test_that("cov_ensemble_atomic errors on bad second argument", {
  p1  <- make_paf(label = "ea_err2a")
  p2  <- make_paf(label = "ea_err2b")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "ens_ea2")
  expect_error(cov_ensemble_atomic(ens, list()),
               "pif_atomic_class")
})

# cov_total_pif ----------------------------------------------------------------

test_that("cov_total_pif runs for two atomics", {
  p1 <- make_paf(label = "tp1")
  p2 <- make_paf(p = 0.4, beta = 1.2, label = "tp2")
  expect_silent(cov_total_pif(p1, p2))
})

test_that("cov_total_pif runs for atomic vs total", {
  p1  <- make_paf(label = "tat1")
  p2  <- make_paf(label = "tat2")
  p3  <- make_paf(p = 0.2, beta = 0.8, label = "tat3")
  tot <- pif_total(p1, p2, weights = c(0.5, 0.5),
                   var_weights = matrix(c(0.1, 0, 0, 0.1), ncol = 2),
                   label = "tot_tat")
  expect_silent(cov_total_pif(p3, tot))
  expect_silent(cov_total_pif(tot, p3))
})

test_that("cov_total_pif runs for two totals", {
  p1  <- make_paf(label = "tt1")
  p2  <- make_paf(label = "tt2")
  tot <- pif_total(p1, p2, weights = c(0.5, 0.5),
                   var_weights = matrix(c(0.1, 0, 0, 0.1), ncol = 2),
                   label = "tot_tt")
  expect_silent(cov_total_pif(tot, tot))
})

test_that("cov_total_pif errors on non-pif_class input", {
  p1 <- make_paf(label = "te1")
  expect_error(cov_total_pif(p1, list()), "pif_class")
})

# Integration: covariance across ensemble hierarchy ----------------------------

test_that("covariance is consistent across hierarchy levels", {
  base1 <- make_paf(label = "h1")
  base2 <- make_paf(p = 0.4, beta = 1.2, label = "h2")
  base3 <- make_paf(p = 0.2, beta = 0.8, label = "h3")
  ens   <- pif_ensemble(base1, base2, weights = c(0.6, 0.4),
                        var_weights = matrix(c(0.1, 0.02, 0.02, 0.1), ncol = 2),
                        label = "h_ens")
  tot   <- pif_total(ens, base3, weights = c(0.7, 0.3),
                     var_weights = matrix(c(0.1, 0, 0, 0.1), ncol = 2),
                     label = "h_tot")

  # Variance of total should be non-negative
  expect_true(variance(tot) >= 0)

  # covariance of total with a base pif
  cm <- covariance(tot, base1, warning = FALSE)
  expect_equal(dim(cm), c(2L, 2L))
  expect_true(isSymmetric(cm))
})
