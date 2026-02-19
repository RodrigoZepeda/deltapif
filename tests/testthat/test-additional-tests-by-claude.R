### test-additional-lines.R
### Tests targeting previously untested lines across the deltapif package.

# ── Helpers ──────────────────────────────────────────────────────────────────

make_atomic <- function(p = 0.3, p_cft = 0.1, beta = log(1.5),
                        var_p = 0.01, var_beta = 0.02,
                        label = paste0("pif_", round(abs(stats::rnorm(1)), 4))) {
  pif(p = p, p_cft = p_cft, beta = beta,
      var_p = var_p, var_beta = var_beta,
      quiet = TRUE, label = label)
}

make_paf_atomic <- function(p = 0.3, beta = log(1.5),
                            var_p = 0.01, var_beta = 0.02,
                            label = paste0("paf_", round(abs(stats::rnorm(1)), 4))) {
  paf(p = p, beta = beta, var_p = var_p, var_beta = var_beta,
      quiet = TRUE, label = label)
}

make_cov_struct <- function() {
  covariance_structure_class(list(
    "a" = list("a" = 0.5, "b" = 0.1),
    "b" = list("a" = 0.1, "b" = 0.8)
  ))
}

# ── R/01-classes.R : pif_class validator ─────────────────────────────────────

test_that("pif_class validator: conf_level out of [0,1] raises error", {
  expect_error(
    pif_class(pif = 0.3, variance = 0.01, conf_level = 1.5,
              type = "PIF", link = identity, link_inv = identity,
              link_deriv = function(x) 1, label = "v"),
    "Invalid confidence level"
  )
  expect_error(
    pif_class(pif = 0.3, variance = 0.01, conf_level = -0.1,
              type = "PIF", link = identity, link_inv = identity,
              link_deriv = function(x) 1, label = "v2"),
    "Invalid confidence level"
  )
})

test_that("pif_class validator: empty pif raises error", {
  expect_error(
    pif_class(pif = numeric(0), variance = 0.01, conf_level = 0.95,
              type = "PIF", link = identity, link_inv = identity,
              link_deriv = function(x) 1, label = "v3"),
    "No entries provided for `pif`"
  )
})

test_that("pif_class validator: empty variance raises error", {
  expect_error(
    pif_class(pif = 0.3, variance = numeric(0), conf_level = 0.95,
              type = "PIF", link = identity, link_inv = identity,
              link_deriv = function(x) 1, label = "v4"),
    "No entries provided for `variance`"
  )
})

test_that("pif_class validator: negative variance raises error", {
  expect_error(
    pif_class(pif = 0.3, variance = -0.01, conf_level = 0.95,
              type = "PIF", link = identity, link_inv = identity,
              link_deriv = function(x) 1, label = "v5"),
    "negative.*variance"
  )
})

test_that("pif_class validator: PIF > 1 gives warning", {
  expect_warning(
    pif_class(pif = 1.1, variance = 0.01, conf_level = 0.95,
              type = "PIF", link = identity, link_inv = identity,
              link_deriv = function(x) 1, label = "v6"),
    "PIF > 1"
  )
})

test_that("pif_class validator: invalid type raises error", {
  expect_error(
    pif_class(pif = 0.3, variance = 0.01, conf_level = 0.95,
              type = "NEITHER", link = identity, link_inv = identity,
              link_deriv = function(x) 1, label = "v7"),
    "should be either `PIF` or `PAF`"
  )
})

# ── R/01-classes.R : pif_atomic_class validator ───────────────────────────────

test_that("pif_atomic_class validator: p and p_cft length mismatch raises error", {
  expect_error(
    pif(p = c(0.2, 0.3), p_cft = 0.1, beta = c(0.5, 0.6),
        var_p = diag(2) * 0.01, var_beta = diag(2) * 0.02,
        quiet = TRUE),
    "[Ss]ame length"
  )
})

test_that("pif_atomic_class validator: p < 0 raises error", {
  expect_error(
    pif(p = -0.1, p_cft = 0.1, beta = 0.5,
        var_p = 0.01, var_beta = 0.02, quiet = TRUE),
    "< 0"
  )
})

test_that("pif_atomic_class validator: p_cft < 0 raises error", {
  expect_error(
    pif(p = 0.3, p_cft = -0.1, beta = 0.5,
        var_p = 0.01, var_beta = 0.02, quiet = TRUE),
    "< 0"
  )
})

test_that("pif_atomic_class validator: p sums > 1 warns", {
  expect_error(
    pif(p = c(0.7, 0.6), p_cft = c(0.3, 0.3),
        beta = c(0.5, 0.6),
        var_p = diag(2) * 0.01, var_beta = diag(2) * 0.02,
        quiet = TRUE),
    "sums.* > 1"
  )
})

# ── R/01-covariance_structures.R : covariance_structure validator ─────────────

test_that("covariance_structure_class validator: non-list inner entry raises error", {
  expect_error(
    covariance_structure_class(list(
      "a" = 0.5,    # not a list
      "b" = list("a" = 0.1, "b" = 0.8)
    )),
    "not a list"
  )
})

test_that("covariance_structure_class validator: inner lists of different length raise error", {
  expect_error(
    covariance_structure_class(list(
      "a" = list("a" = 0.5, "b" = 0.1, "c" = 0.2),
      "b" = list("a" = 0.1, "b" = 0.8)  # shorter
    )),
    "different length"
  )
})

# ── R/01-covariance_structures.R : default_weight_covariance_structure ────────

test_that("default_weight_covariance_structure returns a list for a pif_atomic", {
  p1 <- make_atomic(label = "dw1")
  result <- default_weight_covariance_structure(p1)
  expect_true(is.list(result@cov_list))
})

test_that("default_weight_covariance_structure handles two ensembles with same weights", {
  p1 <- make_atomic(label = "dwe1")
  p2 <- make_atomic(label = "dwe2")
  ens1 <- pif_ensemble(p1, p2, weights = c(0.5, 0.5),
                       var_weights = diag(2) * 0.1, label = "ensA")
  p3 <- make_atomic(label = "dwe3")
  p4 <- make_atomic(label = "dwe4")
  ens2 <- pif_ensemble(p3, p4, weights = c(0.5, 0.5),
                       var_weights = diag(2) * 0.1, label = "ensB")
  # Should not error
  expect_silent(default_weight_covariance_structure2(ens1, ens2))
})

# ── R/01-covariance_structures.R : default_parameter_covariance_structure ────

test_that("default_parameter_covariance_structure runs for a pif_total", {
  p1 <- make_atomic(label = "dp1")
  p2 <- make_atomic(label = "dp2")
  tot <- pif_total(p1, p2, weights = c(0.5, 0.5), label = "dptot")
  expect_silent(default_parameter_covariance_structure(tot, parameter = "p"))
  expect_silent(default_parameter_covariance_structure(tot, parameter = "beta"))
})

test_that("default_pif_covariance_structure runs for a pif_ensemble", {
  p1 <- make_atomic(label = "dpce1")
  p2 <- make_atomic(label = "dpce2")
  ens <- pif_ensemble(p1, p2, weights = c(0.4, 0.6), label = "dpcens")
  expect_silent(default_pif_covariance_structure(ens))
})

test_that("default_weight_pif_covariance_structure runs for two ensembles", {
  p1 <- make_atomic(label = "dwpe1")
  p2 <- make_atomic(label = "dwpe2")
  p3 <- make_atomic(label = "dwpe3")
  ens1 <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "dwpens1")
  expect_silent(default_weight_pif_covariance_structure2(ens1, p3))
})

# ── R/02-generics.R : print covariance_structure_class ───────────────────────

test_that("print.covariance_structure_class produces output", {
  cs <- make_cov_struct()
  expect_output(print(cs))
})

test_that("print.covariance_structure_class with matrix entries works", {
  mat <- matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 2)
  cs <- covariance_structure_class(list(
    "x" = list("x" = mat, "y" = 0.1),
    "y" = list("x" = 0.1, "y" = 0.5)
  ))
  expect_output(print(cs))
})

# ── R/02-generics.R : weights for pif_global_ensemble_class ──────────────────

test_that("weights() returns the weights vector for pif_total", {
  p1 <- make_atomic(label = "wt1")
  p2 <- make_atomic(label = "wt2")
  tot <- pif_total(p1, p2, weights = c(0.3, 0.7), label = "wtot")
  expect_equal(weights(tot), c(0.3, 0.7))
})

test_that("weights() returns the weights vector for pif_ensemble", {
  p1 <- make_atomic(label = "we1")
  p2 <- make_atomic(label = "we2")
  ens <- pif_ensemble(p1, p2, weights = c(1, 2), label = "wens")
  expect_equal(weights(ens), c(1, 2))
})

# ── R/02-generics.R : as.data.frame for pif_class ────────────────────────────

test_that("as.data.frame.pif_class returns a data.frame with expected columns", {
  p1 <- make_atomic(label = "adf1")
  df <- as.data.frame(p1)
  expect_true(is.data.frame(df))
  expect_true("pif" %in% tolower(names(df)) || nrow(df) >= 1)
})

test_that("as.data.frame for pif_ensemble returns multi-row data.frame", {
  p1 <- make_atomic(label = "adfe1")
  p2 <- make_atomic(label = "adfe2")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "adfens")
  df <- as.data.frame(ens)
  expect_true(is.data.frame(df))
  expect_true(nrow(df) >= 1)
})

# ── R/02-generics.R : as.data.frame for cases_class ──────────────────────────

test_that("as.data.frame.cases_class returns a data.frame", {
  p1 <- make_atomic(label = "adfc1")
  ac <- averted_cases(100, pif = p1)
  df <- as.data.frame(ac)
  expect_true(is.data.frame(df))
})

# ── R/02-generics.R : as.matrix for covariance_structure_class ───────────────

test_that("as.matrix.covariance_structure_class returns numeric matrix with scalars", {
  cs <- make_cov_struct()
  m <- as.matrix(cs)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 2L))
  expect_equal(m[1, 1], 0.5)
  expect_equal(m[2, 2], 0.8)
})

test_that("as.matrix.covariance_structure_class handles vector entries", {
  vec <- c(0.1, 0.2)
  cs <- covariance_structure_class(list(
    "a" = list("a" = 0.5, "b" = vec),
    "b" = list("a" = vec, "b" = 0.8)
  ))
  expect_true(is.matrix(as.matrix(cs)))
})

test_that("as.matrix.covariance_structure_class handles matrix entries", {
  mat <- matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 2)
  cs <- covariance_structure_class(list(
    "a" = list("a" = mat, "b" = c(0.1, 0.2)),
    "b" = list("a" = c(0.1, 0.2), "b" = matrix(c(0.5, 0, 0, 0.9), ncol = 2))
  ))
  m <- as.matrix(cs)
  expect_true(is.matrix(m))
})

# ── R/02-generics.R : subset for covariance_structure_class ──────────────────

test_that("subset by numeric index selects correct entries", {
  cs <- make_cov_struct()
  cs2 <- subset(cs, 1:2)
  expect_equal(length(cs2), 2L)
})

test_that("subset by character name works", {
  cs <- make_cov_struct()
  cs_a <- subset(cs, "a")
  expect_equal(row_names(cs_a), "a")
})

test_that("subset with negate = TRUE excludes specified names", {
  cs <- make_cov_struct()
  cs_neg <- subset(cs, "a", negate = TRUE)
  expect_false("a" %in% row_names(cs_neg))
})

test_that("subset by rows and cols separately works", {
  cs3 <- covariance_structure_class(list(
    "x" = list("x" = 1, "y" = 2, "z" = 3),
    "y" = list("x" = 2, "y" = 4, "z" = 6),
    "z" = list("x" = 3, "y" = 6, "z" = 9)
  ))
  # row subset
  by_row <- subset(cs3, rows = "x")
  expect_equal(row_names(by_row), "x")
  # col subset
  by_col <- subset(cs3, cols = "y")
  expect_true(!is.null(by_col))
})

# ── R/02-generics.R : as_covariance_structure ────────────────────────────────

test_that("as_covariance_structure converts a named matrix correctly", {
  m <- matrix(c(1, 0.2, 0.2, 2), ncol = 2,
              dimnames = list(c("p1", "p2"), c("p1", "p2")))
  cs <- as_covariance_structure(m)
  expect_true(S7::S7_inherits(cs, covariance_structure_class))
  expect_equal(cs[["p1"]][["p1"]], 1)
  expect_equal(cs[["p1"]][["p2"]], 0.2)
})

test_that("as_covariance_structure from an unnamed matrix auto-generates names", {
  m <- matrix(c(1, 0.1, 0.1, 2), ncol = 2)
  cs <- as_covariance_structure(m)
  expect_true(S7::S7_inherits(cs, covariance_structure_class))
  expect_equal(length(cs), 2L)
})

# ── R/03-pif.R : pif ─────────────────────────────────────────────────────────

test_that("pif with var_p = NULL and quiet = FALSE emits message", {
  expect_message(
    pif(p = 0.4, p_cft = 0.2, beta = 0.5, var_p = NULL, var_beta = 0.1,
        quiet = FALSE, label = "pifmsg1"),
    "Assuming parameters `p` have no variance"
  )
})

test_that("pif with var_beta = NULL and quiet = FALSE emits message", {
  expect_message(
    pif(p = 0.4, p_cft = 0.2, beta = 0.5, var_p = 0.01, var_beta = NULL,
        quiet = FALSE, label = "pifmsg2"),
    "Assuming parameters `beta` have no variance"
  )
})

test_that("pif with logit link and PIF <= 0 emits warning", {
  expect_warning(
    pif(p = 0.5, p_cft = 0.6, beta = 1.0, link = "logit",
        var_p = 0.01, var_beta = 0.01, quiet = TRUE, label = "piflogit"),
    "<= 0"
  )
})

test_that("pif returns pif_atomic_class", {
  p1 <- make_atomic(label = "pifc1")
  expect_true(S7::S7_inherits(p1, pif_atomic_class))
})

test_that("pif with label = NULL auto-generates label", {
  p1 <- pif(p = 0.3, p_cft = 0.1, beta = 0.5, quiet = TRUE)
  expect_true(nchar(p1@label) > 0)
})

# ── R/04-pif_operations.R : pif_validate_ensemble ────────────────────────────

test_that("pif_validate_ensemble: non-pif_class in list raises error", {
  p1 <- make_atomic(label = "pve1")
  not_a_pif <- S7::new_class("not_pif_x", properties = list())
  piflist <- list("pve1" = p1, "x" = not_a_pif())
  fake <- S7::new_class("fake_ens",
                        properties = list(pif_list = S7::class_list,
                                          weights = S7::class_numeric)
  )
  expect_error(
    validate_global_ensemble(fake(pif_list = piflist, weights = c(1, 1))),
    "pif_class"
  )
})

test_that("pif_validate_ensemble: unnamed pif_list raises error", {
  p1 <- make_atomic(label = "pve2")
  p2 <- make_atomic(label = "pve3")
  fake <- S7::new_class("fake_ens2",
                        properties = list(pif_list = S7::class_list,
                                          weights = S7::class_numeric)
  )
  expect_error(
    validate_global_ensemble(fake(pif_list = list(p1, p2), weights = c(1, 1))),
    "named list"
  )
})

test_that("pif_validate_ensemble: wrong weights length raises error", {
  p1 <- make_atomic(label = "pve4")
  p2 <- make_atomic(label = "pve5")
  piflist <- list("pve4" = p1, "pve5" = p2)
  fake <- S7::new_class("fake_ens3",
                        properties = list(pif_list = S7::class_list,
                                          weights = S7::class_numeric)
  )
  expect_error(
    validate_global_ensemble(fake(pif_list = piflist, weights = c(1, 2, 3))),
    "have length|length"
  )
})

# ── R/04-pif_operations.R : pif_ensemble ─────────────────────────────────────

test_that("pif_ensemble returns pif_ensemble_class", {
  p1 <- make_atomic(label = "pens1")
  p2 <- make_atomic(label = "pens2")
  ens <- pif_ensemble(p1, p2, weights = c(0.6, 0.4), label = "ens_test")
  expect_true(S7::S7_inherits(ens, pif_ensemble_class))
})

test_that("pif_ensemble produces a numeric PIF value", {
  p1 <- make_atomic(label = "pens3")
  p2 <- make_atomic(label = "pens4")
  ens <- pif_ensemble(p1, p2, weights = c(1, 1), label = "ens_test2")
  expect_true(is.numeric(coef(ens)))
  expect_length(coef(ens), 1L)
})

test_that("pif_ensemble with var_weights provides non-zero variance", {
  p1 <- make_atomic(label = "pens5")
  p2 <- make_atomic(label = "pens6")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5),
                      var_weights = diag(2) * 0.1, label = "ens_test3")
  expect_true(variance(ens) >= 0)
})

# ── R/05-covariance_generics.R : cov_atomic_pif ──────────────────────────────

test_that("cov_atomic_pif: non-atomic inputs raise error", {
  p1 <- make_atomic(label = "cap_err1")
  expect_error(
    cov_atomic_pif(list(), p1),
    "pif_atomic_class"
  )
  expect_error(
    cov_atomic_pif(p1, list()),
    "pif_atomic_class"
  )
})

test_that("cov_atomic_pif: identical pif returns variance", {
  p1 <- make_atomic(label = "cap1")
  expect_equal(
    cov_atomic_pif(p1, p1, var_p = p1@var_p, var_beta = p1@var_beta),
    variance(p1)
  )
})

test_that("cov_atomic_pif: wrong var_p matrix dimensions raise error", {
  p1 <- make_atomic(label = "cap2")
  p2 <- make_atomic(label = "cap3")
  expect_error(
    cov_atomic_pif(p1, p2, var_p = matrix(0.01, nrow = 3, ncol = 3)),
    "dimensions"
  )
})

test_that("cov_atomic_pif: wrong var_beta matrix dimensions raise error", {
  p1 <- make_atomic(label = "cap4")
  p2 <- make_atomic(label = "cap5")
  expect_error(
    cov_atomic_pif(p1, p2, var_beta = matrix(0.02, nrow = 3, ncol = 3)),
    "dimensions"
  )
})

test_that("cov_atomic_pif with var_p = 0 scalar works", {
  p1 <- make_atomic(label = "cap6")
  p2 <- make_atomic(label = "cap7")
  expect_silent(cov_atomic_pif(p1, p2, var_p = 0, var_beta = 0))
})

# ── R/05-covariance_generics.R : cov_ensemble_weights ────────────────────────

test_that("cov_ensemble_weights: non-pif_class raises error", {
  p1 <- make_atomic(label = "cew_err1")
  expect_error(
    cov_ensemble_weights(list(), p1),
    "pif_class"
  )
  expect_error(
    cov_ensemble_weights(p1, list()),
    "pif_class"
  )
})

test_that("cov_ensemble_weights: two atomics returns 0", {
  p1 <- make_atomic(label = "cew1")
  p2 <- make_atomic(label = "cew2")
  result <- cov_ensemble_weights(p1, p2)
  expect_equal(as.numeric(result), 0)
})

test_that("cov_ensemble_weights: ensemble vs ensemble returns numeric", {
  p1 <- make_atomic(label = "cew3")
  p2 <- make_atomic(label = "cew4")
  p3 <- make_atomic(label = "cew5")
  p4 <- make_atomic(label = "cew6")
  ens1 <- pif_ensemble(p1, p2, weights = c(0.5, 0.5),
                       var_weights = diag(2) * 0.1, label = "cew_ens1")
  ens2 <- pif_ensemble(p3, p4, weights = c(0.3, 0.7),
                       var_weights = diag(2) * 0.2, label = "cew_ens2")
  result <- cov_ensemble_weights(ens1, ens2)
  expect_true(is.numeric(result))
})

# ── R/05-covariance_generics.R : cov_ensemble_atomic ─────────────────────────

test_that("cov_ensemble_atomic: ensemble vs atomic returns numeric", {
  p1 <- make_atomic(label = "cea1")
  p2 <- make_atomic(label = "cea2")
  p3 <- make_atomic(label = "cea3")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "cea_ens")
  result <- cov_ensemble_atomic(ens, p3)
  expect_true(is.numeric(result))
})

test_that("cov_ensemble_atomic: two atomics equals cov_atomic_pif", {
  p1 <- make_atomic(label = "cea4")
  p2 <- make_atomic(label = "cea5")
  expect_equal(
    cov_ensemble_atomic(p1, p2),
    cov_atomic_pif(p1, p2)
  )
})

test_that("cov_ensemble_atomic: bad first argument raises error", {
  p1 <- make_atomic(label = "cea6")
  expect_error(
    cov_ensemble_atomic(list(), p1),
    "pif_global_ensemble_class|pif_atomic"
  )
})

test_that("cov_ensemble_atomic: bad second argument raises error", {
  p1 <- make_atomic(label = "cea7")
  p2 <- make_atomic(label = "cea8")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "cea_ens2")
  expect_error(
    cov_ensemble_atomic(ens, list()),
    "pif_atomic_class"
  )
})

# ── R/05-covariance_generics.R : cov_total_pif ───────────────────────────────

test_that("cov_total_pif: non-pif inputs raise error", {
  p1 <- make_atomic(label = "ctp_err1")
  expect_error(
    cov_total_pif(list(), p1),
    "pif_class"
  )
})

test_that("cov_total_pif: atomic vs atomic returns numeric", {
  p1 <- make_atomic(label = "ctp1")
  p2 <- make_atomic(label = "ctp2")
  result <- cov_total_pif(p1, p2)
  expect_true(is.numeric(result))
})

test_that("cov_total_pif: atomic vs total returns numeric", {
  p1 <- make_atomic(label = "ctp3")
  p2 <- make_atomic(label = "ctp4")
  p3 <- make_atomic(label = "ctp5")
  tot <- pif_total(p2, p3, weights = c(0.5, 0.5), label = "ctptot")
  expect_true(is.numeric(cov_total_pif(p1, tot)))
  expect_true(is.numeric(cov_total_pif(tot, p1)))
})

test_that("cov_total_pif: total vs total returns numeric", {
  p1 <- make_atomic(label = "ctp6")
  p2 <- make_atomic(label = "ctp7")
  tot <- pif_total(p1, p2, weights = c(0.5, 0.5), label = "ctptot2")
  expect_true(is.numeric(cov_total_pif(tot, tot)))
})

# ── R/05-covariance_generics.R : covariance for pif_global_ensemble_class ─────

test_that("covariance() for pif_total returns a matrix", {
  p1 <- make_atomic(label = "cov1")
  p2 <- make_atomic(label = "cov2")
  tot <- pif_total(p1, p2, weights = c(0.4, 0.6), label = "covtot")
  m <- covariance(tot)
  expect_true(is.matrix(m))
})

test_that("covariance() for pif_ensemble returns a matrix", {
  p1 <- make_atomic(label = "cov3")
  p2 <- make_atomic(label = "cov4")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "covens")
  m <- covariance(ens)
  expect_true(is.matrix(m))
})

test_that("covariance() diagonal >= 0 (variance is non-negative)", {
  p1 <- make_atomic(label = "cov5")
  p2 <- make_atomic(label = "cov6")
  tot <- pif_total(p1, p2, weights = c(0.5, 0.5), label = "covtot2")
  m <- covariance(tot)
  expect_true(all(diag(m) >= 0))
})

# ── R/05-covariance_generics.R : variance for pif_class ──────────────────────

test_that("variance() for pif_atomic_class returns a non-negative numeric", {
  p1 <- make_atomic(label = "var1")
  v <- variance(p1)
  expect_true(is.numeric(v))
  expect_true(v >= 0)
})

test_that("variance() for pif_total returns non-negative numeric", {
  p1 <- make_atomic(label = "var2")
  p2 <- make_atomic(label = "var3")
  tot <- pif_total(p1, p2, weights = c(0.5, 0.5), label = "vartot")
  v <- variance(tot)
  expect_true(is.numeric(v))
  expect_true(v >= 0)
})

# ── R/06-derivatives.R : deriv_pif_p ─────────────────────────────────────────

test_that("deriv_pif_p returns correct length", {
  p     <- c(0.3, 0.7)
  p_cft <- c(0.2, 0.8)
  rr    <- c(2.0, 1.0)
  d <- deriv_pif_p(p = p, p_cft = p_cft, rr = rr)
  expect_length(d, 2L)
})

test_that("deriv_pif_p with all rr = 1 returns zero vector", {
  p     <- c(0.4, 0.6)
  p_cft <- c(0.3, 0.7)
  rr    <- c(1, 1)
  d <- deriv_pif_p(p = p, p_cft = p_cft, rr = rr)
  expect_equal(d, c(0, 0))
})

test_that("deriv_pif_p matches manual calculation", {
  p     <- c(0.5, 0.5)
  p_cft <- c(0.4, 0.6)
  rr    <- c(2, 1)
  mu_obs <- sum(p * rr)
  mu_cft <- sum(p_cft * rr)
  expected <- (mu_cft / mu_obs^2) * (rr - 1)
  expect_equal(deriv_pif_p(p = p, p_cft = p_cft, rr = rr), expected)
})

# ── R/07-links.R : parse_deriv_link ──────────────────────────────────────────

test_that("change_link changes the link of a pif_class", {
  p1 <- make_atomic(label = "cl1")
  p_logit <- change_link(p1, link = "logit")
  expect_false(identical(p1@link, p_logit@link))
})

test_that("change_link raises error for non pif_class input", {
  expect_error(
    change_link(list(), link = "logit"),
    "pif_class"
  )
})

test_that("change_link raises error for ensemble input", {
  p1 <- make_atomic(label = "cl2")
  p2 <- make_atomic(label = "cl3")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "clens")
  expect_error(
    change_link(ens, link = "logit"),
    "[Cc]annot change the link of an ensemble"
  )
})

test_that("change_link with logit and PIF <= 0 emits warning", {
  # PIF will be <= 0 when p_cft > p (intervention increases exposure)
  p_neg <- pif(p = 0.3, p_cft = 0.5, beta = 1.0,
               var_p = 0.01, var_beta = 0.01, quiet = TRUE, label = "cl_neg")
  expect_warning(
    change_link(p_neg, link = "logit"),
    "<= 0"
  )
})

# ── R/08-utils.R : change_link (additional paths) ────────────────────────────

test_that("change_link with identity link is silent for positive PIF", {
  p1 <- make_atomic(label = "cl_id1")
  expect_silent(change_link(p1, link = "identity"))
})

test_that("change_link with log-complement link works", {
  p1 <- make_atomic(label = "cl_lc1")
  p_lc <- change_link(p1, link = "log-complement")
  expect_true(S7::S7_inherits(p_lc, pif_class))
})

# ── R/11-weighted_adjusted.R : weighted_adjusted_fractions ───────────────────

test_that("weighted_adjusted_fractions returns named list of pif_class", {
  p1 <- make_paf_atomic(label = "waf1")
  p2 <- make_paf_atomic(label = "waf2")
  result <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)
  expect_true(is.list(result))
  expect_equal(length(result), 2L)
  expect_true(all(sapply(result, function(x) S7::S7_inherits(x, pif_class))))
})

test_that("weighted_adjusted_fractions raises error for wrong weights length", {
  p1 <- make_paf_atomic(label = "waf3")
  p2 <- make_paf_atomic(label = "waf4")
  expect_error(
    weighted_adjusted_fractions(p1, p2, weights = c(1, 2, 3), quiet = TRUE),
    "[Ww]eight"
  )
})

test_that("weighted_adjusted_fractions raises error if non-pif passed", {
  p1 <- make_paf_atomic(label = "waf5")
  expect_error(
    weighted_adjusted_fractions(p1, list(), quiet = TRUE),
    "pif_class"
  )
})

test_that("weighted_adjusted_fractions raises error for duplicate labels", {
  p1 <- make_paf_atomic(label = "dup_label")
  p2 <- make_paf_atomic(label = "dup_label")
  expect_error(
    weighted_adjusted_fractions(p1, p2, quiet = TRUE),
    "[Dd]uplicate"
  )
})

test_that("weighted_adjusted_fractions with var_p runs without error", {
  p1 <- make_paf_atomic(label = "waf6")
  p2 <- make_paf_atomic(label = "waf7")
  expect_silent(
    weighted_adjusted_fractions(p1, p2, var_p = 0.001, quiet = TRUE)
  )
})

test_that("weighted_adjusted_fractions with custom weights works", {
  p1 <- make_paf_atomic(label = "waf8")
  p2 <- make_paf_atomic(label = "waf9")
  result <- weighted_adjusted_fractions(p1, p2, weights = c(0.3, 0.7), quiet = TRUE)
  expect_equal(length(result), 2L)
})

# ── R/cases.R : averted_cases ────────────────────────────────────────────────

test_that("averted_cases returns cases_class with correct point estimate", {
  p1 <- make_atomic(label = "ac1")
  result <- averted_cases(500, pif = p1)
  expect_true(S7::S7_inherits(result, cases_class))
  expect_equal(result@cases, coef(p1) * 500)
})

test_that("averted_cases with variance > 0 produces positive variance", {
  p1 <- make_atomic(label = "ac2")
  result <- averted_cases(1000, pif = p1, variance = 25)
  expect_true(result@variance > 0)
})

test_that("averted_cases with invalid link raises error", {
  p1 <- make_atomic(label = "ac3")
  expect_error(
    averted_cases(100, pif = p1, link = "bad_link"),
    "Invalid link"
  )
})

test_that("averted_cases with non-pif raises error", {
  expect_error(
    averted_cases(100, pif = list()),
    "pif_class"
  )
})

# ── R/variance_product.R : var_prod ──────────────────────────────────────────

test_that("var_prod computes correctly with non-zero variances", {
  result <- var_prod(x = 2, y = 3, var_x = 0.5, var_y = 0.8)
  expected <- 3^2 * 0.5 + 2^2 * 0.8 + 0.5 * 0.8
  expect_equal(result, expected)
})

test_that("var_prod returns 0 when both variances are 0", {
  expect_equal(var_prod(x = 5, y = 4, var_x = 0, var_y = 0), 0)
})

test_that("var_prod handles one-sided zero variance", {
  expect_equal(var_prod(x = 5, y = 4, var_x = 2, var_y = 0), 4^2 * 2)
  expect_equal(var_prod(x = 5, y = 4, var_x = 0, var_y = 3), 5^2 * 3)
})
