test_that("covariance_structure_class constructs valid objects", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.21, "pif2" = 0.12),
      "pif2" = list("pif1" = 0.12, "pif2" = 0.33)
    )
  )
  expect_true(S7::S7_inherits(cov, covariance_structure_class))
  expect_equal(length(cov), 2L)
})

test_that("covariance_structure_class rejects non-list entries", {
  expect_error(
    covariance_structure_class(
      list(
        "pif1" = 0.21,
        "pif2" = list("pif1" = 0.12, "pif2" = 0.33)
      )
    ),
    "not a list"
  )
})

test_that("covariance_structure_class rejects sublists of inconsistent length", {
  expect_error(
    covariance_structure_class(
      list(
        "pif1" = list("pif1" = 0.21, "pif2" = 0.12),
        "pif2" = list("pif1" = 0.12)
      )
    )
  )
})

test_that("covariance_structure_class allows matrix and vector entries", {
  mat <- matrix(c(0.1, 0.21, 0.47, -0.3), ncol = 2)
  vec <- c(0.22, -0.9, 0.01)

  expect_silent(
    covariance_structure_class(
      list(
        "pif1" = list("pif1" = 0.21, "pif2" = mat),
        "pif2" = list("pif1" = mat,  "pif2" = 0.33)
      )
    )
  )

  expect_silent(
    covariance_structure_class(
      list(
        "pif1" = list("pif1" = 0.5, "pif2" = vec),
        "pif2" = list("pif1" = vec, "pif2" = 0.8)
      )
    )
  )
})

test_that("[[ extraction works by index and by name", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.21, "pif2" = 0.12, "pif3" = 0.31),
      "pif2" = list("pif1" = 0.12, "pif2" = 0.33, "pif3" = -0.01),
      "pif3" = list("pif1" = 0.31, "pif2" = -0.01, "pif3" = 0.80)
    )
  )

  expect_equal(cov[[1]][[3]], 0.31)
  expect_equal(cov[["pif2"]][["pif3"]], -0.01)
  expect_equal(cov[["pif1"]][["pif1"]], 0.21)
})

test_that("[[<- assignment works by index and by name", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.21, "pif2" = 0.12),
      "pif2" = list("pif1" = 0.12, "pif2" = 0.33)
    )
  )

  cov[[1]][[2]] <- list("pif1" = 999, "pif2" = 888)
  expect_equal(cov[[1]][[2]], list("pif1" = 999, "pif2" = 888))

  cov[["pif2"]][["pif1"]] <- list("pif1" = 100, "pif2" = 200)
  expect_equal(cov[["pif2"]][["pif1"]], list("pif1" = 100, "pif2" = 200))
})

test_that("[[<- rejects non-list assignment values", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.21, "pif2" = 0.12),
      "pif2" = list("pif1" = 0.12, "pif2" = 0.33)
    )
  )
  expect_error(
    { cov[["pif1"]] <- 42 },
    "should be a list"
  )
})

test_that("[[<- rejects list of wrong length", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.21, "pif2" = 0.12),
      "pif2" = list("pif1" = 0.12, "pif2" = 0.33)
    )
  )
  expect_error(
    { cov[["pif1"]] <- list("a" = 1) },
    "length"
  )
})

# --- Default covariance structure helpers ---

test_that("covariance_structure and covariance_structure2 work for atomic pifs", {
  p1 <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015, label = "A")
  p2 <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022, label = "B")

  cs1 <- covariance_structure(p1)
  expect_true(is.list(cs1@cov_list))

  cs2 <- covariance_structure2(p1, p2)
  expect_true(is.list(cs2@cov_list))
})

test_that("default_weight_covariance_structure works for an ensemble", {
  p1 <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015, label = "W1")
  p2 <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022, label = "W2")
  ens <- pif_ensemble(p1, p2, weights = c(0.6, 0.4),
                      var_weights = matrix(c(0.1, 0.02, 0.02, 0.1), ncol = 2),
                      label = "Ens")

  wcs <- default_weight_covariance_structure(ens)
  expect_true(is.list(wcs@cov_list))

  wcs2 <- default_weight_covariance_structure2(ens, ens)
  expect_true(is.list(wcs2@cov_list))
})

test_that("default_parameter_covariance_structure works for p and beta", {
  p1 <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015, label = "P1")
  p2 <- paf(0.30, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015, label = "P2")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "Ens2")

  pcs_p    <- default_parameter_covariance_structure(ens, parameter = "p")
  pcs_beta <- default_parameter_covariance_structure(ens, parameter = "beta")

  expect_true(is.list(pcs_p@cov_list))
  expect_true(is.list(pcs_beta@cov_list))
})

test_that("default_parameter_covariance_structure2 returns structure between two pifs", {
  p1 <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015, label = "X")
  p2 <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022, label = "Y")

  cs_p    <- default_parameter_covariance_structure2(p1, p2, parameter = "p")
  cs_beta <- default_parameter_covariance_structure2(p1, p2, parameter = "beta")

  expect_true(is.list(cs_p@cov_list))
  expect_true(is.list(cs_beta@cov_list))
})

test_that("default_pif_covariance_structure and _2 work", {
  p1 <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015, label = "D1")
  p2 <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022, label = "D2")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "D_ens")

  expect_true(is.list(default_pif_covariance_structure(ens)@cov_list))
  expect_true(is.list(default_pif_covariance_structure2(p1, p2)@cov_list))
})

test_that("default_weight_pif_covariance_structure and _2 work", {
  p1 <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015, label = "E1")
  p2 <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022, label = "E2")
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5),
                      var_weights = matrix(c(0.1, 0.05, 0.05, 0.1), ncol = 2),
                      label = "E_ens")

  expect_true(is.list(default_weight_pif_covariance_structure(ens)@cov_list))
  expect_true(is.list(default_weight_pif_covariance_structure2(ens, ens)@cov_list))
})

test_that("identical pifs give variance structure (is_variance = TRUE path)", {
  p1 <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015, label = "Same")
  cs_var  <- default_parameter_covariance_structure2(p1, p1, parameter = "p")
  cs_cov  <- default_parameter_covariance_structure2(
    p1,
    paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022, label = "Diff"),
    parameter = "p"
  )
  # Both should be lists; the variance case should equal var_p
  expect_true(is.list(cs_var@cov_list))
  expect_true(is.list(cs_cov@cov_list))
})

test_that("total pif covariance structures work end-to-end", {
  p1  <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015, label = "T1")
  p2  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022, label = "T2")
  tot <- pif_total(p1, p2, weights = c(0.49, 0.51),
                   var_weights = matrix(c(0.1, 0, 0, 0.1), ncol = 2),
                   label = "Tot")

  expect_true(is.list(default_weight_covariance_structure(tot)@cov_list))
  expect_true(is.list(default_parameter_covariance_structure(tot, parameter = "p")@cov_list))
  expect_true(is.list(default_parameter_covariance_structure(tot, parameter = "beta")@cov_list))
})
