# Helpers ----------------------------------------------------------------------

make_paf2 <- function(p = 0.2, beta = 2.2, var_p = 0.001, var_beta = 0.01,
                      label = paste0("paf_", round(abs(stats::rnorm(1)), 5))) {
  paf(p = p, beta = beta, var_p = var_p, var_beta = var_beta,
      quiet = TRUE, label = label)
}

make_pif2 <- function(p = 0.2, p_cft = 0.1, beta = log(2.2),
                      var_p = 0.001, var_beta = 0.01,
                      label = paste0("pif_", round(abs(stats::rnorm(1)), 5))) {
  pif(p = p, p_cft = p_cft, beta = beta, var_p = var_p, var_beta = var_beta,
      quiet = TRUE, label = label)
}

# weighted_adjusted_fractions — return structure --------------------------------

test_that("weighted_adjusted_fractions returns a named list of pif_class objects", {
  p1  <- make_paf2(label = "Lead")
  p2  <- make_paf2(p = 0.1, beta = 1.2, label = "Radiation")

  adj <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)

  expect_type(adj, "list")
  expect_length(adj, 2L)
  expect_named(adj, c("Lead", "Radiation"))
  expect_true(S7::S7_inherits(adj[["Lead"]],      pif_class))
  expect_true(S7::S7_inherits(adj[["Radiation"]], pif_class))
})

test_that("weighted_adjusted_fractions works with a single fraction", {
  p1  <- make_paf2(label = "Solo")
  adj <- weighted_adjusted_fractions(p1, quiet = TRUE)
  expect_length(adj, 1L)
  expect_named(adj, "Solo")
  expect_true(S7::S7_inherits(adj[["Solo"]], pif_class))
})

test_that("adjusted fractions are pif_class with non-negative variance", {
  p1  <- make_paf2(label = "A")
  p2  <- make_paf2(p = 0.15, label = "B")
  p3  <- make_paf2(p = 0.05, beta = 1.5, label = "C")

  adj <- weighted_adjusted_fractions(p1, p2, p3, quiet = TRUE)

  for (nm in c("A", "B", "C")) {
    expect_true(adj[[nm]]@variance >= 0,
                info = paste("Non-negative variance for", nm))
    expect_true(adj[[nm]]@pif > 0,
                info = paste("Positive PIF for", nm))
    expect_true(adj[[nm]]@pif < 1,
                info = paste("PIF < 1 for", nm))
  }
})

# Proportionality property -----------------------------------------------------

test_that("adjusted fractions are proportional to the original fractions", {
  p1  <- make_paf2(label = "X")
  p2  <- make_paf2(p = 0.1, beta = 1.5, label = "Y")

  adj <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)

  ratio_X <- (adj[["X"]]@pif / coef(p1)) |> as.numeric()
  ratio_Y <- (adj[["Y"]]@pif / coef(p2)) |> as.numeric()

  # Both fractions are scaled by the same factor
  expect_equal(ratio_X, ratio_Y, tolerance = 1e-10)
})

test_that("adjusted fractions sum to the ensemble pif", {
  p1  <- make_paf2(label = "S1")
  p2  <- make_paf2(p = 0.12, beta = 1.8, label = "S2")

  adj     <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)
  pif_ens <- paf_ensemble(p1, p2, quiet = TRUE)

  adj_sum <- sum(sapply(adj, coef))

  expect_equal(adj_sum, coef(pif_ens), tolerance = 1e-10)
})

# Label management -------------------------------------------------------------

test_that("adjusted output labels get '_adj' suffix on pif slot", {
  p1  <- make_paf2(label = "Lead")
  p2  <- make_paf2(label = "Rad")
  adj <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)

  expect_equal(adj[["Lead"]]@label, "Lead_adj")
  expect_equal(adj[["Rad"]]@label,  "Rad_adj")
})

test_that("weighted_adjusted_fractions errors on duplicate labels", {
  p1 <- make_paf2(label = "dup")
  p2 <- make_paf2(label = "dup")

  expect_error(
    weighted_adjusted_fractions(p1, p2, quiet = TRUE),
    "Duplicate labels"
  )
})

test_that("custom label_ensemble and label_sum are accepted silently", {
  p1  <- make_paf2(label = "L1")
  p2  <- make_paf2(label = "L2")

  expect_silent(
    weighted_adjusted_fractions(p1, p2,
                                label_ensemble = "my_ensemble",
                                label_sum      = "my_sum",
                                quiet = TRUE)
  )
})

# Input validation -------------------------------------------------------------

test_that("weighted_adjusted_fractions errors on non-pif_class input", {
  p1 <- make_paf2(label = "ok")
  expect_error(
    weighted_adjusted_fractions(p1, list(), quiet = TRUE),
    "pif_class"
  )
  expect_error(
    weighted_adjusted_fractions(42, quiet = TRUE),
    "pif_class"
  )
})

# Weights ----------------------------------------------------------------------

test_that("custom weights change the adjusted values", {
  p1   <- make_paf2(label = "W1")
  p2   <- make_paf2(p = 0.1, beta = 1.5, label = "W2")

  adj_equal  <- weighted_adjusted_fractions(p1, p2,
                                            weights = c(1, 1),
                                            quiet = TRUE)
  adj_skewed <- weighted_adjusted_fractions(p1, p2,
                                            weights = c(0.9, 0.1),
                                            quiet = TRUE)

  # Different ensemble weights should lead to different adjusted values
  expect_false(isTRUE(all.equal(adj_equal[["W1"]]@pif, adj_skewed[["W1"]]@pif)))
})

test_that("weights with var_weights runs without error", {
  p1 <- make_paf2(label = "Wv1")
  p2 <- make_paf2(label = "Wv2")

  vw <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)
  expect_silent(
    weighted_adjusted_fractions(p1, p2,
                                weights     = c(0.6, 0.4),
                                var_weights = vw,
                                quiet = TRUE)
  )
})

# Correlated prevalences / betas -----------------------------------------------

test_that("var_p argument is accepted and changes variance estimates", {
  p1 <- make_paf2(label = "vp1")
  p2 <- make_paf2(label = "vp2")

  adj_indep <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)
  adj_corr  <- weighted_adjusted_fractions(p1, p2, var_p = 0.0005, quiet = TRUE)

  # Point estimates are unaffected by var_p
  expect_equal(adj_indep[["vp1"]]@pif, adj_corr[["vp1"]]@pif, tolerance = 1e-10)

  # But variances may differ
  # (we just check they are both non-negative)
  expect_true(adj_corr[["vp1"]]@variance >= 0)
  expect_true(adj_corr[["vp2"]]@variance >= 0)
})

test_that("var_beta argument is accepted and changes variance estimates", {
  p1 <- make_paf2(label = "vb1")
  p2 <- make_paf2(label = "vb2")

  adj <- weighted_adjusted_fractions(p1, p2, var_beta = 0.005, quiet = TRUE)
  expect_true(adj[["vb1"]]@variance >= 0)
  expect_true(adj[["vb2"]]@variance >= 0)
})

# conf_level -------------------------------------------------------------------

test_that("conf_level is propagated to the output pif_class objects", {
  p1 <- make_paf2(label = "cl1")
  p2 <- make_paf2(label = "cl2")

  adj90 <- weighted_adjusted_fractions(p1, p2, conf_level = 0.90, quiet = TRUE)
  adj99 <- weighted_adjusted_fractions(p1, p2, conf_level = 0.99, quiet = TRUE)

  expect_equal(adj90[["cl1"]]@conf_level, 0.90)
  expect_equal(adj99[["cl2"]]@conf_level, 0.99)
})

# Three or more fractions ------------------------------------------------------

test_that("weighted_adjusted_fractions handles three fractions correctly", {
  p1  <- make_paf2(label = "T1")
  p2  <- make_paf2(p = 0.15, beta = 1.8, label = "T2")
  p3  <- make_paf2(p = 0.05, beta = 1.4, label = "T3")

  adj <- weighted_adjusted_fractions(p1, p2, p3, quiet = TRUE)

  expect_length(adj, 3L)
  expect_named(adj, c("T1", "T2", "T3"))

  # All positive and < 1
  for (nm in c("T1", "T2", "T3")) {
    expect_true(adj[[nm]]@pif > 0)
    expect_true(adj[[nm]]@pif < 1)
    expect_true(adj[[nm]]@variance >= 0)
  }

  # Sum equals ensemble
  ens     <- paf_ensemble(p1, p2, p3, quiet = TRUE)
  adj_sum <- sum(sapply(adj, coef))
  expect_equal(adj_sum, coef(ens), tolerance = 1e-10)
})

test_that("negative weights fail", {

  p1  <- make_paf2(p = 0.15, beta = 1.8, label = "p1")
  p2  <- make_paf2(p = 0.15, beta = 1.8, label = "p2")
  p3  <- make_paf2(p = 0.15, beta = 1.8, label = "p3")


  expect_error(
    adj <- weighted_adjusted_fractions(p1, p2, p3, quiet = TRUE, weights = rep(0, 3)),
    "strictly positive"
  )

  expect_error(
    adj <- weighted_adjusted_fractions(p1, p2, p3, quiet = TRUE, weights = c(1, -1, 0.4)),
    "strictly positive"
  )

})

test_that("weighted_adjusted_fractions handles same fraction correctly", {
  p1  <- paf(p = 0.15, beta = 1.8, label = "p1", quiet = TRUE, var_p = 0.001)
  p2  <- paf(p = 0.15, beta = 1.8, label = "p2", quiet = TRUE, var_p = 0.001)
  p3  <- paf(p = 0.15, beta = 1.8, label = "p3", quiet = TRUE, var_p = 0.001)

  adj <- weighted_adjusted_fractions(p1, p2, p3, quiet = TRUE, weights = rep(1/coef(p1), 3))

  expect_length(adj, 3L)

  #Now each of the adjusted should be 1/3
  expect_equal(as.numeric(coef(adj$p1)), 1/3)
  expect_equal(as.numeric(coef(adj$p2)), 1/3)
  expect_equal(as.numeric(coef(adj$p3)), 1/3)

  #Modify a bit to check equality
  adj$p1@label <- "p2_adj"
  names(adj$p1@pif) <- "p2"
  expect_equal(adj$p1, adj$p2)


})

# weighted_adjusted_paf --------------------------------------------------------

test_that("weighted_adjusted_paf returns named list of PAF pif_class", {
  p1  <- make_paf2(label = "PafA")
  p2  <- make_paf2(label = "PafB")

  adj <- weighted_adjusted_paf(p1, p2, quiet = TRUE)

  expect_type(adj, "list")
  expect_length(adj, 2L)
  expect_named(adj, c("PafA", "PafB"))
  expect_true(S7::S7_inherits(adj[["PafA"]], pif_class))
})

test_that("weighted_adjusted_paf errors when first argument is a PIF", {
  pi1 <- make_pif2(label = "not_a_paf")
  expect_error(
    weighted_adjusted_paf(pi1, quiet = TRUE),
    "PAF"
  )
})

test_that("weighted_adjusted_paf produces non-negative variances", {
  p1  <- make_paf2(label = "Pa1")
  p2  <- make_paf2(p = 0.1, label = "Pa2")

  adj <- weighted_adjusted_paf(p1, p2, quiet = TRUE)
  expect_true(adj[["Pa1"]]@variance >= 0)
  expect_true(adj[["Pa2"]]@variance >= 0)
})

# weighted_adjusted_pif --------------------------------------------------------

test_that("weighted_adjusted_pif returns named list of PIF pif_class", {
  pi1 <- make_pif2(label = "PifA")
  pi2 <- make_pif2(p = 0.15, p_cft = 0.05, label = "PifB")

  adj <- weighted_adjusted_pif(pi1, pi2, quiet = TRUE)

  expect_type(adj, "list")
  expect_length(adj, 2L)
  expect_named(adj, c("PifA", "PifB"))
  expect_true(S7::S7_inherits(adj[["PifA"]], pif_class))
})

test_that("weighted_adjusted_pif produces non-negative variances", {
  pi1 <- make_pif2(label = "Pv1")
  pi2 <- make_pif2(p = 0.1, p_cft = 0.03, label = "Pv2")

  adj <- weighted_adjusted_pif(pi1, pi2, quiet = TRUE)
  expect_true(adj[["Pv1"]]@variance >= 0)
  expect_true(adj[["Pv2"]]@variance >= 0)
})

test_that("weighted_adjusted_pif adjusted values are proportional to originals", {
  pi1 <- make_pif2(label = "Pr1")
  pi2 <- make_pif2(p = 0.12, p_cft = 0.04, label = "Pr2")

  adj    <- weighted_adjusted_pif(pi1, pi2, quiet = TRUE)
  ratio1 <- (adj[["Pr1"]]@pif / coef(pi1)) |> as.numeric()
  ratio2 <- (adj[["Pr2"]]@pif / coef(pi2)) |> as.numeric()

  expect_equal(ratio1, ratio2, tolerance = 1e-10)
})

test_that("weighted_adjusted_pif sum equals ensemble pif", {
  pi1 <- make_pif2(label = "Ps1")
  pi2 <- make_pif2(p = 0.08, p_cft = 0.02, label = "Ps2")

  adj <- weighted_adjusted_pif(pi1, pi2, quiet = TRUE)
  ens <- pif_ensemble(pi1, pi2, quiet = TRUE)

  adj_sum <- sum(sapply(adj, coef))
  expect_equal(adj_sum, coef(ens), tolerance = 1e-10)
})

# Consistency between wrappers and weighted_adjusted_fractions -----------------

test_that("weighted_adjusted_paf matches weighted_adjusted_fractions for PAFs", {
  p1 <- make_paf2(label = "C1")
  p2 <- make_paf2(label = "C2")

  via_paf <- weighted_adjusted_paf(p1, p2, quiet = TRUE)
  via_gen <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)

  expect_equal(via_paf[["C1"]]@pif,      via_gen[["C1"]]@pif,      tolerance = 1e-12)
  expect_equal(via_paf[["C1"]]@variance, via_gen[["C1"]]@variance,  tolerance = 1e-12)
  expect_equal(via_paf[["C2"]]@pif,      via_gen[["C2"]]@pif,      tolerance = 1e-12)
})

test_that("weighted_adjusted_pif matches weighted_adjusted_fractions for PIFs", {
  pi1 <- make_pif2(label = "Cm1")
  pi2 <- make_pif2(p = 0.1, p_cft = 0.03, label = "Cm2")

  via_pif <- weighted_adjusted_pif(pi1, pi2, quiet = TRUE)
  via_gen <- weighted_adjusted_fractions(pi1, pi2, quiet = TRUE)

  expect_equal(via_pif[["Cm1"]]@pif,      via_gen[["Cm1"]]@pif,     tolerance = 1e-12)
  expect_equal(via_pif[["Cm1"]]@variance, via_gen[["Cm1"]]@variance, tolerance = 1e-12)
})

# Zero variance edge-case ------------------------------------------------------

test_that("weighted_adjusted_fractions handles zero variance fractions", {
  p1 <- paf(0.2, 2.2, var_p = 0, var_beta = 0, quiet = TRUE, label = "Zv1")
  p2 <- paf(0.1, 1.5, var_p = 0, var_beta = 0, quiet = TRUE, label = "Zv2")

  adj <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)
  expect_equal(adj[["Zv1"]]@variance, 0)
  expect_equal(adj[["Zv2"]]@variance, 0)
})
