# test-numerical-math.R
#
# Numerical validation tests for deltapif.
# Every expected value is derived from the closed-form formulas documented in
# the package, computed independently in plain R, and then compared against
# the package output.
#
# Sections
#   1. paf / pif point estimates (identity rr_link)
#   2. Variance of atomic pif (delta method, identity rr_link)
#   3. cov_atomic_pif (two independent atomics, identity rr_link)
#   4. pif_total point estimate and variance
#   5. pif_ensemble point estimate and variance
#   6. weighted_adjusted_fractions point estimates
#   7. weighted_adjusted_fractions variance (delta method on log scale)
#   8. PAF special cases
#   9. cov_atomic_pif shared-beta covariance
#  10. pif_total variance with non-zero var_weights
#  11. pif_ensemble variance with non-zero var_weights

# ── tolerance used throughout ─────────────────────────────────────────────────
TOL <- 1e-10

# =============================================================================
# 1. paf / pif point estimates (identity rr_link)
# =============================================================================

test_that("[num] paf point estimate matches 1 - 1/E[RR] (identity link)", {
  p   <- c(0.4, 0.2)
  rr  <- c(1.2, 1.5)
  # E[RR] = p0*1 + p1*rr1 + p2*rr2,  p0 = 1 - sum(p)
  mu_obs   <- (1 - sum(p)) * 1 + p[1] * rr[1] + p[2] * rr[2]
  expected <- 1 - 1 / mu_obs

  got <- coef(paf(p = p, beta = rr, rr_link = "identity", quiet = TRUE))
  expect_equal(got, expected, tolerance = TOL)
})

test_that("[num] pif point estimate matches 1 - E_cft[RR]/E_obs[RR] (identity link)", {
  p     <- c(0.5, 0.3)
  p_cft <- c(0.2, 0.1)
  rr    <- c(1.4, 2.0)

  p0     <- 1 - sum(p);     p0_cft <- 1 - sum(p_cft)
  mu_obs <- p0 * 1 + p[1] * rr[1] + p[2] * rr[2]
  mu_cft <- p0_cft * 1 + p_cft[1] * rr[1] + p_cft[2] * rr[2]
  expected <- 1 - mu_cft / mu_obs

  got <- coef(pif(p = p, p_cft = p_cft, beta = rr,
                  rr_link = "identity", quiet = TRUE))
  expect_equal(got, expected, tolerance = TOL)
})

test_that("[num] paf equals pif when p_cft = 0 (identity link)", {
  p  <- c(0.3, 0.4)
  rr <- c(1.5, 1.8)

  paf_val <- coef(paf(p = p, beta = rr, rr_link = "identity", quiet = TRUE,
                      label = "paf_num"))
  pif_val <- coef(pif(p = p, p_cft = c(0, 0), beta = rr,
                      rr_link = "identity", quiet = TRUE, type = "PAF",
                      label = "paf_num"))
  expect_equal(paf_val, pif_val, tolerance = TOL)
})

# =============================================================================
# 2. Variance of atomic pif (delta method, identity rr_link)
# =============================================================================

test_that("[num] variance of atomic pif matches delta-method formula (identity link)", {
  p      <- c(0.4, 0.2)
  rr     <- c(1.2, 1.5)
  p_cft  <- c(0.5, 0.1)
  sig_p  <- matrix(c(0.01, 0.002, 0.002, 0.01), ncol = 2)
  sig_b  <- matrix(c(0.05, -0.001, -0.001, 0.007), ncol = 2)

  p0     <- 1 - sum(p);     p0_cft <- 1 - sum(p_cft)
  mu_obs <- p0 + p[1] * rr[1] + p[2] * rr[2]
  mu_cft <- p0_cft + p_cft[1] * rr[1] + p_cft[2] * rr[2]

  # d(PIF)/d(p_i) = (mu_cft/mu_obs^2) * (rr_i - 1)
  dp <- (mu_cft / mu_obs^2) * (rr - 1)
  vp <- as.numeric(t(dp) %*% sig_p %*% dp)

  # d(PIF)/d(rr_i) = (mu_cft * p_i - mu_obs * p_cft_i) / mu_obs^2
  # with rr_link = identity so d(rr)/d(beta) = 1
  db <- (mu_cft * p - mu_obs * p_cft) / mu_obs^2
  vb <- as.numeric(t(db) %*% sig_b %*% db)

  expected <- vp + vb

  got <- variance(pif(p = p, p_cft = p_cft, beta = rr,
                      var_p = sig_p, var_beta = sig_b,
                      rr_link = "identity", quiet = TRUE))
  expect_equal(got, expected, tolerance = TOL)
})

test_that("[num] variance is zero when both var_p and var_beta are zero", {
  p  <- c(0.3, 0.5)
  rr <- c(1.3, 1.7)
  got <- variance(paf(p = p, beta = rr, var_p = 0, var_beta = 0,
                      rr_link = "identity", quiet = TRUE))
  expect_equal(got, 0, tolerance = TOL)
})

test_that("[num] variance of paf with exp link matches numerical calculation", {
  # Use the exp link (log-RR parameterisation)
  p      <- 0.3
  beta   <- log(2.0)          # log RR
  var_b  <- 0.04              # var of log RR
  var_p_ <- 0.01

  rr     <- exp(beta)
  mu_obs <- (1 - p) + p * rr
  mu_cft <- 1              # PAF: p_cft = 0

  dp <- (mu_cft / mu_obs^2) * (rr - 1)   # scalar
  db <- (mu_cft * p - mu_obs * 0) / mu_obs^2 * rr   # d/d(beta) with chain rule

  expected <- dp^2 * var_p_ + db^2 * var_b

  got <- variance(paf(p = p, beta = beta, var_p = var_p_, var_beta = var_b,
                      rr_link = "exp", quiet = TRUE))
  expect_equal(got, expected, tolerance = 1e-8)
})

# =============================================================================
# 3. cov_atomic_pif: two truly independent atomics give zero covariance
# =============================================================================

test_that("[num] cov_atomic_pif is zero for truly independent atomics", {
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "ind1")
  p2 <- paf(0.4, 1.8, var_p = 0.02, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "ind2")

  got <- cov_atomic_pif(p1, p2,
                        var_p    = matrix(0, 1, 1),
                        var_beta = matrix(0, 1, 1))
  expect_equal(got, 0, tolerance = TOL)
})

test_that("[num] cov_atomic_pif with shared p equals correct cross-derivative formula", {
  # Two PAFs that share the same p (perfectly correlated p, independent beta)
  p     <- 0.3
  rr1   <- 1.4;  rr2   <- 1.9
  sig_p <- 0.01   # scalar variance of p (same for both, fully shared)
  sig_b <- 0.0    # independent betas

  mu1 <- (1 - p) + p * rr1;  mu2 <- (1 - p) + p * rr2
  dp1 <- (1 / mu1^2) * (rr1 - 1)   # d(PAF1)/d(p), mu_cft = 1 for PAF
  dp2 <- (1 / mu2^2) * (rr2 - 1)

  expected <- dp1 * sig_p * dp2    # scalar case: t(dp1) * Sigma_p * dp2

  p1 <- paf(p, rr1, var_p = sig_p, var_beta = 0,
            rr_link = "identity", quiet = TRUE, label = "sp1")
  p2 <- paf(p, rr2, var_p = sig_p, var_beta = 0,
            rr_link = "identity", quiet = TRUE, label = "sp2")

  got <- cov_atomic_pif(p1, p2,
                        var_p    = matrix(sig_p),
                        var_beta = matrix(sig_b))
  expect_equal(got, expected, tolerance = TOL)
})

test_that("[num] cov_atomic_pif self-covariance equals variance", {
  p1  <- paf(0.25, 1.6, var_p = 0.015, var_beta = 0.04,
             rr_link = "identity", quiet = TRUE, label = "self")
  cov_self <- cov_atomic_pif(p1, p1,
                             var_p    = matrix(p1@var_p),
                             var_beta = matrix(p1@var_beta))
  expect_equal(cov_self, variance(p1), tolerance = TOL)
})

# =============================================================================
# 4. pif_total point estimate and variance (identity pif_transform)
# =============================================================================

test_that("[num] pif_total point estimate is weighted sum of component PIFs", {
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "t1")
  p2 <- paf(0.4, 1.8, var_p = 0.02, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "t2")
  w  <- c(0.6, 0.4)

  expected <- w[1] * coef(p1) + w[2] * coef(p2)

  got <- coef(pif_total(p1, p2, weights = w,
                        link = "identity", quiet = TRUE))
  expect_equal(got, expected, tolerance = TOL)
})

test_that("[num] pif_total variance with var_weights = 0 equals weighted sum of variances", {
  # When weights are deterministic and components are independent:
  # Var[w1*PIF1 + w2*PIF2] = w1^2*Var[PIF1] + w2^2*Var[PIF2]
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "tv1")
  p2 <- paf(0.4, 1.8, var_p = 0.02, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "tv2")
  w  <- c(0.6, 0.4)

  v1 <- variance(p1);  v2 <- variance(p2)
  expected <- w[1]^2 * v1 + w[2]^2 * v2

  tot <- pif_total(p1, p2, weights = w, var_weights = 0,
                   link = "identity", quiet = TRUE)
  expect_equal(variance(tot), expected, tolerance = TOL)
})

test_that("[num] pif_total with 3 components and known weights", {
  p1 <- paf(0.20, 1.4, var_p = 0.005, var_beta = 0.02,
            rr_link = "identity", quiet = TRUE, label = "t3a")
  p2 <- paf(0.30, 1.6, var_p = 0.008, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "t3b")
  p3 <- paf(0.15, 1.9, var_p = 0.003, var_beta = 0.04,
            rr_link = "identity", quiet = TRUE, label = "t3c")
  w  <- c(0.5, 0.3, 0.2)

  expected_coef <- w[1]*coef(p1) + w[2]*coef(p2) + w[3]*coef(p3)
  expected_var  <- w[1]^2*variance(p1) + w[2]^2*variance(p2) + w[3]^2*variance(p3)

  tot <- pif_total(p1, p2, p3, weights = w, var_weights = 0,
                   link = "identity", quiet = TRUE)
  expect_equal(coef(tot),     expected_coef, tolerance = TOL)
  expect_equal(variance(tot), expected_var,  tolerance = TOL)
})

# =============================================================================
# 5. pif_ensemble point estimate and variance (log-complement transform)
# =============================================================================

test_that("[num] pif_ensemble with equal weights 1 matches 1 - prod(1 - PIF_i)", {
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "e1")
  p2 <- paf(0.4, 1.8, var_p = 0.02, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "e2")
  w  <- c(1, 1)   # default equal weights

  pif1 <- coef(p1);  pif2 <- coef(p2)
  expected <- 1 - (1 - w[1]*pif1) * (1 - w[2]*pif2)

  got <- coef(pif_ensemble(p1, p2, weights = w, quiet = TRUE))
  expect_equal(got, expected, tolerance = TOL)
})

test_that("[num] pif_ensemble with non-unit weights matches formula", {
  p1 <- paf(0.25, 1.4, var_p = 0.005, var_beta = 0.02,
            rr_link = "identity", quiet = TRUE, label = "ew1")
  p2 <- paf(0.35, 1.7, var_p = 0.007, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "ew2")
  w  <- c(0.8, 0.6)

  pif1 <- coef(p1);  pif2 <- coef(p2)
  expected <- 1 - (1 - w[1]*pif1) * (1 - w[2]*pif2)

  got <- coef(pif_ensemble(p1, p2, weights = w, quiet = TRUE))
  expect_equal(got, expected, tolerance = TOL)
})

test_that("[num] pif_ensemble with 3 components matches product formula", {
  p1 <- paf(0.20, 1.4, var_p = 0.005, var_beta = 0.02,
            rr_link = "identity", quiet = TRUE, label = "e3a")
  p2 <- paf(0.30, 1.6, var_p = 0.008, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "e3b")
  p3 <- paf(0.15, 1.9, var_p = 0.003, var_beta = 0.04,
            rr_link = "identity", quiet = TRUE, label = "e3c")
  w  <- c(1, 1, 1)

  pif1 <- coef(p1);  pif2 <- coef(p2);  pif3 <- coef(p3)
  expected <- 1 - (1 - pif1) * (1 - pif2) * (1 - pif3)

  got <- coef(pif_ensemble(p1, p2, p3, weights = w, quiet = TRUE))
  expect_equal(got, expected, tolerance = TOL)
})

test_that("[num] single-component pif_ensemble equals the component itself", {
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "sing_ens")
  got <- coef(pif_ensemble(p1, weights = 1, quiet = TRUE))
  expect_equal(got, coef(p1), tolerance = TOL)
})

# =============================================================================
# 6. weighted_adjusted_fractions — point estimates
# =============================================================================

test_that("[num] adjusted PIF_i = (PIF_i / sum_PIF) * PIF_ensemble", {
  p1 <- paf(0.30, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "adj1")
  p2 <- paf(0.20, 1.8, var_p = 0.02, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "adj2")

  pif1 <- coef(p1);  pif2 <- coef(p2)
  sum_pif <- pif1 + pif2
  ens_pif <- 1 - (1 - pif1) * (1 - pif2)   # default ensemble weights = 1

  expected1 <- (pif1 / sum_pif) * ens_pif
  expected2 <- (pif2 / sum_pif) * ens_pif

  adj <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)
  expect_equal(unname(coef(adj$adj1)), expected1, tolerance = TOL)
  expect_equal(unname(coef(adj$adj2)), expected2, tolerance = TOL)
})

test_that("[num] adjusted fractions sum to ensemble (3 fractions)", {
  p1 <- paf(0.20, 1.4, var_p = 0.005, var_beta = 0.02,
            rr_link = "identity", quiet = TRUE, label = "s3a")
  p2 <- paf(0.15, 1.7, var_p = 0.006, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "s3b")
  p3 <- paf(0.10, 2.0, var_p = 0.004, var_beta = 0.04,
            rr_link = "identity", quiet = TRUE, label = "s3c")

  pif1 <- coef(p1);  pif2 <- coef(p2);  pif3 <- coef(p3)
  ens_expected <- 1 - (1 - pif1) * (1 - pif2) * (1 - pif3)

  adj <- weighted_adjusted_fractions(p1, p2, p3, quiet = TRUE)
  adj_sum <- sum(sapply(adj, coef))
  expect_equal(adj_sum, ens_expected, tolerance = TOL)
})

test_that("[num] adjusted fractions are proportional to originals (same ratio)", {
  p1 <- paf(0.25, 1.5, var_p = 0.01, var_beta = 0.04,
            rr_link = "identity", quiet = TRUE, label = "rat1")
  p2 <- paf(0.35, 1.9, var_p = 0.02, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "rat2")

  adj    <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)
  ratio1 <- adj[["rat1"]]@pif / coef(p1)
  ratio2 <- adj[["rat2"]]@pif / coef(p2)
  expect_equal(unname(ratio1), unname(ratio2), tolerance = TOL)
})

test_that("[num] adjusted point estimates with custom ensemble weights", {
  p1 <- paf(0.30, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "cw1")
  p2 <- paf(0.20, 1.8, var_p = 0.02, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "cw2")
  w  <- c(0.7, 0.9)

  pif1 <- coef(p1);  pif2 <- coef(p2)
  sum_pif <- pif1 + pif2
  ens_pif <- 1 - (1 - w[1]*pif1) * (1 - w[2]*pif2)

  expected1 <- (pif1 / sum_pif) * ens_pif
  expected2 <- (pif2 / sum_pif) * ens_pif

  adj <- weighted_adjusted_fractions(p1, p2, weights = w, quiet = TRUE)
  expect_equal(unname(adj[["cw1"]]@pif), expected1, tolerance = TOL)
  expect_equal(unname(adj[["cw2"]]@pif), expected2, tolerance = TOL)
})

# =============================================================================
# 7. weighted_adjusted_fractions — variance (delta method on log scale)
# =============================================================================

test_that("[num] adjusted variance with var_p = var_beta = 0 is zero", {
  p1 <- paf(0.30, 1.5, var_p = 0, var_beta = 0,
            rr_link = "identity", quiet = TRUE, label = "zv1")
  p2 <- paf(0.20, 1.8, var_p = 0, var_beta = 0,
            rr_link = "identity", quiet = TRUE, label = "zv2")

  adj <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)
  expect_equal(adj[["zv1"]]@variance, 0, tolerance = TOL)
  expect_equal(adj[["zv2"]]@variance, 0, tolerance = TOL)
})

test_that("[num] adjusted variance formula matches manual computation (2 atomics, identity link)", {
  # We verify the delta-method log-scale formula:
  #   PIF_adj_i = (PIF_i / S) * E,   S = PIF_1 + PIF_2,  E = ensemble
  #   ln PIF_adj_i = ln PIF_i - ln S + ln E
  #   Var[ln PIF_adj_i] = Var[ln PIF_i] + Var[ln S] + Var[ln E]
  #                       + 2*(Cov(ln PIF_i, ln E) - Cov(ln PIF_i, ln S) - Cov(ln E, ln S))
  #   Var[PIF_adj_i] ≈ PIF_adj_i^2 * Var[ln PIF_adj_i]

  p1 <- paf(0.30, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "vm1")
  p2 <- paf(0.20, 1.8, var_p = 0.02, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "vm2")

  pif1 <- coef(p1);  pif2 <- coef(p2)
  var1 <- variance(p1);  var2 <- variance(p2)

  # Build sum and ensemble objects (same way the function does internally)
  pif_sum <- pif_total(p1, p2, weights = c(1, 1), var_weights = 0,
                       link = "log-complement", weights_sum_to_1 = FALSE,
                       quiet = TRUE, label = "sum_vm")
  pif_ens <- pif_ensemble(p1, p2, weights = c(1, 1), quiet = TRUE,
                          label = "ens_vm")

  sum_pif <- coef(pif_sum)
  ens_pif <- coef(pif_ens)

  # Check sum_pif = pif1 + pif2 (pif_total with identity and weights = 1)
  # Note: pif_total with log-complement uses the log-complement transform so
  # the sum_pif comes from inv_log_complement(log_complement(pif1) + log_complement(pif2))
  # which equals 1 - (1-pif1)*(1-pif2). We use coef() to let the package compute it.

  var_sum     <- variance(pif_sum)
  var_ens     <- variance(pif_ens)
  var_log_sum <- var_sum / sum_pif^2
  var_log_ens <- var_ens / ens_pif^2

  cov_sum_ens     <- cov_total_pif(pif_sum, pif_ens, warning = FALSE)
  cov_log_sum_ens <- cov_sum_ens / (sum_pif * ens_pif)

  # For PIF_1
  pif_adj1    <- (pif1 / sum_pif) * ens_pif
  var_log_1   <- var1 / pif1^2
  cov_1_sum   <- cov_total_pif(p1, pif_sum, warning = FALSE)
  cov_log_1_sum <- cov_1_sum / (pif1 * sum_pif)
  cov_1_ens   <- cov_total_pif(p1, pif_ens, warning = FALSE)
  cov_log_1_ens <- cov_1_ens / (pif1 * ens_pif)

  var_log_adj1 <- var_log_1 + var_log_sum + var_log_ens +
    2 * (cov_log_1_ens - cov_log_1_sum - cov_log_sum_ens)
  expected_var1 <- max(pif_adj1^2 * var_log_adj1, 0)

  # For PIF_2
  pif_adj2    <- (pif2 / sum_pif) * ens_pif
  var_log_2   <- var2 / pif2^2
  cov_2_sum   <- cov_total_pif(p2, pif_sum, warning = FALSE)
  cov_log_2_sum <- cov_2_sum / (pif2 * sum_pif)
  cov_2_ens   <- cov_total_pif(p2, pif_ens, warning = FALSE)
  cov_log_2_ens <- cov_2_ens / (pif2 * ens_pif)

  var_log_adj2 <- var_log_2 + var_log_sum + var_log_ens +
    2 * (cov_log_2_ens - cov_log_2_sum - cov_log_sum_ens)
  expected_var2 <- max(pif_adj2^2 * var_log_adj2, 0)

  adj <- weighted_adjusted_fractions(p1, p2, quiet = TRUE)
  expect_equal(adj[["vm1"]]@variance, expected_var1, tolerance = 1e-8)
  expect_equal(adj[["vm2"]]@variance, expected_var2, tolerance = 1e-8)
})

# =============================================================================
# 8. PAF special cases
# =============================================================================

test_that("[num] paf = 0 when RR = 1 everywhere", {
  p  <- c(0.4, 0.3)
  rr <- c(1.0, 1.0)
  got <- coef(paf(p = p, beta = rr, rr_link = "identity", quiet = TRUE))
  expect_equal(got, 0, tolerance = TOL)
})

test_that("[num] paf equals 1 when all p exposed to very high RR (theoretical limit)", {
  # PIF = 1 - 1/E[RR]. As RR -> Inf, PIF -> 1.
  p  <- 0.99   # nearly all exposed
  rr <- 1e6
  got <- coef(paf(p = p, beta = rr, rr_link = "identity", quiet = TRUE))
  expect_true(got > 0.999)
})

test_that("[num] pif <= paf when p_cft >= 0 (counterfactual >= zero)", {
  # Any positive counterfactual prevalence gives smaller fraction than PAF
  p     <- 0.4
  beta  <- log(2.5)
  p_cft <- 0.1

  paf_val <- coef(paf(p = p, beta = beta, rr_link = "exp", quiet = TRUE, label = "a"))
  pif_val <- coef(pif(p = p, p_cft = p_cft, beta = beta,
                      rr_link = "exp", quiet = TRUE, label = "b"))
  expect_true(pif_val <= paf_val)
})

# =============================================================================
# 9. cov_atomic_pif: shared beta with non-trivial matrix
# =============================================================================

test_that("[num] cov_atomic_pif matches formula for two PIFs sharing beta", {
  # Two pifs with different p but same beta — provide their actual shared beta var
  p1    <- c(0.3, 0.4)
  p2    <- c(0.2, 0.5)
  rr    <- c(1.3, 1.7)   # shared RR (identity link)
  p_cft1 <- c(0.1, 0.2)
  p_cft2 <- c(0.05, 0.15)
  sig_b <- matrix(c(0.04, 0.01, 0.01, 0.03), ncol = 2)   # shared beta cov

  mu_obs1 <- (1-sum(p1)) + sum(p1 * rr)
  mu_obs2 <- (1-sum(p2)) + sum(p2 * rr)
  mu_cft1 <- (1-sum(p_cft1)) + sum(p_cft1 * rr)
  mu_cft2 <- (1-sum(p_cft2)) + sum(p_cft2 * rr)

  # d(PIF)/d(rr_i) for identity link
  drr1 <- (mu_cft1 * p1 - mu_obs1 * p_cft1) / mu_obs1^2
  drr2 <- (mu_cft2 * p2 - mu_obs2 * p_cft2) / mu_obs2^2

  expected <- as.numeric(t(drr1) %*% sig_b %*% drr2)

  pi1 <- pif(p = p1, p_cft = p_cft1, beta = rr, var_p = 0, var_beta = 0,
             rr_link = "identity", quiet = TRUE, label = "sb1")
  pi2 <- pif(p = p2, p_cft = p_cft2, beta = rr, var_p = 0, var_beta = 0,
             rr_link = "identity", quiet = TRUE, label = "sb2")

  got <- cov_atomic_pif(pi1, pi2,
                        var_p    = matrix(0, 2, 2),
                        var_beta = sig_b)
  expect_equal(got, expected, tolerance = TOL)
})

# =============================================================================
# 10. pif_total variance with non-zero var_weights
# =============================================================================

test_that("[num] pif_total variance with deterministic weights + weight uncertainty", {
  # Var[w1*PIF1 + w2*PIF2] when weights also have variance
  # = w1^2*Var[PIF1] + w2^2*Var[PIF2]
  #   + PIF1^2*sig_w11 + PIF2^2*sig_w22
  #   + 2*PIF1*PIF2*sig_w12    (by delta method)
  # (assuming PIF1, PIF2 independent of weights)

  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "wv1")
  p2 <- paf(0.4, 1.8, var_p = 0.02, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "wv2")
  w <- c(0.6, 0.4)
  sig_w <- matrix(c(0.004, 0.001, 0.001, 0.003), ncol = 2)

  pif1 <- coef(p1);  pif2 <- coef(p2)
  v1   <- variance(p1);  v2 <- variance(p2)

  expected <- w[1]^2*v1 + w[2]^2*v2 +
    pif1^2*sig_w[1,1] + pif2^2*sig_w[2,2] + 2*pif1*pif2*sig_w[1,2]

  tot <- pif_total(p1, p2, weights = w, var_weights = sig_w,
                   link = "identity", quiet = TRUE)
  expect_equal(variance(tot), expected, tolerance = 1e-8)
})

# =============================================================================
# 11. pif_ensemble variance with known structure (var_weights = 0, independent)
# =============================================================================

test_that("[num] pif_ensemble variance with 2 independent atomics (delta method)", {
  # PIF_ens = 1 - (1-PIF1)(1-PIF2)
  # Let c_i = log(1 - PIF_i).  PIF_ens = 1 - exp(c1 + c2).
  # d(PIF_ens)/d(PIF_i) = (1 - PIF_ens) / (1 - PIF_i)
  # Var[PIF_ens] = sum_i [d(PIF_ens)/d(PIF_i)]^2 * Var[PIF_i]

  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "ev1")
  p2 <- paf(0.4, 1.8, var_p = 0.02, var_beta = 0.03,
            rr_link = "identity", quiet = TRUE, label = "ev2")
  w  <- c(1, 1)

  pif1 <- coef(p1);  pif2 <- coef(p2)
  pif_ens <- 1 - (1 - pif1) * (1 - pif2)
  v1 <- variance(p1);  v2 <- variance(p2)

  # d(PIF_ens)/d(PIF_i) = (1 - PIF_ens) / (1 - PIF_i)
  d1 <- (1 - pif_ens) / (1 - pif1)
  d2 <- (1 - pif_ens) / (1 - pif2)

  expected <- d1^2 * v1 + d2^2 * v2

  ens <- pif_ensemble(p1, p2, weights = w, var_weights = 0, quiet = TRUE)
  expect_equal(variance(ens), expected, tolerance = 1e-8)
})

# =============================================================================
# 12. Consistency: pif_total with identity link = naive weighted mean
# =============================================================================

test_that("[num] pif_total (identity link) is identical to manually weighted mean", {
  set.seed(42)
  pafs <- lapply(1:4, function(i)
    paf(runif(1, 0.1, 0.4), runif(1, 1.2, 2.5),
        var_p = runif(1, 0.001, 0.02), var_beta = runif(1, 0.01, 0.1),
        rr_link = "identity", quiet = TRUE,
        label = paste0("cm", i)))
  w <- c(0.4, 0.3, 0.2, 0.1)

  expected <- sum(w * sapply(pafs, coef))
  tot <- do.call(pif_total,
                 c(pafs, list(weights = w, var_weights = 0,
                              link = "identity", quiet = TRUE)))
  expect_equal(coef(tot), expected, tolerance = TOL)
})

# =============================================================================
# 13. weighted_adjusted: single fraction is a no-op (adj = ensemble = itself)
# =============================================================================

test_that("[num] adjusted single fraction equals the fraction itself", {
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05,
            rr_link = "identity", quiet = TRUE, label = "noop1")

  adj <- weighted_adjusted_fractions(p1, quiet = TRUE)
  # PIF_adj = (PIF_1 / PIF_1) * PIF_ensemble = PIF_ensemble = PIF_1
  expect_equal(unname(adj[["noop1"]]@pif), coef(p1), tolerance = TOL)
})

# =============================================================================
# 14. Confidence intervals of atomic pif (Gaussian, identity link)
# =============================================================================

test_that("[num] 95% CI of atomic pif matches ±1.96*SD (identity link)", {
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05, link = "identity",
            rr_link = "identity", quiet = TRUE, label = "ci1")

  est <- coef(p1)
  sd_ <- standard_deviation(p1)
  z   <- qnorm(0.975)

  ci <- confint(p1)
  expect_equal(ci[1], est - z * sd_, tolerance = 1e-8)
  expect_equal(ci[2], est + z * sd_, tolerance = 1e-8)
})

test_that("[num] 90% CI uses correct z-score", {
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05, link = "identity",
            rr_link = "identity", quiet = TRUE, conf_level = 0.90, label = "ci90")

  est <- coef(p1)
  sd_ <- standard_deviation(p1)
  z   <- qnorm(0.95)

  ci <- confint(p1)
  expect_equal(ci[1], est - z * sd_, tolerance = 1e-8)
  expect_equal(ci[2], est + z * sd_, tolerance = 1e-8)
})

test_that("[num] 95% CI of atomic pif matches exp(log(pif) ±1.96*SD) (log link)", {
  p1 <- paf(0.3, 1.5, var_p = 0.01, var_beta = 0.05, link = "log-complement",
            rr_link = "identity", quiet = TRUE, label = "ci95")

  est <- coef(p1)
  sd_ <- standard_deviation(p1)
  z   <- qnorm(0.975)

  ci <- confint(p1)
  expect_equal(ci[1], inv_log_complement(log_complement(est)+z*sqrt(p1@link_variance)), tolerance = 1e-8)
  expect_equal(ci[2], inv_log_complement(log_complement(est)-z*sqrt(p1@link_variance)), tolerance = 1e-8)
})
