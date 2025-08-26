#Collection of checks to make sure the math part of the pif
#is correct

test_that("pif and paf compute to known values", {

  #Compute the population attributable fraction
  p      <- c(0.4, 0.2)
  rr     <- c(1.2, 1.5)
  expect_equal(
    coef(paf(p = p, beta = rr, quiet = TRUE, rr_link = "identity")),
    1 - 1/(p[1]*rr[1] + p[2]*rr[2] + (1 - sum(p))) # 1 - 1 / E[RR] where E[RR] = p_0*1 + p_1+rr_1 + ... + p_n rr_n
  )

  #Compute the potential impact fraction
  p      <- c(0.4, 0.2)
  rr     <- c(1.2, 1.5)
  p_cft  <- c(0.5, 0.1)
  expect_equal(
    coef(pif(p = p, p_cft = p_cft, beta = rr, quiet = TRUE, rr_link = "identity")),
    1 - (p_cft[1]*rr[1] + p_cft[2]*rr[2] + (1 - sum(p_cft)))/(p[1]*rr[1] + p[2]*rr[2] + (1 - sum(p))) # 1 - 1 / E[RR] where E[RR] = p_0*1 + p_1+rr_1 + ... + p_n rr_n
  )

  #Expect the pif to be equal to the paf when p_cft = 0
  expect_equal(
    paf(p = p, beta = rr, quiet = TRUE, label = "0", rr_link = "identity"),
    pif(p = p, p_cft = c(0.0, 0.0), beta = rr, quiet = TRUE, type = "PAF", label = "0", rr_link = "identity")
  )

  #Expect that adding specification for rr = 1 does not change values
  p2      <- c(0.4, 0.4, 0.2)
  rr2     <- c(1.0, 1.2, 1.5)
  p_cft2  <- c(0.4, 0.5, 0.1)
  expect_equal(
    coef(pif(p = p, p_cft = p_cft, beta = rr, quiet = TRUE, rr_link = "identity")),
    coef(pif(p = p2, p_cft = p_cft2, beta = rr2, quiet = TRUE, rr_link = "identity"))
  )

  #Removing the reference category should not change anything
  p      <- c(0.4, 0.2)
  rr     <- c(1, 0.3)
  p_cft  <- c(0.0, 1.0)
  expect_equal(
    coef(pif(p = p, p_cft = p_cft, beta = rr, quiet = TRUE, rr_link = "identity")),
    1 - (0.3 / (0.8*1 + 0.2*0.3))
  )

  #Changing the reference category gives same pif
  p      <- c(0.4, 0.2)
  rr     <- c(1.4, 1.5)
  p_cft  <- c(0.3, 0.2)
  rr2    <- c(1.0/1.5, rr[1]/1.5)
  p_cft2 <- c(1 - sum(p_cft), p_cft[1])
  p2     <- c(1 - sum(p), p[1])
  expect_equal(
    coef(pif(p = p, p_cft = p_cft, beta = rr, quiet = TRUE, rr_link = "identity")),
    coef(pif(p = p2, p_cft = p_cft2, beta = rr2, quiet = TRUE, rr_link = "identity"))
  )

  #Variance computation
  p      <- c(0.4, 0.2)
  rr     <- c(1.2, 1.5)
  p_cft  <- c(0.5, 0.1)
  sigma_p <- matrix(c(0.01, 0.002, 0.002, 0.01), ncol = 2)
  sigma_b <- matrix(c(0.05, -0.001, -0.001, 0.007), ncol = 2)
  mu_obs  <- 1 + as.numeric(p %*% (rr - 1))
  mu_cft  <- 1 + as.numeric(p_cft %*% (rr - 1))

  varp_component <- 0
  for (i in 1:length(p)){
    for(j in 1:length(p)){
      varp_component <- varp_component + sigma_p[i,j]*(rr[i] - 1)*(rr[j] - 1)
    }
  }
  varp_component <- (mu_cft / mu_obs^2)^2*varp_component

  varb_component <- 0
  for (i in 1:length(rr)){
    for(j in 1:length(rr)){
      varb_component <- varb_component +
        sigma_b[i,j]*rr[i]*rr[j]*(mu_cft*p[i] - mu_obs*p_cft[i])*(mu_cft*p[j] - mu_obs*p_cft[j])
    }
  }
  varb_component <- varb_component / mu_obs^4

  variance_pif <- as.numeric(varp_component + varb_component)

  mypif <- pif(p = p, p_cft = p_cft, beta = log(rr), var_p = sigma_p, var_beta = sigma_b,
               rr_link = "exp")

  #Check the p component of the variance
  expect_equal(
    as.numeric(from_parameters_covariance_p_component(p, p, p_cft, p_cft, rr, rr,
                                           sigma_p, upper_bound = FALSE)),
    varp_component
  )

  #Check the beta component of the variance
  expect_equal(
    as.numeric(from_parameters_covariance_beta_component(p, p, p_cft, p_cft, rr, rr, rr, rr,
                                            sigma_b, upper_bound = FALSE)),
    varb_component
  )

  #Check the complete variance
  expect_equal(
    from_parameters_pif_variance(
      p, p_cft, rr, rr,  sigma_p, sigma_b
    ),
    variance_pif
  )

  #Check the complete variance in the pif object
  expect_equal(
    mypif@variance,
    variance_pif
  )

  #Check the complete variance in the pif object
  #FIXME: Error in var_p
  # expect_equal(
  #   variance(mypif),
  #   variance_pif
  # )

})

test_that("nothing changes when constructing a pif object", {

  p      <- c(0.4, 0.2)
  rr     <- c(1.2, 1.5)
  p_cft  <- c(0.5, 0.1)
  sigma_p <- matrix(c(0.01, 0.002, 0.002, 0.01), ncol = 2)
  sigma_b <- matrix(c(0.05, -0.001, -0.001, 0.007), ncol = 2)
  mu_obs  <- as.numeric(1 + p %*% (rr - 1))
  mu_cft  <- as.numeric(1 + p_cft %*% (rr - 1))

  mypif <- pif(p = p, p_cft = p_cft, beta = log(rr), var_p = sigma_p, var_beta = sigma_b,
               rr_link = "exp", link = "logit", conf_level = 0.9)

  expect_equal(mypif@p, p)
  expect_equal(mypif@p_cft, p_cft)
  expect_equal(mypif@rr, rr)
  expect_equal(mypif@rr_link_deriv_vals, rr)
  expect_equal(mypif@var_p, sigma_p)
  expect_equal(mypif@var_beta, sigma_b)
  expect_equal(mypif@mu_obs, mu_obs)
  expect_equal(mypif@mu_cft, mu_cft)
  expect_equal(mypif@rr_link, exp)
  expect_equal(mypif@rr_link_deriv, exp)
  expect_equal(mypif@link, logit)
  expect_equal(mypif@link_inv, inv_logit)
  expect_equal(mypif@link_deriv, deriv_logit)
  expect_equal(mypif@conf_level, 0.9)



})

test_that("covariance calculations", {

  pif1 <- pif(p = 0.52, p_cft = 0.1, beta = 1.4, var_p = 0, var_beta = 0.1, rr_link = "identity")
  pif2 <- pif(p = 0.55, p_cft = 0.12, beta = 1.3, var_p = 0, var_beta = 0.2, rr_link = "identity")

  #For the covariance by hand
  beta_1 <- deriv_pif_beta(0.52, 0.1, 1.4, 1)
  beta_2 <- deriv_pif_beta(0.55, 0.12, 1.3, 1)

  p_1 <- deriv_pif_p(0.52, 0.1, 1.4)
  p_2 <- deriv_pif_p(0.55, 0.12, 1.3)

  expect_equal(
    cov_atomic_pif(pif1, pif2, var_beta = 0.3, var_p = 0.2),
    beta_1*0.3*beta_2 + p_1*0.2*p_2
  )



})

#FIXME: varp error

# test_that("covariance of a thing with itself is variance", {
#
#   p      <- c(0.4, 0.2)
#   rr     <- c(1.2, 1.5)
#   p_cft  <- c(0.5, 0.1)
#   pif1   <- pif(p = p, p_cft = p_cft, beta = rr, quiet = TRUE, var_p = matrix(c(0.01, 0.0002, 0.0002, 0.01), ncol = 2))
#   pif2   <- pif(p = p, p_cft = c(0.1, 0.2), beta = rr, quiet = TRUE, var_p = matrix(c(0.01, 0.0002, 0.0002, 0.01), ncol = 2))
#
#   # #Covariance of a thing with itself should be the variance
#   # expect_equal(
#   #   matrix(rep(variance(pif1), 4), ncol = 2),
#   #   covariance(pif1, pif1)
#   # )
#   #
#   # #Correlation of a thing with itself should be 1
#   # expect_equal(
#   #   matrix(rep(1, 4), ncol = 2),
#   #   correlation(pif1, pif1, uncorrelated_beta = FALSE, uncorrelated_p = FALSE)
#   # )
#   #
#   # #Covariance matrix should have the variance of pif1 and pif2 in the entries
#   # #of the diagonal
#   # expect_equal(
#   #   c(variance(pif1), variance(pif2)),
#   #   diag(covariance(pif1, pif2, uncorrelated_beta = FALSE, uncorrelated_p = FALSE))
#   # )
#   #
#   # #Correlation matrix should have entries 1
#   # expect_equal(
#   #   rep(1, 2),
#   #   diag(correlation(pif1, pif2, uncorrelated_beta = FALSE, uncorrelated_p = FALSE))
#   # )
#   #
#   # expect_equal(
#   #   from_parameters_pif_variance(0.2, 0.1, 1.2, 1.2, 0.2*(1.2-1), 0.1*(1.2-1), 0.01, 0.01),
#   #   from_parameters_pif_covariance(0.2, 0.2, 0.1, 0.1, 1.2, 1.2, 1.2, 1.2,
#   #                                  0.2*(1.2-1), 0.2*(1.2-1), 0.1*(1.2-1),0.1*(1.2-1), 0.01, 0.01)
#   # )
#   #
#   # expect_true(
#   #   #Covariance of a thing with itself should be the variance
#   #   isSymmetric(
#   #     covariance(pif1, pif2, uncorrelated_beta = FALSE, uncorrelated_p = FALSE)
#   #   )
#   # )
#
#
# })
#
