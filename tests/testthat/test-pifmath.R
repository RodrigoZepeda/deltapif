#Collection of checks to make sure the math part of the pif
#is correct

test_that("pif and paf compute to known values", {

  #Compute the population attributable fraction
  p      <- c(0.4, 0.2)
  rr     <- c(1.2, 1.5)
  expect_equal(
    coef(paf(p = p, beta = rr, quiet = TRUE)),
    1 - 1/(p[1]*rr[1] + p[2]*rr[2] + (1 - sum(p))) # 1 - 1 / E[RR] where E[RR] = p_0*1 + p_1+rr_1 + ... + p_n rr_n
  )

  #Compute the potential impact fraction
  p      <- c(0.4, 0.2)
  rr     <- c(1.2, 1.5)
  p_cft  <- c(0.5, 0.1)
  expect_equal(
    coef(pif(p = p, p_cft = p_cft, beta = rr, quiet = TRUE)),
    1 - (p_cft[1]*rr[1] + p_cft[2]*rr[2] + (1 - sum(p_cft)))/(p[1]*rr[1] + p[2]*rr[2] + (1 - sum(p))) # 1 - 1 / E[RR] where E[RR] = p_0*1 + p_1+rr_1 + ... + p_n rr_n
  )

  #Expect the pif to be equal to the paf when p_cft = 0
  expect_equal(
    paf(p = p, beta = rr, quiet = TRUE),
    pif(p = p, p_cft = c(0.0, 0.0), beta = rr, quiet = TRUE, type = "PAF")
  )

  #Expect that adding specification for rr = 1 does not change values
  p2      <- c(0.4, 0.4, 0.2)
  rr2     <- c(1.0, 1.2, 1.5)
  p_cft2  <- c(0.4, 0.5, 0.1)
  expect_equal(
    coef(pif(p = p, p_cft = p_cft, beta = rr, quiet = TRUE)),
    coef(pif(p = p2, p_cft = p_cft2, beta = rr2, quiet = TRUE))
  )

  #TODO: Add computations by hand of the variance


})


test_that("covariance of a thing with itself is variance", {

  p      <- c(0.4, 0.2)
  rr     <- c(1.2, 1.5)
  p_cft  <- c(0.5, 0.1)
  pif1   <- pif(p = p, p_cft = p_cft, beta = rr, quiet = TRUE, var_p = matrix(c(0.01, 0.0002, 0.0002, 0.01), ncol = 2))
  pif2   <- pif(p = p, p_cft = c(0.1, 0.2), beta = rr, quiet = TRUE, var_p = matrix(c(0.01, 0.0002, 0.0002, 0.01), ncol = 2))

  #Covariance of a thing with itself should be the variance
  expect_equal(
    matrix(rep(variance(pif1), 4), ncol = 2),
    covariance(pif1, pif1)
  )

  #Correlation of a thing with itself should be 1
  expect_equal(
    matrix(rep(1, 4), ncol = 2),
    correlation(pif1, pif1)
  )

  #Covariance matrix should have the variance of pif1 and pif2 in the entries
  #of the diagonal
  expect_equal(
    c(variance(pif1), variance(pif2)),
    diag(covariance(pif1, pif2))
  )

  #Correlation matrix should have entries 1
  expect_equal(
    rep(1, 2),
    diag(correlation(pif1, pif2))
  )

  expect_equal(
    from_parameters_pif_variance(0.2, 0.1, 1.2, 1.2, 0.2*(1.2-1), 0.1*(1.2-1), 0.01, 0.01),
    from_parameters_pif_covariance(0.2, 0.2, 0.1, 0.1, 1.2, 1.2, 1.2, 1.2,
                                   0.2*(1.2-1), 0.2*(1.2-1), 0.1*(1.2-1),0.1*(1.2-1), 0.01, 0.01)
  )



})
