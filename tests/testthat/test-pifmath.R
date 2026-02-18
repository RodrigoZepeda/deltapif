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
               rr_link = "exp", quiet = TRUE)

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
  expect_equal(
    variance(mypif),
    variance_pif
  )

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
               rr_link = "exp", link = "logit", conf_level = 0.9, quiet = TRUE)

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

  pif1 <- pif(p = 0.52, p_cft = 0.1, beta = 1.4, var_p = 0, var_beta = 0.1, rr_link = "identity", quiet = TRUE)
  pif2 <- pif(p = 0.55, p_cft = 0.12, beta = 1.3, var_p = 0, var_beta = 0.2, rr_link = "identity", quiet = TRUE)

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
    from_parameters_pif_variance(p = 0.2, p_cft = 0.1, rr = 1.2,
                                 rr_link_deriv_vals = 1.2,
                                 var_p = 0.2*(1.2-1), var_beta = 0.1*(1.2-1)),
    from_parameters_pif_covariance(p1 = 0.2, p2 = 0.2, p1_cft = 0.1,
                                   p2_cft = 0.1, rr1 = 1.2, rr2 = 1.2,
                                   rr_link_deriv_vals1 = 1.2,
                                   rr_link_deriv_vals2 = 1.2,
                                   var_p = 0.2*(1.2-1),
                                   var_beta = 0.1*(1.2-1))
  )

  expect_true(
    #Covariance of a thing with itself should be the variance
    isSymmetric(
      covariance(pif1, pif2)
    )
  )


})

test_that("computation of ensemble variance", {

  pif1 <- paf(0.2, 0.1, var_p = 0.01, quiet = TRUE)
  pif2 <- paf(0.2, 0.3, var_p = 0.01, quiet = TRUE)

  var_weights = matrix(c(0.1, 0.2, 0.2, 0.3), ncol = 2, nrow = 2)
  weights = c(0.4, 0.6)
  result  <- pif_ensemble(pif1, pif2, weights = weights, var_weights = var_weights, quiet = TRUE)
  cov_mat <- covariance(pif1, pif2)
  val <- 0
  for (i in 1:2){
    for (j in 1:2){
      val <- val + (
        weights[i] %*% cov_mat[i,j] %*% weights[j] +
          result@coefs[i]* var_weights[i,j] * result@coefs[j]
      ) / ((1 - weights[i]*result@coefs[i])*(1 - weights[j]*result@coefs[j]))
    }
  }
  val <- val * (1 - result@pif)^2
  expect_equal(result@variance, as.numeric(val))


})

test_that("computation of the variance", {

  pvar_men <- matrix(
    c(
      0.005, -0.002, 0, -0, 0, 0, -0, -0, 0, -0, -0, -0.002, 0.002, -0, 0, -0, -0,
      0, -0, 0, -0, -0, 0, -0, 0, 0, -0, -0, -0, -0, 0, -0, 0, -0, 0, 0, 0, -0, 0,
      -0, -0, 0, -0, -0, 0, -0, -0, -0, 0, -0, 0, 0, -0, -0, 0, 0, -0, -0, 0, -0, 0,
      0, -0, 0, -0, -0, -0, 0, -0, -0, 0, 0, 0.001, 0, -0, -0, -0, -0, -0, -0, -0,
      0, -0, 0, 0.001, -0, 0, 0, 0, 0, 0, 0, -0, 0, -0, -0, 0, 0, -0, -0, -0, -0,
      -0, -0, -0, -0, 0, 0, 0.001, -0, -0, -0, 0, -0, 0, -0, -0, 0, -0, -0, 0.001
    ),
    nrow = 11L,
    ncol = 11L,
    dimnames = list(
      c(
        ">0-<1 g", "1-<2g", "2-<5g", "5-<10g", "10-<15g", "15-<20g", "20-<30g",
        "30-<40g", "40-<50g", "50-<60g", "\U{2265}60g"
      ),
      c(
        ">0-<1 g", "1-<2g", "2-<5g", "5-<10g", "10-<15g", "15-<20g", "20-<30g",
        "30-<40g", "40-<50g", "50-<60g", "\U{2265}60g"
      )
    )
  )

  pvar_women <- matrix(
    c(
      0, -0, -0, 0, 0, -0, 0, -0, -0, -0, 0, -0, 0.002, 0, 0, 0, 0, 0, -0, -0, 0,
      -0, -0, 0, 0, -0, -0, 0, 0, 0, -0, -0, -0, 0, 0, -0, 0, -0, 0, 0, -0, -0, 0,
      -0, 0, 0, -0, -0, 0, 0, -0, 0, -0, 0, -0, -0, 0, 0, 0, 0, 0, -0, -0, 0, -0,
      -0, 0, 0, 0, 0, -0, -0, 0.001, 0, -0, 0, 0, -0, -0, 0, -0, 0, -0, 0, 0, 0, -0,
      -0, -0, -0, -0, -0, -0, 0, -0, 0, 0, -0, 0, -0, 0, -0, 0, 0, -0, 0, -0, -0, 0,
      -0, 0, -0, -0, -0, -0, -0, 0, -0, 0, -0, 0
    ),
    nrow = 11L,
    ncol = 11L,
    dimnames = list(
      c(
        ">0-<1 g", "1-<2g", "2-<5g", "5-<10g", "10-<15g", "15-<20g", "20-<30g",
        "30-<40g", "40-<50g", "50-<60g", "\U{2265}60g"
      ),
      c(
        ">0-<1 g", "1-<2g", "2-<5g", "5-<10g", "10-<15g", "15-<20g", "20-<30g",
        "30-<40g", "40-<50g", "50-<60g", "\U{2265}60g"
      )
    )
  )

  beta_var_men <- matrix(
    c(
      8.687334105923459e-07, 1.7374668211846918e-06, 4.235075376637686e-06,
      7.710009019007069e-06, 1.3682551216829448e-05, 1.8894951680383518e-05,
      2.6604960699390592e-05, 3.7681311684443006e-05, 4.7997520935227106e-05,
      5.9073871920279516e-05, 8.481009920907776e-05, 1.7374668211846918e-06,
      3.4749336423693837e-06, 8.470150753275372e-06, 1.5420018038014138e-05,
      2.7365102433658896e-05, 3.7789903360767037e-05, 5.3209921398781185e-05,
      7.536262336888601e-05, 9.599504187045421e-05, 0.00011814774384055903,
      0.00016962019841815551, 4.235075376637686e-06, 8.470150753275372e-06,
      2.064599246110872e-05, 3.758629396765947e-05, 6.670243718204356e-05,
      9.211288944186967e-05, 0.00012969918340952915, 0.00018369639446165966,
      0.00023398791455923217, 0.00028798512561136266, 0.0004134492336442541,
      7.71000901900707e-06, 1.542001803801414e-05, 3.758629396765947e-05,
      6.842633004368774e-05, 0.00012143264204936135, 0.00016769269616340375,
      0.00023611902620709153, 0.00033442164119943166, 0.0004259779983001406,
      0.0005242806132924807, 0.0007526896304805652, 1.3682551216829448e-05,
      2.7365102433658896e-05, 6.670243718204356e-05, 0.00012143264204936135,
      0.0002155001816650638, 0.00029759548896604045, 0.00041902813101540187,
      0.0005934806590299774, 0.000755960954729827, 0.0009304134827444024,
      0.0013357590625429749, 1.889495168038352e-05, 3.778990336076704e-05,
      9.211288944186967e-05, 0.00016769269616340375, 0.00029759548896604045,
      0.0004109651990483415, 0.0005786578952117453, 0.0008195685291366352,
      0.0010439460803411895, 0.0012848567142660792, 0.0018446196577974411,
      2.6604960699390596e-05, 5.320992139878119e-05, 0.00012969918340952915,
      0.0002361190262070915, 0.00041902813101540187, 0.0005786578952117453,
      0.000814776921418837, 0.001153990170336067, 0.0014699240786413303,
      0.0018091373275585603, 0.0025973092882780064, 3.7681311684443006e-05,
      7.536262336888601e-05, 0.00018369639446165966, 0.00033442164119943166,
      0.0005934806590299774, 0.0008195685291366353, 0.0011539901703360672,
      0.0016344268943127155, 0.002081892470565476, 0.002562329194542124,
      0.003678638053193748, 4.799752093522711e-05, 9.599504187045423e-05,
      0.00023398791455923217, 0.0004259779983001406, 0.0007559609547298271,
      0.0010439460803411897, 0.0014699240786413305, 0.002081892470565476,
      0.002651863031671298, 0.0032638314235954435, 0.0046857579813015466,
      5.9073871920279516e-05, 0.00011814774384055903, 0.00028798512561136266,
      0.0005242806132924807, 0.0009304134827444024, 0.0012848567142660792,
      0.0018091373275585603, 0.002562329194542124, 0.0032638314235954435,
      0.004017023290579007, 0.005767086746217288, 8.481009920907777e-05,
      0.00016962019841815554, 0.0004134492336442541, 0.0007526896304805652,
      0.0013357590625429749, 0.0018446196577974411, 0.002597309288278007,
      0.003678638053193748, 0.0046857579813015466, 0.005767086746217288,
      0.008279585935286217
    ),
    nrow = 11L,
    ncol = 11L
  )

  beta_var_women <- matrix(
    c(
      8.687334105923459e-07, 1.7374668211846918e-06, 4.235075376637686e-06,
      7.710009019007069e-06, 1.3682551216829448e-05, 1.8894951680383518e-05,
      2.6604960699390592e-05, 3.7681311684443006e-05, 4.7997520935227106e-05,
      5.9073871920279516e-05, 8.481009920907776e-05, 1.7374668211846918e-06,
      3.4749336423693837e-06, 8.470150753275372e-06, 1.5420018038014138e-05,
      2.7365102433658896e-05, 3.7789903360767037e-05, 5.3209921398781185e-05,
      7.536262336888601e-05, 9.599504187045421e-05, 0.00011814774384055903,
      0.00016962019841815551, 4.235075376637686e-06, 8.470150753275372e-06,
      2.064599246110872e-05, 3.758629396765947e-05, 6.670243718204356e-05,
      9.211288944186967e-05, 0.00012969918340952915, 0.00018369639446165966,
      0.00023398791455923217, 0.00028798512561136266, 0.0004134492336442541,
      7.71000901900707e-06, 1.542001803801414e-05, 3.758629396765947e-05,
      6.842633004368774e-05, 0.00012143264204936135, 0.00016769269616340375,
      0.00023611902620709153, 0.00033442164119943166, 0.0004259779983001406,
      0.0005242806132924807, 0.0007526896304805652, 1.3682551216829448e-05,
      2.7365102433658896e-05, 6.670243718204356e-05, 0.00012143264204936135,
      0.0002155001816650638, 0.00029759548896604045, 0.00041902813101540187,
      0.0005934806590299774, 0.000755960954729827, 0.0009304134827444024,
      0.0013357590625429749, 1.889495168038352e-05, 3.778990336076704e-05,
      9.211288944186967e-05, 0.00016769269616340375, 0.00029759548896604045,
      0.0004109651990483415, 0.0005786578952117453, 0.0008195685291366352,
      0.0010439460803411895, 0.0012848567142660792, 0.0018446196577974411,
      2.6604960699390596e-05, 5.320992139878119e-05, 0.00012969918340952915,
      0.0002361190262070915, 0.00041902813101540187, 0.0005786578952117453,
      0.000814776921418837, 0.001153990170336067, 0.0014699240786413303,
      0.0018091373275585603, 0.0025973092882780064, 3.7681311684443006e-05,
      7.536262336888601e-05, 0.00018369639446165966, 0.00033442164119943166,
      0.0005934806590299774, 0.0008195685291366353, 0.0011539901703360672,
      0.0016344268943127155, 0.002081892470565476, 0.002562329194542124,
      0.003678638053193748, 4.799752093522711e-05, 9.599504187045423e-05,
      0.00023398791455923217, 0.0004259779983001406, 0.0007559609547298271,
      0.0010439460803411897, 0.0014699240786413305, 0.002081892470565476,
      0.002651863031671298, 0.0032638314235954435, 0.0046857579813015466,
      5.9073871920279516e-05, 0.00011814774384055903, 0.00028798512561136266,
      0.0005242806132924807, 0.0009304134827444024, 0.0012848567142660792,
      0.0018091373275585603, 0.002562329194542124, 0.0032638314235954435,
      0.004017023290579007, 0.005767086746217288, 8.481009920907777e-05,
      0.00016962019841815554, 0.0004134492336442541, 0.0007526896304805652,
      0.0013357590625429749, 0.0018446196577974411, 0.002597309288278007,
      0.003678638053193748, 0.0046857579813015466, 0.005767086746217288,
      0.008279585935286217
    ),
    nrow = 11L,
    ncol = 11L
  )
  p_men <- c(
    0.013999999999999999, 0.022000000000000002, 0.081, 0.081, 0.069, 0.046, 0.073,
    0.042, 0.032, 0.02, 0.048
  )
  p_women <- c(
    0.02, 0.054000000000000006, 0.091, 0.081, 0.061, 0.026000000000000002, 0.034,
    0.015, 0.006999999999999999, 0.003, 0.005
  )
  beta_men <- c(
    0.007624814384345996, 0.015249628768691991, 0.030499257537383983,
    0.0676702276610707, 0.12009082655344941, 0.16583971285952537,
    0.2335099405205961, 0.33072632392100754, 0.4212709947351162,
    0.5184873781355276, 0.8120427319328485
  )
  beta_women <- c(
    0.007624814384345996, 0.015249628768691991, 0.037170970123686725,
    0.0676702276610707, 0.12009082655344941, 0.16583971285952537,
    0.2335099405205961, 0.33072632392100754, 0.4212709947351162,
    0.5184873781355276, 0.7443725042717777
  )

  paf_men <- paf(p = p_men, beta = beta_men,
                 var_p = pvar_men,
                 var_beta = beta_var_men,
                 label = "Men", quiet = TRUE)

  paf_women <- paf(p = p_women, beta = beta_women,
                   var_p = pvar_women,
                   var_beta = beta_var_women,
                   label = "Women", quiet = TRUE)

  total <- paf_total(paf_men, paf_women, weights = c(0.5, 0.5), link = "identity",
                     var_weights = matrix(c(0.01, 0.0, 0.0, 0.021), ncol = 2))
  expect_equal(
    total@variance,
    as.numeric(
      (
        c(0.5, 0.5) %*% covariance(paf_men, paf_women) %*% c(0.5, 0.5) +
          c(coef(paf_men), coef(paf_women)) %*%  matrix(c(0.01, 0.0, 0.0, 0.021), ncol = 2) %*% c(coef(paf_men), coef(paf_women))
      )
    )
  )

  #Variance identical regardless of transformation
  total2 <- paf_total(paf_men, paf_women, weights = c(0.5, 0.5), link = "log-complement",
                     var_weights = matrix(c(0.01, 0.0, 0.0, 0.021), ncol = 2))
  expect_equal(total@variance, total2@variance)



})
