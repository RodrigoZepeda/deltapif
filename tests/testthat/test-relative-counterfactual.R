#Numerical gradient used to check the delta-method variances
num_grad <- function(f, x, h = 1e-6) {
  sapply(seq_along(x), function(i) {
    e <- replace(numeric(length(x)), i, h)
    (f(x + e) - f(x - e)) / (2 * h)
  })
}

test_that("A relative reduction of a binary exposure equals r times the PAF", {
  sd_beta <- (log(2.20) - log(1.15)) / (2 * qnorm(0.975))
  paf_obj <- paf(p = 0.085, beta = log(1.59), var_beta = sd_beta^2,
                 var_p = 0.01^2, quiet = TRUE)
  pif_obj <- pif(p = 0.085, p_cft = function(p) 0.85 * p, beta = log(1.59),
                 var_beta = sd_beta^2, var_p = 0.01^2, quiet = TRUE)

  expect_equal(coef(pif_obj), 0.15 * coef(paf_obj))
  expect_equal(variance(pif_obj), 0.15^2 * variance(paf_obj), tolerance = 1e-6)

  #A fixed target with the same value has a larger variance
  pif_fixed <- pif(p = 0.085, p_cft = 0.85 * 0.085, beta = log(1.59),
                   var_beta = sd_beta^2, var_p = 0.01^2, quiet = TRUE)
  expect_equal(coef(pif_fixed), coef(pif_obj))
  expect_gt(variance(pif_fixed), variance(pif_obj))
})

test_that("The variance of a categorical relative counterfactual matches the delta method", {
  p_obs <- c(0.2, 0.1, 0.05)
  beta  <- log(c(1.2, 1.5, 2.5))
  var_p <- matrix(c(1e-4, -2e-5, -1e-5, -2e-5, 5e-5, -5e-6, -1e-5, -5e-6, 2e-5), ncol = 3)
  var_b <- diag(c(0.01, 0.02, 0.03))
  g     <- function(p) c(p[1], p[2] + 0.5 * p[3], 0.5 * p[3])

  pif_obj <- pif(p = p_obs, p_cft = g, beta = beta, var_p = var_p,
                 var_beta = var_b, quiet = TRUE)

  grad_p <- num_grad(function(p) pif_fun(p, g(p), exp(beta)), p_obs)
  grad_b <- num_grad(function(b) pif_fun(p_obs, g(p_obs), exp(b)), beta)
  expected <- as.numeric(t(grad_p) %*% var_p %*% grad_p + t(grad_b) %*% var_b %*% grad_b)

  expect_equal(coef(pif_obj), pif_fun(p_obs, g(p_obs), exp(beta)))
  expect_equal(variance(pif_obj), expected, tolerance = 1e-6)
})

test_that("Analytic and numerical jacobians of the counterfactual agree", {
  g     <- function(p) c(0.5 * p[1], p[2] + 0.5 * p[1])
  g_jac <- function(p) matrix(c(0.5, 0.5, 0, 1), ncol = 2)

  args <- list(p = c(0.3, 0.2), beta = log(c(1.4, 2)), var_p = diag(c(1e-3, 1e-3)),
               var_beta = diag(c(0.01, 0.02)), quiet = TRUE)

  numeric_jac  <- do.call(pif, c(args, list(p_cft = g)))
  function_jac <- do.call(pif, c(args, list(p_cft = g, p_cft_deriv = g_jac)))
  matrix_jac   <- do.call(pif, c(args, list(p_cft = g, p_cft_deriv = g_jac(c(0.3, 0.2)))))

  expect_equal(variance(numeric_jac), variance(function_jac), tolerance = 1e-6)
  expect_equal(variance(matrix_jac), variance(function_jac))
  expect_equal(numeric_jac@p_cft_jacobian, g_jac(c(0.3, 0.2)), tolerance = 1e-6)
})

test_that("The covariance of relative counterfactuals sharing parameters matches the delta method", {
  p_obs <- c(0.2, 0.1)
  beta  <- log(c(1.3, 1.8))
  var_p <- diag(c(1e-4, 5e-5))
  var_b <- diag(c(0.01, 0.02))
  g1    <- function(p) 0.5 * p
  g2    <- function(p) c(p[1] + p[2], 0)

  pif1 <- pif(p = p_obs, p_cft = g1, beta = beta, var_p = var_p, var_beta = var_b,
              quiet = TRUE, label = "half")
  pif2 <- pif(p = p_obs, p_cft = g2, beta = beta, var_p = var_p, var_beta = var_b,
              quiet = TRUE, label = "shift")

  grad <- function(g) {
    list(p = num_grad(function(p) pif_fun(p, g(p), exp(beta)), p_obs),
         b = num_grad(function(b) pif_fun(p_obs, g(p_obs), exp(b)), beta))
  }
  d1 <- grad(g1)
  d2 <- grad(g2)
  expected <- as.numeric(t(d1$p) %*% var_p %*% d2$p + t(d1$b) %*% var_b %*% d2$b)

  expect_equal(covariance(pif1, pif2)[1, 2], expected, tolerance = 1e-6)
})

test_that("Fixed counterfactuals keep a zero jacobian and reject p_cft_deriv", {
  pif_obj <- pif(p = c(0.3, 0.2), p_cft = c(0.1, 0.1), beta = c(0.2, 0.4),
                 var_p = 0, var_beta = 0, quiet = TRUE)
  expect_equal(pif_obj@p_cft_jacobian, matrix(0, 2, 2))

  expect_error(
    pif(p = 0.3, p_cft = 0.1, beta = 0.2, p_cft_deriv = matrix(1), quiet = TRUE),
    "p_cft_deriv"
  )
})
