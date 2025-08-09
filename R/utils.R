#' Additional helper functions
#'
#' @inheritParams pif
#' @param rr Values for the relative risk at `beta` for each of the exposure categories.
#' @param rr_link_deriv_vals Values for the relative risk's derivative at `beta` for each of the exposure categories.
#' @param link_vals The link function evaluated at the potential impact fraction
#' @param link_deriv_vals The link function's derivative evaluated at the potential impact fraction
#' @param mu_obs The average relative risk for the observed prevalence
#' @param mu_cft The average relative risk for the counterfactual prevalence
#' @param upper_bound Boolean indicating whether we are calculating an
#' upper bound for the variance (when covariances are unknown) or
#' if we are calculating the variance (when covariances are known).
#' @name helpers
#' @keywords internal

#' Compute the mean under the observed prevalence
#' @rdname helpers
mu_obs_fun <- function(p, rr) {
  as.numeric(1 + p %*% (rr - 1))
}

#' Compute the mean under the observed prevalence
#' @rdname helpers
mu_cft_fun <- function(p_cft, rr) {
  mu_obs_fun(p_cft, rr)
}

#' Compute the potential impact fraction
#' @rdname helpers
pif_fun <- function(p, p_cft, rr) {
  # Calculate the mean p and rr
  mu_obs <- mu_obs_fun(p, rr)
  mu_cft <- mu_cft_fun(p_cft, rr)

  # Calculate the potential impact fraction
  pif_fun2(mu_obs, mu_cft)
}

#' Compute the potential impact fraction
#' @rdname helpers
pif_fun2 <- function(mu_obs, mu_cft) {
  # Calculate the potential impact fraction
  1 - (mu_cft / mu_obs)
}

#' Partial derivative of PIF respect to p
#' @rdname helpers
deriv_pif_p <- function(p, p_cft, rr, mu_obs = NULL, mu_cft = NULL) {
  if (is.null(mu_obs)) {
    mu_obs <- mu_obs_fun(p, rr)
  }

  if (is.null(mu_cft)) {
    mu_cft <- mu_cft_fun(p_cft, rr)
  }

  # The derivative of pif
  (mu_cft / (mu_obs)^2) * (rr - 1)
}

#' Partial derivative of PIF respect to beta
#' @rdname helpers
deriv_pif_beta <- function(p, p_cft, rr, rr_link_deriv_vals, mu_obs = NULL,
                           mu_cft = NULL) {
  if (is.null(mu_obs)) {
    mu_obs <- mu_obs_fun(p, rr)
  }

  if (is.null(mu_cft)) {
    mu_cft <- mu_cft_fun(p_cft, rr)
  }

  # The derivative of pif
  ((mu_obs * p_cft - mu_cft * p) / (mu_obs^2)) * rr_link_deriv_vals
}

#' Covariance component with respect to p
#' @rdname helpers
covariance_p_component <- function(p1, p2, p1_cft, p2_cft, rr1, rr2, mu_obs1,
                                   mu_obs2, mu_cft1, mu_cft2, sigma_p,
                                   upper_bound = FALSE) {
  # Get the derivatives
  vcp_deriv1 <- deriv_pif_p(
    p = p1, p_cft = p1_cft, rr = rr1, mu_obs = mu_obs1,
    mu_cft = mu_cft1
  )
  vcp_deriv2 <- deriv_pif_p(
    p = p2, p_cft = p2_cft, rr = rr2, mu_obs = mu_obs2,
    mu_cft = mu_cft2
  )

  # Whether we are working with an upper bound for the variance
  if (upper_bound) {
    vcp_deriv1 <- abs(vcp_deriv1)
    vcp_deriv2 <- abs(vcp_deriv2)
  }

  t(vcp_deriv1) %*% sigma_p %*% vcp_deriv2
}

#' Covariance component with respect to beta
#' @rdname helpers
covariance_beta_component <- function(p1, p2, p1_cft, p2_cft, rr1, rr2,
                                      rr_link_deriv_vals1, rr_link_deriv_vals2, mu_obs1,
                                      mu_obs2, mu_cft1, mu_cft2, sigma_beta,
                                      upper_bound = FALSE) {
  # Get the derivatives
  vcbeta_deriv1 <- deriv_pif_beta(
    p = p1, p_cft = p1_cft, rr = rr1,
    rr_link_deriv_vals = rr_link_deriv_vals1, mu_obs = mu_obs1,
    mu_cft = mu_cft1
  )

  vcbeta_deriv2 <- deriv_pif_beta(
    p = p2, p_cft = p2_cft, rr = rr2,
    rr_link_deriv_vals = rr_link_deriv_vals2, mu_obs = mu_obs2,
    mu_cft = mu_cft2
  )

  # Whether we are working with an upper bound for the variance
  if (upper_bound) {
    vcbeta_deriv1 <- abs(vcbeta_deriv1)
    vcbeta_deriv2 <- abs(vcbeta_deriv2)
  }

  t(vcbeta_deriv1) %*% sigma_beta %*% vcbeta_deriv2
}

#' Covariance between two impact fractions
#' @rdname helpers
pif_covariance <- function(p1, p2, p1_cft, p2_cft, rr1, rr2,
                           rr_link_deriv_vals1, rr_link_deriv_vals2, mu_obs1,
                           mu_obs2, mu_cft1, mu_cft2, sigma_p, sigma_beta,
                           link_deriv_1, link_deriv_2,
                           upper_bound_p = FALSE,
                           upper_bound_beta = FALSE) {
  p_component <- covariance_p_component(
    p1 = p1,
    p2 = p2,
    p1_cft = p1_cft,
    p2_cft = p2_cft,
    rr1 = rr1, rr2 = rr2,
    mu_obs1 = mu_obs1,
    mu_obs2 = mu_obs2,
    mu_cft1 = mu_cft1,
    mu_cft2 = mu_cft2,
    sigma_p = sigma_p,
    upper_bound = upper_bound_p
  )

  beta_component <- covariance_beta_component(
    p1 = p1,
    p2 = p2,
    p1_cft = p1_cft,
    p2_cft = p2_cft,
    rr1 = rr1, rr2 = rr2,
    rr_link_deriv_vals1 = rr_link_deriv_vals1,
    rr_link_deriv_vals2 = rr_link_deriv_vals2,
    mu_obs1 = mu_obs1,
    mu_obs2 = mu_obs2,
    mu_cft1 = mu_cft1,
    mu_cft2 = mu_cft2,
    sigma_beta = sigma_beta,
    upper_bound = upper_bound_beta
  )

  as.numeric(link_deriv_1 * link_deriv_2 * (p_component + beta_component))
}

#' Variance for a potential impact fraction
#' @rdname helpers
pif_variance <- function(p, p_cft, rr, rr_link_deriv_vals, mu_obs, mu_cft, sigma_p,
                         sigma_beta, link_deriv_vals,
                         upper_bound_p = FALSE,
                         upper_bound_beta = FALSE) {
  pif_covariance(
    p1 = p,
    p2 = p,
    p1_cft = p_cft,
    p2_cft = p_cft,
    rr1 = rr,
    rr2 = rr,
    rr_link_deriv_vals1 = rr_link_deriv_vals,
    rr_link_deriv_vals2 = rr_link_deriv_vals,
    mu_obs1 = mu_obs,
    mu_obs2 = mu_obs,
    mu_cft1 = mu_cft,
    mu_cft2 = mu_cft,
    sigma_p = sigma_p,
    sigma_beta = sigma_beta,
    link_deriv_1 = link_deriv_vals,
    link_deriv_2 = link_deriv_vals,
    upper_bound_p = upper_bound_p,
    upper_bound_beta = upper_bound_beta
  )
}

#' Confidence interval for a potential impact fraction
#' @rdname helpers
pif_ci <- function(link_vals, variance, alpha_confint, link_inv) {
  b1 <- link_inv(link_vals - stats::qnorm(alpha_confint) * sqrt(variance))
  b2 <- link_inv(link_vals + stats::qnorm(alpha_confint) * sqrt(variance))
  sort(c(b1, b2))
}

#' Function to parse the link from a word to a function
#' @param link_name The name of the link or a function
#' @rdname helpers
parse_link <- function(link_name) {
  if (is.function(link_name)) {
    return(link_name)
  }

  if (link_name == "identity") {
    return(identity)
  } else if (link_name == "logit") {
    return(logit)
  } else if (link_name == "log-complement") {
    return(log_complement)
  } else if (link_name == "hawkins") {
    return(hawkins)
  } else if (link_name == "exponential") {
    return(exp)
  } else {
    cli::cli_abort("Cannot find link {.val {link_name}}. Please specify the function using
                   {.code rr_link}")
  }
}

#' Function to parse the inverse link from a word to a function
#' @param link_name The name of the link or a function for the inverse of the link
#' @rdname helpers
parse_inv_link <- function(link_name) {
  if (is.function(link_name)) {
    return(link_name)
  }

  if (link_name == "identity") {
    return(identity)
  } else if (link_name == "logit") {
    return(inv_logit)
  } else if (link_name == "log-complement") {
    return(inv_log_complement)
  } else if (link_name == "hawkins") {
    return(inv_hawkins)
  } else if (link_name == "exponential") {
    return(log)
  } else {
    cli::cli_abort("Cannot find link {.val {link_name}}. Please specify the function using
                   {.code rr_link}")
  }
}

#' Get the covariances between two pif classes (atomic pifs)
#'
#' @param pif1 A `pif_class` object
#' @param pif2 A `pif_class` object
#' @param sigma_p Covariance matrix for the prevalences in both `pif1` and `pif2`.
#' @param sigma_beta Covariance matrix for the parameter `beta` in both `pif1`
#' and `pif2`.
#' @param independent_p If `pif1` and `pif2` share the same prevalence data. Either
#' `TRUE`, `FALSE` or `guess` (default).
#' @param independent_beta If `pif1` and `pif2` share the same `beta` parameter. Either
#' `TRUE`, `FALSE` or `guess` (default).
#'
#' @section Covariance matrices:
#' By default if no `sigma_p` is specified this assumes the covariances
#' between the parameters `p` of `pif1` and `pif2` are uncorrelated.
#' However, if `pif1` and `pif2` share the same prevalence estimates
#' (i.e. share the same `p`s, then the user should set `independent_p = TRUE`
#' to account for that correlation).
#'
#' The same thing happens with `sigma_beta`. If no `sigma_beta` is specified
#' then the assumption is that the `beta` parameters from `pif1` and
#' from `pif2` are uncorrelated unless `independent_beta` is set to `TRUE`.
#'
#' If the user provides a covariance matrix for `sigma_p` then `independent_p`
#' is disregarded. Similarly, if the user provides `sigma_beta` then
#' `independent_beta` is ignored.
#'
#' @keywords internal
pif_class_covariance <- function(pif1, pif2, sigma_p = NULL, sigma_beta = NULL,
                                 independent_p = "guess", independent_beta = "guess",
                                 quiet = FALSE) {
  # Check they are pif objects
  if (!S7::S7_inherits(pif1, pif_class) || !S7::S7_inherits(pif2, pif_class)) {
    cli::cli_abort(
      paste0(
        "One of the arguments passed to covariance estimation is not a ",
        "`PIFCI::pif_class` object."
      )
    )
  }

  # Check that both p's have the same length
  if (length(pif1@p) != length(pif2@p)) {
    cli::cli_abort(
      paste0(
        "pif1 has {.code length(pif1@p) = {length(pif1@p)}} exposure",
        "categories and pif2 has {.code length(pif2@p) = {length(pif2@p)}}.",
        "They don't seem to represent the same population. They need to ",
        "be of the same length. Are you missing some categories with 0% ",
        "prevalence perhaps?"
      )
    )
  }

  # Check that both p's have the same length
  if (length(pif1@p_cft) != length(pif2@p_cft)) {
    cli::cli_abort(
      paste0(
        "pif1 has {.code length(pif1@p_cft) = {length(pif1@p_cft)}} counterfactual exposure",
        "categories and pif2 has {.code length(pif2@p_cft) = {length(pif2@p_cft)}}.",
        "They don't seem to represent the same population. They need to ",
        "be of the same length. Are you missing some categories with 0% ",
        "prevalence perhaps?"
      )
    )
  }

  # Check that the ps don't appear similar
  if (independent_p == "guess" && is.null(sigma_p) && all(pif1@p == pif2@p) && all(pif1@sigma_p == pif2@sigma_p)) {
    independent_p <- FALSE
    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "The prevalence parameters p for the potential impact fractions appear ",
          "to be the same. If they are, set `independent_p = FALSE`. Otherwise ",
          "set `independent_p = TRUE"
        )
      )
    }
  }

  if (independent_beta == "guess" && is.null(sigma_beta) && all(pif1@beta == pif2@beta) && all(pif1@sigma_beta == pif2@sigma_beta)) {
    independent_beta <- FALSE
    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "The relative risk beta parameters for the potential impact fractions appear ",
          "to be the same. If they are, set `independent_beta = FALSE`. Otherwise ",
          "set `independent_beta = TRUE"
        )
      )
    }
  }

  # If p's are the same set sigma_p as the covariance, otherwise assume independence
  if (is.null(sigma_p) && (independent_p == "guess" || independent_p)) {
    sigma_p <- matrix(0, ncol = length(pif1@p), nrow = length(pif1@p))
  } else if (is.null(sigma_p)) {
    sigma_p <- pif1@sigma_p
  }

  # If beta's are the same set sigma_theta as the covariance, otherwise assume independence
  if (is.null(sigma_beta) && (independent_beta == "guess" || independent_beta)) {
    sigma_beta <- matrix(0, ncol = length(pif1@beta), nrow = length(pif1@beta))
  } else if (is.null(sigma_beta)) {
    sigma_beta <- pif1@sigma_beta
  }

  pif_covariance(
    p1 = pif1@p,
    p2 = pif2@p,
    p1_cft = pif1@p_cft,
    p2_cft = pif2@p_cft,
    rr1 = pif1@rr,
    rr2 = pif2@rr,
    mu_obs1 = pif1@mu_obs,
    mu_obs2 = pif2@mu_obs,
    mu_cft1 = pif1@mu_cft,
    mu_cft2 = pif2@mu_cft,
    sigma_p = sigma_p,
    sigma_beta = sigma_beta,
    rr_link_deriv_vals1 = pif1@rr_link_deriv_vals,
    rr_link_deriv_vals2 = pif2@rr_link_deriv_vals,
    link_deriv_1 = pif1@link_deriv_vals,
    link_deriv_2 = pif2@link_deriv_vals
  )
}

#' Apply a function to a property of a pif_class or to the first
#' pif_class available in a pif_list_class
#'
#' @param Either a `pif_class` or a `pif_list_class`
#'
#' @keywords internal
pif_list_class_apply_1st <- function(x, fun, property){

  if (S7::S7_inherits(x, pif_class)){
    return(
      fun(S7::prop(x, property))
    )
  }

  if (S7::S7_inherits(x, pif_list_class)){
    return(
      pif_list_class_apply_1st(x[[1]], fun, property)
    )
  }

  cli::cli_abort("Invalid class for `x` should be a `pif_class` or `pif_list_class`")
}


#' Create a list of `pif_list_class` or `pif_list` objects:
#'
#' @param x A `pif_list_class` or `pif_list` objects
#' @param ... Additional `pif_list_class` or `pif_list` objects
#'
#' @return A `pif_list_class` containing all the objects
#' @keywords internal
pif_list <- function(x, ...){
  list_of_pifs <- append(list(x), list(...))
  # if (is.null(q)){
  #   q <- rep(1, length(list_of_pifs))
  # }
  pif_list_class(
    list_of_pifs#,
    #q = q
  )
}

#' Get the type of the fraction
#'
#' Obtain whether a fraction is a potential impact fraction or a
#' population attributable fraction
#'
#' @param x A `pif_class` object
#'
#' @return Either `PIF` or `PAF` depending on the object
fraction_type <- function(x){
  x@type
}
