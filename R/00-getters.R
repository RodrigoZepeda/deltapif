#' Getters
#'
#' Collection of functions for the potential impact fraction class
#' so that they work as getters of the properties. Each
#' function is constructed as `get_property`
#'
#' @param self A `pif_class` or `pif_atomic_class` object created with S7.
#'
#' @name getters
#'
#' @keywords internal
NULL

#' Get the alpha / 2 value for the confidence interval
#' @rdname getters
get_alpha_confint <- function(self) {
  (1 - self@conf_level) / 2
}

#' Get the relative risks
#' @rdname getters
get_rr <- function(self) {
  rrvals <- rep(NA_real_, length(self@beta))
  for (i in 1:length(self@beta)) {
    rrvals[i] <- self@rr_link(self@beta[i])
  }
  return(rrvals)
}

#' Get the mean of the relative risk under the observed prevalence
#' @rdname getters
get_mu_obs <- function(self) {
  mu_obs_fun(self@p, self@rr)
}

#' Get the mean of the relative risk under the counterfactual prevalence
#' @rdname getters
get_mu_cft <- function(self) {
  mu_cft_fun(self@p_cft, self@rr)
}

#' Get the potential impact fraction
#' @rdname getters
get_pif <- function(self) {
  pif_fun2(self@mu_obs, self@mu_cft)
}

#' Get the potential impact fraction transformed by the link function
#' @rdname getters
get_link_vals <- function(self) {
  pif_vals <- rep(NA_real_, length(self@pif))
  for (i in 1:length(self@pif)) {
    pif_vals[i] <- self@link(self@pif[i])
  }
  return(pif_vals)
}

#' Get the derivative of the link function evaluated at the pif
#' @rdname getters
get_link_deriv_vals <- function(self) {
  deriv_vals <- rep(NA_real_, length(self@pif))
  for (j in 1:length(deriv_vals)) {
    deriv_vals[j] <- self@link_deriv(self@pif[j])
  }
  return(deriv_vals)
}

#' Get the derivative of the rr-link function evaluated at beta
#' @rdname getters
get_rr_link_deriv_vals <- function(self) {
  deriv_vals <- rep(NA_real_, length(self@pif))
  for (j in 1:length(deriv_vals)) {
    deriv_vals[j] <- self@rr_link_deriv(self@beta[j])
  }
  return(deriv_vals)
}

#' Get the variance
#' @rdname getters
get_variance_atomic <- function(self) {
  pif_variance(
    p = self@p, p_cft = self@p_cft, rr = self@rr,
    rr_link_deriv_vals = self@rr_link_deriv_vals, mu_obs = self@mu_obs,
    mu_cft = self@mu_cft, sigma_p = self@sigma_p,
    sigma_beta = self@sigma_beta,
    link_deriv_vals = self@link_deriv_vals,
    upper_bound_p = self@sigma_p_upper_bound,
    upper_bound_beta = self@sigma_beta_upper_bound
  )
}

#' Get the confidence interval
#' @rdname getters
get_ci <- function(self) {
  pif_ci(
    link_vals = self@link_vals, variance = self@variance,
    alpha_confint = get_alpha_confint(self), link_inv = self@link_inv
  )
}

#' Get the coefficients
#' @rdname getters
get_total_coefs <- function(self){
  sapply(self@pif_list, coef)
}

#' Get the coefficients
#' @rdname getters
get_total_pif <- function(self){
  t(self@weights) %*% self@coefs
}

#' Get types
#' @rdname getters
get_total_type <- function(self){
  pif_types <- sapply(self@pif_list, fraction_type)
  ifelse(any(pif_types == "PIF"), "PIF", "PAF")
}

#' Get the covariance matrix between elements of the list
#' @rdname getters
get_total_covariance <- function(self){
  diag(0, nrow = length(self@pif_list))
}

#' Get the covariance matrix between elements of the list
#' @rdname getters
get_total_variance <- function(self){
  t(self@weights) %*% self@covariance %*% self@weights +
    t(self@coefs) %*% self@sigma_weights %*% self@coefs +
    sum(diag(self@covariance*self@sigma_weights))
}
