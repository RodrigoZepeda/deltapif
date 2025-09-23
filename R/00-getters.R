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

#' Get the cases transformed by the link function
#' @rdname getters
get_link_vals_cases <- function(self) {
  self@link(self@cases)
}

#' Get the derivative of the link function evaluated at the pif
#' @rdname getters
get_link_deriv_vals_cases <- function(self) {
  self@link_deriv(self@cases)
}

#' Get the cases = paf * overall_cases
#' @rdname getters
get_cases <- function(self){
  self@overall_cases * self@pif_obj@pif
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
  deriv_vals <- rep(NA_real_, length(self@beta))
  for (j in 1:length(deriv_vals)) {
    deriv_vals[j] <- self@rr_link_deriv(self@beta[j])
  }
  return(deriv_vals)
}

#' Get the link_variance
#' @rdname getters
get_link_variance <- function(self){
  (self@link_deriv_vals)^2 * self@variance
}

#' Get the variance of the cases
#' @rdname getters
get_variance_cases <- function(self){
  var_prod(x = self@overall_cases, y = self@pif_obj@pif,
           var_x = self@variance_cases, var_y = self@pif_obj@variance)
}

#' Get the link_variance
#' @rdname getters
get_variance_atomic <- function(self) {
  from_parameters_pif_variance(
    p = self@p, p_cft = self@p_cft, rr = self@rr,
    rr_link_deriv_vals = self@rr_link_deriv_vals,
    var_p = self@var_p,
    var_beta = self@var_beta,
    upper_bound_p = self@upper_bound_p,
    upper_bound_beta = self@upper_bound_beta
  )
}

#' Get the confidence interval
#' @rdname getters
get_ci <- function(self) {
  pif_atomic_ci(
    link_vals = self@link_vals, link_variance = self@link_variance,
    conf_level = self@conf_level, link_inv = self@link_inv
  )
}

#' Get each of the PIF_i
#' @rdname getters
get_ensemble_coefs <- function(self){
  sapply(self@pif_list, coef)
}

#' Get the sum of g(q_i PIF_i)
#' @rdname getters
get_sum_transformed_weighted_coefs <- function(self){
  sum(sapply(self@coefs*self@weights, self@pif_transform))
}

#' Get g^{-1} sum(g(q_i PIF_i))
#' @rdname getters
get_global_ensemble_pif <- function(self){
  self@pif_inverse_transform(self@sum_transformed_weighted_coefs)
}

#' Get the ensemble coefficients
#' @rdname getters
get_ensemble_pif <- function(self){
  1 - exp(sum(log(1 - self@weights*self@coefs)))
}

#' Get the coefficients
#' @rdname getters
get_total_pif <- function(self){
  as.numeric(t(self@weights) %*% self@coefs)
}

#' Get types
#' @rdname getters
get_ensemble_type <- function(self){
  pif_types <- sapply(self@pif_list, fraction_type)
  ifelse(any(pif_types == "PIF"), "PIF", "PAF")
}

#' Get the covariance of the summands of pif total
#' @rdname getters
get_covariance <- function(self){
  npifs <- length(self@pif_list)
  if (npifs > 1){
    cov_mat <- do.call(covariance,
                       list(x = self@pif_list[[1]], ... = self@pif_list[[2:npifs]]))
    return(cov_mat)
  } else {
    return(
      variance(self@pif_list[[1]])
    )
  }
}

#' Get the variance of pif total
#' @rdname getters
get_variance <- function(self){
  variance(self)
}
