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

#' Get the link_variance
#' @rdname getters
get_link_variance <- function(self){
  (self@link_deriv_vals)^2 * self@variance
}

#' Get the link_variance
#' @rdname getters
get_variance_atomic <- function(self) {
  from_parameters_pif_variance(
    p = self@p, p_cft = self@p_cft, rr = self@rr,
    rr_link_deriv_vals = self@rr_link_deriv_vals,
    mu_obs = self@mu_obs,
    mu_cft = self@mu_cft, sigma_p = self@sigma_p,
    sigma_beta = self@sigma_beta,
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

#' Get the coefficients
#' @rdname getters
get_total_coefs <- function(self){
  sapply(self@pif_list, coef)
}

#' Get the coefficients
#' @rdname getters
get_ensemble_coefs <- function(self){
  get_total_coefs(self)
}

#' Get the ensemble coefficients
get_ensemble_pif <- function(self){
  exp(sum(log(self@coefs)))
}

#' Get the coefficients
#' @rdname getters
get_total_pif <- function(self){
  as.numeric(t(self@weights) %*% self@coefs)
}

#' Get types
#' @rdname getters
get_total_type <- function(self){
  pif_types <- sapply(self@pif_list, fraction_type)
  ifelse(any(pif_types == "PIF"), "PIF", "PAF")
}

#' Get types
#' @rdname getters
get_ensemble_type <- function(self){
  get_total_type(self)
}

#' Get the covariance of the summands of pif total
#' @rdname getters
get_covariance_total <- function(self){
  npifs <- length(self@pif_list)
  if (npifs > 1){
    cov_mat <- matrix(0, ncol = npifs, nrow = npifs)
    for (i in 1:(npifs - 1)){
      for (j in (i + 1):npifs){
        cov_mat[i,j] <- cov_total_pif(self@pif_list[[i]], self@pif_list[[j]])
      }
    }
    cov_mat <- cov_mat + t(cov_mat) + diag(sapply(self@pif_list, var))
    return(cov_mat)
  } else {
    return(
      var(self@pif_list[[1]])
    )
  }
}

#' Get the covariance of the summands of pif total
#' @rdname getters
get_ensemble_covariance <- function(self){
  npifs <- length(self@pif_list)
  if (npifs > 1){
    cov_mat <- matrix(0, ncol = npifs, nrow = npifs)
    for (i in 1:(npifs - 1)){
      for (j in (i + 1):npifs){
        #In ensemble the covariance is given by ln(1 - pif[i])'*ln(1 - pif[j])'*cov(pif[i],pif[j])
        g_i_prime <- link_deriv_vals(self@pif_list[[i]])
        g_j_prime <- link_deriv_vals(self@pif_list[[j]])
        cov_mat[i,j] <- g_i_prime*g_j_prime*cov_total_pif(self@pif_list[[i]], self@pif_list[[j]])
      }
    }
    cov_mat <- cov_mat + t(cov_mat) + diag(sapply(self@pif_list, var)*(sapply(self@pif_list, link_deriv_vals)^2))
    return(cov_mat)
  } else {
    return(
      var(self@pif_list[[1]])*(link_deriv_vals(self@pif_list[[1]])^2)
    )
  }
}

#' Get the variance of pif total
#' @rdname getters
get_variance_total <- function(self){
  as.numeric(
    t(self@weights) %*% self@covariance %*% self@weights +
    t(self@coefs) %*% self@sigma_weights %*% self@coefs +
    sum(diag(self@covariance*self@sigma_weights))
  )
}

#' Get the variance of pif total
#' @rdname getters
get_ensemble_variance <- function(self){
  as.numeric(
    sum(self@covariance)
  )
}
