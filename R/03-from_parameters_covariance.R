#' Covariance component with respect to p
#'
#' Calculates the covariance `p` or `beta` component between two
#' potential impact fractions or the complete covariance.
#'
#' @param p1 Observed proportion exposed in the first potential impact fraction
#' @param p2 Observed proportion exposed in the first potential impact fraction
#' @param p1_cft Counterfactual proportion exposed in the first potential
#' impact fraction
#' @param p2_cft Counterfactual proportion exposed in the first potential
#' impact fraction
#' @param rr1 Relative risk in the first potential impact fraction
#' @param rr2 Relative risk in the second potential impact fraction
#' @param mu_obs1 Average relative risk in the population of the first
#' potential impact fraction.
#' @param mu_obs2 Average relative risk in the population of the second
#' potential impact fraction.
#' @param mu_cft1 Average relative risk in the counterfactual population of the
#' first potential impact fraction.
#' @param mu_cft2 Average relative risk in the counterfactual of the second
#' potential impact fraction.
#' @param var_p Covariance matrix with the entry `var_p[i,j]` corresponding
#' to the covariance between `p1[i]` and `p2[j]`.
#' @param upper_bound Whether the variance should be calculated or an upper
#' bound for the variance (assuming perfect correlation).
#' @param rr_link_deriv_vals1 Values for the derivative of the relative
#' risk function evaluated at `theta` (the relative risk parameter).
#' @param rr_link_deriv_vals2 Values for the derivative of the relative
#' risk function evaluated at `theta` (the relative risk parameter).
#'
#' @inheritParams derivatives
#' @inheritParams pifpaf
#' @inheritParams classes
#'
#' @return The covariance component for either `p` or `beta`
#'
#' @section Formulas:
#' The following represents the covariance components.
#' For `p`:
#'
#' \deqn{
#' \text{CC}_{p} = \dfrac{
#'  \partial \textrm{PIF}_i}{\partial p_i}^{\top}
#'  \textrm{covariance}\big( \hat{p}_i, \hat{p}_j\big)
#'  \dfrac{\partial \textrm{PIF}_j}{\partial p_j}
#' }
#'
#' and for `beta`:
#'
#' \deqn{
#' \text{CC}_{\theta} =  \dfrac{
#' \partial \textrm{PIF}_i}{\partial \theta_i}^{\top}
#' \textrm{covariance}\big( \hat{\theta}_i, \hat{\theta}_j\big)
#' \dfrac{\partial \textrm{PIF}_j}{\partial \theta_j}
#' }
#'
#' @name covariance_from_parameters
#'
#' @note This estimation does not require a `pif_class` object
#' but instead is built from parameters.
#'
#' @seealso [derivatives] for the definition of the derivatives involved
#' in the formulas.
NULL

#' @rdname covariance_from_parameters
#' @export
from_parameters_covariance_p_component <- function(
    p1, p2, p1_cft, p2_cft, rr1, rr2, mu_obs1, mu_obs2, mu_cft1, mu_cft2,
    var_p, upper_bound) {

  # Get the derivatives
  vcp_deriv1 <- deriv_pif_p(
    p = p1, p_cft = p1_cft, rr = rr1, mu_obs = mu_obs1,
    mu_cft = mu_cft1
  )
  vcp_deriv2 <- deriv_pif_p(
    p = p2, p_cft = p2_cft, rr = rr2, mu_obs = mu_obs2,
    mu_cft = mu_cft2
  )

  # Whether we are working with an upper bound for the link_variance
  if (upper_bound) {
    vcp_deriv1 <- abs(vcp_deriv1)
    vcp_deriv2 <- abs(vcp_deriv2)
  }

  t(vcp_deriv1) %*% var_p %*% vcp_deriv2
}

#' @rdname covariance_from_parameters
#' @export
from_parameters_covariance_beta_component <- function(
    p1, p2, p1_cft, p2_cft, rr1, rr2, rr_link_deriv_vals1, rr_link_deriv_vals2,
    mu_obs1, mu_obs2, mu_cft1, mu_cft2, var_beta, upper_bound) {
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

  # Whether we are working with an upper bound for the link_variance
  if (upper_bound) {
    vcbeta_deriv1 <- abs(vcbeta_deriv1)
    vcbeta_deriv2 <- abs(vcbeta_deriv2)
  }

  t(vcbeta_deriv1) %*% var_beta %*% vcbeta_deriv2
}

#' @rdname covariance_from_parameters
#' @export
from_parameters_pif_covariance <- function(
    p1, p2, p1_cft, p2_cft, rr1, rr2, rr_link_deriv_vals1, rr_link_deriv_vals2,
    mu_obs1, mu_obs2, mu_cft1, mu_cft2, var_p, var_beta,
    upper_bound_p = FALSE, upper_bound_beta = FALSE) {

  p_component <- from_parameters_covariance_p_component(
    p1 = p1,
    p2 = p2,
    p1_cft = p1_cft,
    p2_cft = p2_cft,
    rr1 = rr1, rr2 = rr2,
    mu_obs1 = mu_obs1,
    mu_obs2 = mu_obs2,
    mu_cft1 = mu_cft1,
    mu_cft2 = mu_cft2,
    var_p = var_p,
    upper_bound = upper_bound_p
  )

  beta_component <- from_parameters_covariance_beta_component(
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
    var_beta = var_beta,
    upper_bound = upper_bound_beta
  )

  as.numeric((p_component + beta_component))
}

#' @rdname covariance_from_parameters
#' @export
from_parameters_pif_variance <- function(
    p, p_cft, rr, rr_link_deriv_vals, mu_obs, mu_cft, var_p,
    var_beta, upper_bound_p = FALSE, upper_bound_beta = FALSE) {

  from_parameters_pif_covariance(
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
    var_p = var_p,
    var_beta = var_beta,
    upper_bound_p = upper_bound_p,
    upper_bound_beta = upper_bound_beta
  )
}





