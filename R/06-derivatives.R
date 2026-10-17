#' Partial derivatives of PIF
#'
#' Calculates the partial derivatives of a potential impact fraction
#' with respect to the parameters `p` or `beta`.
#'
#' @inheritParams pifpaf
#' @param rr The relative risk for each of the exposure levels.
#' @param rr_link_deriv_vals The derivative of the relative risk function `g`
#' with respect to the parameter `beta` evaluated at `beta`.
#' @param mu_obs The average value of the relative risk in the observed population.
#' @param mu_cft The average value of the counterfactual relative risk in the population.
#' @param p_cft_jacobian (Optional) Jacobian of the counterfactual prevalence
#' `p_cft` with respect to `p` when the counterfactual is a function of the
#' observed prevalence (for example a relative reduction `p_cft = 0.85 * p`).
#' Entry `[j, i]` is the derivative of `p_cft[j]` with respect to `p[i]`.
#' `NULL` (default) assumes `p_cft` is a fixed target.
#'
#' @return The partial derivative (usually a vector)
#'
#' @section Formulas:
#' The partial derivative of `PIF` with respect to `p` is:
#' \deqn{
#' \dfrac{\partial \textrm{PIF}}{\partial p} =
#' \dfrac{\mu^{\text{cft}}}{\big(\mu^{\text{obs}}\big)^2} \cdot \big( \text{RR}(\beta) - 1\big)
#' }
#' If the counterfactual is a function of the prevalence, \eqn{p_{*} = g(p)},
#' the total derivative adds the term
#' \eqn{- J_g^{\top} \big( \text{RR}(\beta) - 1\big) / \mu^{\text{obs}}}
#' where \eqn{J_g} is the Jacobian of \eqn{g} evaluated at \eqn{p}.
#' The partial derivative of `PIF` with respect to `beta` is:
#' \deqn{
#' \dfrac{\partial \textrm{PIF}}{\partial \beta} =
#' \Bigg(
#'  \dfrac{
#'    \mu^{\text{cft}} \cdot p - \mu^{\text{obs}} \cdot p_{*}
#'  }{
#'    \Big( \mu^{\text{obs}}\Big)^2
#'  }
#' \Bigg)\odot \text{RR}^{\prime}(\beta)
#' }
#' with \eqn{\odot} representing the Hadamard (elementwise) product.
#'
#' @note As `p` and `beta` are usually vectors these are vector-valued
#' derivatives.
#'
#' @name derivatives
NULL

#' @rdname derivatives
#' @export
deriv_pif_p <- function(p, p_cft, rr, mu_obs = NULL, mu_cft = NULL,
                        p_cft_jacobian = NULL) {

  if (any(p < 0)){
    cli::cli_abort(
      "Invalid probability p < 0"
    )
  }

  if (any(p_cft < 0)){
    cli::cli_abort(
      "Invalid counterfactual probability p_cft < 0"
    )
  }

  if (is.null(mu_obs)) {
    mu_obs <- mu_obs_fun(p, rr)
  }

  if (is.null(mu_cft)) {
    mu_cft <- mu_cft_fun(p_cft, rr)
  }

  if (is.nan(mu_obs) || is.na(mu_obs) || is.infinite(mu_obs) || mu_obs == 0){
    cli::cli_abort(
      "Invalid value for `mu_obs` {mu_obs}"
    )
  }

  if (is.nan(mu_cft) || is.na(mu_cft) || is.infinite(mu_cft)){
    cli::cli_abort(
      "Invalid value for `mu_cft` {mu_cft}"
    )
  }

  # The derivative of pif
  deriv_p <- (mu_cft / (mu_obs)^2) * (rr - 1)

  # Add the change in the counterfactual if it is a function of p
  if (!is.null(p_cft_jacobian)){
    deriv_p <- deriv_p - as.numeric(t(p_cft_jacobian) %*% (rr - 1)) / mu_obs
  }

  deriv_p
}

#' Partial derivative of PIF respect to beta
#' @rdname derivatives
deriv_pif_beta <- function(p, p_cft, rr, rr_link_deriv_vals, mu_obs = NULL,
                           mu_cft = NULL) {

  if (any(p < 0)){
    cli::cli_abort(
      "Invalid probability p < 0"
    )
  }

  if (any(p_cft < 0)){
    cli::cli_abort(
      "Invalid counterfactual probability p_cft < 0"
    )
  }


  if (is.null(mu_obs)) {
    mu_obs <- mu_obs_fun(p, rr)
  }

  if (is.null(mu_cft)) {
    mu_cft <- mu_cft_fun(p_cft, rr)
  }

  if (length(rr) != length(rr_link_deriv_vals)){
    cli::cli_abort(
      "rr and rr_link_deriv_vals have different lengths"
    )
  }

  if (any(is.na(rr_link_deriv_vals))){
    cli::cli_abort(
      "Missing values in `rr_link_deriv_vals`"
    )
  }

  # The derivative of pif
  ((mu_cft * p - mu_obs * p_cft) / (mu_obs^2)) * rr_link_deriv_vals
}
