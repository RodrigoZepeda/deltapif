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
#'
#' @return The partial derivative (usually a vector)
#'
#' @section Formulas:
#' The partial derivative of `PIF` with respect to `p` is:
#' \deqn{
#' \dfrac{\partial \textrm{PIF}}{\partial p} =
#' \dfrac{\mu^{\text{cft}}}{\big(\mu^{\text{obs}}\big)^2} \cdot \big( \text{RR}(\beta) - 1)
#' }
#' The partial derivative of `PIF` with respect to `beta` is:
#' \deqn{
#' \dfrac{\partial \textrm{PIF}}{\partial \beta} =
#' \Bigg(\dfrac{
#'  \mu^{\text{obs}} \cdot p_{*} - \mu^{\text{cft}} \cdot p
#' }{
#'  \Big( \mu^{\text{obs}}\Big)^2
#' }\Bigg)\odot \text{RR}'(\beta)
#'
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
deriv_pif_p <- function(p, p_cft, rr, mu_obs = NULL, mu_cft = NULL) {
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
  (mu_cft / (mu_obs)^2) * (rr - 1)
}

#' Partial derivative of PIF respect to beta
#' @rdname derivatives
deriv_pif_beta <- function(p, p_cft, rr, rr_link_deriv_vals, mu_obs = NULL,
                           mu_cft = NULL) {
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
  ((mu_obs * p_cft - mu_cft * p) / (mu_obs^2)) * rr_link_deriv_vals
}
