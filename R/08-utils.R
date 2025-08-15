#' Compute the mean under the observed prevalence
#'
#' @inheritParams derivatives
#'
#' @return The mean relative risk under the observed prevalence
#'
#' @keywords internal
mu_obs_fun <- function(p, rr) {
  if (any(is.na(p)) || any(is.na(rr))){
    cli::cli_abort(
      "One of the following parameters: `p`, `p_cft`, `rr` has missing values"
    )
  }

  if (length(p) < 1 || length(rr) < 1){
    cli::cli_abort(
      "Not enough values passed to `p`, `p_cft` or `rr` (length < 0)!"
    )
  }

  if (sum(p) > 1){
    cli::cli_abort(
      "Probabilities cannot sum > 1. They sum to {sum(p)}"
    )
  }

  if (length(p) != length(rr)){
    cli::cli_abort(
      "p and r have different lengths"
    )
  }
  if (any(p < 0) || any(p > 1)){
    cli::cli_abort(
      "Invalid probability has value {p[c(which(p < 0), which(p > 1))]} outside [0,1]."
    )
  }
  if (any(rr <= 0)){
    cli::cli_abort(
      "Invalid relative risk has value {rr[which(rr <= 0)]} which is <= 0"
    )
  }

  as.numeric(1 + p %*% (rr - 1))
}

#' Compute the mean under the counterfactual prevalence
#'
#' @inheritParams derivatives
#'
#' @return The mean relative risk under the counterfactual prevalence
#'
#' @keywords internal
mu_cft_fun <- function(p_cft, rr) {
  mu_obs_fun(p_cft, rr)
}

#' Compute the potential impact fraction
#'
#' @inheritParams derivatives
#'
#' @return A potential impact fraction (numeric)
#'
#' @keywords internal
pif_fun <- function(p, p_cft, rr) {
  # Calculate the mean p and rr
  mu_obs <- mu_obs_fun(p, rr)
  mu_cft <- mu_cft_fun(p_cft, rr)

  # Calculate the potential impact fraction
  pif_fun2(mu_obs, mu_cft)
}

#' Compute the potential impact fraction
#'
#' @inheritParams derivatives
#'
#' @return A potential impact fraction (numeric)
#'
#' @keywords internal
pif_fun2 <- function(mu_obs, mu_cft) {
  # Calculate the potential impact fraction
  if (is.na(mu_obs) || is.na(mu_cft)){
    cli::cli_abort(
      "Missing input for `mu_obs` or `mu_cft`"
    )
  }
  # Calculate the potential impact fraction
  if (length(mu_obs) != 1 || length(mu_cft) != 1){
    cli::cli_abort(
      "Variables `mu_obs` and `mu_cft` have to be of length 1."
    )
  }
  1 - (mu_cft / mu_obs)
}

#' Confidence interval for a potential impact fraction
#'
#' @param link_vals Values of the link function evaluated at `pif` (i.e. `link(pif)`)
#' @param link_variance Link_variance estimate for the linked potential impact
#' fraction (i.e. for `link(pif)`)
#' @param conf_level Confidence level for the interval
#' @param link_inv Inverse of the link function used to compute `link_vals` and
#' `link_variance`.
#'
#' @return A vector with the lower and upper bounds
#'
#' @keywords internal
pif_atomic_ci <- function(link_vals, link_variance, conf_level, link_inv) {
  if (is.na(link_vals) || is.na(link_variance)){
    cli::cli_abort(
      "Missing values in confidence interval"
    )
  }

  if (length(link_vals) != 1 || length(link_variance) != 1){
    cli::cli_abort(
      "Variables `link_vals` and `link_variance` should be of length 1."
    )
  }

  b1 <- link_inv(link_vals - stats::qnorm((1 - conf_level)/2) * sqrt(link_variance))
  b2 <- link_inv(link_vals + stats::qnorm((1 - conf_level)/2) * sqrt(link_variance))
  sort(c(b1, b2))
}

#' Apply a function to a property of a pif_class or to the first
#' pif_class available in a pif_total_class
#'
#' @param x Either a `pif_class` or a `pif_total_class`
#' @param fun A function to apply to the first element of a `pif_total_class`
#' @param property A property of interest from the `pif_class` to extract
#'
#' @keywords internal
pif_class_apply_1st <- function(x, fun, property){

  if (S7::S7_inherits(x, pif_class)){
    return(
      fun(S7::prop(x, property))
    )
  }

  if (S7::S7_inherits(x, pif_total_class)){
    return(
      pif_class_apply_1st(x@pif_list[[1]], fun, property)
    )
  }

  cli::cli_abort(
    "Invalid class for `x` should be a `pif_class` or `pif_class_apply_1st`"
  )
}

#' Get the values of the derivative of link
#'
#' Obtain the values of the derivative of the link function of
#' a potential impact fraction (PIF) or a population attributable
#' fraction (PAF) at `pif`
#'
#' @param x A `pif_class` object
#'
#' @return A number indicating the derivative of `link(pif)`
#'
#' @examples
#' #A potential impact fraction
#' pif1 <- pif(p = 0.2, p_cft = 0.1, beta = 1.2, quiet = TRUE)
#' link_deriv_vals(pif1)
#'
#' @export
link_deriv_vals <- function(x){
  x@link_deriv_vals
}

#' Get the type of the fraction
#'
#' Obtain whether a fraction is a potential impact fraction (PIF) or a
#' population attributable fraction (PAF)
#'
#' @param x A `pif_class` object
#'
#' @return A character either `PIF` or `PAF` depending on the object
#'
#' @examples
#' #A potential impact fraction
#' pif1 <- pif(p = 0.2, p_cft = 0.1, beta = 1.2, quiet = TRUE)
#' fraction_type(pif1)
#'
#' #A population attributable fraction
#' paf1 <- paf(p = 0.2, beta = 1.2, quiet = TRUE)
#' fraction_type(paf1)
#'
#' @export
fraction_type <- function(x){
  if (S7::S7_inherits(x, pif_class)){
    return(x@type)
  }
  cli::cli_abort(
    "Invalid class for `x` should be a `pif_class`"
  )
}

#' Change the link for the fraction's variance calculation
#'
#' Change the link function for the potential impact fraction
#' or population attributable fraction to a different link.
#'
#' @param x A `pif_class` object
#' @inheritParams pifpaf
#'
#' @return A `pif_class` object with a different `link`.
#'
#' @examples
#' #A potential impact fraction
#' pif1 <- pif(p = 0.2, p_cft = 0.1, beta = 1.2, var_p = 0.01,
#'   var_beta = 0.2)
#' pif1
#'
#' #Now change the pif to logit to control the negatives
#' pif1_logit <- change_link(pif1, link = "logit")
#' pif1_logit
#'
#' @export
change_link <- function(x, link = "identity", link_inv = NULL, link_deriv = NULL){

  if (!S7::S7_inherits(x, pif_class)){
    cli::cli_abort(
      "Invalid class for `x` should be a `pif_class`"
    )
  }

  if (S7::S7_inherits(x, pif_ensemble_class)){
    cli::cli_abort(
      "Cannot change the link of an ensemble (`pif_ensemble_class`)"
    )
  }

  # Check that is link is given as character then link inv and link_deriv are not specified
  link_name <- link
  check_links(link = link, link_inv = link_inv, link_deriv = link_deriv)

  # Get the inverse function
  if (is.null(link_inv) & is.character(link)) {
    link_inv <- parse_inv_link(link)
  }

  if (is.null(link_deriv) & is.character(link)){
    link_deriv <- parse_deriv_link(link)
  }

  # Get the functions
  link    <- parse_link(link)

  if (!is.function(link_deriv) && is.null(link_deriv)) {
    link_deriv <- Deriv::Deriv(link)
  }


  x@link       <- link
  x@link_inv   <- link_inv
  x@link_deriv <- link_deriv

  if (is.character(link_name) && link_name == "logit" && coef(x) <= 0){
    cli::cli_alert_danger(
      paste0(
        "Value for {fraction_type(x)} = {round(coef(x),2)} <= 0. ",
        "Change link to a different value as `logit` is only ",
        "valid for strictly positive {fraction_type(x)}s."
      )
    )
  }

  return(x)
}
