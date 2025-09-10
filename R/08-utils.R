#' Compute the mean under the observed prevalence
#'
#' Calculates the mean average relative risk of a population
#' under the observed prevalence:
#' \deqn{
#' E[\text{RR}] = \sum\limits_{i=1}^{n} p_i \cdot \text{RR}_i
#' }
#'
#' @inheritParams derivatives
#'
#' @return The mean relative risk under the observed prevalence
#'
#' @examples
#' \dontrun{
#' #Consider a popualtion with 3 exposure categories with
#' #relative risks 1.0, 1.1, 1.3 and prevalences 0.6, 0.3, 0.1
#' #Notice that the reference relative risk (1) is not added
#' pval  <- c(0.3, 0.1) #Prevalences
#' rrval <- c(1.1, 1.3) #Risks
#' mu_obs_fun(p = pval, rr = rrval) #Average relative risk
#' }
#'
#' @keywords internal
#' @seealso [mu_cft_fun()]
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
#' Calculates the mean average relative risk of a population
#' under the counterfactual prevalence:
#' \deqn{
#' E[\text{RR}] = \sum\limits_{i=1}^{n} p_i^{\text{cft}} \cdot \text{RR}_i
#' }
#'
#' @inheritParams derivatives
#'
#' @return The mean relative risk under the counterfactual prevalence
#'
#' @examples
#' \dontrun{
#' #Consider a popualtion with 3 exposure categories with
#' #relative risks 1.0, 1.1, 1.3 and prevalences 0.6, 0.3, 0.1
#' #Notice that the reference relative risk (1) is not added
#' pval  <- c(0.4, 0.0) #Prevalences on counterfactual scenarios
#' rrval <- c(1.1, 1.3) #Risks
#' mu_cft_fun(p_cft = pval, rr = rrval) #Average relative risk
#' }
#'
#' @keywords internal
#' @seealso [mu_obs_fun()]
mu_cft_fun <- function(p_cft, rr) {
  mu_obs_fun(p_cft, rr)
}

#' Compute the potential impact fraction
#'
#' Calculates the potential impact fraction following the formula:
#' \deqn{
#' \text{PIF} = \dfrac{
#' \sum\limits_{i=1}^n p_i \cdot \text{RR}_i - \sum\limits_{i=1}^n p_i^{\text{cft}} \cdot \text{RR}_i
#' }{
#' \sum\limits_{i=1}^n p_i \cdot \text{RR}_i
#' }
#' }
#'
#' @inheritParams derivatives
#'
#' @return A potential impact fraction (numeric)
#'
#' @examples
#' \dontrun{
#' #This example comes from Levin 1953
#' #Relative risk of lung cancer given smoking was 3.6
#' #Proportion of individuals smoking where 49.9%
#' #Counterfactual is 0% smoking
#' pif_fun(0.499, 0.0, 3.6)
#'
#' #You can also use multiple exposure categories. As an example
#' #These are the relative risks for age-group 25-59 for each BMI level:
#' rr <- c("<18.5" = 1.38, "25 to <30" = 0.83,
#'   "30 to <35" = 1.20, ">=35" = 1.83)
#'
#' #While the prevalences are:
#' p <- c("<18.5" = 1.9, "25 to <30" = 34.8,
#'   "30 to <35" = 17.3, ">=35" = 13.3) / 100
#'
#' #A counterfactual of reducing the prevalence of >=35 by half and
#' #having them be on the "30 to <35" category instead
#' p_cft <- c("<18.5" = 1.9, "25 to <30" = 34.8,
#'   "30 to <35" = 17.3 + 13.3/2, ">=35" = 13.3/2) / 100
#'
#' pif_fun(p = p, p_cft = p_cft, rr = rr)
#' }
#' @keywords internal
#' @seealso [mu_cft_fun()], [mu_obs_fun()], and [pif_fun2()] to calculate from
#' the average relative risks. For the main function in the package see [pif()]
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
#' @examples
#' \dontrun{
#' #This example comes from Levin 1953
#' #Relative risk of lung cancer given smoking was 3.6
#' #Proportion of individuals smoking where 49.9%
#' #Counterfactual is 0% smoking
#' mu_obs <- mu_obs_fun(p = 0.499, rr = 3.6)
#' mu_cft <- mu_cft_fun(p_cft = 0.0, rr = 3.6)
#' pif_fun2(mu_obs, mu_cft)
#'
#' #You can also use multiple exposure categories. As an example
#' #These are the relative risks for age-group 25-59 for each BMI level:
#' rr <- c("<18.5" = 1.38, "25 to <30" = 0.83,
#'   "30 to <35" = 1.20, ">=35" = 1.83)
#'
#' #While the prevalences are:
#' p <- c("<18.5" = 1.9, "25 to <30" = 34.8,
#'   "30 to <35" = 17.3, ">=35" = 13.3) / 100
#'
#' #A counterfactual of reducing the prevalence of >=35 by half and
#' #having them be on the "30 to <35" category instead
#' p_cft <- c("<18.5" = 1.9, "25 to <30" = 34.8,
#'   "30 to <35" = 17.3 + 13.3/2, ">=35" = 13.3/2) / 100
#'
#' mu_obs <- mu_obs_fun(p = p, rr = rr)
#' mu_cft <- mu_cft_fun(p_cft = p_cft, rr = rr)
#' pif_fun2(mu_obs, mu_cft)
#' }
#' @keywords internal
#' @seealso [mu_cft_fun()], [mu_obs_fun()], and [pif_fun()] to calculate from
#' the prevalence and relative risks. For the main function in the package
#' see [pif()]
pif_fun2 <- function(mu_obs, mu_cft) {

  if (length(mu_obs) != 1 || length(mu_cft) != 1){
    cli::cli_abort(
      "Invalid length > 1 for either mu_obs or mu_cft"
    )
  }
  # Calculate the potential impact fraction
  if (is.na(mu_obs) || is.na(mu_cft)){
    cli::cli_abort(
      "Missing input for `mu_obs` or `mu_cft`"
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
#' @return A vector with the lower and upper bounds of the confidence interval
#'
#' @keywords internal
#' @seealso [confint]
pif_atomic_ci <- function(link_vals, link_variance, conf_level, link_inv) {

  if (length(link_vals) != 1 || length(link_variance) != 1){
    cli::cli_abort(
      "Variables `link_vals` and `link_variance` should be of length 1."
    )
  }

  if (is.na(link_vals) || is.na(link_variance)){
    cli::cli_abort(
      "Missing values in confidence interval"
    )
  }

  if (conf_level > 1 || conf_level < 0){
    cli::cli_abort(
      "Invalid confidence level not in (0,1)"
    )
  }

  b1 <- link_inv(link_vals - stats::qnorm((1 - conf_level)/2) * sqrt(link_variance))
  b2 <- link_inv(link_vals + stats::qnorm((1 - conf_level)/2) * sqrt(link_variance))
  sort(c(b1, b2))
}


#' Get the values of the derivative of link
#'
#' Obtain the values of the derivative of the link function of
#' a potential impact fraction (PIF) or a population attributable
#' fraction (PAF) at `pif`
#'
#' @param x A `pif_class` object
#'
#' @return A number indicating the derivative of `link()` evaluated at `pif`.
#'
#' @examples
#' \dontrun{
#' #Create a pif object
#' pif_obj <- pif(p = 0.499, beta = log(3.6), p_cft = 0.499/2, var_p = 0.001,
#'   var_beta = 0.1, link = "logit", quiet = TRUE)
#'
#' #Obtain the value of the link at the derivative
#' link_deriv_vals(pif_obj)
#'
#' #This is the same as the derivative of logit evaluated at 0.2823
#' deriv_logit(coef(pif_obj))
#' }
#'
#'
#' @keywords internal
link_deriv_vals <- function(x){
  if (S7::S7_inherits(x, pif_class)){
    return(x@link_deriv_vals)
  }
  cli::cli_abort(
    "Invalid class for `x` should be a `pif_class`"
  )
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

#' Change the link
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
    cli::cli_warn(
      paste0(
        "Value for {fraction_type(x)} = {round(coef(x),2)} <= 0. ",
        "Change link to a different value as `logit` is only ",
        "valid for strictly positive {fraction_type(x)}s."
      )
    )
  }

  return(x)
}

#' Names of a PIF's components
#'
#' Returns a character vector of the names of all of the fractions
#' that make up a `pif` object. This is particularly useful for `pif_total`
#' and `pif_ensemble` to obtain the names that build them up.
#'
#' @param pif A potential impact fraction of either `pif_atomic_class`
#' or `pif_global_ensemble_class`
#'
#'
#' @return A character vector with the names of all the fractions
#' that make up the `pif`.
#'
#' @examples
#' paf_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001,
#'         label = "Women lead")
#' paf_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001,
#'         label = "Women radiation")
#' paf_women      <- paf_ensemble(paf_lead_women, paf_rad_women,
#'         label = "Women")
#' paf_lead_men   <- paf(0.30, 2.2, quiet = TRUE, var_p = 0.001,
#'         label = "Men lead")
#' paf_rad_men    <- paf(0.10, 1.2, quiet = TRUE, var_p = 0.001,
#'         label = "Men radiation")
#' paf_men        <- paf_ensemble(paf_lead_men, paf_rad_men,
#'         label = "Men")
#' paf_tot        <- paf_total(paf_men, paf_women, weights = c(0.49, 0.51),
#'         label = "Population")
#'
#' #For a single PIF return the names
#' flatten_names(paf_lead_women)
#'
#' #For an ensemble return the ones that make them up
#' flatten_names(paf_women)
#'
#' #For totals return the ones that make them up
#' flatten_names(paf_tot)
#'
#' @export
flatten_names <- function(pif){

  if (S7::S7_inherits(pif, pif_atomic_class)){
    return(pif@label)
  } else if (S7::S7_inherits(pif, pif_global_ensemble_class)){
    labs <- c(pif@label, unlist(sapply(pif@pif_list, flatten_names)))
    names(labs) <- NULL
    return(labs)
  } else {
    cli::cli_abort(
      paste0(
        "The `flatten_names` function is only available for ",
        "`pif_atomic_class` and `pif_global_ensemble_class` objects."
      )
    )
  }
}
